#!/usr/bin/env bash
# Annotate assembled receptor sequences against the pinned IMGT germline genes.
# TRUST4 assembles candidate receptors from RNA-seq reads. Its AIRR table names
# each reported sequence and its read support; neither count measures cells.
# parse_trust4_airr.R first removes records with a TRUST4-imputed CDR3, the
# junction region used to distinguish rearrangements.
#
# TRUST4's locus selects the receptor family: antibody chains IGH/IGK/IGL query
# the IgBLAST `ig` databases; T-cell receptor chains TRA/TRB/TRD/TRG query `tr`.
# IgBLAST assigns variable (V), diversity (D), joining (J) and constant (C)
# genes and reports junctions, alignments and productivity. Sequence IDs and read
# support stay linked to the TRUST4 observations. MakeDb converts the results
# to AIRR tables. CreateGermlines reconstructs IG germlines for later comparison
# of observed V-region bases with their inferred unmutated sequence.
#
# A V call naming several genes is an equal-score alignment tie, not several
# receptors. The params.igblast.num_alignments_* settings control how many
# candidate alignments IgBLAST can report. igblastn runs directly because
# Change-O's AssignGenes wrapper does not expose those limits. Other flags
# follow Change-O 1.3.4 runIgBLASTN. build_reannotated_airr.R records whether
# a best-score list reaches its limit before classifying V/J assignments;
# filter_airr.R applies common record QC after this wrapper finishes.
#
# References:
#   Change-O IgBLAST workflow, whose command this reproduces:
#     https://changeo.readthedocs.io/en/stable/examples/igblast.html
#   igblastn options and their defaults: `igblastn -help` in the pinned image.
#   IgBLAST 1.22.0 igblastn (Ye et al. 2013, Nucleic Acids Res,
#     doi:10.1093/nar/gkt382).
#   Change-O 1.3.4 MakeDb and CreateGermlines (Gupta et al. 2015,
#     Bioinformatics, doi:10.1093/bioinformatics/btv359).
#
# The reannotate_igblast rule supplies paths, the library label, thread count
# and alignment limits in the argument order below. This wrapper fixes human
# IMGT annotation, fmt7 fields, MakeDb conversion options and IG dmask germlines.
# The rule creates and removes the scratch directory. Empty receptor families
# skip the tool calls; build_reannotated_airr.R still writes the output schema.
set -euo pipefail

if [[ $# -ne 19 ]]; then
    echo "usage: $0 <trust4_airr> <out_igblast_airr> <scratch_dir> <sample_label> <nproc> <igblast_db> <germline_dir> <trust4_cdr3_out> <out_qc_audit> <out_makedb_fail> <out_makedb_log> <num_alignments_v> <num_alignments_j> <num_alignments_d> <num_alignments_c> <out_makedb_pass_gz> <auxiliary_data> <reference_metadata> <out_igblast_hits_gz>" >&2
    exit 2
fi

TRUST4_AIRR="$1"
OUT_AIRR="$2"
IGB="$3"
SAMPLE="$4"
NPROC="$5"
IGBLAST_DB="$6"
GERMLINE_DIR="$7"
CDR3_OUT="$8"
OUT_QC_AUDIT="$9"
OUT_MAKEDB_FAIL="${10}"
OUT_MAKEDB_LOG="${11}"
NUM_ALIGNMENTS_V="${12}"
NUM_ALIGNMENTS_J="${13}"
NUM_ALIGNMENTS_D="${14}"
NUM_ALIGNMENTS_C="${15}"
OUT_MAKEDB_PASS="${16}"
AUXILIARY_DATA="${17}"
REFERENCE_METADATA="${18}"
OUT_IGBLAST_HITS="${19}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The canonical CIGAR adapter belongs to this producer convention. A changed
# analysis image needs its annotation regressions reviewed before using it.
python - <<'PY'
import changeo
if changeo.__version__ != "1.3.4":
    raise RuntimeError("CIGAR normalization requires the pinned Change-O 1.3.4 producer")
PY

# MakeDb keeps equally best-scoring candidates. Lists reaching the configured
# cap cannot establish complete resolution: an unreported tied allele could
# change the assignment. Lower-scoring hits at the cap do not create that
# uncertainty. Pipeline defaults are V 50 and J 20 in params.igblast.
# Raising a cap can change the first tied allele used for germline reconstruction.

# The pipeline default keeps D at IgBLAST's default of 3. A present D call must
# agree with the assigned locus, as V and J calls must. A missing D call passes
# that criterion. Reporting more D alignments can change d_call and which
# records pass QC, so the default stays at IgBLAST's own value.

# C needs the full tied constant-gene set to resolve isotypes at class and
# subclass levels. Different IGHG genes can resolve to IgG while leaving the
# subclass unresolved. IgBLAST's default of 2 truncates any C tie longer than
# two alleles, so a record can read as confidently subclassed while the evidence
# ties several IgG genes. IGHG4A is resolved as IGHG4 downstream.
#
# 125 is how many alleles the pinned human ig C database holds, so no tie can be
# truncated; tr holds 23, and asking for more alignments than exist is harmless.
# The extra fmt7 output lives in the scratch directory and is deleted with it.

# Normalise TRUST4 rows and drop imputed CDR3s (cdr3.out score 0.01) before IgBLAST.
# The parse stage also opens the QC audit; its rows are handed to
# build_reannotated_airr.R below, which extends them with the MakeDb outcome and
# writes the retained shard.
PARSED_TRUST4_AIRR="$IGB/trust4_sequences.tsv"
PARSED_QC_AUDIT="$IGB/qc_audit_parse.tsv"
Rscript "$SCRIPT_DIR/parse_trust4_airr.R" \
    --input_file="$TRUST4_AIRR" \
    --output_file="$PARSED_TRUST4_AIRR" \
    --cdr3_file="$CDR3_OUT" \
    --qc_audit_file="$PARSED_QC_AUDIT" \
    --library_id="$SAMPLE"

IG_FASTA="$IGB/ig.fasta"
TR_FASTA="$IGB/tr.fasta"
: > "$IG_FASTA"
: > "$TR_FASTA"

echo "=== all-locus IgBLAST reannotation ==="
echo "Sample:       $SAMPLE"
echo "TRUST4 AIRR:  $TRUST4_AIRR"
echo "IgBLAST DB:   $IGBLAST_DB"
echo "Germline dir: $GERMLINE_DIR"
echo "Auxiliary:    $AUXILIARY_DATA"

# Split the parsed records into the two IgBLAST query sets. These two patterns
# are the gate that decides which records IgBLAST ever sees: a record on any
# other locus is written to neither FASTA and is left unannotated. That loss is
# not silent (parse_trust4_airr.R counts it in the `igblast_locus_routing` row
# of the QC audit, against the same seven loci named there), so the two lists
# must stay in step. Columns are looked up by header name so the FASTA cannot be
# built from the wrong field if the parsed schema is reordered.
awk -F'\t' -v ig="$IG_FASTA" -v tr="$TR_FASTA" '
  NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    sid = h["sequence_id"]; seq = h["sequence"]; locus = h["locus"]
    if (!sid || !seq || !locus) {
      print "FATAL: sequence_id, sequence, or locus column missing from TRUST4 AIRR" > "/dev/stderr"
      exit 2
    }
    next
  }
  $sid != "" && $seq != "" && $locus ~ /^IG[HKL]$/ {
    printf(">%s\n%s\n", $sid, $seq) >> ig
  }
  $sid != "" && $seq != "" && $locus ~ /^TR[ABDG]$/ {
    printf(">%s\n%s\n", $sid, $seq) >> tr
  }
' "$PARSED_TRUST4_AIRR"

run_family() {
    local family="$1"
    local fasta="$2"
    local nrec seqtype
    # awk rather than `grep -c`, which exits 1 on an empty FASTA and would need
    # its exit status discarded, which would mask a genuine read failure.
    nrec=$(awk '/^>/ { n++ } END { print n + 0 }' "$fasta")
    echo "$family FASTA records: $nrec"
    if [[ "$nrec" -eq 0 ]]; then
        return
    fi

    # IgBLAST's own name for the receptor type, which differs from the database
    # prefix: the ig databases are searched as "Ig", the tr ones as "TCR".
    seqtype=Ig
    [[ "$family" == tr ]] && seqtype=TCR

    # IGDATA is how igblastn finds the internal IMGT domain-boundary data that
    # -domain_system imgt needs; the germline and auxiliary paths below are
    # explicit so a missing database fails here rather than silently falling
    # back. The auxiliary data file is the one obtain_imgt_reference wrote:
    # the image's NCBI file plus entries for the J alleles it lacks, which
    # workflow/lib/germline_reference.py documents. IgBLAST reports no CDR3
    # for a query whose first-listed top J hit has no entry.
    # `7 std qseq sseq btop` is the format MakeDb parses and the format
    # workflow/lib/reannotated_airr.R reads the C-region alignment out of.
    IGDATA="$IGBLAST_DB" igblastn \
        -query "$fasta" \
        -out "$IGB/${family}_igblast.fmt7" \
        -num_threads "$NPROC" \
        -ig_seqtype "$seqtype" \
        -organism human \
        -auxiliary_data "$AUXILIARY_DATA" \
        -germline_db_V "$IGBLAST_DB/database/imgt_human_${family}_v" \
        -germline_db_D "$IGBLAST_DB/database/imgt_human_${family}_d" \
        -germline_db_J "$IGBLAST_DB/database/imgt_human_${family}_j" \
        -c_region_db "$IGBLAST_DB/database/imgt_human_${family}_c" \
        -outfmt '7 std qseq sseq btop' \
        -domain_system imgt \
        -num_alignments_V "$NUM_ALIGNMENTS_V" \
        -num_alignments_J "$NUM_ALIGNMENTS_J" \
        -num_alignments_D "$NUM_ALIGNMENTS_D" \
        -num_alignments_C "$NUM_ALIGNMENTS_C"

    # --nproc is deliberately 1: Change-O parallelises MakeDb across input
    # files, not across the records of one file, so a larger value would only
    # reserve cores it cannot use. igblastn above is the stage that scales.
    # --log keeps the per-record parser decision, which is the only place the
    # REASON a record was rejected is written down; --failed keeps the records
    # themselves. Both are copied out of this scratch directory below.
    MakeDb.py igblast \
        -i "$IGB/${family}_igblast.fmt7" \
        -s "$fasta" \
        -r "$GERMLINE_DIR" \
        --regions default --extended --asis-id --format airr --failed \
        --nproc 1 \
        --log "$IGB/${family}_makedb.log" \
        --outdir "$IGB" --outname "$family"

    if [[ -s "$IGB/${family}_db-fail.tsv" ]]; then
        echo "$family MakeDb failures: $(tail -n +2 "$IGB/${family}_db-fail.tsv" | wc -l)"
    fi
}

run_family ig "$IG_FASTA"
run_family tr "$TR_FASTA"

# Retain the native hit coordinates and aligned bases before scratch cleanup.
{
    for family in ig tr; do
        if [[ -f "$IGB/${family}_igblast.fmt7" ]]; then
            cat "$IGB/${family}_igblast.fmt7"
        fi
    done
} | gzip -n > "${OUT_IGBLAST_HITS}.tmp"
mv "${OUT_IGBLAST_HITS}.tmp" "$OUT_IGBLAST_HITS"


# Retain unmodified MakeDb-pass rows before canonical CIGAR conversion. The
# source tool owns these columns. This evidence survives scratch cleanup and
# stays outside the public release. Both families use the same AIRR writer.
PASS_HEADER=""
: > "$IGB/makedb_pass.tsv"
for family in ig tr; do
    pass="$IGB/${family}_db-pass.tsv"
    [[ -s "$pass" ]] || continue
    header=$(head -n 1 "$pass")
    if [[ -z "$PASS_HEADER" ]]; then
        PASS_HEADER="$header"
        cat "$pass" > "$IGB/makedb_pass.tsv"
    elif [[ "$header" != "$PASS_HEADER" ]]; then
        echo "FATAL: ig and tr MakeDb pass tables have different columns" >&2
        exit 2
    else
        tail -n +2 "$pass" >> "$IGB/makedb_pass.tsv"
    fi
done
if [[ ! -s "$IGB/makedb_pass.tsv" ]]; then
    printf 'sequence_id\tv_cigar\td_cigar\tj_cigar\n' > "$IGB/makedb_pass.tsv"
fi
gzip -c "$IGB/makedb_pass.tsv" > "${OUT_MAKEDB_PASS}.tmp"
mv "${OUT_MAKEDB_PASS}.tmp" "$OUT_MAKEDB_PASS"

# Retain MakeDb's own account of what it rejected. The QC audit records how many
# records MakeDb could not annotate, per locus, but not why, and the rest of this
# scratch directory is deleted when the rule finishes, so without this a
# rejected contig could only be re-examined by re-running IgBLAST.
#
# Each block is one input record: its V/D/J calls and productivity as parsed, plus
# an `ERROR>` line where MakeDb had one. A rejection without an `ERROR>` line
# gives no reason in this log, so the failure table below keeps the parsed fields
# of every rejected record. The log is gzipped because it holds one block per
# input record.
: > "$IGB/makedb_combined.log"
for family in ig tr; do
    if [[ -s "$IGB/${family}_makedb.log" ]]; then
        echo "# receptor_family: ${family}" >> "$IGB/makedb_combined.log"
        cat "$IGB/${family}_makedb.log" >> "$IGB/makedb_combined.log"
    fi
done
gzip -c "$IGB/makedb_combined.log" > "$OUT_MAKEDB_LOG"

# The rejected records themselves, with whatever MakeDb did manage to parse,
# tagged by the database they were run against. MakeDb writes both families
# through the same AIRR writer, so the two headers are expected to match; if a
# future version makes them diverge, stop rather than concatenate mismatched
# columns under one header.
FAIL_HEADER=""
: > "$OUT_MAKEDB_FAIL"
for family in ig tr; do
    fail="$IGB/${family}_db-fail.tsv"
    [[ -s "$fail" ]] || continue
    header=$(head -n 1 "$fail")
    if [[ -z "$FAIL_HEADER" ]]; then
        FAIL_HEADER="$header"
        printf 'receptor_family\t%s\n' "$header" > "$OUT_MAKEDB_FAIL"
    elif [[ "$header" != "$FAIL_HEADER" ]]; then
        echo "FATAL: ig and tr MakeDb failure tables have different columns" >&2
        exit 2
    fi
    awk -F'\t' -v fam="$family" 'NR > 1 { printf("%s\t%s\n", fam, $0) }' "$fail" \
        >> "$OUT_MAKEDB_FAIL"
done
# A library where nothing failed still owes the rule its (header-only) file.
if [[ ! -s "$OUT_MAKEDB_FAIL" ]]; then
    printf 'receptor_family\tsequence_id\n' > "$OUT_MAKEDB_FAIL"
fi

IG_DBPASS=""
TR_DBPASS=""
IG_GERM=""
IG_FMT7=""
TR_FMT7=""
[[ -s "$IGB/ig_db-pass.tsv" ]] && IG_DBPASS="$IGB/ig_db-pass.tsv"
[[ -s "$IGB/tr_db-pass.tsv" ]] && TR_DBPASS="$IGB/tr_db-pass.tsv"
[[ -s "$IGB/ig_igblast.fmt7" ]] && IG_FMT7="$IGB/ig_igblast.fmt7"
[[ -s "$IGB/tr_igblast.fmt7" ]] && TR_FMT7="$IGB/tr_igblast.fmt7"

# MakeDb opens each output table on its first record, so ig_db-pass.tsv exists
# only when at least one IG record was annotated; no emptiness check is needed.
if [[ -n "$IG_DBPASS" ]]; then
    # --failed so a record whose germline could not be reconstructed is named
    # here rather than only showing up downstream as an SHM rate that is
    # silently absent. build_reannotated_airr.R reports the attached count.
    CreateGermlines.py \
        -d "$IG_DBPASS" \
        -g dmask \
        -r "$GERMLINE_DIR" \
        --format airr --failed \
        --outdir "$IGB" --outname ig
    [[ -s "$IGB/ig_germ-pass.tsv" ]] && IG_GERM="$IGB/ig_germ-pass.tsv"
    if [[ -s "$IGB/ig_germ-fail.tsv" ]]; then
        echo "IG CreateGermlines failures: $(tail -n +2 "$IGB/ig_germ-fail.tsv" | wc -l)"
    fi
fi

# Each of the five IgBLAST/Change-O products is passed even when the step that
# would have written it had nothing to run on, in which case the value is empty.
Rscript "$SCRIPT_DIR/build_reannotated_airr.R" \
    --ig_db_file="$IG_DBPASS" \
    --tr_db_file="$TR_DBPASS" \
    --ig_germ_file="$IG_GERM" \
    --ig_fmt7_file="$IG_FMT7" \
    --tr_fmt7_file="$TR_FMT7" \
    --trust4_file="$PARSED_TRUST4_AIRR" \
    --reference_metadata_file="$REFERENCE_METADATA" \
    --num_alignments_v="$NUM_ALIGNMENTS_V" --num_alignments_d="$NUM_ALIGNMENTS_D" \
    --num_alignments_j="$NUM_ALIGNMENTS_J" --num_alignments_c="$NUM_ALIGNMENTS_C" \
    --output_file="$OUT_AIRR" \
    --library_id="$SAMPLE" \
    --qc_audit_in_file="$PARSED_QC_AUDIT" \
    --qc_audit_file="$OUT_QC_AUDIT"

echo "=== all-locus IgBLAST reannotation complete ==="
