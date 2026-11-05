#!/usr/bin/env bash
# count_seq_depth.sh: raw sequencing depth of the FASTQ that enter TRUST4.
#
# The workflow discovers every matched R1/R2 FASTQ for a library (primary,
# rerun, or top-up) and passes them to TRUST4 as repeated -1/-2 arguments. The
# resulting aggregate input depth is not recorded anywhere else in the pipeline.
# TRUST4's own logs report candidate and assembled reads, not the input supplied
# to it.
# This helper counts reads and bases directly from the raw gzipped FASTQ using
# host `zcat` and `awk`. The `count_sequencing_depth` rule invokes it with no
# container or conda directive, so it runs on the host.
#
# Usage:
#   count_seq_depth.sh OUTPUT_TSV SAMPLE --r1 R1_FILE... --r2 R2_FILE...
#
# The R1/R2 file lists are the exact per-library lists Snakemake resolves via
# find_all_fastq_files() and passes to run_trust4 as -1/-2, so the depth
# reported here is precisely what entered assembly.
#
# Output (one data row, tab-separated):
#   gussid  n_fastq_pairs  total_read_pairs  total_reads  total_bases
# where total_read_pairs = R1 read count, total_reads = R1 + R2 read count,
# and total_bases = R1 + R2 sequenced bases. Each matched file pair must contain
# the same number of records, and a mismatch fails this job. The depth job runs
# independently of the assembly job, so it does not gate assembly. Because the
# check has to pass here, total_reads is exactly twice total_read_pairs.

set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "ERROR: usage: $(basename "$0") OUTPUT_TSV SAMPLE --r1 R1_FILE... --r2 R2_FILE..." >&2
    exit 2
fi

out="$1"
sample="$2"
shift 2

# Parse the two file lists. --r1/--r2 flip which array subsequent paths land in.
r1=()
r2=()
cur=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1) cur="r1" ;;
        --r2) cur="r2" ;;
        *)
            case "$cur" in
                r1) r1+=("$1") ;;
                r2) r2+=("$1") ;;
                *)  echo "ERROR: file argument '$1' given before --r1/--r2" >&2; exit 2 ;;
            esac
            ;;
    esac
    shift
done

if [[ ${#r1[@]} -eq 0 || ${#r2[@]} -eq 0 ]]; then
    echo "ERROR: need at least one R1 and one R2 file (got ${#r1[@]} R1, ${#r2[@]} R2)" >&2
    exit 2
fi
if [[ ${#r1[@]} -ne ${#r2[@]} ]]; then
    echo "ERROR: unequal R1/R2 file counts for ${sample}: ${#r1[@]} R1, ${#r2[@]} R2" >&2
    exit 2
fi

# Validate unwrapped four-line FASTQ while counting in one decompression pass.
# Headers start with @, separators start with +, and sequence/quality lengths
# agree. LF and CRLF are accepted, including a final line without a newline.
# Empty files have zero depth. Each file must end at a complete record boundary.
# Quality-score ranges and mate identifiers remain the assembler's concern.
count_stream() {
    zcat -f -- "$1" \
        | awk -v input="$1" -v library="$sample" '
            function invalid(reason) {
                printf "ERROR: invalid FASTQ for %s in %s at line %d: %s\n", \
                    library, input, NR, reason > "/dev/stderr"
                failed = 1
                exit 1
            }
            {
                sub(/\r$/, "")
                if (NR % 4 == 1 && $0 !~ /^@./) invalid("expected @ header")
                if (NR % 4 == 2) sequence_length = length($0)
                if (NR % 4 == 3 && $0 !~ /^\+/) invalid("expected + separator")
                if (NR % 4 == 0) {
                    if (length($0) != sequence_length) invalid("sequence/quality lengths differ")
                    reads++
                    bases += sequence_length
                }
            }
            END {
                if (failed) exit 1
                if (NR % 4 != 0) invalid("incomplete four-line record")
                printf "%.0f\t%.0f\n", reads, bases
            }'
}

echo "count_seq_depth: sample=${sample}  R1 files=${#r1[@]}  R2 files=${#r2[@]}"
echo "  R1: ${r1[*]}"
echo "  R2: ${r2[*]}"

# Count mates pair by pair so mismatches in different lanes cannot cancel out
# in the library totals. Command substitution under `set -o pipefail` aborts if
# zcat fails (for example, on a truncated or corrupt member).
r1_reads=0
r2_reads=0
r1_bases=0
r2_bases=0
for index in "${!r1[@]}"; do
    r1_out="$(count_stream "${r1[$index]}")"
    r2_out="$(count_stream "${r2[$index]}")"
    r1_pair_reads="${r1_out%%$'\t'*}"
    r1_pair_bases="${r1_out##*$'\t'}"
    r2_pair_reads="${r2_out%%$'\t'*}"
    r2_pair_bases="${r2_out##*$'\t'}"

    if [[ "$r1_pair_reads" -ne "$r2_pair_reads" ]]; then
        printf 'ERROR: R1/R2 read-count mismatch for %s pair %d: %s=%d, %s=%d\n' \
            "$sample" "$((index + 1))" \
            "${r1[$index]}" "$r1_pair_reads" \
            "${r2[$index]}" "$r2_pair_reads" >&2
        exit 1
    fi

    echo "  pair $((index + 1)): reads=${r1_pair_reads} R1_bases=${r1_pair_bases} R2_bases=${r2_pair_bases}"
    r1_reads=$(( r1_reads + r1_pair_reads ))
    r2_reads=$(( r2_reads + r2_pair_reads ))
    r1_bases=$(( r1_bases + r1_pair_bases ))
    r2_bases=$(( r2_bases + r2_pair_bases ))
done

n_pairs="${#r1[@]}"
total_read_pairs="$r1_reads"
total_reads=$(( r1_reads + r2_reads ))
total_bases=$(( r1_bases + r2_bases ))

echo "  R1 reads=${r1_reads} bases=${r1_bases}"
echo "  R2 reads=${r2_reads} bases=${r2_bases}"

mkdir -p "$(dirname "$out")"
temporary="$(mktemp "${out}.tmp.XXXXXX")"
trap 'rm -f "$temporary"' EXIT
{
    printf 'gussid\tn_fastq_pairs\ttotal_read_pairs\ttotal_reads\ttotal_bases\n'
    printf '%s\t%d\t%d\t%d\t%d\n' "$sample" "$n_pairs" "$total_read_pairs" "$total_reads" "$total_bases"
} > "$temporary"
chmod "$(printf '%03o' "$((0666 & ~$(umask)))")" "$temporary"
mv -f "$temporary" "$out"

echo "count_seq_depth: wrote ${out} (read_pairs=${total_read_pairs}, total_reads=${total_reads}, total_bases=${total_bases})"
