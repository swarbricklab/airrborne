# Public test data

This guide describes the public test inputs and how to review changes to them. The four synthetic libraries represent two donors. Their receptor sequences come from public germline references and seeded simulation, with constructed metadata and read representation. All sequences, abundance vectors and metadata are public or synthetic.

Normal testing needs no sequence generator. From the repository root, validate the inputs with:

```bash
python test/system/check_inputs.py
```

The [test guide](../README.md) covers setup and acceptance commands. The `post-trust4` mode starts from `post_trust4/`. The `fastq` mode starts from `reads/` and adds FASTQ discovery, depth counting and assembly. Both modes write fresh outputs in separate run directories.

## Data catalogue

The catalogue lists the inputs to each stage and the evidence behind the expected results. Normal tests use these committed files. Scripts in `test/data_builders/` create candidates for deliberate data changes, and independent reference calculations sit beside the component integration tests.

| File or directory | Origin and purpose | Stage and validation consumer |
|---|---|---|
| `samplesheet.csv` | Constructed library, tissue, donor and read-directory identities, with library display metadata. | FASTQ discovery, scope preparation and report preparation read the system configuration. Input validation checks library membership. |
| `membership.csv` | Constructed project membership and one `fail_qc` exclusion. | Scope preparation; the system checker independently reconstructs included libraries. |
| `tissue_metadata.csv` | Constructed tissue display annotations. | Report preparation; report component tests check scope isolation. |
| `reads/` | `generate_receptors.py` writes 150-base paired reads from declared transcript fragments. | FASTQ discovery, depth counting and TRUST4 assembly; input validation checks every mate against its fragment. |
| `post_trust4/` | `assemble_reads.py` retains authentic native AIRR, CDR3 rows, per-base votes, producer receipts and measured FASTQ depth. | Post-TRUST4 staging, parsing, annotation and SHM; input validation checks native links, counts and producer output hashes. |
| `construction/construction.json` | The generator records seed, tool version, reference digest, read layout and library representation. | Input validation checks membership and the assembly reference; the annotation reference calculation checks the run's reference receipt. |
| `construction/transcripts.json` | The generator records alleles, inserts, trimming, ancestors, mutations and final transcript joins. Positions are zero-based with exclusive ends. | Input validation reconstructs ancestors and reads, checks mutation coordinates and compares native junction sets. |
| `construction/reference_alleles.json` | The generator selects public reference bases with matched GenAIRR anchors and source aliases. | Input validation reconstructs the ancestors used for FASTQ generation. |
| `construction/family_0.tsv`, `construction/family_1.tsv` | GenAIRR lineage node sequences and ancestry. | Input validation checks selected descendants and roots; reviewed clone partitions refer to these constructed families. |
| `construction/fragments.tsv.gz` | The generator records one row per fragment: transcript, coordinates, replicate and lane. | Input validation connects every read pair to its transcript and retained sequencing depth. |
| `construction/native_recovery.tsv` | Observed native record identifiers, junction sources and recovered support from the frozen assembly. | Input validation compares the native outputs; reviewed record expectations anchor downstream decisions. |
| `annotation/segments.fasta`, `annotation/joins.tsv` | Public allele subset and declared seven-locus V-J joins, described in [the annotation guide](annotation/README.md). | Actual pinned annotation tools, common QC and focused SHM/isotype cases. These inputs have no FASTQ assembly claim. |
| `expected/expected_records.tsv` | Reviewed annotations, QC decisions, qualified C calls, representatives and donor/scope clone partitions for every native row. | The record oracle checks parsing, annotation, common QC, support transfer and clone membership. |
| `expected/expected_values.tsv` | Reviewed cells with explicit types and a `review_basis` for each value. | The system checker checks retained intermediate and scope tables with the existing numerical tolerances. |
| `expected/overlap_reference.tsv` | Five constructed repertoire pairs and frozen immunarch answers, retaining the measured version and image identity. | Overlap unit tests compare Jaccard and Morisita-Horn values independently of production code. |
| `checksums.json` | SHA-256 identities of every non-documentation data file except itself. | Input validation requires exact file coverage and matching bytes. |

The checksum catalogue records file identity. The construction checks and reviewed scientific expectations establish the relationships those bytes protect.

## What the inputs exercise

Both donors have two libraries, named `SYN_D1_T1`, `SYN_D1_T2`, `SYN_D2_T1` and `SYN_D2_T2`. `project_a` contains all four. `project_b` excludes `SYN_D1_T1` through its membership `fail_qc` value; the cohort still contains that library. Each library has two matched FASTQ lanes.

The 11 transcript segments include two related IGH families, a separate IGH rearrangement, three TRB junctions, one short constant segment and one IGH descendant with a stop codon. Related IGH sequences have different mutations and constant calls. Sharing and unequal support vary by library, so clone membership, transferred support, isotypes and paired overlap have meaningful inputs.

The short constant segment provides 30 reference bases before the shared J-to-C nucleotide is removed. Its qualification is tested against the unchanged 40-informative-position setting. Native recovered support exceeds the configured IGH 1000 and TRB 100 targets in each library.

These are receptor transcript segments. They include the full selected variable sequence and either 270 or 30 bases of the first constant exon, with the shared J-to-C nucleotide included once. Their extent supplies evidence for annotation, whole-V SHM and isotype qualification.

Leaders, full constant isoforms and UTRs are outside these cases. Reads have no sequencing errors or background transcriptome. IGK, IGL, TRA, TRD and TRG retain focused contract coverage in the existing suites; this connected fixture claims IGH and TRB coverage.

`SYN_D1_T2` has authentic alternative-CDR3 outputs within one IGH backbone and one TRB backbone. Its late IGH junction occurs on two assembled backbones. These outputs exercise the pipeline's handling of alternative CDR3 variants and their support. A junction match in `construction/native_recovery.tsv` identifies its construction source; it does not assert that the entire assembled sequence equals that transcript. The producer files retain TRUST4's actual output checksums, software version, reference and settings.

Representation counts are the repetitions of each tiled fragment. They differ from FASTQ fragment counts and from TRUST4 `consensus_count`. `consensus_count` is recovered read support. None of these quantities represents cells or molecules.

The native recovery table records the observed assembly results. Reviewed pipeline results live in [`expected_values.tsv`](expected/expected_values.tsv) and [`expected_records.tsv`](expected/expected_records.tsv). Independent metric checks reconstruct the released calculations.

The input checker reconstructs ancestors from the reference alleles, trimming and junction inserts, checks recorded mutations, compares every mate to its fragment, validates retained depth and checks native file links. It requires exact agreement between the constructed and recovered junction sets for each library. Its corruption tests detect a wrong mate name, a shortened vote array, a wrong mutation coordinate and a broken backbone link:

```bash
python -m unittest test.components.inputs.unit.test_data_integrity
```

## Read the expected values

[expected_values.tsv](expected/expected_values.tsv) names individual output cells and the evidence used to review them:

| Column | Meaning |
|---|---|
| `table` | Output table path relative to the run's results directory. |
| `selector` | JSON column/value pairs that identify exactly one row; values compare as strings. |
| `column` | The output column to check. |
| `expected` | The reviewed answer for that cell. |
| `type` | `string`, `integer` or `number`, selecting the comparison rule. |
| `review_basis` | How the answer was reviewed; this text explains the expectation. |

The [comparison code](../support/expectations.py) defines numerical tolerances. A mismatch reports the selected cell and its observed and expected values.

## Reviewed record outcomes

[`expected_records.tsv`](expected/expected_records.tsv) specifies every native record's expected annotation, QC, C-call qualification, representative and IGH clone partition. The native IDs, support and junction sources come from `construction/native_recovery.tsv`; those are observed assembly outcomes, independently checked against the frozen files.

V/J genes follow the selected public alleles in `construction/transcripts.json`. Annotation must return every record with its native locus and support. The common-QC failure list contains only the deliberately introduced IGH stop codon. Every other criterion passes these particular records.

The short C segment stays in common QC but cannot meet the 40-position evidence floor. Other IGH records have sufficient evidence for their named C gene. `SYN_D1_T2/assemble0_1` carries the late junction on the early IgG backbone, so its qualified gene is IGHG1; a junction-source label alone does not set its isotype. That BCR alternative transfers support 131 to representative `assemble0_0` (support 714), once, yielding 845. Both TRB alternatives remain representatives.

Clone partition labels describe constructed sequence relationships within each donor. `IGH_root`, `IGH_early` and `IGH_late` belong together, and `IGH_other` and `IGH_other_late` form another family. The deliberately distinct junction and short-C rearrangement remain separate.

The fixture's automatic fallback is 0.1, checked in `expected/expected_values.tsv`. Assertions compare pairwise membership within donor and scope, so renumbering clone IDs is harmless but merging or splitting these families fails. Review the expectations explicitly if sequences or scientific settings change; do not regenerate them from pipeline output.

[Reviewed record checks](../components/record_qc/integration/reference_records.py) use these reviewed decisions to check all parsing, annotation and common-QC audit rows at library and scope level. This fixture has no imputed calls, exact duplicates, routing loss or annotation loss; separate parser and real-tool annotation tests exercise those cases. The oracle does not duplicate the production biological predicates.

## Reproduce a candidate

Curation uses Python 3.11, GenAIRR 2.3.1, Biopython 1.85 and PyYAML 6.0.3, pinned in [environment.yaml](../data_builders/environment.yaml). Create that separate environment and put it on `PATH`; the workflow and test environments stay unchanged:

```bash
conda env create -p .snakemake/fixture-generator -f test/data_builders/environment.yaml
export PATH="$PWD/.snakemake/fixture-generator/bin:$PATH"
```

The reference comes from the pinned Immcantation image through the production reference builder. Its SHA-256 is `085536aa1a10e8492de129f3420f36713d71992b3e48d426e82fa166fd8ea50e`.

From the repository root with the public images available:

```bash
analysis_image=immcantation/suite@sha256:44e228a45b047bdb352363ab20a670336b003e98f3d1bac505f78af41a91ab67
workflow/scripts/container_exec.sh "docker://$analysis_image" "$analysis_image" \
  python3 workflow/scripts/obtain_imgt_reference.py \
  --germline-dir /usr/local/share/germlines/imgt/human/vdj \
  --constant-dir /usr/local/share/germlines/imgt/human/constant \
  --auxiliary-data /usr/local/share/igblast/optional_file/human_gl.aux \
  --output test/resources/trust4/human_IMGT+C.fa \
  --metadata-output test/resources/trust4/human_IMGT+C.fa.metadata.tsv \
  --auxiliary-output test/resources/trust4/human_gl.aux
python test/data_builders/generate_receptors.py \
  --reference test/resources/trust4/human_IMGT+C.fa \
  --output experiments/synthetic_candidate
python test/data_builders/assemble_reads.py \
  --reference test/resources/trust4/human_IMGT+C.fa \
  --output experiments/synthetic_assembly_candidate
```

Both output directories must be new:

- `generate_receptors.py` writes `construction/` and `reads/` in its candidate directory for review.
- `assemble_reads.py` assembles the committed reads and writes `libraries/` in a separate candidate directory. It uses the production container wrapper, provenance and depth scripts.

Review sequence changes before replacing committed reads. Review native outputs separately from expected pipeline results, then copy accepted native `libraries/` under `post_trust4/`. Both commands leave the committed inputs and expected results unchanged.

The generator uses seed 1337, `run_records()` with exposed provenance, and `clonal_lineage()` with two roots and four generations.

The adapter matches GenAIRR's curated V/D/J bases exactly to the pinned pipeline reference and keeps their validated anchors. It selects two genes per segment. The TRB records use TRBD1 with TRBJ1 and TRBC1. Full translated-sequence checks reject unintended stops; one declared stop-containing record remains as a filtering case.

## Promote reviewed inputs

Generate sequences and assemble committed reads in separate new candidate directories. Review transcript changes before replacing reads. Review native outputs and their actual producer history separately from pipeline expectations. Preserve the native output and receipt bytes when copying an accepted bundle.

Update the native recovery ledger from observed native records. Review the record and numerical expectations using the declared sequence relationships, independent calculations and their existing tolerances. Do not generate expectations from pipeline output. Update `checksums.json` only after those data changes have been reviewed. Include every non-Markdown file under `test/data/` except the catalogue itself, using its relative path and SHA-256 digest.

After reviewing the changed data and its catalogue entries, write the checksums:

```bash
python - <<'PY'
import hashlib
import json
from pathlib import Path

root = Path("test/data")
catalogue = root / "checksums.json"
checksums = {
    path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in sorted(root.rglob("*"))
    if path.is_file() and path.suffix != ".md" and path != catalogue
}
catalogue.write_text(json.dumps(checksums, indent=2, sort_keys=True) + "\n")
PY
```

Run input validation, affected component tests and both fresh system modes.

The separate empty-library scenario remains a constructed post-TRUST4 case; `test/data_builders/empty_libraries.py` creates it in a new run directory and `test/system/check_empty_libraries.py` checks its completed outputs.

The permuted case is built the same way by `test/data_builders/permuted_indices.py` and checked by `test/system/check_permuted.py`. Generated runs, references and image caches do not belong under `test/data/`.

## Sources and attribution

The variable and constant reference bases derive from IMGT, the international ImMunoGeneTics information system, Montpellier, France, through the human reference tree in Immcantation Suite 4.8.0. The snapshot identity is the pinned image digest above and the combined-reference SHA-256. The included allele subset and synthetic descendants contain normalisation, trimming, recombination, mutation and transcript joins as recorded in the evidence files. IMGT supplies its data and metadata under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); its [terms of use](https://www.imgt.org/) distinguish those data from IMGT tools (checked 24 September 2026). Cite the IMGT resource when reusing these bases.

GenAIRR, by Thomas Konstantinovsky and Ayelet Peres, supplies the recombination and lineage simulation. The installed release is 2.3.1. The inspected source is revision `1ededa676e8904430f7ea0bdcf03092de4888a41`, including its [recombination API](https://github.com/MuteJester/GenAIRR/blob/1ededa676e8904430f7ea0bdcf03092de4888a41/src/GenAIRR/_experiment/recombination.py) and [clonal lineage API](https://github.com/MuteJester/GenAIRR/blob/1ededa676e8904430f7ea0bdcf03092de4888a41/src/GenAIRR/_experiment/clonal.py). The software has a [GPL-3.0-or-later licence](https://github.com/MuteJester/GenAIRR/blob/1ededa676e8904430f7ea0bdcf03092de4888a41/LICENSE). Its source aliases identify the matching bundled allele from which the anchor was obtained. All retained allele bases match the attributed pipeline reference.

Native inputs were produced by TRUST4 v1.1.9-r619 in this pinned OCI image:

```text
ghcr.io/swarbricklab/trust4@sha256:100b5553815b377837e2864e20ef8b15cc7a84544bc2a5b805fe5c64c277c09f
```

These frozen bundles retain their original producer history byte for byte. Their receipts name the original `test/synthetic/reads/` locations. Those inputs now live under `test/data/reads/` with the same library directories and filenames. The receipts describe the original production and are not rewritten when data moves. `checksums.json` covers their unchanged bytes.

Full FASTQ runs use the assembly image in `test/system/config.yaml` and compare against the same reviewed outcomes. Changing that image does not regenerate the frozen bundles. The frozen producer image and the image configured for a new FASTQ run have separate identities.

The command follows the pinned [TRUST4 documentation](https://github.com/liulab-dfci/TRUST4/blob/v1.1.9/README.md). See Song et al., [Nature Methods 2021](https://doi.org/10.1038/s41592-021-01142-2) for the assembly method. Tests protect this pipeline's use of the native files and its downstream decisions. The panel does not measure tool performance.
