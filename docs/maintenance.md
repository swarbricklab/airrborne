# Maintain the pipeline

This guide covers changes to the pipeline: finding analytical decisions, updating dependencies, investigating reference problems and validating deployments. For running an analysis and reading its results, start with the [running guide](running.md) and [interpretation guide](interpretation.md).

The [test guide](../test/README.md) gives the acceptance commands, and the [coverage map](../test/COVERAGE.md) links scientific requirements to tests.

## Follow a library through the code

Bulk RNA-seq contains reads from many cells mixed together. TRUST4 assembles receptor sequences from those reads, and `consensus_count` records their recovered read support. Cell identities and receptor-chain pairings are unavailable from this analysis.

Start with the rule in `workflow/rules/` to find the input files, settings and command. `Snakefile` includes the rule files in dependency order and omits `assembly.smk` for analysis-only runs. Follow each rule's entry script in `workflow/scripts/` to the functions in `workflow/lib/` that make the decisions. Some tool choices are in the rule or shell script itself, as listed below.

| Rule file | Responsibility |
|---|---|
| [common.smk](../workflow/rules/common.smk) | Shared code and resource paths, wildcard constraints and the library processing-code list. |
| [inputs.smk](../workflow/rules/inputs.smk) | Scope membership, analysis image, curated reference and reference audit. |
| [assembly.smk](../workflow/rules/assembly.smk) | FASTQ pairing, TRUST4 assembly and independent sequencing-depth counts. |
| [library.smk](../workflow/rules/library.smk) | Reannotation, common QC, grouping preparation and library measurements. |
| [igh.smk](../workflow/rules/igh.smk) | Donor pools, scope threshold, clone assignment, SHM summaries, diversity and selection. |
| [scope.smk](../workflow/rules/scope.smk) | Scope collectors, headline and donor overlap calculations. |
| [provenance.smk](../workflow/rules/provenance.smk) | Tool versions, source evidence and effective scientific settings. |
| [exports.smk](../workflow/rules/exports.smk) | Shared library AIRR exports, scope metadata linkage and transactional publication. |
| [report.smk](../workflow/rules/report.smk) | Optional display metadata, report preparation and rendering. |

1. **Read the assembly.** [parse_trust4_airr.R](../workflow/scripts/parse_trust4_airr.R) reads TRUST4's tables. [trust4_parsing.R](../workflow/lib/trust4_parsing.R) identifies imputed CDR3 calls for removal by the script. One assembled sequence can have several reported CDR3 alternatives, each with its own AIRR row. [assembly_backbone.R](../workflow/lib/assembly_backbone.R) later decides which alternatives take part in grouping.
2. **Annotate the sequences.** [run_igblast_all_loci.sh](../workflow/scripts/run_igblast_all_loci.sh) runs IgBLAST and Change-O to assign receptor genes and describe how the sequence aligns to the reference. [reannotated_airr.R](../workflow/lib/reannotated_airr.R) joins those annotations to the original sequence IDs and read support. It also records whether a gene assignment can be used for grouping.
3. **Apply common QC.** [record_qc.R](../workflow/lib/record_qc.R) decides which annotated records pass the checks shared by all receptor loci. `filter_airr.R` writes the passing table. A record that passes can still lack the evidence needed by one later analysis.
4. **Prepare and group the records.** Non-IGH clonotypes require an exact nucleotide junction and resolved V/J assignments. IGH clone assignment links similar junctions at a chosen distance threshold. A chain of links can join sequences farther apart than that threshold. Both paths keep links to the original records and account for their read support.
5. **Calculate measurements.** Isotype, somatic hypermutation (SHM), selection and diversity apply their own requirements. The report displays their results. Change the function that calculates a result or decides its availability, then follow its output to the report.

The [glossary](glossary.md) explains the table's terms. The [QC matrix](qc_matrix.md) lists the applied checks and their recorded outcomes. The table below points to the code that makes those choices.

## Locations of analytical decisions

Paths in this table are relative to `workflow/`. Function names are search terms within the linked file. The named settings are under `params` in `config/airrborne/config.yaml`; their valid values are in `workflow/schemas/config.schema.yaml`.

| Change or review | Where the decision is made | Where to follow it |
|---|---|---|
| Common record QC | [lib/record_qc.R](../workflow/lib/record_qc.R), `airr_qc_flags()` | [scripts/filter_airr.R](../workflow/scripts/filter_airr.R) applies the flags. `airr_qc.max_alignment_n_fraction` controls the allowance for unknown aligned bases. |
| Multiple V or J assignments | [lib/gene_calls.R](../workflow/lib/gene_calls.R), `vj_tie_kind()` and `vj_canonical_token()` | [lib/reannotated_airr.R](../workflow/lib/reannotated_airr.R) checks whether IgBLAST may have omitted equally good matches and attaches the resolved assignments. [lib/germline_reference.py](../workflow/lib/germline_reference.py), `build_reference_metadata()`, establishes which complete reference sequences are identical. |
| Several CDR3 alternatives from one assembly, or repeated full sequences | [lib/assembly_backbone.R](../workflow/lib/assembly_backbone.R), `prepare_assembly_mapping()` | The IGH and non-IGH preparation functions below supply eligibility. This shared function chooses representative rows, keeps their source links and adds each eligible source's read support once. |
| Non-IGH clonotype definition | [lib/rearrangement_groups.R](../workflow/lib/rearrangement_groups.R), `build_rearrangement_group_table()` | [lib/rearrangement_keys.R](../workflow/lib/rearrangement_keys.R), `exact_rearrangement_group_key()`, builds the value used to identify matching junction and V/J assignments. [scripts/build_clonotypes.R](../workflow/scripts/build_clonotypes.R) writes the groups and membership. |
| Which IGH sequences can enter clone inference | [lib/igh_grouping.R](../workflow/lib/igh_grouping.R), `prepare_igh_grouping()` | [lib/igh_clones.R](../workflow/lib/igh_clones.R), `clone_junction_eligible()`, checks the junction alphabet and `lineage.max_n` allowance. [scripts/prepare_igh_grouping.R](../workflow/scripts/prepare_igh_grouping.R) writes the shared input for fitting and assignment. |
| Fit an IGH threshold from a donor's sequences | [lib/donor_distance_threshold.R](../workflow/lib/donor_distance_threshold.R), `fit_donor_distance_threshold()` | [scripts/fit_donor_threshold.R](../workflow/scripts/fit_donor_threshold.R) passes `lineage.min_fit_distances`, `fit_edge`, `fit_specificity` and `fit_subsample`. The function defines the distance calculation and SHazaM fit. |
| Choose the threshold used across a scope | [lib/scope_clone_threshold.R](../workflow/lib/scope_clone_threshold.R), `choose_scope_clone_threshold()` | [scripts/resolve_scope_clone_threshold.R](../workflow/scripts/resolve_scope_clone_threshold.R) reads donor fits and `lineage.clonal_threshold`. This is where a configured number, the median usable fit or the fallback is selected. |
| Assign IGH clones | [scripts/assign_scope_donor_clones.R](../workflow/scripts/assign_scope_donor_clones.R), the `hierarchicalClones()` call | The call chooses how junction differences join sequences into clones. [lib/igh_clones.R](../workflow/lib/igh_clones.R) builds membership and clone IDs. Review the distance calculation in donor fitting alongside any assignment change. |
| Accept a constant-region call or name its isotype | [lib/isotype_calls.R](../workflow/lib/isotype_calls.R), `c_call_qc_flags()`, and [lib/isotype_resolution.R](../workflow/lib/isotype_resolution.R), `resolve_isotype()` | The first file checks evidence and blanks failed calls through `mask_low_evidence_c_calls()`. The second resolves gene names to isotype labels. [lib/isotype_metrics.R](../workflow/lib/isotype_metrics.R) decides which read support contributes to each reported fraction. |
| Require enough evidence for SHM or selection | [lib/shm_summaries.R](../workflow/lib/shm_summaries.R), `shm_region_eligible()` and `choose_shm_evidence()`; [lib/baseline_selection.R](../workflow/lib/baseline_selection.R), `selection_record_eligible()` | [lib/igh_features.R](../workflow/lib/igh_features.R) and [lib/shm_regions.R](../workflow/lib/shm_regions.R) measure SHM evidence. Selection prepares its own clone sequences in `collapse_clones_to_effective_sequences()`. Their requirements are separate from common QC and clone assignment. |

### Nucleotide and amino-acid clonotypes

A change from nucleotide to amino-acid grouping starts in `build_rearrangement_group_table()` and `exact_rearrangement_group_key()`. `junction` contains nucleotides. `junction_aa` contains their amino-acid translation. Different nucleotide sequences can encode the same amino-acid sequence, so this change can merge groups that are separate today.

Review these parts together before implementing that change:

1. **Membership.** Decide which loci change and whether V/J assignments still need to match. Define what to do with missing, partial or ambiguous amino-acid sequences. The present eligibility checks expect nucleotides.
2. **Reported sequences and support.** Decide how a group describes its member nucleotide sequences. Review representative selection and support totals in `rearrangement_groups.R`; keep links back to all contributing records.
3. **Comparisons.** Trace `clonotype_key` from the grouped tables to [scripts/calculate_donor_overlap.R](../workflow/scripts/calculate_donor_overlap.R) and its `group_key` input to [lib/overlap_metrics.R](../workflow/lib/overlap_metrics.R). The same key function also records exact junction and V/J matches for IGH in `build_clone_membership_frame()`. Changing that field does not change the SCOPer clone definition.
4. **Checks and descriptions.** Start with [clonotype definitions](../test/components/non_igh_clonotypes/unit/test-clonotype_table.R), [retained-sequence grouping](../test/components/grouping/unit/test-retained_sequences.R) and the overlap tests. Add cases where different nucleotide junctions share one translation and review the expected membership, support and comparison results. Update the dictionary, QC matrix, interpretation guide and generated methods together.

### Tool calls and their settings

| Tool | Command or function call | Preparation and settings |
|---|---|---|
| TRUST4 | [rules/assembly.smk](../workflow/rules/assembly.smk), `run_trust4` | The rule passes `params.trust4` settings and the curated reference. Native output parsing is in `trust4_parsing.R`. |
| IgBLAST, Change-O MakeDb and CreateGermlines | [scripts/run_igblast_all_loci.sh](../workflow/scripts/run_igblast_all_loci.sh) | [rules/library.smk](../workflow/rules/library.smk), `reannotate_igblast`, passes database paths and alignment limits from `params.igblast`. The wrapper sets the other tool arguments. |
| SHazaM threshold fitting | [lib/donor_distance_threshold.R](../workflow/lib/donor_distance_threshold.R), `fit_donor_distance_threshold()` | `distToNearest()` measures junction differences; `findThreshold()` fits a cutoff. The module also controls the fitting seed. Review its package-internal seed replacement when upgrading SHazaM. |
| SCOPer clone assignment | [scripts/assign_scope_donor_clones.R](../workflow/scripts/assign_scope_donor_clones.R), `hierarchicalClones()` | Prepared sequences are grouped separately for each donor within each scope, using that scope's saved threshold. |

`deps.trust4_container` and `deps.immcantation_docker` select the tool images. The [annotation integration test](../test/components/annotation/integration/test-annotation.R) exercises the annotation tools. The [threshold and assignment test](../test/components/igh_clones/integration/test-threshold_assignment.R) follows fitting through scope threshold selection to clone assignment.

Keep changes to scientific choices separate from documentation edits. Follow the [test guide](../test/README.md#check-a-complete-change) for acceptance, and review changed expectations against the intended method. Even comments in tracked workflow source can cause jobs to rerun and processing IDs to change.

## Scientific policy rationale

The table below explains the scientific policies implemented in the pipeline and their trade-offs. The linked code and tests show how each policy is applied. Changing a policy requires an independently justified scientific comparison.

| Policy and problem | Choice and alternative | Trade-off, outputs and evidence |
|---|---|---|
| Several IG junction variants come from one contig. | Choose the dominant common-QC variant before grouping eligibility. Do not replace it with a minor eligible sibling. | Prevents eligibility from choosing the preferred junction; can lose support when the dominant variant fails. Changes grouping, clones and all their summaries. See `assembly_backbone.R` and grouping retained-sequence tests. |
| Several source observations can represent one retained sequence. | Transfer eligible support once through contig and exact-sequence representative links; retain source rows. Counting every variant independently is the alternative. | Conserves recovered support without counting every source as an inference member. Light-chain transferred support need not confirm the selected junction independently. See `rearrangement_groups.R`, membership tables and grouping tests. |
| V/J calls can span distinct genes. | Resolve a composite only from identical complete reference sequences and complete candidate lists. Choosing one top name is the alternative. | Reduces unsupported assignments but excludes unresolved support. Changes clone partitions, clonotypes and gene usage. See `gene_calls.R`, annotation and gene-resolution tests. |
| Donors differ in their fitting evidence. | Use the median usable donor threshold, then 0.1 if none qualifies; retain all fits. Pooling distances or fitting every donor independently would define different clones. | One scope threshold makes its assignment rule explicit, but the median and fallback remain modelling policies. See `donor_distance_threshold.R` and threshold-assignment integration tests. |
| Identical member sequences can have several partial evidence sources. | Require one source to pass a regional SHM gate; do not pool votes or partial coverage. | Keeps coverage tied to one assembly, at the cost of unavailable SHM. Changes coverage and SHM summaries while preserving membership. See `shm_summaries.R` and SHM evidence tests. |
| SHM and selection answer different questions. | SHM gates supported regional rates; selection retains partial aligned members, collapses them and masks incomplete effective codons. | Selection can use partial evidence but is exploratory; its eligibility does not certify full coverage. See `baseline_selection.R`, `shm_regions.R` and their independent numerical tests. |

The [decision map](#follow-a-library-through-the-code) locates the owners. For a proposed change, record the problem, alternative policies, scientific comparison, affected outputs and retained tests with the change. Do not infer biological validation from passing the computational regression suite.

## Output compatibility and documentation ownership

`OUTPUT_SCHEMA_VERSION` in `workflow/lib/release_bundle.py` identifies the published interface and is copied into provenance and the completion receipt. Use these rules when changing it:

| Change | Version policy |
|---|---|
| Denominator, row population, identifier construction, blank-versus-zero meaning, existing column name/type or enum meaning | Major version. Existing analyses can change meaning or fail. Treat enum additions as breaking for closed-category readers. |
| Scientific method or analytical default that changes estimates or eligibility | Major version even when the table shape is unchanged. Record the method and setting change separately. |
| Optional additional table or column with existing meanings preserved | Minor version, with a migration note; exact-header readers still need updating. |
| Prose, diagnostics or validation that leaves the scientific interface and values unchanged | Patch version when a release is issued. A documentation-only source commit does not rewrite an existing freeze. |

Structural compatibility does not establish scientific comparability. Compare source revision, effective settings, reference identities, assembly history and population as well as the schema version. IDs only identify rows within their stated release and scope context.

Schema 1.0.0 to 2.0.0 changes AIRR orientation semantics. The stored sequence and its coordinates stay in V orientation. Preserve the former `rev_comp` in `changeo_rev_comp` and set standard `rev_comp` to F. Existing analyses that reverse the stored sequence when the old flag is T need correction. Produce a new validated release with fresh provenance and hashes; do not modify a frozen release in place.

A hypothetical change from measured-member to grouped-member SHM denominators would also need a major version, even with the same name and type. That denominator change is not part of 2.0.0.

| Fact | Authoritative owner | Repeated views |
|---|---|---|
| Defaults and bounds | Config template and configuration schema, checked against rule inputs. | Running guide and effective config. |
| Columns, types, keys and categories | Table writer plus `docs/data_dictionary.csv`, checked against outputs. | Release file catalogue and report labels. |
| Scientific rationale | The policy explanations above and interpretation guide. | Short README limitations. |
| Reuse commands | Executable scripts in `examples/` with import/join tests. | Running guide and portable README snippets. |
| Actual choices | Effective config, source evidence, provenance and operative threshold. | Generated methods and report. |

The release README repeats selected facts so it remains useful away from the pipeline checkout. Generate those facts from the authoritative values and check them with the existing contract tests. Keep independent scientific calculations independent of production estimators.

## Update the analysis image

The analysis image contains IgBLAST, Change-O, SHazaM and SCOPer. An upgrade is a scientific compatibility change until comparisons establish otherwise.

1. Record the old and candidate digests and query the installed versions. Update `deps.immcantation_docker` in the public and test configurations, then update deployment configurations deliberately. Preserve old images, references, receipts and native assembly histories.
2. Prepare fresh reference paths with the candidate image. Run the [reference audit](#audit-the-germline-reference), including IgBLAST database agreement, completeness metadata, supplemented J auxiliary entries and remaining gaps. Do not overwrite the old reference snapshot for comparison.
3. Check the annotation panel and coordinate adapters with forward/reversed queries, insertions/deletions and reverse-D hits. Check native-to-standard conversion, selected reference bases, germline reconstruction availability and the explicit Change-O restore path.
4. Compare donor distances, mixture fits, operative thresholds and membership by sequence. Then compare SHM measured/eligible populations, mutation counts, effective observed/germline selection sequences and numerical estimates. Use the annotation, grouping, IGH-clone, SHM, selection and release components.
5. Run fresh post-TRUST4, FASTQ, empty and permuted acceptance with the intended source revision. Inspect copied offline reports. Review changed expected records and cells against independent evidence before updating them.
6. Decide the output version and migration requirements above. Update fixed methods, generated descriptions and citations where behaviour changed. Record acceptance SHAs, tool/reference identities and retained comparison outputs before publishing the new release.

### SHazaM seed replacement

SHazaM 1.3.2 resets the random seed with `set.seed(NULL)` inside `gmmFit` and `rocSpace`, after the pipeline has seeded the fit. The wrapper `override_shazam_gmm_seed` in `donor_distance_threshold.R` replaces exactly one such reset in each function with seed 1337 in the running R session. This makes fits repeatable without changing the mixture family or cutoff rule. It fails if the expected package code shape changes. Threshold component tests check repeatability and the replacement; they do not establish biological fit quality.

On upgrade, inspect both functions in the candidate package. Remove the replacement if the package preserves the caller's seed or provides an explicit seed argument, then rerun the independent threshold and assignment checks. Do not silently patch a changed function or weaken the shape check to make an upgrade pass.

## Retain the files needed for reuse

| Task | Required material | Included in portable `release/`? |
|---|---|---|
| Browse results offline | Entire release, including report assets and detail files. | Yes. |
| Analyse tables with project metadata | Release plus the authoritative metadata and its join-key definition. | Metadata is separate. |
| Reuse AIRR with native Change-O tools | Release, matching checkout, pinned analysis image and matching IMGT-gapped reference. | Checkout, image and reference are separate. |
| Rerun after TRUST4 | Retained native AIRR, CDR3, final.out, depth files and producer receipts, plus manifests, settings, code and reference/software. | No; use the [retained-input list](running.md#reuse-retained-assemblies). |
| Audit SHM vote attribution | Release membership/coverage, source IGH features, consensus-vote tables and native final.out. | Position votes and source features are separate. |
| Reproduce assembly | Original paired FASTQs, samplesheet, original assembly image/reference/settings and producer receipts. | No. |
| Reconstruct modified pipeline source | Base Git revision, submodules, binary working/index patches and retained untracked source files. | Source evidence is separate. |

For a run with local code changes, retain `shared/provenance/source_evidence/<state_sha256>/` and the base Git objects. `write_source_evidence` saves binary patches and untracked source there. The hash identifies this evidence; reconstructing the source requires the files themselves.

Automatic capture covers the paths listed in `SOURCE_PATHS` in `provenance.py`. Save changes outside those paths and dirty submodule patches separately. For symlinked untracked source, also archive the target file.

For each acceptance run, record `git rev-parse HEAD`, submodule revisions, host, date, exact commands and validation receipts. The historical Gadi table below combines several runs and lacks a separate pipeline SHA for every row. Commit `b78741cc8f7afee9d2b4f7bc31729e43665c264c` introduced that account; use each execution's provenance for its exact source identity.

## Update the assembly image

1. Select a published [BioContainers tag](https://quay.io/repository/biocontainers/trust4?tab=tags) and inspect its architecture, digest and installed `run-trust4` version.
2. Set `deps.trust4_container` to the same digest in `config/airrborne/config.yaml` and `test/system/config.yaml`. Update parent analysis configurations separately.
3. Run both fresh system modes and the [acceptance checks](../test/README.md#check-a-complete-change). Keep frozen native bundles and their producer records unchanged. Review scientific differences before changing expectations.
4. Run `setup.sh` again, then use a normal FASTQ run to produce new assemblies. Under Apptainer or Singularity, setup pulls the new image as a SIF into `deps.resources`, named by the new digest, with its receipt beside it. Keep historical caches. `--post-trust4` retains the original assemblies and their producer versions.

### Update the release environment

Edit [envs/release.yaml](../envs/release.yaml), then regenerate the Linux x86-64 pins. From the clone's root, the following command uses Docker on Linux or macOS. Start Docker first.

```bash
release_pin="$(mktemp)"
docker run --rm --platform linux/amd64 -v "$PWD:/work" -w /work \
  condaforge/miniforge3 bash -c 'conda env create -q -p /tmp/release \
  -f envs/release.yaml >&2 && conda list --explicit --md5 -p /tmp/release' \
  > "$release_pin" && mv "$release_pin" envs/release.linux-64.pin.txt
```

The existing pin file is replaced only after a successful solve and export. Update the pin-generation date in the YAML header and review the package changes. Snakemake hashes the YAML and platform pin file, so changing either invalidates existing workflow environments. Run `bash test/test_env.sh` to refresh the supplemental environment from its input signature, then run the [acceptance checks](../test/README.md#check-a-complete-change).

Snakemake 9.27.1 falls back to solving the YAML if pin installation fails. Resolve that warning before claiming an environment used the recorded pins.

## Audit the germline reference

Publication requires the curated reference, supplemented auxiliary data and allele metadata sidecar to match the pinned analysis image. The audit rebuilds these files and compares IgBLAST's databases with the reference and source FASTAs. To run the audit from the clone's root with the default results path:

```bash
./run_module.sh -- results/airrborne/shared/provenance/imgt_reference.json
```

The receipt records these fields:

| Field | Content |
|---|---|
| `reference`, `image`, `imgt_snapshot` | File and image identities, IMGT source and download date. |
| `germline_files`, `curation` | Source FASTAs, allele counts and repeated-name choices. |
| `reference_metadata`, `identity_classes` | Metadata-sidecar identity and complete-reference sequence classes with their exact allele membership. |
| `igblast_databases` | Allele and sequence comparisons for each database. |
| `auxiliary_data` | Original and supplemented entries, copied allele pairs and remaining gaps. |
| `agreement` | Reference, metadata and auxiliary rebuilds, and database agreement. |

A disagreement stops the reference audit without writing a receipt. Read the differences named in its log. If prepared files have been altered, keep the reference named by `refs.imgt_reference_name`, its `.metadata.tsv` sidecar and `human_gl.aux` for comparison.

Choose a new `deps.resources` directory and rerun the workflow to rebuild the files from the pinned image. Retain the previous reference and receipts with their original results. The [reference explanation](interpretation.md#germline-reference) describes which missing auxiliary entries the pipeline can supplement.

## Validate a Gadi deployment

Deployment validation covers job submission, container execution and output checks. A dry run checks the planned jobs, and a run inside an allocation checks workflow execution. To validate PBS submission as well, complete the profile and accounting checks below. The [acceptance record](#gadi-acceptance-record) lists completed runs.

### Prepare in a copyq job

Compute nodes have no network access, and login nodes limit each process, so prepare the software in a `copyq` job, following the [Gadi section of the running guide](running.md#gadi). A checkout without a deployment configuration uses the public test configuration:

```bash
git submodule update --init --recursive
qsub -q copyq -P "$PROJECT" -W block=true -j oe -o setup.log \
  -l ncpus=1,mem=32GB,walltime=04:00:00,storage=gdata/$PROJECT+scratch/$PROJECT,wd <<'EOF'
./setup.sh --config test/system/config.yaml --tests --no-browsers
EOF
```

Setup downloads both images to `deps.resources`, records their identities beside them and creates the release environment. It writes the resources path to `.snakemake/env.sh` for the tests.

The release environment uses `.snakemake/conda` when `conda-prefix` is uncommented in `profiles/workflow/config.yaml`. Otherwise it uses the shared prefix in `profiles/global/nci/config.yaml`. Check that this path, `PROJECT`, `storage` and container binds match the deployment. Keep the checkout, inputs, caches and outputs on project storage accessible to compute nodes.

Assembly needs the prepared TRUST4 SIF. A `docker://` URI cannot replace it on these nodes: with a populated cache, SingularityCE 3.11.3 still contacted quay.io and failed on a Gadi compute node on 1 October 2026.

Post-TRUST4 execution does not obtain the assembly image. Its version collection records unavailable versions when the configured OCI URI cannot be resolved. Check that behaviour during the first live deployment.

Record the pipeline and shared-profile revisions, effective configuration, image identities and reference digests. Validate the deployment configuration and project membership against their schemas and scope manifests.

### Run the component tests in an allocation

Component R tests run their containers from the launching process, so run them inside an allocation:

```bash
qsub -I -q normal -P "$PROJECT" \
  -l ncpus=4,mem=16GB,walltime=01:00:00,storage=gdata/$PROJECT+scratch/$PROJECT,wd
source .snakemake/env.sh
python test/run_tests.py
```

### Submit the system cases through PBS

From a login node, the runner submits each rule as a PBS job through `profiles/global/nci` and checks the results when the workflow ends:

```bash
source .snakemake/env.sh
./test/system/run.sh post-trust4
./test/system/run.sh fastq
./test/system/run.sh post-trust4 --case empty
./test/system/run.sh post-trust4 --case permuted
```

Each case stages a fresh directory under `.snakemake/test-runs/`, runs the workflow and checks the released tables, reviewed record decisions and independent calculations. Fresh acceptance requires a directory that has not already run the workflow.

Cases running at the same time from one checkout share the prepared images, receipts and reference. The cases also run inside one allocation with `-- --profile none` after [preparing local execution](running.md#run-on-a-cluster).

Keep job logs, benchmarks, PBS identifiers and release validation results. Check requested memory and walltime against PBS accounting before scaling up. Wrapper-only measurements omit container and child-process use.

The shared adapter maps the `threads` resource to PBS `ncpus`, while commands use Snakemake's native `threads`. Keep `set-threads` and `set-resources.threads` matched for each multithreaded rule. Confirm the command thread count, core allocation and `ncpus` agree for each such rule and for a single-core job. Check memory, walltime, jobfs, queue and storage mappings separately.

Each submission attempt receives a unique log path. The status adapter reads the complete signed PBS exit code; a nonzero code means failure. Queued, running and held jobs remain pending. Check these behaviours through actual submission and accounting. Changes to the shared adapter need validation with its other users.

Browser checks need system libraries that setup cannot install without root access. Run them on a workstation or in CI, following the [test guide](../test/README.md#check-a-complete-change).

### Gadi acceptance record

The runs below passed on Gadi on 29 and 30 October 2026 from the pipeline commits that introduced this record, with shared profile revision `2c7a5a6`, SingularityCE 3.11.3, Snakemake 9.27.1 and the pinned analysis and assembly images.

| Run | Where | Outcome |
|---|---|---|
| `setup.sh --config test/system/config.yaml --tests --no-browsers` | `copyq`, fresh clone, 1 core, 32 GB requested | Both images, their receipts, the release environment and the test environment prepared in 69 minutes using 8.5 GB of disk. |
| `python test/run_tests.py` | Compute node, 8 cores, 32 GB requested | 170 Python cases, every R suite and the JavaScript suite passed. |
| `test/system/run.sh post-trust4` | The same allocation, `--profile none -- --cores 8 --resources mem_mb=30000` | Passed. The workflow took 6 minutes, and a repeated launch had nothing to do. |
| `test/system/run.sh fastq` and `post-trust4 --case empty` at the same time | The same allocation, 4 cores and 15000 MB each | Both passed from one checkout. The workflows took 9 minutes and 5 minutes. |
| `test/system/run.sh post-trust4 -- --profile none` | A later 8-core, 32 GB allocation | The launcher set `--cores 8` and `mem_mb=32768` from the allocation. Passed in 3 minutes with the reference already prepared. |
| `test/system/run.sh post-trust4` and `fastq --profile profiles/global/nci` at the same time | Login node, PBS submission | 164 and 172 submitted jobs finished in 71 and 81 minutes with no failures or retries, and a repeat of the FASTQ case took 77 minutes. Release, reviewed-record and independent-calculation checks passed for every run. |
| `test/system/run.sh post-trust4` | Login node of the fresh clone prepared above, PBS submission | 165 submitted jobs finished in 80 minutes, and the runner's checks passed. |
| `python test/run_tests.py` with the resources directory named by `.snakemake/env.sh` | Compute node, 4 cores, 16 GB requested | 170 Python cases and every R suite passed in 10 minutes. |
| `test/system/run.sh post-trust4` after `setup.sh --config test/system/config.yaml --tests` | Login node, PBS submission | 165 submitted jobs finished in 52 minutes, and the runner's checks passed. |
| `setup.sh --tests --no-browsers` with a deployment configuration naming `resources/airrborne/`, then `test/system/run.sh post-trust4` | Fresh clone: `copyq`, then a login node with PBS submission | Setup reused the pulled images in 1 minute. The check ran with those images, created no `test/resources/`, and 165 submitted jobs finished in 50 minutes with the runner's checks passed. |
| `test/system/run.sh post-trust4 -- --profile none` after the launcher gained Slurm sizing and `profiles/local` lost its `jobs` line | A 4-core, 16 GB allocation | The launcher set `--cores 4` and `mem_mb=16384` from the allocation. Passed; the job used 7 minutes of walltime and 2.47 GB. |

The prepared files measured 2.5 GB for the analysis SIF, 40 MB for the TRUST4 SIF, 415 MB for the release environment and 645 MB for the test environment. The allocation job used 2.9 GB and 27 minutes of walltime for all of its runs.

PBS accounting for the submitted jobs confirmed the mapping: single-core rules received `ncpus=1`, and `reannotate_igblast` received `ncpus=8`, 32000 MB and 12 hours. Every job used under 1.5 GB and under one minute on the public inputs, so the profile allocations remain provisional for production data.

On an earlier run that day, one FASTQ attempt stopped in `build_clonotypes` when SingularityCE reported `failed to mount squashfs filesystem: file exists` while eight containers started together on one node, and the rerun passed without change. The local profile now retries a failed job once for this reason.

## Redraw the rule graph

`docs/graphs/workflow_dag.svg` is the rule graph of a run that starts from FASTQs. Redraw it after a change to rules or their inputs. The [renderer](graphs/render_rulegraph.py) wraps labels, assigns rules to the stages in the [README overview](../README.md#analysis-overview) and removes redundant arrows with Graphviz `tred`. Stage names and colours match the Mermaid diagram. When adding or renaming a rule, update the renderer's `STAGES` mapping; an unmapped rule stops diagram generation. Every rule and dependency path is retained except the aggregate target. An omitted arrow has another directed path between the same rules; this display does not enumerate every direct file dependency. The renderer changes only the diagram, not the workflow. Native SVG titles contain the exact rule names. Open the SVG itself in a browser to access tooltips; image previews may not expose them.

The public test configuration supplies FASTQs in any clone:

```bash
conda install -n snakemake_9.27.1 -c conda-forge graphviz
./run_module.sh -c test/system/config.yaml --profile none -d -n
```

- `-d -n` writes the graph over the committed file without executing workflow jobs. `-d` alone writes it and then runs the workflow.
- The `conda install` line adds Graphviz `dot` and `tred` when the host lacks them.
- A `--post-trust4` launch with `-d` leaves out `count_sequencing_depth` and `run_trust4`, so redraw the graph in FASTQ mode before committing it.

## Prepare a source release

Run the [release acceptance checks](../test/README.md#check-a-complete-change) for the intended source revision, including the manual CI run. Check the author order and identifiers in `CITATION.cff`. Add the chosen source version, release date and assigned DOI when available; keep these consistent with the tag and deposition metadata.

Test setup from a fresh recursive clone of the intended release tag. Build the deposited archive from tracked source and inspect its contents. Git source archives omit Git metadata and submodule contents, so [execution requires a clone](running.md#prepare-the-environment).

Record the pipeline revision and the profile revision from `git ls-tree HEAD profiles/global` in the deposition metadata. Include the pinned profile source when archiving the Gadi deployment files.

### Register on WorkflowHub

WorkflowHub supports [Snakemake workflows](https://about.workflowhub.eu/docs/registering_workflows/supported-workflow-types/). Use its [Git repository import](https://about.workflowhub.eu/docs/registering_workflows/adding-files/) to register `https://github.com/swarbricklab/airrborne`, select the tested release tag and choose `Snakefile` as the main workflow file. Select Snakemake as the workflow type if it is not detected. The repository supplies these files:

| File | Use during registration |
|---|---|
| `Snakefile` | Main executable workflow, with supporting rules and scripts in the repository. |
| `README.md` | Description, setup instructions and a synthetic example. |
| `LICENSE` | MIT licence for the pipeline source. |
| `CITATION.cff` | Ordered creators, affiliations and ORCIDs. |
| `docs/graphs/workflow_dag.svg` | Optional workflow diagram. |
| `test/system/config.yaml` and `test/data/` | Example configuration and public synthetic inputs. |

The Git import route does not require a committed `ro-crate-metadata.json` or an abstract CWL file. WorkflowHub [generates RO-Crate metadata](https://about.workflowhub.eu/docs/registering_workflows/complete-workflow-metadata/) during registration. Abstract CWL and a diagram are optional attachments; the synthetic example supports reproducibility but is not a registration prerequisite.

Complete the registration form with a descriptive title, description, maturity, the submitting team, MIT licence, sharing permissions, tags and creators in citation order. Check each creator's ORCID and affiliation. Review the imported description and ensure that it directs users to clone the repository recursively: a downloaded workflow file or source archive alone does not provide the Git metadata and submodule checkout required by this pipeline.

For a [WorkflowHub DOI](https://about.workflowhub.eu/docs/citable/), review the metadata, freeze the workflow version, make it public and then generate the DOI. Keep the source tag and profile revision recorded with that entry. Link any separate software archive and paper from the workflow record. Deposit study results separately from the reusable workflow source.
