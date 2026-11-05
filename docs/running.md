# Run the pipeline

Follow this guide to run airrborne on your own data. For a first run with synthetic inputs, use the [README quick start](../README.md#quick-start). The [glossary](glossary.md#software-and-execution) explains Snakemake, container images, profiles and bind paths.

| Task | Section |
|---|---|
| Describe the libraries and FASTQs | [Inputs](#prepare-inputs) |
| Change analysis settings | [Analysis settings](#configure-the-analysis) |
| Install the software | [Environment](#prepare-the-environment) |
| Run locally or reuse assemblies | [Launch](#launch) |
| Add annotations to the report | [Report metadata](#add-report-metadata) |
| Submit cluster jobs | [Run on a cluster](#run-on-a-cluster), including [Gadi](#gadi) |
| Find results or diagnose a failure | [Outputs](#outputs) and [troubleshooting](#troubleshooting) |
| Find the table that answers a question | [Interpretation guide](interpretation.md#find-the-table-for-a-question) |

## Prepare inputs

### Samplesheet

Replace the rows of [config/airrborne/samplesheet.csv](../config/airrborne/samplesheet.csv) with one row per sequencing library. The required columns are `library_id`, `tissue_id`, `donor_id` and `path`:

```csv
library_id,tissue_id,donor_id,path,sequencing_batch
library_a,tissue_1,donor_1,batch_1/library_a,batch_1
library_b,tissue_2,donor_1,batch_1/library_b,batch_1
library_c,tissue_3,donor_2,batch_2/library_c,batch_2
library_d,tissue_4,donor_2,batch_2/library_d,batch_2
```

- Identifiers are strings, so `0460` and `NA` stay labels. Each starts with a letter or digit and contains only letters, digits, `.`, `_` and `-`.
- A `donor_id` names one person across the whole samplesheet. Each scope infers IGH clones jointly from that donor's libraries. When two studies reuse a donor code, prefix it with the study, so `donor_1` becomes `project_a_donor_1`.
- Each tissue belongs to one donor, and each library to one tissue.
- `path` is a directory below `deps.data`. Each library needs paired FASTQs with the [filenames below](#fastq-layout).
- Other columns are optional [report metadata](#add-report-metadata).
- Every row is analysed. If the study requires one technical replicate per tissue, select it before running. Libraries from the same tissue remain separate observations. Put all lanes and sequencing deliveries of one library in the same directory.

**From a sample spreadsheet.** Map each person to a `donor_id`, each biopsy or sample to a `tissue_id`, and each sequencing library to a `library_id`. With one library per sample, `tissue_id` may equal `library_id`.

A `tissue_id` identifies one specimen. A donor's tumour biopsy and blood sample need different IDs. Put descriptions such as tumour or blood in a separate column in the [tissue metadata file](#add-report-metadata).

The `path` column names a library directory, so omit columns containing individual R1 and R2 paths. For FASTQs under `/data/rnaseq/<library>/`, set `deps.data: /data/rnaseq/` and write:

```csv
library_id,tissue_id,donor_id,path
P01_T1,P01_T1,P01,P01_T1
P01_T2,P01_T2,P01,P01_T2
P02_T1,P02_T1,P02,P02_T1
```

Use de-identified codes for `donor_id` and `tissue_id`, because they appear in the released tables. The files in `deps.data` are read-only inputs to the pipeline.

**From one directory of FASTQs.** When every library's files sit in one directory, make a directory of links with one subdirectory per library and point `deps.data` at it. For Illumina names such as `P01_T1_S7_L001_R1_001.fastq.gz`, with or without the `_L001` lane part:

```bash
mkdir -p /data/links
for r1 in /data/raw/*_R1_001.fastq.gz; do
  library=$(basename "$r1" | sed -E 's/_S[0-9]+(_L[0-9]+)?_R1_001\.fastq\.gz$//')
  mkdir -p "/data/links/$library"
  ln -s "$r1" "${r1/_R1_001/_R2_001}" "/data/links/$library/"
done
ls /data/links
```

The command removes the sample number, lane and read suffix to find the library ID. It then puts links to all lanes of that library in one directory. Check the directories listed by `ls`. For a different naming scheme, adjust the `sed` expression to extract the library ID.

For a library delivered in several folders, link all its FASTQs into the same library directory. Give each link a unique name while keeping a supported read-pair pattern. The containers also need access to the original files: add their directory, here `/data/raw`, to `AIRRBORNE_CONTAINER_BINDS` unless it is already mounted. See [container selection](#container-selection).

**From SRA downloads.** An SRA experiment (`SRX`) describes a combination of sample, library and sequencing settings. Its runs (`SRR`) identify submitted read files, as the [SRA metadata guide](https://www.ncbi.nlm.nih.gov/sra/docs/submitmeta/) explains. Check that the runs describe the same physical library before combining them in one pipeline directory. When each library has one run, the run accession can name the library:

```bash
mkdir -p /data/links
for r1 in /data/raw/*_1.fastq.gz; do
  library=$(basename "$r1" _1.fastq.gz)
  mkdir -p "/data/links/$library"
  ln -s "$r1" "${r1%_1.fastq.gz}_2.fastq.gz" "/data/links/$library/"
done
```

**Start with a few libraries.** For a first run, include only a few libraries in the samplesheet. After checking their results, add the remaining libraries and rerun. The pipeline reuses the completed assemblies and annotations, processes the new libraries and repeats scope-level steps such as clone inference and report creation. In the dry run, `run_trust4` lists only the new libraries.

Check `metrics/headline.tsv` in the first run's release to assess receptor recovery:

- `igh_yield_per_million_read_pairs` gives TRUST4-reported full-CDR3 support per million input read pairs, before quality control. There is a corresponding column for each of the other six loci.
- `igh_rarefaction_status` and `trb_rarefaction_status` show whether diversity can be compared at a common read-support depth. `available` means the library reached the target (IGH 1000 or TRB 400) and met any configured diversity minimum. Other values explain why the estimate is unavailable.

Real-data recovery benchmarks by sequencing depth or tissue are not recorded in this guide. Use the first few libraries to assess recovery for your dataset, and check their [resource use](#resources-on-a-cluster) before planning the full run.

### FASTQ layout

Each library directory holds gzipped paired-end reads named with one of these patterns:

```text
<prefix>_R1.fastq.gz      <prefix>_R2.fastq.gz
<prefix>_R1_001.fastq.gz  <prefix>_R2_001.fastq.gz
<prefix>_R1.fq.gz         <prefix>_R2.fq.gz
<prefix>_1.fastq.gz       <prefix>_2.fastq.gz
<prefix>_1.fq.gz          <prefix>_2.fq.gz
```

The last two patterns support SRA and ENA download names. The `_R1` and `_R2` forms also accept an underscore suffix when it matches in both mates. For example, `<prefix>_R1_1.fastq.gz` pairs with `<prefix>_R2_1.fastq.gz`.

Put every lane of a library in the same directory. The pipeline matches read pairs by prefix, suffix and extension, then assembles and counts them together. For example, `P01_T1_S7_L001_R1_001.fastq.gz` and `P01_T1_S7_L002_R1_001.fastq.gz`, with their R2 mates, form one library.

Files with other naming patterns are ignored. A recognised R1 or R2 file without its mate stops the run. Use symbolic links to give files supported names if needed. The pipeline uses the reads as supplied, without trimming or filtering before assembly.

FASTQ paths, `deps.resources` and `outs.results_dir` may contain only ASCII letters, digits or one of `_./:+-`. TRUST4 builds shell commands from these paths.

### Project membership

An optional membership CSV, named by `membership_csv`, declares projects:

```csv
project,library_id,usage
project_a,library_a,analysed
project_a,library_b,analysed
project_a,library_c,analysed
project_b,library_c,analysed
project_b,library_d,analysed
project_b,library_e,fail_qc
```

- `usage` is `analysed` or `fail_qc`. A `fail_qc` row leaves that library out of that project.
- A `fail_qc` row can name a library the samplesheet does not carry, which records why the library is absent. Only `analysed` rows are checked against the samplesheet.
- A library can belong to several projects.

[config/airrborne/membership.csv](../config/airrborne/membership.csv) is a template.

### Scopes

A scope is a set of libraries analysed together and published in one `release/` folder. The reserved `cohort` scope contains every library in the samplesheet. Without `membership_csv`, it is the only scope. Each project in the membership file adds a separate scope alongside `cohort`.

- A run builds every declared scope by default.
- `--scope NAME` builds one scope and processes only its libraries. The others, including `cohort`, produce nothing in that run.
- Each scope lists its included libraries in `manifest/scope_manifest.tsv`.

The pipeline infers IGH clones separately for each donor within each scope, using one threshold for the scope. A library included in two scopes can therefore have different clone assignments in each. Both scopes reuse the same library annotation and library-level measurements.

To analyse two studies separately and together, declare each study as a project. Each study then has its own release, and `cohort` analyses both at one threshold. Compare the studies within the `cohort` release, as the [interpretation guide](interpretation.md#igh-clones) explains.

For comparisons such as responders versus non-responders or tumour versus blood, keep both groups in one scope. Add the group labels to [tissue metadata](#add-report-metadata) so the report can group and colour results by them. For statistical analysis, join the same metadata to the released tables.

## Configure the analysis

Edit [config/airrborne/config.yaml](../config/airrborne/config.yaml) to set the input and output paths. Its comments give each key's pipeline default and distinguish it from the tool's default. Numerical QC and evidence settings live under `params`, and thread, memory and walltime requests live in `profiles/`.

Every run checks the configuration against `workflow/schemas/config.schema.yaml`. Unknown keys and missing required keys stop the run. The schema supplies defaults for optional settings where defined.

### Diversity and rarefaction

- `params.rarefaction.targets` sets the depth at which diversity is compared. The production targets are IGH 1000 and TRB 400 [support units](glossary.md#assembly-and-annotation).
- Omitting a locus from the targets removes its target value and curve. Its diversity row and support vector remain.
- `params.diversity.min_consensus_count` defaults to zero at every locus. A positive value makes diversity, dominance and target-depth values unavailable below it. Support totals, richness and curves stay available.
- Neither setting filters input libraries or sets a report-wide sample-size gate.

The [interpretation guide](interpretation.md#diversity-and-dominance) explains the resulting statuses.

### Pair comparisons

`params.overlap.loci` selects the loci for within-donor pair comparisons, which need at least two libraries from one donor in the scope. Their support minimum is separate from the diversity minima and rarefaction targets:

```yaml
params:
  overlap:
    loci: [IGH, TRB]
    min_consensus_count: {IGH: 0, TRB: 0}
```

- `loci` selects from IGH, TRA, TRB, TRD and TRG and cannot be empty. IGK and IGL have clonotype summaries but are outside the overlap interface.
- Minimum values are nonnegative integers. Zero, including for an omitted locus, applies no minimum.
- A positive minimum must be met at both endpoints for default report inclusion. Pairs below it stay in `donor_overlap_pairs.tsv` with their `pair_status`.

Changing the pair minimum recalculates overlap and refreshes the release. The [interpretation guide](interpretation.md#support-and-report-inclusion) explains how the report uses it.

## Choose settings from the evidence

Use the defaults in `config/airrborne/config.yaml` as a starting point, then review the evidence below before changing them. The synthetic tests check the implementation. The suitability of a cutoff depends on the dataset and the comparison being made. Each release records the settings and methods used.

| Choice | Evidence to inspect | Consequence and basis |
|---|---|---|
| IGH clone threshold | Donor fits, component separation, fallback status and member co-assignment under alternative thresholds. | Changes clones and all clone summaries. The donor median, 500-distance minimum and 0.1 fallback are pipeline policies. |
| Rarefaction target | Grouped support distribution, target status and which libraries remain comparable. | Changes comparison depth and availability. The default comparison depths are IGH 1000 and TRB 400. Assess suitability using the support available in the libraries being compared. |
| SHM evidence rules | Complete coverage audit, reconstruction availability and source vote-mapping diagnostics. | Changes the eligible member population. The vote and occupancy floors are pipeline evidence policies. |
| C-call evidence minimum | Callable class/subclass fractions and informative C-alignment positions. | Changes the composition denominator. The 40-position default is a pipeline policy. |
| Pair-support minimum | Both endpoint totals and sensitivity of the comparison. | Changes default report inclusion; defined below-minimum values stay released. Zero disables the positive minimum. |
| Selection consensus frequency | Member coverage, effective consensus and usable clone densities. | Changes observed effective sequences. The default 0.6 follows SHazaM's canonical frequency choice; the pipeline adds incomplete-codon masking. |

Keep alternative analyses in separate results directories. Report changed settings and eligibility populations alongside the result. The [threshold guidance](interpretation.md#assess-a-threshold-fit) explains why a computational fit status alone is insufficient evidence.

## Platform and runtime evidence

Check the testing status of your platform before a production run. The [Gadi record](maintenance.md#gadi-acceptance-record) gives the hosts, dates and resource measurements for completed tests.

| Platform and runtime | Evidence status | Route |
|---|---|---|
| Linux x86-64, Gadi SingularityCE 3.11.3, local allocation or PBS | Tested in the recorded Gadi runs. | Use the bundled Gadi profile or local execution inside an allocation. |
| Linux x86-64, Apptainer or Docker | Implemented and expected; runtime-specific acceptance is not recorded in this guide. | Install the runtime and run the public fixture before production. |
| Intel macOS, Docker Desktop | Expected but unverified. | Use the macOS installers and local profile. |
| Apple Silicon macOS, Docker Desktop with amd64 emulation | Expected but unverified; native arm64 analysis images are not supplied. | Use the matching host installer; Docker runs the pinned amd64 image. |
| Windows through WSL2 with a Linux runtime | Unverified; no acceptance claim. | Validate the Linux workflow inside WSL2 before use. |
| Native Windows shell | Unsupported. | The launcher requires a Unix shell and tools. |
| Other clusters, including Slurm | Deployment-dependent and unverified here. | Validate submission, allocation and runtime separately. |

## Prepare the environment

Use a Git clone, and run commands from its root directory unless a section says otherwise. Downloaded ZIP and tar archives lack the Git metadata needed to record the pipeline version. To use a specific release, check out its tag and run `git submodule update --init --recursive`.

Setup needs network access to create environments and download container images. On a cluster, read [run on a cluster](#run-on-a-cluster) before starting setup.

Install conda and a container runtime. [Miniforge](https://github.com/conda-forge/miniforge) installs conda without administrator rights. Install it on storage with room for environments, and put it on `PATH`. This example is for Linux x86-64:

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /data/projects/study/miniforge3
source /data/projects/study/miniforge3/etc/profile.d/conda.sh
```

For macOS, use the matching Intel or Apple Silicon installer in the [Miniforge platform instructions](https://github.com/conda-forge/miniforge#install). Install and start [Docker Desktop for Mac](https://docs.docker.com/desktop/setup/install/mac-install/). The pipeline requests Linux amd64 images, using emulation on Apple Silicon. These macOS paths are unverified by the acceptance record below.

Repeat the `source` line for the actual installation in each new shell and job script. Then prepare the software:

```bash
./setup.sh
```

Setup requires `conda` on `PATH`. It creates the Snakemake environment, downloads both pinned container images and prepares the environment used for exports and reports. This last step checks the configured workflow, so the samplesheet and its inputs must already exist. For retained assemblies, add `--post-trust4` instead.

If setup is interrupted, run it again. It reuses the environments and images already prepared. Also rerun setup after changing a pinned image or switching to another configuration. `./setup.sh --help` lists the options.

Use the table below to plan storage. On Gadi, setup including the test environment took 69 minutes and used 8.5 GB of disk. The [acceptance record](maintenance.md#gadi-acceptance-record) gives the dates and conditions of those measurements.

| What | Where | Size |
|---|---|---|
| `snakemake_9.27.1` environment | The directory conda uses for environments, shown by `conda info --envs`. | |
| Downloaded conda packages | The conda package cache, `pkgs_dirs` in `.condarc`. | |
| Both images as SIFs, with their receipts | `deps.resources` | 2.5 GB for the analysis image and 40 MB for the TRUST4 image. |
| Image layers while a SIF is built | `TMPDIR`, or `APPTAINER_TMPDIR` or `SINGULARITY_TMPDIR`. The pull bypasses the runtime's cache. | About 2.7 GB for the analysis image, measured on GitHub Actions (`ubuntu-24.04`) on 24 September 2026. |
| Release environment | `.snakemake/conda` in the launch directory for local execution. With a site profile, the `conda-prefix` the profiles set, as the [Gadi section](#gadi) describes. | 415 MB |
| Test environment, with `--tests` | `.snakemake/test-env` in the clone. | 645 MB |

A Miniforge installed as above keeps its environments and packages under its own directory. When conda comes from a site module or the home directory, and the home directory has a small quota, add these lines to `~/.condarc` before setup:

```yaml
envs_dirs:
  - /data/projects/study/conda/envs
pkgs_dirs:
  - /data/projects/study/conda/pkgs
```

Point `TMPDIR` at storage with room for the image layers.

Add `--tests` to prepare the [test environment](../test/README.md). Setup writes `.snakemake/env.sh`; source this file in the shell where tests will run. It activates Snakemake and selects the runtime, test interpreter, configuration and resources prepared by setup.

Normal pipeline launches do not need `.snakemake/env.sh`. Once the test environment exists, later setup runs continue to include it.

The launcher uses Snakemake 9.27.1 when it is already on `PATH`. Otherwise it activates the `snakemake_9.27.1` Conda environment and checks the version. Start Docker before launching if it is the selected runtime.

| Requirement | Version or source |
|---|---|
| Snakemake | 9.27.1, from [the environment definition](../envs/snakemake_9.27.1.yaml). |
| Container runtime | Apptainer, Singularity or Docker on Linux; Docker Desktop on macOS. |
| Assembly and analysis images | OCI digests pinned in [the configuration](../config/airrborne/config.yaml). |
| Release environment | [envs/release.yaml](../envs/release.yaml), with explicit package pins for Linux x86-64 (`linux-64`). Setup creates it for the host's execution method, and a launch creates it when it is absent; other platforms solve the YAML. |

### Container selection

The pipeline uses the first available runtime in this order: Apptainer, Singularity, Docker. Set `AIRRBORNE_CONTAINER_RUNTIME` to choose one explicitly.

With Apptainer or Singularity, setup stores both container images as SIF files in `deps.resources`. It records each image's identity in a file beside it. Runs sharing this directory reuse those images, so compute jobs need no network access. If an image is missing at launch, the main Snakemake process downloads it.

Docker stores the images in its daemon. Both images use `linux/amd64`, which requires emulation on Apple Silicon.

The container can access the launch directory automatically. For data, resources or outputs elsewhere, set `AIRRBORNE_CONTAINER_BINDS` to a comma-separated list of existing absolute directory paths. For symbolic links, include the directories containing the original files.

## Launch

Run these commands from the clone's root. The first previews the jobs without running them; the second starts the analysis. Relative paths in the configuration are resolved from this directory:

```bash
./run_module.sh -n
./run_module.sh
```

On hosts outside NCI Gadi, jobs run on the launching machine with `profiles/local`. The profile uses at most 8 cores and gives assembly and annotation up to 8 threads each, bounded by the core limit. For example, to limit the whole workflow to two cores:

```bash
./run_module.sh -- --cores 2
```

Set `--resources mem_mb=N` to limit the total requested memory of jobs running at once, where `N` is the memory available to the run in megabytes. Snakemake uses the per-rule requests in `profiles/local/config.yaml` to decide which jobs can run together. These requests do not reserve memory or enforce a limit on a running job.

If a job exceeds the host or allocation's available memory, the system ends it. Raise that rule's `mem_mb` request so fewer jobs run together, or use a host or allocation with more memory.

| Option | Effect |
|---|---|
| `-c`, `--config FILE` | Configuration file. Default: `config/airrborne/config.yaml`. |
| `-n`, `--dry-run` | Dry run. |
| `--post-trust4` | Start from retained assemblies and depth records. |
| `--scope NAME` | Build one declared scope. |
| `--profile DIR` | Select a site profile for execution, graphs and unlock. Use `none` for local execution, including on Gadi. Inside a PBS or Slurm job, local execution takes `--cores` and the `mem_mb` limit from the allocation. |
| `--workflow-profile DIR` | Replace the bundled per-rule resources, `profiles/workflow` with a site profile or `profiles/local` otherwise, for execution, graphs and unlock. |
| `-u`, `--unlock` | Unlock the launch directory, then exit. |
| `-d`, `--dag` | Maintainers only: redraw the committed rule graph, as the [maintenance guide](maintenance.md#redraw-the-rule-graph) describes. |
| `-h`, `--help` | Print the options and notes, then exit. |
| `-- ARGS` | Pass further arguments to Snakemake. A run normally needs only `--cores N` and `--resources mem_mb=N`. Graph and unlock commands do not receive them. |

Pass the launcher options `--post-trust4` and `--scope` on the command line.

The dry run prints a table of rules and job counts. It reports missing samplesheet columns, FASTQs or mates before starting any jobs. A successful analysis ends with `Snakemake workflow execution finished.` and writes `release/.complete.json` for each selected scope.

After a failure or interruption, rerun the same command to continue from the completed jobs. For a long run over SSH, use a `tmux` or `screen` session so a dropped connection does not end the launcher.

### Reuse retained assemblies

Use this mode when an earlier run's TRUST4 assemblies and depth counts already exist under `outs.results_dir`:

```bash
./run_module.sh --post-trust4 -n
./run_module.sh --post-trust4
```

Add `--scope NAME` to select an existing scope. Without project membership, the only scope is `cohort`.

Each selected library needs these files under `libraries/<library_id>/`:

```text
trust4/<library_id>_airr.tsv
trust4/<library_id>_cdr3.out
trust4/<library_id>_final.out
metrics/<library_id>_seq_depth.tsv
```

FASTQ discovery, depth counting and assembly are skipped. Reference preparation still runs when needed. A missing retained file stops the run. Use a normal FASTQ run after changing reads or assembly settings.

Keep `trust4/<library_id>_assembly_provenance.json` beside reused outputs when available. Without it, the release records their producer history as unknown.

## Add report metadata

Add metadata to show annotations and filters in the report. Without metadata, the report still supports donor and tissue selection. Metadata does not change the scientific tables. To use it in downstream analysis, [join it to those tables](interpretation.md#project-metadata) by library or tissue ID.

### Declare the files

For a tissue table keyed by `tissue_id`, only the file is required:

```yaml
scope_metadata:
  project_a:
    tissue:
      file: config/airrborne/tissue_metadata.csv
```

Without `membership_csv`, the scope to name is `cohort`.

Use `columns` to select which annotations appear. Add a `library` entry for annotations keyed by `library_id`:

```yaml
scope_metadata:
  project_a:
    tissue:
      file: config/airrborne/tissue_metadata.csv
      columns: [timepoint, timepoint_order, response]
    library:
      file: config/airrborne/samplesheet.csv
      columns: [sequencing_batch]
```

- Files are CSV or TSV. Without `columns`, the report uses every non-identifier column.
- The report shows the original headers. Include units in the header where needed. No column definitions or type declarations are required.
- The samplesheet is read as report metadata only when a scope names it.
- Columns selected from the tissue and library files must have distinct names.
- Core identifiers validate the join and are never imported as annotations. Rows for libraries or tissues outside the scope stay out of the report.
- Paths are relative to the launch directory, and every scope name must be a declared scope.

[config/airrborne/tissue_metadata.csv](../config/airrborne/tissue_metadata.csv) is a template.

### Column types

The report infers each column's type from the full source column:

- A column with at least one finite number, and otherwise only numbers, empty cells or `NA`, is numeric and gets continuous colours.
- Mixed text and numbers, and integer codes with leading zeroes, stay categorical.
- An entirely empty column stays categorical and unavailable.
- Empty cells are missing. `NA` is missing in a numeric column and a literal label in a categorical column. Zero is a value.
- Units are not inferred.

### Timepoints

Put timepoints in the tissue metadata file, under the exact headers `timepoint` and `timepoint_order`:

```csv
tissue_id,donor_id,timepoint,timepoint_order,response
tissue_1,donor_1,baseline,1,stable
tissue_2,donor_1,follow_up,2,stable
tissue_3,donor_2,baseline,1,progressed
tissue_4,donor_2,follow_up,2,progressed
```

- `timepoint` supplies tissue-level labels, including numeric-looking ones. It can group, colour and filter plots, like any categorical column such as `response`.
- `timepoint_order` must hold finite numbers, empty cells or `NA`. It orders visits and supports continuous colouring and range filters. It never means elapsed time, and it is not a grouping field.
- With both fields, the report offers reference-visit comparisons when each label has one order and different labels have different orders. With `timepoint` alone, labelled panels and pair filters remain available.
- Report preparation rejects conflicting label and order mappings within the scope and names the tissues and values. One label can repeat on several tissues with the same order. Tissues outside the scope cause no conflict.
- A missing label or order stays missing and is not filled from another tissue.
- When `columns` restricts the file, select both fields for ordered comparisons.

The report never reads a visit label or order as a treatment phase.

### Join rules

- The join keeps every library and every repeated tissue. A missing metadata row stays missing.
- Duplicate keys, missing selected columns, malformed rows and mismatched donor or tissue identifiers stop report preparation.
- Annotations stay on their source rows. The pipeline does not decide which fields are donor characteristics or spread values across tissues.

### Refresh a report

After a metadata change, rerun the command that produced the run:

```bash
./run_module.sh
```

Changing metadata values or declarations rebuilds the affected reports and releases without rerunning the scientific analysis. A shared metadata file refreshes every report that reads it. If the FASTQs are no longer available, add `--post-trust4` to use the retained assemblies.

## Deploy in an analysis repository

Add the pipeline as a submodule of the repository holding the data and configuration. Run from that repository's root:

```bash
git submodule add https://github.com/swarbricklab/airrborne.git modules/airrborne
git submodule update --init --recursive
cp -r modules/airrborne/config/airrborne config/airrborne
```

Edit the copied `config/airrborne/config.yaml` for the analysis repository's data, resources and outputs. Setup and the launcher read that path from the launch directory by default, so run from the repository's root:

```bash
./modules/airrborne/setup.sh
./modules/airrborne/run_module.sh -n
./modules/airrborne/run_module.sh
```

Commit the submodule reference to record the pipeline revision used by the analysis repository. Results and logs are written relative to that repository's root. The pipeline still reads its bundled scripts and profiles from the submodule.

## Run on a cluster

Choose one of two ways to run on a cluster: run the whole workflow inside one allocation, or use a site profile to submit each rule as a separate scheduler job. Without a site profile, the launcher runs jobs on the machine where it starts.

Prepare the software on a host with network access. Make conda and the container runtime available in job scripts as well as at login. The launcher needs either `conda` on `PATH` or Snakemake 9.27.1 already active.

If no runtime is available on `PATH`, setup and the launcher try to load a module named `singularity`, or the value of `AIRRBORNE_CONTAINER_RUNTIME`. If the site's module has another name, load it before setup. The generated `.snakemake/env.sh` repeats that load for tests.

**Prepare the software.** On Gadi, use a `copyq` job as described [below](#gadi). At another site, ask which host permits environment creation and image downloads. This may be a data-transfer node, an interactive job or a partition with network access. Allow enough time for the downloads and installation; the recorded first setup on Gadi took 69 minutes, including test software.

**Inside one allocation.** Request an allocation using the site's normal procedure, load the container runtime and set `AIRRBORNE_CONTAINER_BINDS` for data outside the clone. Then use the [launch commands](#launch). This method needs no scheduler profile.

Inside a PBS or Slurm job, the launcher reads the allocation's CPU and memory limits and prints them. For PBS, it reads `PBS_NCPUS` and `PBS_VMEM`. For Slurm, it reads `SLURM_CPUS_PER_TASK` and `SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` multiplied by the CPU count. To set the limits explicitly, pass values after `--`:

```bash
./run_module.sh -- --cores 8 --resources mem_mb=30000
```

`--cores` limits the CPUs used by the workflow, and `mem_mb` limits the total memory requested by concurrent jobs in megabytes. Keep both within the allocation. On Gadi, add `--profile none` to select this execution method; see the [Gadi instructions](#gadi) for environment preparation.

Use this Slurm script as a starting point and replace the paths with your own. Allocation detection has been tested with PBS on Gadi, but not on a Slurm host. Check the CPU and memory limits printed in the job output:

```bash
#!/bin/bash
#SBATCH --job-name=airrborne
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=airrborne-%j.out
# Put conda on PATH as at login, with this line or the site's conda module.
source /path/to/conda/etc/profile.d/conda.sh
cd /path/to/airrborne
module load apptainer
export AIRRBORNE_CONTAINER_RUNTIME=apptainer
export AIRRBORNE_CONTAINER_BINDS=/data/projects/study
./run_module.sh
```

- Add the account and partition options the site requires, such as `#SBATCH --account` and `#SBATCH --partition`.
- Set `--time` within the partition's limit. Submit the script again if the limit ends the run: it continues from the finished jobs.
- Request memory with either `--mem` or `--mem-per-cpu`.
- Snakemake's progress and any error appear in the `--output` file. Each job's own log is under `outs.log_dir`.

The [resource notes](#resources-on-a-cluster) below give the per-rule bounds.

**One scheduler job per rule.** The bundled `profiles/global/nci` directory submits to PBS Pro on [Gadi](#gadi). For another scheduler:

1. Install the site's Snakemake executor plugin in the environment that runs Snakemake 9.27.1 and write a site profile that selects it. The [executor catalogue](https://snakemake.github.io/snakemake-plugin-catalog/) documents plugins, including Slurm and cluster-generic.
2. Pass the site directory with the launcher's `--profile DIR` option, before `--`. The launcher combines it with `profiles/workflow` for per-rule resources. Pass `--workflow-profile DIR` to replace those settings. Each profile directory needs a `config.yaml`. Keep `software-deployment-method: [conda]` for release and report rules.
3. Load the container runtime on compute nodes and pass the runtime selection and required bind paths into their job environments. Keep the checkout, inputs, resources, outputs and conda cache visible at the same paths.
4. Review CPU, memory, walltime and storage requests against the site's scheduler mappings. On nodes without network access, run setup on a host with access first.
5. Run both fresh system modes from the [test guide](../test/README.md) and check scheduler accounting before scaling up.

Jobs submitted through a site profile other than the bundled NCI directory receive no module loading or bind setup from the launcher. The profile or the job environment must supply the runtime and binds.

### Resources on a cluster

Run [a few libraries](#samplesheet) first to estimate the resources needed for the full dataset. No real-data benchmark for memory, runtime or output size is recorded here.

For jobs submitted separately, `profiles/workflow/config.yaml` requests up to 8 cores, 96 GB and 16 hours per library for TRUST4, and 8 cores, 32 GB and 12 hours for annotation. These values set the resources requested from the scheduler. `profiles/local/config.yaml` has smaller requests for scheduling jobs within one host or allocation.

After this first run, inspect the scheduler accounting summaries in `logs/joblogs/` and the benchmarks in `logs/airrborne/benchmarks/`. These record the resources the jobs actually used.

### Gadi

On Gadi login and `copyq` nodes, the launcher selects the bundled PBS profile automatically. It loads Singularity and submits each workflow rule as a PBS job. Prepare the software in a `copyq` job, then launch the workflow from a login node. Compute nodes have no network access, and login-node limits can interrupt setup.

The automatic selection checks for `nci` in the hostname. To select the profile explicitly, pass `--profile profiles/global/nci`. Unless `AIRRBORNE_CONTAINER_BINDS` is set, the launcher mounts the existing `/g/data/$PROJECT` and `/scratch/$PROJECT` directories into containers.

1. **Clone and set the project** (login node). Put the clone under `/scratch/$PROJECT` or `/g/data/$PROJECT`, where compute jobs can read it. Clone with `--recursive`, or run `git submodule update --init --recursive` in an existing clone.

   `PROJECT` selects the project charged for jobs. Check it with `echo "$PROJECT"` and export the correct value in the shell that will launch the workflow. If conda is not installed, follow [prepare the environment](#prepare-the-environment), using `/g/data/$PROJECT` for Miniforge to avoid the small home-directory quota.
2. **Set the conda location** (login node). The shared PBS profile uses a conda directory in the maintainers' project, a56. For another project, uncomment this line in `profiles/workflow/config.yaml` before setup:

   ```yaml
   conda-prefix: .snakemake/conda
   ```

   This creates the release environment inside the clone and overrides the shared profile's `conda-prefix`. Leave the shared setting unchanged.

   Jobs can access your project's `/g/data` and `/scratch` directories by default. If inputs or software are in another project, edit `STORAGE` in `profiles/global/nci/config.yaml` to include every required project, for example `gdata/vb89+gdata/qz73`.
3. **Configure** (login node). Prepare the [inputs](#prepare-inputs) and set the paths in `config/airrborne/config.yaml`. Point `deps.data` at the existing FASTQ directories on project storage. Do this before setup, because setup reads the samplesheet and checks its inputs.
4. **Prepare the software** (`copyq`). Include `--tests` to prepare for the deployment check in the next step:

   ```bash
   qsub -q copyq -P "$PROJECT" -W block=true -j oe -o setup.log \
     -l ncpus=1,mem=32GB,walltime=04:00:00,storage=gdata/$PROJECT+scratch/$PROJECT,wd <<'EOF'
   # Add any source or module load line that conda needs at login.
   ./setup.sh --tests --no-browsers
   EOF
   ```

   Add any conda setup command that you normally type at login to the job script. Login scripts run automatically, but commands entered by hand do not.

   `-W block=true` keeps the terminal waiting until the job ends. If the connection drops, the job continues; check it with `qstat -u "$USER"`. Successful setup ends with `Setup complete` in `setup.log`.

   Set `storage` to include every project holding the clone, inputs or conda cache. Add `--post-trust4` when reusing assemblies without FASTQs. If an image download runs out of space, follow the [temporary-storage troubleshooting advice](#troubleshooting).
5. **Check the deployment** (login node). This replaces the README quick start on Gadi. It runs the public synthetic example through PBS and compares the results with the expected answers:

   ```bash
   source .snakemake/env.sh
   ./test/system/run.sh post-trust4
   ```

   Success prints `System test passed: post-trust4 (public)`. To check assembly from FASTQs as well, run `./test/system/run.sh fastq`. The [acceptance record](maintenance.md#gadi-acceptance-record) gives measured durations for setup and both tests.
6. **Launch** (login node). Preview the plan, then run inside `tmux` or `screen`, because the main process stays up for the whole run:

   ```bash
   ./run_module.sh -n
   ./run_module.sh
   ```

   Watch the jobs with `qstat -u "$USER"`. Job output and PBS accounting summaries are under `logs/joblogs/<rule>/`. Rule logs are under `logs/airrborne/`.

   Keep the launcher running while the submitted jobs finish. A `tmux` session stays on the login node where it started, so record that node's name for reconnecting. If the launcher stops, rerun the same command to continue from completed jobs.
7. **Read the results** (workstation). Copy the whole `scopes/<scope>/release/` folder under `outs.results_dir`, such as with `rsync -a` or `scp -r`, and open its `report/index.html` in a browser. To view the report without copying it, serve it from the login node, as [open the report](#open-the-report) describes.

Gadi charges each job for the walltime it used at a rate set by the larger of its core count and its share of the node's memory; the [resource notes](#resources-on-a-cluster) give the per-rule bounds, and `profiles/global/nci/queues.yaml` holds the queue rates.

To run inside one allocation instead, create the environment for local execution once from a login node with `./run_module.sh --profile none -- --conda-create-envs-only`, or switch on the `conda-prefix` line in `profiles/workflow/config.yaml`, then launch with `--profile none` inside the allocation. Maintainers record validated Gadi runs in the [maintenance guide](maintenance.md#gadi-acceptance-record).

## Outputs

`outs.results_dir` holds these directories:

```text
libraries/<library_id>/    shared assembly, annotation and library measurements
scopes/<scope>/           working tables for one scope
scopes/<scope>/release/   portable results for that scope
shared/provenance/        reference audit and tool-version records
```

Every scope's `release/` has the same layout:

```text
report/index.html     the offline report
metrics/              per-library and per-pair tables, such as headline.tsv,
                      diversity.tsv, isotype_usage.tsv and donor_overlap_pairs.tsv
metrics/detail/       per-clone, per-member and per-clonotype tables, such as
                      clone_members.tsv and shared_groups.tsv
igh_clone_threshold/  the operative IGH threshold and each donor's fit status
airr/                 one annotation table per library and the repertoire metadata
manifest/             the libraries of the scope and the size and hash of each file
README.md, data_dictionary.csv, provenance.json, effective_config.yaml
```

The release README lists the files and explains how to join them. The [interpretation guide](interpretation.md#find-the-table-for-a-question) points to tables for common questions.

After changing an analysis setting, rerun the workflow before using the release. Publication stops if libraries in a scope were processed with different code or settings. If publication fails, fix the reported problem and relaunch.

The release contains the tables needed for downstream analysis. To rerun from assemblies, also keep the native assembly and depth files listed under [reuse retained assemblies](#reuse-retained-assemblies).

### Open the report

Open `report/index.html` in a web browser. It needs no installation, server or network connection. Keep the surrounding files together, because the page loads data and assets from them. The `report/details/` directory contains a detail file for each library.

- On the machine that holds the results, open `report/index.html` in a browser, by double-clicking it or with the browser's File > Open.
- On a remote host, copy the whole `release/` folder to a workstation with `rsync -a` or `scp -r` and open it there, or serve the folder in place:

  ```bash
  cd scopes/<scope>/release/report
  python3 -m http.server 8765 --bind 127.0.0.1
  ```

  Start the server on the remote host. On the workstation, run `ssh -L 8765:localhost:8765 <host>`, then open `http://localhost:8765` in a browser. Stop the remote server with Ctrl-C when finished.

  Binding to `127.0.0.1` prevents direct connections through other network interfaces, but other users of the same remote host can still access the server. VS Code Remote-SSH can forward the port automatically.

### Copy and verify a portable release

Copy the entire `release/` folder when moving or sharing results. Published tables can be hard links to the working results, so save edited versions separately.

The example below uses `rsync -aL` to make an independent copy in a new temporary directory. The `-L` option also copies the original files if a data-management checkout uses symbolic links. It then checks every copied file against the release's recorded sizes and SHA-256 hashes:

```bash
source_release=results/airrborne/scopes/cohort/release
portable_release=$(mktemp -d)
rsync -aL "$source_release/" "$portable_release/"
cd "$portable_release"
python3 - <<'PYCODE'
import csv
import hashlib
import json
from pathlib import Path

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

manifest = Path("manifest/release_files.tsv")
receipt = json.loads(Path(".complete.json").read_text())
assert sha256(manifest) == receipt["release_files_sha256"], "Manifest hash mismatch"
with manifest.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        path = Path(row["path"])
        assert path.stat().st_size == int(row["n_bytes"]), path
        assert sha256(path) == row["sha256"], path
print("Verified copied release sizes and SHA-256 hashes")
PYCODE
```

Replace `source_release` with the actual release path before running the example. Open the copied `report/index.html` in a browser and follow a library to its detail view to check that its files load. To view it on another machine, copy the entire release there first.

AIRR-to-Change-O conversion also requires the matching pipeline checkout and reference. See the [reuse guidance](interpretation.md#airr-files).

## Troubleshooting

- **A rule fails.** Read its log under `outs.log_dir` and check the named inputs. Empty repertoires have valid outputs; missing or malformed inputs stop the run.
- **Snakemake reports a locked directory.** Check that no workflow is still running in that directory. Then run the launcher with `-u` and relaunch the original command.
- **A container cannot read an input.** Add its containing directory to `AIRRBORNE_CONTAINER_BINDS`, as described under [container selection](#container-selection).
- **A run was interrupted.** Relaunch the same command. Finished jobs are kept and the rest run again.
- **Setup was ended by a login-node limit.** Run it again; it reuses the environment and images it finished. If the limit ends it again, run it where a long process is allowed, such as a `copyq` job on Gadi.
- **A job ran out of memory.** The host or scheduler ends a job that exceeds the available memory. Raise that rule's `mem_mb` in `profiles/local/config.yaml` for local execution, or in `profiles/workflow/config.yaml` for submitted jobs, and relaunch.
- **A test reports a missing Python module.** Run setup once with `--tests`, then `source .snakemake/env.sh` in the shell that runs the tests.
- **An R test cannot start.** R tests use the configured analysis image. Rerun `./setup.sh --config test/system/config.yaml --tests`, source `.snakemake/env.sh` and check that the selected container runtime is available.
- **The wrong runtime runs.** Set `AIRRBORNE_CONTAINER_RUNTIME` to `apptainer`, `singularity` or `docker`.
- **The image pull runs out of space.** Point `TMPDIR`, or `APPTAINER_TMPDIR` under Apptainer and `SINGULARITY_TMPDIR` under Singularity, at storage with room for the downloaded layers and SIF. The [setup table](#prepare-the-environment) lists their sizes.
