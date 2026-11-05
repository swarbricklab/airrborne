# Tests

Use this guide to check changes to the pipeline. To check a first installation, follow the [README quick start](../README.md#quick-start) on a workstation or the [Gadi deployment check](../docs/running.md#gadi).

Component tests check individual pipeline decisions. System tests run the workflow on public synthetic data and compare the results with reviewed expectations. The [glossary](../docs/glossary.md) explains AIRR terms; the [data guide](data/README.md) describes the inputs and expected results. The [coverage map](COVERAGE.md) links scientific decisions to their tests.

## Run your first test

On Linux or macOS, install `git`, `conda` and a container runtime: Apptainer, Singularity or Docker on Linux, or Docker Desktop on macOS. Start Docker if that is the selected runtime. Setup needs network access. R and its packages are supplied by the pinned analysis image.

Run from the repository root. If the root README quick start has already passed, continue with [component tests](#run-a-component). On a cluster, follow [test on a cluster](#test-on-a-cluster).

```bash
git submodule update --init --recursive
./setup.sh --config test/system/config.yaml --tests --no-browsers
source .snakemake/env.sh
./test/system/run.sh post-trust4 --work-dir .snakemake/quickstart -- --cores 2
```

This test starts from the TRUST4 assemblies included in the repository: four libraries from two donors. A successful run prints `System test passed: post-trust4 (public)`. Open the report in a web browser:

```text
.snakemake/quickstart/results/airrborne/scopes/project_a/release/report/index.html
```

The work directory must be absent or empty. The runner retains its inputs, results and logs. Omit `--work-dir` on subsequent runs to get a fresh directory automatically. Browser automation is part of the complete checks below.

## Run a component

Use the prepared environment in each new terminal:

```bash
source .snakemake/env.sh
python test/run_tests.py --component shm
python test/run_tests.py --component annotation --level integration
```

Component names match the directories in [components/](components/) and appear in `python test/run_tests.py --help`. Omitting `--level` runs both levels:

- `unit` checks functions and scientific boundaries, including pinned packages.
- `integration` checks scripts, external tools and file handoffs.

Some components have only one runnable level. Files named `reference_calculations.py` or `reference_records.py` support the [system result checker](system/check_results.py); they do not run as component tests. System runs do not repeat the component suites.

## Check a complete change

For a complete check of a pipeline change, run the component tests and system cases below. Prepare Chromium for the optional offline-report check:

```bash
source .snakemake/env.sh
bash test/test_env.sh
```

On Linux, install its system libraries if needed:

```bash
python -m playwright install-deps chromium
```

Run the components and each system case. System commands choose fresh output directories under `.snakemake/test-runs/` and print the report paths:

```bash
python test/run_tests.py
./test/system/run.sh post-trust4 --check-report -- --cores 2
./test/system/run.sh fastq -- --cores 2
./test/system/run.sh post-trust4 --case empty -- --cores 2
./test/system/run.sh post-trust4 --case permuted -- --cores 2
```

| System case | What it checks |
|---|---|
| `post-trust4` | Annotation through publication, starting from native assemblies and depth records. |
| `fastq` | The same decisions, plus read discovery, depth counting and TRUST4 assembly. |
| `post-trust4 --case empty` | Empty and all-failing libraries, positive analysis minima and an empty published scope. |
| `post-trust4 --case permuted` | Two copies of one library that differ only in TRUST4 CDR3 variant numbering, with one constructed support tie, give identical released values. |

Both public modes check the released tables, AIRR links, reviewed record decisions and independent metric calculations. The empty and permuted cases have their own assertions.

`--check-report` copies the release and opens it offline in Chromium. A successful first run can count as post-TRUST4 acceptance for the same source and inputs. If it omitted browser testing, run the [copied-report check](components/report/integration/check_offline_report.py) with `--release` pointing to that run's release folder.

Options after `--` go to `run_module.sh`. The launcher forwards Snakemake arguments such as `--cores 2`, which limits local CPU use. The system runner supplies the test configuration, starting mode and complete set of scopes. To run a selected scope, use the production launcher.

Pull requests to `main` run component tests and post-TRUST4 acceptance. The [CI path filter](../.github/scripts/fastq_required.sh) selects changes that also need FASTQ acceptance.

Before a release, manually run **public tests** in GitHub Actions for the intended commit. The `contracts`, `post-trust4` and `fastq` jobs must pass. CI and Gadi/PBS acceptance each require checks in their own environment.

## Investigate a failure

| Failure | Where to look |
|---|---|
| Missing command, package or container | Source `.snakemake/env.sh` in the current terminal. If setup is incomplete, rerun the setup command above. |
| Component assertion | Read the reported test name and expected/observed values. Rerun its component and level. |
| Workflow job | Open `workflow.log` in the printed work directory, then the job log named by the error. |
| Validation after the workflow | Keep the terminal error and inspect its named output table, record or column. These validation messages are not in `workflow.log`. |
| Work directory already exists | Omit `--work-dir` for a fresh run. Keep the previous output for comparison. |

For a reviewed-value mismatch, use the table, row selector, column and `review_basis` in [expected_values.tsv](data/expected/expected_values.tsv). Record outcomes are in [expected_records.tsv](data/expected/expected_records.tsv). Investigate the decision before changing its expected answer. Change expectations only for a reviewed scientific change, following the [data promotion instructions](data/README.md#promote-reviewed-inputs).

## Add or change a test

Extend a nearby test when it already owns the behaviour. Otherwise, put a new file directly in `components/<component>/unit/` or `integration/`. The runner recognises these filenames; it does not search deeper subdirectories:

| Language | Filename |
|---|---|
| Python | `test_*.py` |
| R | `test-*.R` |
| JavaScript | `*.test.js` |
| Shell | `test_*.sh` |

Keep constructed inputs small and explain the biological case and expected answer at the top. Use independent arithmetic or reviewed evidence for the answer, with explicit numerical tolerances where needed. Run the component command while editing, then the complete checks appropriate to the change.

The [diversity unit tests](components/diversity/unit/test-observed_diversity.R) show an R example; [input depth tests](components/inputs/unit/test_input_depth.py) show Python. R tests call `load_test_modules()` to name the production modules they need. The runner supplies that helper and the pinned R environment.

| Location | What belongs there |
|---|---|
| `components/` | Scientific tests and independent result calculations, grouped by pipeline component. |
| `data/` | Committed inputs, construction evidence and reviewed expectations. |
| `data_builders/` | Scripts that generate new data candidates or the empty-library case. |
| `system/` | Workflow staging, execution and result checks. |
| `support/` | Shared execution and comparison helpers. |

Normal testing needs no generator environment. Follow the [data guide](data/README.md) to inspect or deliberately change committed inputs.

`test/test_env.sh` adds test tools to the release-derived Python environment and reuses matching definitions. Use `bash test/test_env.sh --rebuild` if that environment needs replacement. Workflow environments, images and retained runs have separate locations.

## Test on a cluster

Maintainers validating a deployment prepare the software on a host with network access, run the component tests inside a compute allocation, and run the system cases through the site's scheduler from a submission host. The [Gadi validation](../docs/maintenance.md#validate-a-gadi-deployment) gives the commands for the bundled PBS profile. Browser checks belong on a workstation or in CI.
