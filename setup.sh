#!/usr/bin/env bash
# Prepare the pinned launcher, both configured images and the release
# environment, so a later launch needs no network access. --tests also prepares
# the separate test environment. Configuration paths use the caller's
# directory, as they do in run_module.sh.
set -eo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config="config/airrborne/config.yaml"
tests=false
browser_options=()
mode_options=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)
            [[ $# -ge 2 && -n "$2" ]] || { echo "--config needs a file." >&2; exit 2; }
            config="$2"; shift 2 ;;
        --post-trust4) mode_options=(--post-trust4); shift ;;
        --tests) tests=true; shift ;;
        --no-browsers) browser_options=(--no-browsers); shift ;;
        -h|--help)
            cat <<'HELP'
Usage: ./setup.sh [-c|--config FILE] [--post-trust4] [--tests] [--no-browsers]

Prepare everything a launch needs, so the launch itself needs no network
access: the snakemake_9.27.1 conda environment, both configured container
images with their identity receipts under deps.resources, and the release
environment for the way run_module.sh executes on this host.

Run it from the launch directory, the clone's root or an analysis repository's
root, where the configuration path resolves. Run it with conda on PATH, on a
host with network access that allows a long process, such as a Gadi copyq job.
docs/running.md lists its duration and disk use. It reads the configured
samplesheet and the inputs it names, so they must already exist.

It selects the runtime named by AIRRBORNE_CONTAINER_RUNTIME, or else the first
of Apptainer, Singularity and Docker on PATH. When the selected runtime is not
on PATH, it loads the environment module of the same name, singularity when
none is selected. Image layers are unpacked under TMPDIR, or APPTAINER_TMPDIR
or SINGULARITY_TMPDIR. Running setup again reuses an existing snakemake_9.27.1
environment and the images already pulled.

  -c, --config FILE  Configuration file (default: config/airrborne/config.yaml)
  --post-trust4      Resolve the workflow from retained assemblies instead of
                     FASTQs, as run_module.sh --post-trust4 does
  --tests            Also prepare the test environment and write
                     .snakemake/env.sh, which a shell sources before running tests
  --no-browsers      With --tests, skip the Chromium download that only the
                     report browser checks use
HELP
            exit 0 ;;
        *) echo "Unknown setup option: $1" >&2; exit 2 ;;
    esac
done
[[ -f "$config" ]] || { echo "Configuration file is absent: $config" >&2; exit 2; }
if ! conda run --no-capture-output -n snakemake_9.27.1 snakemake --version >/dev/null 2>&1; then
    conda env create -f "${root}/envs/snakemake_9.27.1.yaml"
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snakemake_9.27.1
[[ "$(snakemake --version)" == 9.27.1 ]] || { echo "Snakemake 9.27.1 is required." >&2; exit 1; }
config="$(python -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' "$config")"
# The resources directory, resolved as the Snakefile resolves it, then the
# identity receipts the workflow's rules declare. A --post-trust4 workflow has
# no assembly rule, so it obtains only the analysis image.
image_targets="$(python - "$root" "$config" "${#mode_options[@]}" <<'PY'
from pathlib import Path
import sys
import yaml
sys.path.insert(0, sys.argv[1])
from workflow.lib.provenance import image_receipt_path
deps = yaml.safe_load(Path(sys.argv[2]).read_text())["deps"]
resources = Path(deps["resources"]).resolve()
print(resources)
names = ("immcantation_docker",) if sys.argv[3] != "0" else ("immcantation_docker", "trust4_container")
for name in names:
    print(image_receipt_path(resources, deps[name].removeprefix("docker://")))
PY
)"
targets=()
while IFS= read -r value; do targets+=("$value"); done <<< "$image_targets"
resources="${targets[0]}"
targets=("${targets[@]:1}")
runtime="${AIRRBORNE_CONTAINER_RUNTIME:-}"
if [[ -z "$runtime" && "$(uname -s)" == Darwin ]]; then
    runtime=docker
fi
if [[ -z "$runtime" ]]; then
    for candidate in apptainer singularity docker; do
        if command -v "$candidate" >/dev/null 2>&1; then runtime="$candidate"; break; fi
    done
fi
# A cluster can supply the runtime as an environment module named after it.
# env.sh repeats a load made here, so a new shell gets the same runtime.
loaded_module=""
if [[ "$runtime" != docker ]] && ! command -v "${runtime:-singularity}" >/dev/null 2>&1 \
        && command -v module >/dev/null 2>&1 && module load "${runtime:-singularity}" 2>/dev/null; then
    runtime="${runtime:-singularity}"
    loaded_module="$runtime"
fi
case "$runtime" in
    apptainer|singularity)
        command -v "$runtime" >/dev/null 2>&1 || { echo "$runtime is not on PATH." >&2; exit 2; } ;;
    docker) docker info >/dev/null ;;
    *)
        echo "No container runtime is on PATH. Install Apptainer, Singularity or Docker," \
             "or select one with AIRRBORNE_CONTAINER_RUNTIME." >&2
        exit 2 ;;
esac
export AIRRBORNE_CONTAINER_RUNTIME="$runtime"
# Obtain the images through the workflow's own local rules, so Snakemake
# records how each receipt was produced and can rerun the rule when its code,
# settings or runtime change. run_module.sh always adds the whole-workflow
# target in --post-trust4 mode, so this call names the receipts directly.
analysis_only=False
[[ ${#mode_options[@]} -eq 0 ]] || analysis_only=True
snakemake --snakefile "${root}/Snakefile" --configfile "$config" \
    --workflow-profile "${root}/profiles/local" --profile none --cores 1 \
    --config "analysis_only=${analysis_only}" -- "${targets[@]}"
# Create the release environment that the AIRR export, release and report
# rules declare, for the execution method run_module.sh selects on this host.
# Snakemake resolves the configured workflow first, so the samplesheet and the
# inputs it names must exist. A host that selects a site profile can prepare
# local execution as well with: run_module.sh --profile none -- --conda-create-envs-only
echo "Creating the release environment..."
"${root}/run_module.sh" -c "$config" "${mode_options[@]}" -- --conda-create-envs-only
if [[ "$tests" == true ]]; then
    bash "${root}/test/test_env.sh" "${browser_options[@]}"
fi
# The test environment stays selected once it exists, whichever configuration
# a later setup prepares. The tests reuse that configuration's images.
test_env="${root}/.snakemake/test-env"
{
    printf 'source %q\n' "$(conda info --base)/etc/profile.d/conda.sh"
    printf 'conda activate snakemake_9.27.1\n'
    if [[ -n "$loaded_module" ]]; then
        printf 'module load %q\n' "$loaded_module"
    fi
    printf 'export AIRRBORNE_CONTAINER_RUNTIME=%q\n' "$runtime"
    # Leave binds in the caller's environment so the production launcher's
    # automatic Gadi bindings still apply when no explicit binds were supplied.
    if [[ -x "${test_env}/bin/python" ]]; then
        mkdir -p "${root}/.snakemake/test-tmp"
        printf 'export PATH=%q:"$PATH"\n' "${test_env}/bin"
        printf 'export AIRRBORNE_TEST_CONFIG=%q\n' "$config"
        printf 'export AIRRBORNE_TEST_RESOURCES=%q\n' "$resources"
        printf 'export AIRRBORNE_ROOT=%q\n' "$root"
        printf 'export TMPDIR=%q\n' "${root}/.snakemake/test-tmp"
        printf 'export RSCRIPT=%q\n' "${root}/test/support/rscript.sh"
    fi
} > "${root}/.snakemake/env.sh"
printf 'Setup complete. Before running tests: source %q\n' "${root}/.snakemake/env.sh"
