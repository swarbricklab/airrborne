#!/usr/bin/env bash
# Run the committed synthetic inputs through the production launcher and check
# the scientific results. Both modes retain their fresh outputs for inspection.
set -eo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$root"
usage() {
    cat <<'HELP'
Usage: test/system/run.sh post-trust4|fastq [--work-dir DIR] [--case empty|permuted]
                          [--check-report] [-- LAUNCHER_OPTIONS...]

post-trust4 starts from committed native assemblies and depth measurements.
fastq starts from paired reads in a separate fresh output directory.
--case empty selects the constructed empty/all-failing post-TRUST4 scenario.
--case permuted selects two copies of one library that differ only in TRUST4
  CDR3 variant numbering, and checks that their results agree.
--check-report opens a copied public release offline in Chromium.
Execution options after -- go to run_module.sh, for example:
  test/system/run.sh post-trust4 -- --profile none -- --cores 2

Prepare tests with root setup.sh --config test/system/config.yaml --tests,
then source .snakemake/env.sh. Component tests run through test/run_tests.py.
HELP
}
mode="${1:-}"
case "$mode" in
    -h|--help) usage; exit 0 ;;
    post-trust4|fastq) shift ;;
    *) usage >&2; exit 2 ;;
esac
work=""
scenario=public
browser=false
launcher_options=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --work-dir|--case)
            [[ $# -ge 2 && -n "$2" ]] || { echo "$1 needs a value." >&2; exit 2; }
            if [[ "$1" == --work-dir ]]; then work="$2"; else scenario="$2"; fi
            shift 2 ;;
        --check-report) browser=true; shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; launcher_options=("$@"); break ;;
        *) echo "Unknown system-test option: $1" >&2; exit 2 ;;
    esac
done
# The system test owns its input, starting mode and complete scope population.
# A launcher override must not run another configuration outside the fresh tree.
for option in "${launcher_options[@]}"; do
    case "$option" in
        -c|--config|--config=*|--configfile|--configfile=*|--configfiles|--configfiles=*|\
        --scope|--scope=*|--post-trust4|-s|--snakefile|--snakefile=*|\
        --directory|--directory=*|-d|--dag)
            echo "The system test owns configuration, mode, scopes and workflow paths: $option" >&2
            exit 2 ;;
    esac
done
case "$scenario" in
    public) ;;
    empty|permuted)
        [[ "$mode" == post-trust4 && "$browser" == false ]] || {
            echo "The $scenario case requires post-trust4 without --check-report." >&2; exit 2;
        } ;;
    *) echo "The constructed --case values are empty and permuted." >&2; exit 2 ;;
esac
python "${root}/test/system/check_inputs.py"
if [[ -z "$work" ]]; then
    mkdir -p "${root}/.snakemake/test-runs"
    work="$(mktemp -d "${root}/.snakemake/test-runs/${mode}-${scenario}.XXXXXX")"
fi
work="$(python - "$work" <<'PY'
from pathlib import Path
import sys
path = Path(sys.argv[1]).resolve()
if path.exists() and (not path.is_dir() or any(path.iterdir())):
    raise SystemExit(f"Fresh system tests need an absent or empty directory: {path}")
path.mkdir(parents=True, exist_ok=True)
print(path)
PY
)"
# Reuse the images, receipts and reference that setup.sh prepared, which
# .snakemake/env.sh records; otherwise the test configuration's own directory
# is used. Directories outside the checkout are mounted into the containers.
resources_option=()
if [[ -n "${AIRRBORNE_TEST_RESOURCES:-}" ]]; then
    resources_option=(--resources "${AIRRBORNE_TEST_RESOURCES}")
fi
for outside in "$work" "${AIRRBORNE_TEST_RESOURCES:-}"; do
    case "$outside" in
        ""|"$root"/*) ;;
        *) export AIRRBORNE_CONTAINER_BINDS="${outside}${AIRRBORNE_CONTAINER_BINDS:+,${AIRRBORNE_CONTAINER_BINDS}}" ;;
    esac
done
trap 'test_status=$?; if [[ "$test_status" -ne 0 ]]; then printf "System test failed. Inspect %s\n" "$work" >&2; fi; exit "$test_status"' EXIT
config="${work}/config.yaml"
if [[ "$scenario" == empty ]]; then
    python "${root}/test/data_builders/empty_libraries.py" --work-dir "$work" "${resources_option[@]}"
elif [[ "$scenario" == permuted ]]; then
    python "${root}/test/data_builders/permuted_indices.py" --work-dir "$work" "${resources_option[@]}"
else
    seed=()
    if [[ "$mode" == post-trust4 ]]; then seed=(--seed-results "${root}/test/data/post_trust4"); fi
    python "${root}/test/system/stage_inputs.py" --config "${root}/test/system/config.yaml" \
        --work-dir "$work" --out-config "$config" "${seed[@]}" "${resources_option[@]}"
fi
mode_options=()
if [[ "$mode" == post-trust4 ]]; then mode_options=(--post-trust4); fi
# run_module.sh owns profile selection, Conda activation and workflow execution.
"${root}/run_module.sh" --config "$config" "${mode_options[@]}" "${launcher_options[@]}" \
    2>&1 | tee "${work}/workflow.log"
results="${work}/results/airrborne"
python "${root}/test/components/release/integration/check_tables.py" \
    --results "$results" --dictionary "${root}/docs/data_dictionary.csv"
bundles=("${results}"/scopes/*/release)
[[ -d "${bundles[0]}" ]] || { echo "The workflow did not publish a release." >&2; exit 1; }
for bundle in "${bundles[@]}"; do
    python "${root}/test/components/release/integration/check_release.py" --bundle "$bundle" --results "$results"
done
expectations=()
if [[ "$scenario" == public ]]; then
    expectations=(--expected-values "${root}/test/data/expected/expected_values.tsv"
                  --expected-records "${root}/test/data/expected/expected_records.tsv")
fi
python "${root}/test/system/check_results.py" --results "$results" --config "$config" "${expectations[@]}"
if [[ "$scenario" == empty ]]; then
    python "${root}/test/system/check_empty_libraries.py" --work-dir "$work"
elif [[ "$scenario" == permuted ]]; then
    python "${root}/test/system/check_permuted.py" --work-dir "$work"
fi
if [[ "$browser" == true ]]; then
    python "${root}/test/components/report/integration/check_offline_report.py" \
        --release "${results}/scopes/project_a/release"
fi
printf '\nSystem test passed: %s (%s)\nResults: %s\nWorkflow log: %s/workflow.log\n' "$mode" "$scenario" "$results" "$work"
for bundle in "${bundles[@]}"; do printf 'Report: %s/report/index.html\n' "$bundle"; done
