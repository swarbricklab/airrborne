#!/usr/bin/env bash
# Run every R test in the configured pinned analysis image.
# AIRRBORNE_TEST_CONFIG selects the image prepared by root setup.sh.
# The checkout root reaches the container explicitly, including with cleanenv.
set -euo pipefail
fixture_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fixture_config="${AIRRBORNE_TEST_CONFIG:-${fixture_root}/test/system/config.yaml}"
case ",${AIRRBORNE_CONTAINER_BINDS:-}," in
    *,"$fixture_root",*) ;;
    *) export AIRRBORNE_CONTAINER_BINDS="${fixture_root}${AIRRBORNE_CONTAINER_BINDS:+,${AIRRBORNE_CONTAINER_BINDS}}" ;;
esac
image_arguments="$(bash "${fixture_root}/test/support/analysis_image.sh" "${fixture_config}")"
images=()
while IFS= read -r value; do images+=("$value"); done <<< "${image_arguments}"
mkdir -p "${fixture_root}/.snakemake/test-tmp"
exec "${fixture_root}/workflow/scripts/container_exec.sh" "${images[0]}" "${images[1]}" \
    env "AIRRBORNE_ROOT=${fixture_root}" "TMPDIR=${fixture_root}/.snakemake/test-tmp" Rscript "$@"
