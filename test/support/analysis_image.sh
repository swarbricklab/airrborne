#!/usr/bin/env bash
# Print the two image arguments container_exec.sh takes for the analysis image:
# the local SIF the workflow obtains for Apptainer and Singularity, then the
# pinned reference Docker runs. The SIF is under AIRRBORNE_TEST_RESOURCES, the
# directory setup.sh prepared, or else under the CONFIG's deps.resources, read
# relative to the caller's working directory as Snakemake reads it.
set -euo pipefail
if [[ $# -ne 1 ]]; then
    echo "usage: $0 CONFIG_YAML" >&2
    exit 2
fi
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
python - "${root}" "$1" <<'PY'
import os
import sys
from pathlib import Path
import yaml
root, config = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(root))
from workflow.lib.provenance import image_sif_path
deps = yaml.safe_load(config.read_text())["deps"]
resources = Path(os.environ.get("AIRRBORNE_TEST_RESOURCES") or deps["resources"]).resolve()
print(image_sif_path(resources, deps["immcantation_docker"]))
print(deps["immcantation_docker"])
PY
