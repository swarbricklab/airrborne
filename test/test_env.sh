#!/usr/bin/env bash
# Prepare the release-derived Python environment for scientific tests.
# Reuse matching definitions; --rebuild replaces only this environment.
# Root setup.sh owns activation. --no-browsers omits Chromium downloads.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
prefix="${root}/.snakemake/test-env"
browsers=true
rebuild=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-browsers) browsers=false ;;
        --rebuild) rebuild=true ;;
        -h|--help) echo "Usage: test/test_env.sh [--no-browsers] [--rebuild]"; exit 0 ;;
        *) echo "Unknown test environment option: $1" >&2; exit 2 ;;
    esac
    shift
done
platform="$(conda info | awk '$1 == "platform" {print $3}')"
[[ -n "$platform" ]] || { echo "conda info does not report a platform." >&2; exit 1; }
pin="${root}/envs/release.${platform}.pin.txt"
packages=(jsonschema=4.26.0 nodejs=24)
playwright_version=1.63.0
base_python="$(conda info --base)/bin/python"
fingerprint="$("$base_python" - "$root" "$platform" "$pin" "${packages[@]}" "$playwright_version" <<'PY'
import hashlib
from pathlib import Path
import sys
root, platform, pin, *packages = sys.argv[1:]
digest = hashlib.sha256((platform + "\n" + "\n".join(packages)).encode())
digest.update((Path(root) / "envs/release.yaml").read_bytes())
if Path(pin).is_file():
    digest.update(Path(pin).read_bytes())
print(digest.hexdigest())
PY
)"
if [[ "$rebuild" == false && -x "${prefix}/bin/python" && -f "${prefix}/.inputs.sha256" && "$(cat "${prefix}/.inputs.sha256")" == "$fingerprint" ]]; then
    echo "Reusing test environment: $prefix"
else
    # Workflow environments, downloads and run trees have separate owners.
    rm -rf "$prefix"
    solver_options=()
    if [[ -f "$pin" ]]; then
        conda create -y -q -p "$prefix" --file "$pin"
    else
        conda env create -q -p "$prefix" -f "${root}/envs/release.yaml"
        # Conda 26.1.1 libmamba misreads the channel-qualified AIRR history entry.
        solver_options=(--solver classic)
    fi
    CONDA_AGGRESSIVE_UPDATE_PACKAGES='' conda install -y -q -p "$prefix" \
        "${solver_options[@]}" --freeze-installed --override-channels \
        -c conda-forge -c bioconda "${packages[@]}"
    "${prefix}/bin/python" -m pip install --quiet "playwright==${playwright_version}"
    printf '%s\n' "$fingerprint" > "${prefix}/.inputs.sha256"
fi
"${prefix}/bin/python" -c 'import airr, jsonschema, playwright, yaml'
"${prefix}/bin/node" --version
if [[ "$browsers" == true ]]; then
    "${prefix}/bin/python" -m playwright install chromium
fi
printf 'Test environment ready: %s\n' "$prefix"
