#!/usr/bin/env bash
# Run a command in a pinned analysis or assembly image with Apptainer,
# Singularity or Docker.
#
# Usage: container_exec.sh <sif_path> <docker_image> <cmd> [args...]
# Apptainer and Singularity run <sif_path>. Docker runs <docker_image> as the
# calling user on linux/amd64. The run_trust4 rule passes ASSEMBLY_SIF as
# <sif_path> and ASSEMBLY_IMAGE_REFERENCE as <docker_image> after
# obtain_assembly_image prepares the image for the selected runtime.
#
# AIRRBORNE_CONTAINER_RUNTIME names the runtime. Otherwise the first of apptainer,
# singularity and docker on PATH is used, and the script exits with status 127
# when none is. The launch directory is mounted at its own path, as is each
# entry of the comma-separated AIRRBORNE_CONTAINER_BINDS, which must be an existing
# absolute path.
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "usage: $0 <sif_path> <docker_image> <cmd> [args...]" >&2
    exit 2
fi

sif="$1"; image="$2"; shift 2

# Each path is mounted once. Singularity warns when a path repeats, which
# happens when AIRRBORNE_CONTAINER_BINDS names the launch directory.
bind_paths=("$PWD")
if [[ -n "${AIRRBORNE_CONTAINER_BINDS:-}" ]]; then
    IFS=',' read -r -a extra_binds <<< "${AIRRBORNE_CONTAINER_BINDS}"
    for bind_path in "${extra_binds[@]}"; do
        if [[ "${bind_path}" != /* || ! -e "${bind_path}" ]]; then
            echo "ERROR: AIRRBORNE_CONTAINER_BINDS entries must be existing absolute paths: ${bind_path}" >&2
            exit 2
        fi
        for known in "${bind_paths[@]}"; do
            [[ "${bind_path}" == "${known}" ]] && continue 2
        done
        bind_paths+=("${bind_path}")
    done
fi

singularity_bind_args=()
docker_bind_args=()
for bind_path in "${bind_paths[@]}"; do
    singularity_bind_args+=(--bind "${bind_path}:${bind_path}")
    docker_bind_args+=(-v "${bind_path}:${bind_path}")
done

runtime="${AIRRBORNE_CONTAINER_RUNTIME:-}"
if [[ -z "$runtime" ]]; then
    for candidate in apptainer singularity docker; do
        if command -v "$candidate" >/dev/null 2>&1; then
            runtime="$candidate"
            break
        fi
    done
fi
case "$runtime" in
    apptainer|singularity)
        exec "$runtime" exec --cleanenv "${singularity_bind_args[@]}" "$sif" "$@"
        ;;
    docker)
        exec docker run --rm --platform linux/amd64 \
            --user "$(id -u):$(id -g)" \
            "${docker_bind_args[@]}" -w "$PWD" "$image" "$@"
        ;;
    *)
        echo "ERROR: select an installed apptainer, singularity or docker runtime" >&2
        exit 127
        ;;
esac
