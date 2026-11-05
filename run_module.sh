#!/bin/bash

# Launch Snakefile with caller-relative configuration and bundled profiles.
# A site profile or the host selects local execution or cluster submission.
# --post-trust4 leaves FASTQ discovery and assembly out of the DAG.
# Snakefile owns scope validation.

set -eo pipefail  # Stop when a command or a pipeline fails.

# Paths and launch defaults.

# Resolve bundled files relative to this script. Snakemake keeps the caller's
# working directory, so relative configuration and data paths use that directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="${SCRIPT_DIR}"
SNAKEFILE_PATH="${WORKFLOW_ROOT}/Snakefile"
# Resolve this default relative to the caller's working directory. The parent
# workflow uses config/airrborne/config.yaml; override it with -c.
CONFIG_FILE_PATH="config/airrborne/config.yaml"
DRY_RUN=false
GENERATE_DAG=false
UNLOCK=false
POST_TRUST4=false
# Analysis-only mode excludes FASTQ and assembly rules when the DAG is loaded.
ANALYSIS_ONLY_ARG=False

POST_TRUST4_TARGET_RULE="post_trust4_all"
SCOPE=""
SITE_PROFILE=""
WORKFLOW_PROFILE=""
# Core and memory limits taken from a scheduler allocation for local execution.
ALLOCATION_ARGS=()

# Command-line reference.

show_help() {
    cat << 'EOF'
Launcher for the airrborne pipeline
===================================

USAGE:
    ./run_module.sh [OPTIONS] [-- EXTRA_SNAKEMAKE_ARGS...]

OPTIONS:
    -c, --config FILE       Read configuration relative to the launch directory
                            (default: config/airrborne/config.yaml)
    -n, --dry-run           List planned jobs without executing them
    -d, --dag               Maintainers only: redraw the committed
                            docs/graphs/workflow_dag.svg (needs Graphviz dot/tred),
                            then go on to run
    -u, --unlock            Unlock the launch directory, then exit
    --post-trust4           Start from retained assemblies and depth records; skip
                            FASTQ discovery, depth counting and assembly
    --scope NAME            Build only this declared scope
    --profile DIR           Use this site profile for execution, graphs and unlock
                            (use none for local execution, including on Gadi;
                            profiles/global/nci selects the Gadi setup explicitly)
    --workflow-profile DIR  Replace the bundled per-rule resource profile
    -h, --help              Show this help message

EXAMPLES:
    # Run from this repository's root
    ./run_module.sh

    # Preview planned jobs
    ./run_module.sh -n

    # Unlock the launch directory after checking that no run is active
    ./run_module.sh -u

    # Use the parent analysis repository's configuration from its root
    ./modules/airrborne/run_module.sh

    # Reuse retained assemblies and depth records
    ./run_module.sh --post-trust4

    # Pass extra arguments to Snakemake
    ./run_module.sh -- --cores 2

NOTES:
    - Uses Snakemake 9.27.1 on PATH, or activates the snakemake_9.27.1 conda
      environment with the conda on PATH
    - Runs every job on this host with profiles/local unless a site profile
      submits jobs to a scheduler. Arguments after -- go to Snakemake; the two
      that size a run are --cores N and --resources mem_mb=N
    - Inside a PBS or Slurm job, local execution takes --cores and the mem_mb
      limit from the allocation (PBS_NCPUS and PBS_VMEM, or SLURM_CPUS_PER_TASK
      and SLURM_MEM_PER_NODE or SLURM_MEM_PER_CPU) and prints them; values
      after -- replace them
    - Runs containers through workflow/scripts/container_exec.sh with Apptainer,
      Singularity or Docker. AIRRBORNE_CONTAINER_RUNTIME selects one, and
      AIRRBORNE_CONTAINER_BINDS names data directories outside the launch directory
    - --post-trust4 needs the retained TRUST4 assemblies and sequencing-depth
      files of the selected libraries; missing files stop the run
    - Extra Snakemake arguments apply only to the workflow invocation; graph
      and unlock commands do not receive them

    On NCI Gadi:
    - Hosts whose names contain nci submit each rule as a PBS job through
      profiles/global/nci with the profiles/workflow resources, and write each
      job's output under logs/joblogs/ in the launch directory. --profile none
      runs on the host instead, and --profile profiles/global/nci selects the
      Gadi setup on any host
    - Inside a PBS job on Gadi, add --profile none so the run uses the
      allocation instead of submitting further jobs
    - The Gadi setup loads the singularity module and binds existing
      /g/data/$PROJECT and /scratch/$PROJECT unless AIRRBORNE_CONTAINER_BINDS is set
    - profiles/global is a Git submodule: git submodule update --init --recursive

EOF
}

# Parse launcher options and retain Snakemake arguments.

# Unknown options also pass through. After --, every argument goes to Snakemake
# without further launcher parsing, including names that match launcher options.
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -c|--config)
            if [[ $# -lt 2 || -z "$2" || "$2" == -* ]]; then
                echo "ERROR: $1 requires a configuration file, such as config/airrborne/config.yaml." >&2
                exit 2
            fi
            CONFIG_FILE_PATH="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -d|--dag)
            GENERATE_DAG=true
            shift
            ;;
        -u|--unlock)
            UNLOCK=true
            shift
            ;;
        --scope)
            if [[ $# -lt 2 || -z "$2" || "$2" == -* ]]; then
                echo "ERROR: --scope requires a declared scope name. Omit it to build every scope." >&2
                exit 2
            fi
            SCOPE="$2"
            shift 2
            ;;
        --profile|--workflow-profile)
            if [[ $# -lt 2 || -z "$2" || "$2" == -* ]]; then
                echo "ERROR: $1 requires a profile directory or none." >&2
                exit 2
            fi
            if [[ "$1" == --profile ]]; then
                SITE_PROFILE="$2"
            else
                WORKFLOW_PROFILE="$2"
            fi
            shift 2
            ;;
        --post-trust4)
            POST_TRUST4=true
            ANALYSIS_ONLY_ARG=True
            shift
            ;;
        --)
            shift
            EXTRA_ARGS+=("$@")
            break
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

# Show launch settings.

echo "============================================================================"
echo "Launcher for the airrborne pipeline"
echo "============================================================================"
echo "Host:       $(hostname)"
echo "Snakefile:  ${SNAKEFILE_PATH}"
echo "Config:     ${CONFIG_FILE_PATH}"
echo "Dry-run:    ${DRY_RUN}"
echo "Unlock:     ${UNLOCK}"
echo "Post-TRUST4: ${POST_TRUST4}"
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    echo "Extra args: ${EXTRA_ARGS[*]}"
fi
echo "============================================================================"
echo ""

# Check launch files before activating the environment.

if [[ ! -f "${CONFIG_FILE_PATH}" ]]; then
    echo "ERROR: Configuration file not found: ${CONFIG_FILE_PATH}" >&2
    echo "Paths resolve from the launch directory. Copy config/airrborne/ from the pipeline or pass -c FILE." >&2
    exit 1
fi

if [[ ! -f "${SNAKEFILE_PATH}" ]]; then
    echo "ERROR: Snakefile not found: ${SNAKEFILE_PATH}"
    exit 1
fi

# The NCI executor configuration is supplied by the profiles/global submodule.
# A host whose name contains nci selects it unless --profile is given, and
# naming that bundled directory with --profile selects the same Gadi setup.
# Stop before environment setup or submission if a fresh clone has not
# initialised it.
USE_NCI_PROFILE=false
NCI_PROFILE_DIR="${WORKFLOW_ROOT}/profiles/global/nci"
if [[ -z "${SITE_PROFILE}" && "$(hostname)" == *"nci"* ]]; then
    USE_NCI_PROFILE=true
elif [[ -n "${SITE_PROFILE}" && -d "${SITE_PROFILE}" && -d "${NCI_PROFILE_DIR}" \
        && "$(cd "${SITE_PROFILE}" && pwd -P)" == "$(cd "${NCI_PROFILE_DIR}" && pwd -P)" ]]; then
    USE_NCI_PROFILE=true
fi
if [[ "${USE_NCI_PROFILE}" == true && ! -d "${NCI_PROFILE_DIR}" ]]; then
    echo "ERROR: NCI profile not found: ${WORKFLOW_ROOT}/profiles/global/nci" >&2
    echo "Initialise repository submodules, then retry:" >&2
    echo "  git submodule update --init --recursive" >&2
    exit 1
fi

# Use the pinned Snakemake on PATH, including caller-managed environments.
# Fall back to the bundled environment name when that version is unavailable.

if [[ "$(snakemake --version 2>/dev/null || true)" != "9.27.1" ]]; then
    echo "Activating snakemake_9.27.1 conda environment..."
    eval "$(conda shell.bash hook)"
    conda activate snakemake_9.27.1
fi
if [[ "$(snakemake --version)" != "9.27.1" ]]; then
    echo "ERROR: Snakemake 9.27.1 is required. Create envs/snakemake_9.27.1.yaml or activate an environment with that version." >&2
    exit 1
fi
# Host rules use the interpreter beside the selected Snakemake executable.
export PATH="$(dirname "$(command -v snakemake)"):$PATH"
echo ""

# Select profiles and prepare the execution environment.

# Keep each option and path as a separate argument, including checkout paths
# containing spaces. All three Snakemake commands share these profile options.
PROFILE_ARGS=()

if [[ "${USE_NCI_PROFILE}" == true ]]; then
    echo "Using the bundled NCI Gadi profile"

    global_profile_path="${NCI_PROFILE_DIR}"
    workflow_profile_path="${WORKFLOW_ROOT}/profiles/workflow"

    PROFILE_ARGS=(--profile "${global_profile_path}"
                  --workflow-profile "${workflow_profile_path}")

    echo "Loading singularity module..."
    module load singularity

    # Jobs inherit this environment through qsub -V. The wrapper mounts the
    # launch directory itself; inputs and outputs elsewhere on the project's
    # storage need these binds. An explicit setting is left as it is.
    if [[ -z "${AIRRBORNE_CONTAINER_BINDS:-}" && -n "${PROJECT:-}" ]]; then
        binds=()
        for bind_path in "/g/data/${PROJECT}" "/scratch/${PROJECT}"; do
            [[ -d "${bind_path}" ]] && binds+=("${bind_path}")
        done
        if [[ ${#binds[@]} -gt 0 ]]; then
            AIRRBORNE_CONTAINER_BINDS="$(IFS=,; echo "${binds[*]}")"
            export AIRRBORNE_CONTAINER_BINDS
            echo "Container binds: ${AIRRBORNE_CONTAINER_BINDS}"
        fi
    fi

    # Create a directory for PBS job logs. pbs_submit.py writes these to a
    # cwd-relative logs/joblogs/, so create it under the caller's working
    # directory.
    echo "Creating job logs directory..."
    mkdir -p logs/joblogs

    echo "HPC configuration complete"
    echo ""
elif [[ -n "${SITE_PROFILE}" && "${SITE_PROFILE}" != none ]]; then
    PROFILE_ARGS=(--profile "${SITE_PROFILE}"
                  --workflow-profile "${WORKFLOW_ROOT}/profiles/workflow")
    echo "Using site profile: ${SITE_PROFILE}"
else
    # The local workflow profile selects Conda and supplies resource requests.
    PROFILE_ARGS=(--workflow-profile "${WORKFLOW_ROOT}/profiles/local")
    if [[ "${SITE_PROFILE}" == none ]]; then
        PROFILE_ARGS+=(--profile none)
    fi

    echo "Running on this host with the local workflow profile: ${WORKFLOW_ROOT}/profiles/local"
    echo "Containers run through workflow/scripts/container_exec.sh: apptainer, singularity or docker must be on PATH."
    # A cluster can supply the runtime as an environment module named after
    # it, as Gadi does. setup.sh loads the same module.
    runtime="${AIRRBORNE_CONTAINER_RUNTIME:-}"
    if [[ -z "${runtime}" ]]; then
        for candidate in apptainer singularity docker; do
            if command -v "${candidate}" >/dev/null 2>&1; then runtime="${candidate}"; break; fi
        done
    fi
    if [[ "${runtime}" != docker ]] && ! command -v "${runtime:-singularity}" >/dev/null 2>&1 \
            && command -v module >/dev/null 2>&1; then
        echo "Loading the ${runtime:-singularity} module..."
        module load "${runtime:-singularity}" || true
    fi
    # A scheduler job states its allocation, so local execution inside one is
    # sized to it. Explicit --cores and --resources after -- come later and win.
    # PBS gives the core count and the memory in bytes. Slurm gives the cores
    # per task or per node, and the memory per node or per core in megabytes.
    allocation_cores="${PBS_NCPUS:-${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-}}}"
    allocation_mem_mb=""
    if [[ -n "${PBS_VMEM:-}" ]]; then
        allocation_mem_mb=$(( PBS_VMEM / 1048576 ))
    elif [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        allocation_mem_mb="${SLURM_MEM_PER_NODE}"
    elif [[ -n "${SLURM_MEM_PER_CPU:-}" && -n "${allocation_cores}" ]]; then
        allocation_mem_mb=$(( SLURM_MEM_PER_CPU * allocation_cores ))
    fi
    if [[ -n "${allocation_cores}" ]]; then
        ALLOCATION_ARGS+=(--cores "${allocation_cores}")
    fi
    if [[ -n "${allocation_mem_mb}" ]]; then
        ALLOCATION_ARGS+=(--resources "mem_mb=${allocation_mem_mb}")
    fi
    if [[ ${#ALLOCATION_ARGS[@]} -gt 0 ]]; then
        echo "Limits from the scheduler allocation: ${ALLOCATION_ARGS[*]}"
    fi
    echo ""
fi
if [[ -n "${WORKFLOW_PROFILE}" ]]; then
    PROFILE_ARGS+=(--workflow-profile "${WORKFLOW_PROFILE}")
fi

# Unlock and exit before graph generation or execution.

if [[ "${UNLOCK}" == true ]]; then
    echo "Unlocking working directory..."
    snakemake "${PROFILE_ARGS[@]}" \
        --snakefile "${SNAKEFILE_PATH}" \
        --configfile "${CONFIG_FILE_PATH}" \
        --unlock
    echo "Working directory unlocked."
    exit 0
fi

# Write the rule graph before the workflow invocation.

if [[ "${GENERATE_DAG}" == true ]]; then
    if command -v dot >/dev/null 2>&1 && command -v tred >/dev/null 2>&1; then
        echo "Generating workflow DAG visualization..."
        mkdir -p "${WORKFLOW_ROOT}/docs/graphs"
        DAG_TMP="$(mktemp "${WORKFLOW_ROOT}/docs/graphs/workflow_dag.XXXXXX.svg")"

        DAG_CMD=(
            snakemake
            "${PROFILE_ARGS[@]}"
            --snakefile "${SNAKEFILE_PATH}"
            --configfile "${CONFIG_FILE_PATH}"
            --rulegraph dot
        )

        if [[ "${POST_TRUST4}" == true ]]; then
            DAG_CMD+=("${POST_TRUST4_TARGET_RULE}")
            echo "Generating DAG for post-TRUST4 rules only (target: ${POST_TRUST4_TARGET_RULE})"
        fi

        DAG_CMD+=(--config "analysis_only=${ANALYSIS_ONLY_ARG}")
        if [[ -n "${SCOPE}" ]]; then
            DAG_CMD+=("scope=${SCOPE}")
        fi

        if ! "${DAG_CMD[@]}" \
            | python "${WORKFLOW_ROOT}/docs/graphs/render_rulegraph.py" \
            > "${DAG_TMP}"; then
            rm -f "${DAG_TMP}"
            echo "ERROR: DAG generation failed."
            exit 1
        fi
        mv "${DAG_TMP}" "${WORKFLOW_ROOT}/docs/graphs/workflow_dag.svg"

        echo "DAG generated: ${WORKFLOW_ROOT}/docs/graphs/workflow_dag.svg"
        echo ""
    else
        echo "ERROR: Graph generation requires Graphviz dot and tred on PATH." >&2
        echo "Install it with: conda install -n snakemake_9.27.1 -c conda-forge graphviz" >&2
        exit 1
    fi
fi

# Build and execute the workflow command.

echo "Starting Snakemake workflow..."
echo ""

SNAKEMAKE_CMD=(
    snakemake
    "${PROFILE_ARGS[@]}"
    --snakefile "${SNAKEFILE_PATH}"
    --configfile "${CONFIG_FILE_PATH}"
    --rerun-incomplete
    --show-failed-logs
)

if [[ "${DRY_RUN}" == true ]]; then
    SNAKEMAKE_CMD+=(-n)
    echo "DRY-RUN MODE: Snakemake will list planned jobs without executing them."
    echo ""
fi

# Reuse retained assembly and depth products without resolving FASTQ inputs.
if [[ "${POST_TRUST4}" == true ]]; then
    SNAKEMAKE_CMD+=("${POST_TRUST4_TARGET_RULE}")
fi

# Place allocation limits and then forwarded arguments after the profile
# defaults, so command-line overrides win.
if [[ ${#ALLOCATION_ARGS[@]} -gt 0 ]]; then
    SNAKEMAKE_CMD+=("${ALLOCATION_ARGS[@]}")
fi
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    SNAKEMAKE_CMD+=("${EXTRA_ARGS[@]}")
fi

if [[ "${POST_TRUST4}" == true ]]; then
    echo "POST-TRUST4 MODE: Running downstream rules only (target: ${POST_TRUST4_TARGET_RULE})."
    echo ""
fi

# Append launcher settings last, so scope and assembly mode match these options.
SNAKEMAKE_CMD+=(--config "analysis_only=${ANALYSIS_ONLY_ARG}")
if [[ -n "${SCOPE}" ]]; then
    SNAKEMAKE_CMD+=("scope=${SCOPE}")
fi

echo "Command: ${SNAKEMAKE_CMD[@]}"
echo ""
"${SNAKEMAKE_CMD[@]}"

# Report successful completion.

if [[ "${DRY_RUN}" == true ]]; then
    echo ""
    echo "============================================================================"
    echo "Dry-run completed successfully."
    echo "No workflow jobs were executed."
    echo "============================================================================"
else
    echo ""
    echo "============================================================================"
    echo "Snakemake workflow execution finished."
    echo "============================================================================"
fi
