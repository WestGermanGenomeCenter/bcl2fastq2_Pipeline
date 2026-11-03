#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<EOF
Usage:
  $(basename "$0") \\
      --watch-dir DIR \\
      --pipeline FILE \\
      --snakemake-dir DIR \\
      --execution-script FILE \\
      --env-script FILE \\
      [--session-prefix PREFIX]

Required arguments:
  --watch-dir          Directory containing sequencing runs
  --pipeline           Path to runPipeline.sh (currently unused, kept for compatibility)
  --snakemake-dir      Snakemake working directory
  --execution-script   Script submitted via qsub
  --env-script         Script that sets PATH/modules/conda environment

Optional arguments:
  --session-prefix     Prefix used to identify running jobs
                       (default: demx_pipe_auto)
  -h, --help           Show this help

Example:
  $(basename "$0") \
      --watch-dir /gpfs/project/projects/bmfz_gtl/devices/illumina/nextseq2000 \
      --pipeline /gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/runPipeline.sh \
      --snakemake-dir /gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026 \
      --execution-script /gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/autostart/execute_pipeline.sh \
      --env-script /gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/autostart/pipeline_environment.sh
EOF
}

# ------------------------------------------------------------------
# Defaults
# ------------------------------------------------------------------

SESSION_PREFIX="demx_pipe_auto"
LOCKFILE=".pipeline_launched"

WATCH_DIR=""
PIPELINE=""
SNAKEMAKE_DIR=""
EXECUTION_SCRIPT=""
ENV_SCRIPT=""

# ------------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --watch-dir)
            WATCH_DIR="$2"
            shift 2
            ;;
        --pipeline)
            PIPELINE="$2"
            shift 2
            ;;
        --snakemake-dir)
            SNAKEMAKE_DIR="$2"
            shift 2
            ;;
        --execution-script)
            EXECUTION_SCRIPT="$2"
            shift 2
            ;;
        --env-script)
            ENV_SCRIPT="$2"
            shift 2
            ;;
        --session-prefix)
            SESSION_PREFIX="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo
            usage
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------------
# Validate required arguments
# ------------------------------------------------------------------

required=(
    WATCH_DIR
    PIPELINE
    SNAKEMAKE_DIR
    EXECUTION_SCRIPT
    ENV_SCRIPT
)

for var in "${required[@]}"; do
    if [[ -z "${!var}" ]]; then
        echo "Error: Missing required argument for $var"
        echo
        usage
        exit 1
    fi
done

[[ -d "$WATCH_DIR" ]] || { echo "WATCH_DIR does not exist."; exit 1; }
[[ -d "$SNAKEMAKE_DIR" ]] || { echo "SNAKEMAKE_DIR does not exist."; exit 1; }
[[ -f "$EXECUTION_SCRIPT" ]] || { echo "EXECUTION_SCRIPT not found."; exit 1; }
[[ -f "$ENV_SCRIPT" ]] || { echo "ENV_SCRIPT not found."; exit 1; }

# ------------------------------------------------------------------
# Load environment
# ------------------------------------------------------------------

source "$ENV_SCRIPT"

# ------------------------------------------------------------------
# Prevent multiple auto-started jobs
# ------------------------------------------------------------------

if qstat | grep -q "$SESSION_PREFIX"; then
    active=$(qstat | grep "$SESSION_PREFIX" | awk '{print $1}')
    echo "[$(date)] Auto pipeline session still active: $active - exiting"
    exit 0
fi

# ------------------------------------------------------------------
# Scan for new runs
# ------------------------------------------------------------------

for folder in "$WATCH_DIR"/*/; do

    [[ -d "$folder" ]] || continue

    if [[ -f "$folder/$LOCKFILE" ]]; then
        echo "[$(date)] WARNING: $folder contains $LOCKFILE - skipping"
        continue
    fi

    if [[ ! -f "$folder/CopyComplete.txt" ]]; then
        echo "[$(date)] WARNING: $folder missing CopyComplete.txt - skipping"
        continue
    fi

    if [[ ! -f "$folder/config.yaml" ]]; then
        echo "[$(date)] WARNING: $folder has CopyComplete.txt but no config.yaml - skipping"
        continue
    fi

    run_name=$(basename "$folder")

    echo "[$(date)] Found run: $run_name"

    cp "$folder/config.yaml" "$SNAKEMAKE_DIR/config.yaml"

    echo "[$(date)] Copied config.yaml -> $SNAKEMAKE_DIR/config.yaml"

    {
        echo "launched=$(date)"
        echo "session=$SESSION_PREFIX"
    } > "$folder/$LOCKFILE"

    echo "[$(date)] STARTING pipeline."

    sleep 30

    qsub -q default "$EXECUTION_SCRIPT"

done

echo "[$(date)] Finished scanning."