#!/bin/bash
#
# this script can be added to your crontab like this:
#crontab -l
#30 2 * * * cd /gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/autostart/ && bash autostart_cronjob_bclconvert.sh ./logfile_cronjob_pipeline.log 2>&1
#
#
WATCH_DIR="/gpfs/project/projects/bmfz_gtl/devices/illumina/nextseq2000"
PIPELINE="/gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/runPipeline.sh"
SNAKEMAKE_DIR="/gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026"
LOCKFILE=".pipeline_launched"
AUTO_SESSION_PREFIX="demx_pipe_auto"   # unique prefix — only auto-launched sessions use this, only one demuxing at a time allowed
execution_script=/gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/autostart/execute_pipeline.sh
# ── Check if any auto-launched pipeline session is still running ─────────────
# qstat does need to be added to the env that is very minimal: 
# remove the line below or edit it if your username is different /paths are wrong
PATH=/home/daric102/.local/bin:/home/daric102/conda/envs/smk9/bin:/home/daric102/conda/condabin:/usr/lib64/qt-3.3/bin:/usr/local/openssl/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/pbs/bin:/home/daric102/bin

if qstat | grep -q "$AUTO_SESSION_PREFIX"; then
    active=$(qstat | grep "$AUTO_SESSION_PREFIX" | awk '{print $1}')
    echo "[$(date)] Auto pipeline session still active: $active — exiting"
    exit 0
fi

# ── Scan for new runs ─────────────────────────────────────────────────────────
for folder in "$WATCH_DIR"/*/; do
    [[ -d "$folder" ]] || continue

    [[ -f "$folder/$LOCKFILE" ]] && continue

    shopt -s nullglob
    samplesheets=("$folder"/SampleSheet*.csv)
    [[ ${#samplesheets[@]} -eq 0 ]] && continue

    [[ -f "$folder/config.yaml" ]] || {
        echo "[$(date)] WARNING: $folder has SampleSheet but no config.yaml — skipping"
        continue
    }

    run_name=$(basename "$folder")
    session_name="${AUTO_SESSION_PREFIX}"

    echo "[$(date)] Found run: $run_name"

    cp "$folder/config.yaml" "$SNAKEMAKE_DIR/config.yaml"
    echo "[$(date)] Copied config.yaml -> $SNAKEMAKE_DIR/config.yaml"

    echo "launched=$(date)"       > "$folder/$LOCKFILE"
    echo "session=$session_name" >> "$folder/$LOCKFILE"
# split this into a typical pipeline execution script
    qsub -q default $execution_script 


    done

echo "[$(date)] finished scanning."