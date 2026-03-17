#!/bin/bash
#
# this script can be added to your crontab like this:
#crontab -l
#30 2 * * * cd /gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/autostart/ && bash autostart_cronjob_bclconvert.sh ./logfile_cronjob_pipeline.log 2>&1
#
#
WATCH_DIR="/gpfs/project/projects/bmfz_gtl/devices/illumina/novaseqx"
PIPELINE="/gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/runPipeline.sh"
SNAKEMAKE_DIR="/gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026"
LOCKFILE=".pipeline_launched"
AUTO_SESSION_PREFIX="demx_pipe_auto"   # unique prefix — only auto-launched sessions use this, only one demuxing at a time allowed
execution_script=/gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/autostart/execute_pipeline.sh
# ── Check if any auto-launched pipeline session is still running ─────────────
# cron is usually starting scripts in a very minimal env. We here need to change that for qstat to work. So i added my complete PATH to make sure all things are working as in my typical env.
# remove the line below or edit it if your username is different /paths are wrong
PATH=/home/daric102/.local/bin:/home/daric102/conda/envs/smk9/bin:/home/daric102/conda/condabin:/usr/lib64/qt-3.3/bin:/usr/local/openssl/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/pbs/bin:/home/daric102/bin

if qstat | grep -q "$AUTO_SESSION_PREFIX"; then
    active=$(qstat | grep "$AUTO_SESSION_PREFIX" | awk '{print $1}')
    echo "[$(date)] Auto pipeline session still active: $active — exiting"
    exit 0
fi

# ── Scan for new runs ─────────────────────────────────────────────────────────
for folder in "$WATCH_DIR"/*/; do
# skip if not a folder
    [[ -d "$folder" ]] || continue

   # [[ -f "$folder/$LOCKFILE" ]] && continue

# skip if lockfile is inside the folder
     if [[ -f "$folder/$LOCKFILE" ]]; then
        echo "[$(date)] WARNING: $folder does include Lockfile: $LOCKFILE - skipping"
        continue
     fi



 # Skip if CopyComplete.txt doesn't exist
     if [[ ! -f "$folder/CopyComplete.txt" ]]; then
        echo "[$(date)] WARNING: $folder does not include CopyComplete.txt - skipping"
        continue
     fi

# then check for config.yaml
    [[ -f "$folder"/config.yaml ]] || {
        echo "[$(date)] WARNING: $folder has CopyComplete.txt but no config.yaml — skipping"
        continue
    }
# if all checks passed (not doing the continue once) then: gather infos, create lockfile, copy config.yaml, qsub
    run_name=$(basename "$folder")
    session_name="${AUTO_SESSION_PREFIX}"

    echo "[$(date)] Found run: $run_name"

    cp "$folder/config.yaml" "$SNAKEMAKE_DIR/config.yaml"
    echo "[$(date)] Copied config.yaml -> $SNAKEMAKE_DIR/config.yaml"

    echo "launched=$(date)"	  > "$folder/$LOCKFILE"
    echo "session=$session_name" >> "$folder/$LOCKFILE"
    echo "[$(date)] STARTING: no job already running, and $folder meets all criteria: no lockfile present, CopyComplete.txt and config.yaml present. - starting"
    wait 30
    qsub -q default $execution_script


    done

echo "[$(date)] finished scanning."


