#!/usr/bin/env bash



# get the outputfolder out of the config.yaml file, rename the dag and report
source ~/conda/etc/profile.d/conda.sh
conda activate smk8
out="$(grep OutputFolder config.yaml -A1 | tail -n 1 | awk '{print $1}' | sed 's/"//g' )"
mkdir -p $out
time_exec="`date +"%Y_%m_%d_%I_%M_%p"`"


# print before execution what is on and what is off
echo "Pre-run: options enabled: "
echo "(config.yaml option set to True)"
echo "==================================="
grep True config.yaml | awk '{print $2":",$1}' | sed 's/_active//g;s/_on//g;s/://g'
echo ""
echo ""
echo "Pre-run: options disabled: "
echo "(config.yaml option set to False)"
echo "==================================="
grep False config.yaml | awk '{print $2":",$1}' | sed 's/_active//g;s/_on//g;s/://g'
echo ""
echo ""


if ls $out/*config.yaml 1> /dev/null 2>&1; then
    echo "Found files from previous execution, moving them to $out/logs/previous_executions"
    mkdir -p $out/logs/previous_executions
    mv -f $out/*_report*.html $out/logs/previous_executions/.
    mv -f $out/*_config.yaml $out/logs/previous_executions/.
    mv -f $out/*_rulegraph*.pdf $out/logs/previous_executions/.
    echo "Files from old execution moved."
else
    echo "No files from a previous execution found. starting..."
fi 





# actual execution
snakemake -s rules/convert2fastq.smk --forceall --rulegraph | dot -Tpdf > $out/convert2fastq_rulegraph.$time_exec.pdf
snakemake -s rules/convert2fastq.smk --profile pbs
snakemake -s rules/convert2fastq.smk --report $out/convert2fastq_report.$time_exec.html
snakemake -s rules/fastqProcessing.smk --forceall --rulegraph | dot -Tpdf > $out/fastqProcessing_rulegraph.$time_exec.pdf
snakemake -s rules/fastqProcessing.smk --profile pbs
snakemake -s rules/fastqProcessing.smk --report $out/fastqProcessing_report.$time_exec.html

# create checksums of all files created
echo "Completed the run."
echo "Creating Filelist of $out"
find $out -type f -exec ls  -alth --time-style=long-iso {} \; | sort > $out/filelist_project_$out.$time_exec.sha256
echo "Last Task: creating checksums:"
echo "Creating checksumfile $out/checksums_project_$out.$time_exec.sha256 ..."
find $out -type f -exec sha256sum {} \; | sort > $out/checksums_project_$out.$time_exec.sha256
echo "Done."