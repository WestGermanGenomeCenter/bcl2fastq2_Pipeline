#!/usr/bin/env bash
# get the outputfolder out of the config.yaml file, rename the dag and report
# check if conda exec is in ~, then load this to be able to load the smk8 env before executing the snake
if [ -f ~/conda/etc/profile.d/conda.sh ]; then
    source ~/conda/etc/profile.d/conda.sh
else
    echo "File not found: ~/conda/etc/profile.d/conda.sh"
    echo "trying to activate the conda env without loading conda"
fi

conda activate smk8
out="$(grep OutputFolder config.yaml -A1 | tail -n 1 | awk '{print $1}' | sed 's/"//g' )"
mkdir -p $out
time_exec="`date +"%Y_%m_%d_%I_%M_%p"`"
nice snakemake -s rules/convert2fastq.smk --forceall --rulegraph | dot -Tpdf > $out/convert2fastq_rulegraph.$time_exec.pdf
nice snakemake -s rules/convert2fastq.smk --cores all --sdm conda --jobs 7
nice snakemake -s rules/convert2fastq.smk --report $out/convert2fastq_report.$time_exec.html
nice snakemake -s rules/fastqProcessing.smk --forceall --rulegraph | dot -Tpdf > $out/fastqProcessing_rulegraph.$time_exec.pdf
nice snakemake -s rules/fastqProcessing.smk --cores all --sdm conda --jobs 4
nice snakemake -s rules/fastqProcessing.smk --report $out/fastqProcessing_report.$time_exec.html
