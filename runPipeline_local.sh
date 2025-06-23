#!/usr/bin/env bash
# get the outputfolder out of the config.yaml file, rename the dag and report
source ~/conda/etc/profile.d/conda.sh
conda activate smk8
out="$(grep OutputFolder config.yaml -A1 | tail -n 1 | awk '{print $1}' | sed 's/"//g' )"
mkdir -p $out
time_exec="`date +"%Y_%m_%d_%I_%M_%p"`"
nice snakemake -s rules/convert2fastq.smk --forceall --rulegraph | dot -Tpdf > $out/convert2fastq_rulegraph.$time_exec.pdf
nice snakemake -s rules/convert2fastq.smk --cores all
nice snakemake -s rules/convert2fastq.smk --report $out/convert2fastq_report.$time_exec.html
nice snakemake -s rules/fastqProcessing.smk --forceall --rulegraph | dot -Tpdf > $out/fastqProcessing_rulegraph.$time_exec.pdf
nice snakemake -s rules/fastqProcessing.smk --cores all
nice snakemake -s rules/fastqProcessing.smk --report $out/fastqProcessing_report.$time_exec.html
