#!/usr/bin/env bash
# get the outputfolder out of the config.yaml file, rename the dag and report
source ~/conda/etc/profile.d/conda.sh
conda activate smk8
out="$(grep OutputFolder config.yaml -A1 | tail -n 1 | awk '{print $1}' | sed 's/"//g' )"
mkdir -p $out
snakemake -s rules/convert2fastq.smk --forceall --rulegraph | dot -Tpdf > $out/convert2fastq_rulegraph.pdf
snakemake -s rules/convert2fastq.smk --profile pbs
snakemake -s rules/convert2fastq.smk --report $out/convert2fastq_report.html
snakemake -s rules/fastqProcessing.smk --forceall --rulegraph | dot -Tpdf > $out/fastqProcessing_rulegraph.pdf
snakemake -s rules/fastqProcessing.smk --profile pbs
snakemake -s rules/fastqProcessing.smk --report $out/fastqProcessing_report.html
