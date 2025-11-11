snakemake -s rules/convert2fastq.smk --profile pbs
snakemake -s rules/fastqProcessing.smk --profile pbs
