#!/bin/bash
# properties = {properties}
module load Python/3.12.0
module load Snakemake/8.28.0
source /software/conda/3/etc/profile.d/conda.sh
conda activate bcl2fastq2
umask 000
{exec_job} 

