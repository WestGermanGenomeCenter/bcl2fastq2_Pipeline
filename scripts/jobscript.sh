#!/bin/bash
# properties = {properties}
module load Python/3.8.3
module load Snakemake/7.8.5
#module load Snakemake/6.4.0
module load Miniconda/3_snakemake
source /software/conda/3/etc/profile.d/conda.sh
conda activate bcl2fastq2
umask 000
{exec_job} 

