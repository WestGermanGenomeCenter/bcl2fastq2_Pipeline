#!/bin/bash
yes | conda env create -f environment.yml
#yes | conda create --name bcl2fastq2
conda activate bcl2fastq2
yes | conda install -c conda-forge mamba
yes | mamba install snakemake=7.8.5
yes | mamba install pandas
