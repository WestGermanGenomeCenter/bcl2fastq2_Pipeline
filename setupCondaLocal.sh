#!/bin/bash
yes | conda create --name bcl2fastq2
conda activate bcl2fastq2
yes | conda install -c conda-forge mamba=0.17.0
yes | mamba install snakemake-minimal=6.4.0
yes | mamba install pandas
