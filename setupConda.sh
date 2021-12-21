#!/bin/bash
module load Miniconda/3_snakemake
source /software/conda/3/etc/profile.d/conda.sh
yes | conda create -c http://conda.repo.test.hhu.de/bioconda -c http://conda.repo.test.hhu.de/main -c http://conda.repo.test.hhu.de/conda-forge --override-channels --name bcl2fastq2
conda activate bcl2fastq2
yes | conda install -c http://conda.repo.test.hhu.de/conda-forge --override-channels mamba=0.17.0
conda config --add channels http://conda.repo.test.hhu.de/main
conda config --add channels http://conda.repo.test.hhu.de/conda-forge
conda config --add channels http://conda.repo.test.hhu.de/bioconda
yes | mamba install -c http://conda.repo.test.hhu.de/bioconda -c http://conda.repo.test.hhu.de/main -c http://conda.repo.test.hhu.de/conda-forge --override-channels snakemake-minimal=6.4.0
yes | mamba install -c http://conda.repo.test.hhu.de/bioconda -c http://conda.repo.test.hhu.de/main -c http://conda.repo.test.hhu.de/conda-forge --override-channels pandas
