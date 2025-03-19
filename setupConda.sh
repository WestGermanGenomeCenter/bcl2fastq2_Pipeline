#!/bin/bash
module load Miniconda/3_snakemake
source /software/conda/3/etc/profile.d/conda.sh
yes | conda create -c http://conda.repo.test.hhu.de/bioconda -c http://conda.repo.test.hhu.de/main -c http://conda.repo.test.hhu.de/conda-forge --override-channels --name bcl2fastq2 python=3.8.3
conda activate bcl2fastq2
# conda env create -f environment.yml # todo - this could simplify the setup 
yes | conda install -c http://conda.repo.test.hhu.de/conda-forge --override-channels mamba
conda config --add channels http://conda.repo.test.hhu.de/main
conda config --add channels http://conda.repo.test.hhu.de/conda-forge
conda config --add channels http://conda.repo.test.hhu.de/bioconda
yes | mamba install -c http://conda.repo.test.hhu.de/bioconda -c http://conda.repo.test.hhu.de/main -c http://conda.repo.test.hhu.de/conda-forge --override-channels snakemake=7.8.5
yes | mamba install -c http://conda.repo.test.hhu.de/bioconda -c http://conda.repo.test.hhu.de/main -c http://conda.repo.test.hhu.de/conda-forge --override-channels pandas
