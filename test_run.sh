#!/usr/bin/env bash
source /software/conda/3/etc/profile.d/conda.sh
conda activate bcl2fastq2
umask 000 # give rights to overwrite clusterLogs after each run
module load Snakemake/6.4.0

snakemake  -s rules/convert2fastq.smk --conda-frontend conda --resources --use-envmodules --use-conda --cluster-config cluster.json  --cluster "qsub -m n -A {cluster.account} -l select={cluster.nodes}:ncpus={threads}:mem={cluster.mem} -l walltime={cluster.time} -r n -o {cluster.output} -e {cluster.error}" -j 5000 --latency-wait 300 -p  --cluster-status "python snakemake-utils/statuscommand.py" --jobscript scripts/jobscript.sh --max-status-checks-per-second 1 --latency-wait 600 -n -p -r
snakemake  -s rules/fastqProcessing.smk --conda-frontend conda --resources --use-envmodules --use-conda --cluster-config cluster.json  --cluster "qsub -m n -A {cluster.account} -l select={cluster.nodes}:ncpus={threads}:mem={cluster.mem} -l walltime={cluster.time} -r n -o {cluster.output} -e {cluster.error}" -j 5000 --latency-wait 300 -p  --cluster-status "python snakemake-utils/statuscommand.py" --jobscript scripts/jobscript.sh --max-status-checks-per-second 1 --latency-wait 600 -n -p -r
