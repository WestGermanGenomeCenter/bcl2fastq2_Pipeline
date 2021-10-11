#!/usr/bin/env bash
module load Snakemake/6.4.0

[ ! -d "cluserLogs" ] && mkdir -p "clusterLogs"

while getopts hq: flag
do
    case "${flag}" in
        q) queue=${OPTARG};;
        h) hiseq=1;;
    esac
done

s_param=$1
if [ ! -z "$1" -a "$1" = "-q" ]; then
    s_param=""
    if [ ! -z "$3" ]; then
        s_param=$3
    fi
fi

cluster_param=cluster.json

if [ "$hiseq" = "1" ]; then
    cluster_param=cluster_hiseq.json
fi

if [ ! -z "$queue" -a "$queue" != " " ]; then
        echo "Pipeline running with queue: $queue";
        snakemake $s_param -s rules/bcl2fastq.smk --use-envmodules --cluster-config $cluster_param --cluster "qsub -A {cluster.account} -q ${queue} -l select={cluster.nodes}:ncpus={config["bcl2fastq"]["threads"]}:mem={cluster.mem} -l walltime={cluster.time} -r n -o {cluster.output} -e {cluster.error}" -j 99 --latency-wait 300 -p --keep-going --cluster-status "python snakemake-utils/statuscommand.py" --jobscript jobscript.sh
        snakemake $s_param -s rules/fastqProcessing.smk --use-envmodules --cluster-config $cluster_param --cluster "qsub -A {cluster.account} -q ${queue} -l select={cluster.nodes}:ncpus={config["bcl2fastq"]["fastQC_threads"]}:mem={cluster.mem} -l walltime={cluster.time} -r n -o {cluster.output} -e {cluster.error}" -j 5000 --latency-wait 1800 -p --keep-going --cluster-status "python snakemake-utils/statuscommand.py" --jobscript jobscript.sh
else
        echo "Pipeline running with available queue";
        snakemake $s_param -s rules/bcl2fastq.smk --use-envmodules --cluster-config $cluster_param --cluster "qsub -A {cluster.account} -l select={cluster.nodes}:ncpus={config["bcl2fastq"]["threads"]}:mem={cluster.mem} -l walltime={cluster.time} -r n -o {cluster.output} -e {cluster.error}" -j 99 --latency-wait 300 -p --keep-going --cluster-status "python snakemake-utils/statuscommand.py" --jobscript scripts/jobscript.sh
        snakemake $s_param -s rules/fastqProcessing.smk --use-envmodules --cluster-config $cluster_param --cluster "qsub -A {cluster.account} -l select={cluster.nodes}:ncpus={config["bcl2fastq"]["fastQC_threads"]}:mem={cluster.mem} -l walltime={cluster.time} -r n -o {cluster.output} -e {cluster.error}" -j 5000 --latency-wait 1800 -p --keep-going --cluster-status "python snakemake-utils/statuscommand.py" --jobscript scripts/jobscript.sh
fi

#chmod ago+rwx -R clusterLogs/ || true
#chmod ago+rwx -R .snakemake/ || true