#!/usr/bin/env bash
source /software/conda/3/etc/profile.d/conda.sh
conda activate bcl2fastq2
umask 000 # give rights to overwrite clusterLogs after each run
module load Snakemake/6.4.0
# check for clusterLogs folder if it doesn't exist create it
[ ! -d "cluserLogs" ] && mkdir -p "clusterLogs"

# globals
cluster_param=cluster.json
snakemake_params=""
local=false

# validate provided input
options=$(getopt -o q:hl -- "$@")
eval set -- "$options"

while true ; do
    case "$1" in
        -q)
            case "$2" in
                "") QUEUE=" " ; shift 2 ;;
                *) QUEUE="-q $2" ; shift 2 ;;
            esac ;;
        -h) cluster_param=cluster_hiseq.json ; shift ;;
        -l) local=true ; shift ;;
        --) shift; snakemake_params=$@ ; shift ${#@} ; break ;;
        *) snakemake_params=$@ ;  shift ${#@} ;;
    esac
done
    
if $local ; then
    # run pipeline locally in two steps
    set +x
    # overwrite hpc envs with local envs
    cp envs/local/* envs/
    snakemake $snakemake_params -s rules/convert2fastq.smk --conda-frontend mamba --use-conda -p 
    snakemake $snakemake_params -s rules/fastqProcessing.smk --resources --conda-frontend mamba --use-conda -p
    # overwrite local envs with hpc envs again
    cp envs/hpc/* envs/

else
    # print information for user
    if [ ! -z "$QUEUE" -a "$QUEUE" != " " ]; then
            echo "Pipeline running with queue: $QUEUE";
    else
            echo "Pipeline running with available queue";
    fi

    echo "Snakemake parameter: $snakemake_params"

    # run pipeline on hpc in two steps
    snakemake $snakemake_params -s rules/convert2fastq.smk  --conda-frontend mamba --use-conda --cluster-config $cluster_param --cluster "qsub -m n -A {cluster.account} ${QUEUE} -N {cluster.name} -l select={cluster.nodes}:ncpus={config["bcl2fastq"]["threads"]}:mem={cluster.mem} -l walltime={cluster.time} -r n -o {cluster.output} -e {cluster.error}" -j 99 --latency-wait 300 -p --cluster-status "python snakemake-utils/statuscommand.py" --jobscript scripts/jobscript.sh --max-status-checks-per-second 1 --latency-wait 600 --jobname {cluster.name}.{jobid}
    snakemake $snakemake_params -s rules/fastqProcessing.smk --conda-frontend mamba --resources --use-envmodules --use-conda --cluster-config $cluster_param --cluster "qsub -m n -A {cluster.account} ${QUEUE} -N {cluster.name} -l select={cluster.nodes}:ncpus={threads}:mem={cluster.mem} -l walltime={cluster.time} -r n -o {cluster.output} -e {cluster.error}" -j 5000 --latency-wait 300 -p  --cluster-status "python snakemake-utils/statuscommand.py" --jobscript scripts/jobscript.sh --max-status-checks-per-second 1 --latency-wait 600 --jobname {cluster.name}.{jobid}
fi

umask 022
