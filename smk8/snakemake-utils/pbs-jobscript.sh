#!/bin/sh
# properties = {properties}

export SINGULARITY_CACHEDIR=/tmp/.singularity-$(id -u)

echo "Will execute the following jobscript: "
cat $0

# Will be inserted by pbs-submit.py
# <modules>

{exec_job}
