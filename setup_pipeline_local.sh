
# setup the smk8 compatible workflow
# sometimes you need to 'module load condaforge' before the next command works
conda env create -n smk8 -f smk8_environment.yml
conda activate smk8



# make sure to create a snakemake profile and add it to the runPipeline_local.sh scrip with --profile your_profile_name