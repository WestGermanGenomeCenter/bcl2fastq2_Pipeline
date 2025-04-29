
# setup the smk8 compatible workflow
# sometimes you need to 'module load condaforge' before the next command works
conda env create -n smk8 -f smk8_environment.yml
conda activate smk8
mkdir -p ~/.config
mkdir -p ~/.config/snakemake
cp -r /gpfs/project/projects/bmfz_gtl/software/snakemake_profile/pbs ~/.config/snakemake/.

