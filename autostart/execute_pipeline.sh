#PBS -l select=1:ncpus=2:mem=6gb                           
#PBS -l walltime=38:59:59
#PBS -A circs
#PBS -N demx_pipe_auto
#PBS -q BioJob

module load intel/xe2019
module load Java 

conda activate smk9


# PIPELINE="/gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/runPipeline.sh"
cd /gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/
wait 60;
bash runPipeline.sh


echo "`date +"%d.%m.%Y-%T"`" >> $LOGFILE

qstat -f $PBS_JOBID
