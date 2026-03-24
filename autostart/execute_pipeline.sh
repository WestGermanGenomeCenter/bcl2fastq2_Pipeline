#PBS -l select=1:ncpus=2:mem=6gb                           
#PBS -l walltime=38:59:59
#PBS -A circs
#PBS -N demx_pipe_auto
#PBS -q BioJob

module load intel/xe2019
module load Java 

conda activate smk9
echo "`date +"%d.%m.%Y-%T"` starting the pipeline run" >> $LOGFILE


# PIPELINE="/gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/runPipeline.sh"
cd /gpfs/project/projects/bmfz_gtl/software/software_tests/update_2026/

bash runPipeline.sh


echo "`date +"%d.%m.%Y-%T"` ended the pipeline run" >> $LOGFILE

qstat -f $PBS_JOBID
