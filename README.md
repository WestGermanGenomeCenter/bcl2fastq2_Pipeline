# demultiplexing Pipeline for Illumina Machines with extensive QC


This Pipeline uses Snakemake 8 on a PBS Pro HPC Cluster to Demultiplex Illumina Sequencing Runs, either with bcl2fastq2 or the newer bclconvert.


##  deployment
if you are HILBERT HPC user:

 `bash setup_pipeline.sh`

if not: 
 - create your snakemake profile if needed
 - create a conda env with `conda env create -f smk8_environment.yml`
 - `conda activate smk8 ` <- this env should then include snakemake 8.x.x


then:
 - edit the config.yaml file. This is where all the settings are for this Pipeline. Only use tools that you need!

then:

`bash runPipeline.sh`


# Examples
The Pipeline is split into two main Parts:
- the demultiplexing (rules/convert2fastq.smk)
- all other things (fastqProcessing.smk + all other .smk files in rules/)

For example DAGs and Reports, check the examples/ folder inside this repo

# General remarks: The Pipeline 
- is setup to always use the config.yaml in the root folder
- expects bcl2fastq2/bclconvert to be somewhere as executable (needs to be set in config.yaml)
- always does some basic steps like sha256sums / fastqc + multiqc reports
- can deal with .fastq files as input (if demultiplexing already happened), for this enable the skip_demux part in the config.yaml
- can deal with paired and and single read sequencing setups. Diamond does not work with pe setups, and the naming scheme is strict (see config.yaml)
- if qc steps are enabled that require trimming of adapters, but trimming is disabled the pipeline will trimm regardless. This is to ensure the proper input quality of data for each step.
- The mayor focus here is demultiplexing of MiSeq/Nextseq2000/NovaseqX and do RNASeq analysis until transcript count tables are created
- can deal with UMIs, even in paired end setups (check config.yaml)
- collects qc of raw sequencing output in a multiqc file, and output of all analysis steps in a multiqc_report_complete file.
- general deletion of data does not happen. You need to cleanup on your own.
- expects a RunFolder as input that includes a SampleSheet. The SampleSheet needs to have a number (ProjectID) in the Filename. The RunFolder is where the demultiplexing then starts. The Samplesheet cannot be outside the root folder of the RunFolder created by Illuminas machines. See config.yaml for an example.
- SampleSheet and OutputFolder are the only 2 required inputs in the config.yaml, all other steps can be disabled.
- currently uses RSeQC as post-mapping QC for RNASeq runs, but only the for us useful modules of RSeQC.
- includes many QC steps to identify contamination, from most reliable to least: Kraken , BioBloom, BLAST, Diamond. Please be aware that either results are dependent on Quality of the input and the Database used.
- automatically scales up the ressources for each job after a failed attempt. after 3 attempts it gives up - you can configure this in your snakemake profile.
- many results of analyses are shown in the multiqc report complete file, but a lot of tools also put out more data then multiqc catches. Check the output folder after each run for more!