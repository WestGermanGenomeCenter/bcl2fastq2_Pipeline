from snakemake.utils import min_version, validate
import os
import sys
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum
import pandas as pd
import re
# just to keep the rules more readable
dirname=os.path.dirname

min_version("6.4.0") # use-envmodules are not working in versions of 5.10.0 and below for clusters
p = os.path.abspath(".")
#Snakemake configs and setup
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
samplesheet = config["bcl2fastq"]["SampleSheet"]
outputfolder = config["bcl2fastq"]["OutputFolder"]
fastqre   = re.compile(r'\.fastq.gz$')

def validateBefore(outputfolder):
    success = validateSamplesheet(samplesheet)
    outputfolder = validateOutput(outputfolder)
    return True if success and outputfolder is not None else False

validRun = validateBefore(outputfolder)

# Check for projectNum in samplesheet name
projectNum=validateProjectNum(samplesheet)
if not projectNum:
    print("There was no project number to be found in the samplesheet and/or the name of the samplesheet")

# Set expected pipeline output
def getOutput():
    all = list()
    all.extend(expand("{out}/fastq_infiles_list.tx",out=outputfolder))
    all.extend(expand("{out}/untrimmed_fastq/{num}_sha256sums_fastqfiles.sha256",out=outputfolder,num=projectNum))# projectNum+"_sha256sums_fastqfiles.txt"
    all.extend(expand("{out}/{num}_software_environment.tsv",out=outputfolder,num=projectNum))
    all.extend(expand("{out}/{num}_config.yaml",out=outputfolder,num=projectNum))
    if config["interop_plots"]["interop_plots_active"]:
        # done_flag=outputfolder+"/interop_plots/interop_plots_done.flag"
        all.extend(expand("{out}/interop_plots/interop_plots_done.flag",out=outputfolder))
    if config["prevent_revcomp"]["prevent_revcomp_active"]:
        all.extend(expand("{out}/revcomp_prevented.flag",out=outputfolder))
# env_file=config["bcl2fastq"]["OutputFolder"]+"/"+projectNum+"_software_environment.tsv"
    if not validRun:
        all = list()
    return all

def getParentDir(wildcards):
    return dirname(str(config["bcl2fastq"]["SampleSheet"]))

localrules: all, create_fastq_list

rule all:
    input:
        getOutput() if projectNum else ""

# this rule has no parameter/flag/file that it needs to be executed before demultiplexing
# this works just because the ressources for this job are tiny, such that it will 
# all the time start before the demultiplexing starts
rule prevent_revcomp:
    input:
        sample_sheet=config["bcl2fastq"]["SampleSheet"]
    params:
        folder_of_runxml=getParentDir,
        script_to_fix="scripts/prevent_revcomp_nexts2.sh",
        out = config["bcl2fastq"]["OutputFolder"]
    output:
        config["bcl2fastq"]["OutputFolder"]+"/revcomp_prevented.flag" 
        #backup_copy=getParentDir+"backup_original_RunInfo.xml"
    shell:
        """
        mkdir -p {params.out}
        bash {params.script_to_fix} {params.folder_of_runxml}
        touch {output}
        """

rule bcl2fastq2:
    input:
        input = config["bcl2fastq"]["SampleSheet"]
    params:
        barcode_mismatches = config["bcl2fastq"]["barcode_mismatch"],
        threads = config["bcl2fastq"]["threads"],
        infolder = getParentDir,
        additionalOptions=[" "+config["bcl2fastq"]["options"],""][len(config["bcl2fastq"]["options"])>0],
        # add param for bcl convert location
        bcl_convert_path = config["bcl2fastq"]["bcl_convert_path"],
        out = config["bcl2fastq"]["OutputFolder"],
        outfastqs = config["bcl2fastq"]["OutputFolder"] + "/*.fastq.gz", # we want to exclude the undetermined maybe here, maybe after the multiqc
        output_dir = config["bcl2fastq"]["OutputFolder"],
        undetermined = config["bcl2fastq"]["OutputFolder"]  + "/untrimmed_fastq/" + "Undetermined*.fastq.gz",
        out_fastqs_dir = config["bcl2fastq"]["OutputFolder"] + "/untrimmed_fastq/",
        Logs_path= config["bcl2fastq"]["OutputFolder"] + "/Logs",
        bcl_log_path=config["bcl2fastq"]["OutputFolder"] + "/logs/bcl_convert_logs"

    log:
        config["bcl2fastq"]["OutputFolder"]+"/logs/e_bcl.log"
    output:
        bcl2fastq2Output = config["bcl2fastq"]["OutputFolder"]+"/Reports/Demultiplex_Stats.csv"
    threads: config["bcl2fastq"]["threads"]
    message:
        "Running bcl convert"
    envmodules:
        config["bcl2fastq"]["bcl2fastq2_version"]
        # change the cmd params for bcl convert
        # then make a new runPipeline.sh with a -n option
        # or make a config params to switch between bcl convert and bcl2fastq?????
        # barcode mismatches no longer part of command line arguments, now part of the samplesheet
    shell:
        """
        mkdir -p {params.out}
        chmod ago+rwx -R {params.out}
        {params.bcl_convert_path} --bcl-input-directory {params.infolder} --sample-sheet {input[0]} {params.additionalOptions} --no-lane-splitting true --output-directory {params.output_dir} --force --bcl-num-decompression-threads {params.threads} --bcl-num-conversion-threads {params.threads} --bcl-num-compression-threads {params.threads} --bcl-num-parallel-tiles {params.threads}
        cp {input[0]} {params.out}
        mkdir -p {params.output_dir} && mkdir -p {params.out_fastqs_dir} && mv {params.outfastqs} {params.out_fastqs_dir}
        mv {params.undetermined} {params.output_dir}
        chmod 775 -R {params.output_dir}
        mv {params.Logs_path} {params.bcl_log_path}
        
        """

rule create_fastq_list:
    input:
        rules.bcl2fastq2.output
    params:
        #out=config["bcl2fastq"]["OutputFolder"]
        out_fastqs_dir = config["bcl2fastq"]["OutputFolder"] + "/untrimmed_fastq/"
    output:
        fastq_list_file=config["bcl2fastq"]["OutputFolder"]+"/fastq_infiles_list.tx"
    shell:
        " cd {params.out_fastqs_dir} && shopt -s globstar; ls -d ** | grep '.fastq.gz' | grep -v 'Undetermined' >{output}"

rule create_checksums:
    input:
        fastq_list_file=config["bcl2fastq"]["OutputFolder"]+"/fastq_infiles_list.tx"
    params:
        out=config["bcl2fastq"]["OutputFolder"]+"/untrimmed_fastq/"
    output:
        checksum_file=config["bcl2fastq"]["OutputFolder"]+"/untrimmed_fastq/"+projectNum+"_sha256sums_fastqfiles.sha256"
    shell:
        "cd {params.out} && shopt -s globstar; ls -d ** | grep '.fastq.gz' | grep -v 'Undetermined' | xargs sha256sum >>{output}"


rule copy_software_env:
   input:
   	versions_file=config["software_versions_file"]
   params:
   	out=config["bcl2fastq"]["OutputFolder"] # not used
   output:
   	env_file=config["bcl2fastq"]["OutputFolder"]+"/"+projectNum+"_software_environment.tsv",
   	config_file=config["bcl2fastq"]["OutputFolder"]+"/"+projectNum+"_config.yaml"
   shell:
   	"cp {input.versions_file} {output.env_file} && cp config.yaml {output.config_file}"


rule interop_plots:
    input:
        bcl2fastq2Output = config["bcl2fastq"]["OutputFolder"]+"/Reports/Demultiplex_Stats.csv"
    output:
        done_flag=outputfolder+"/interop_plots/interop_plots_done.flag"
    params:
        interop_dir=config["interop_plots"]["interop_binaries_dir"],
        log_dir=config["bcl2fastq"]["OutputFolder"]+"/logs/interop/",
        output_dir=config["bcl2fastq"]["OutputFolder"]+"/interop_plots/",
        script_interop=config["interop_plots"]["interop_script_path"]+"/interop_plots.sh",

        raw_files_place=getParentDir
        # or just simply executing a bash script again
        
    conda:
        p+"/envs/interop.yaml"
    log:
        outputfolder+"/logs/interop/"+projectNum+"interop_plots.log"
    shell:
        " mkdir -p {params.log_dir} &&  mkdir -p {params.output_dir} && cd {params.output_dir} &&  bash {params.script_interop} {params.interop_dir} {params.raw_files_place} {params.output_dir}  && touch {output}"

onsuccess:
    print("Workflow finished without errors")

onerror:
    print("An error occurred during runtime")

onstart:
    print("Setting up and running pipeline")
