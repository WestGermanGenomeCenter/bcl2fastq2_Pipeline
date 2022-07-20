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
    if not validRun:
        all = list()
    return all

def getParentDir(wildcards):
    return dirname(str(config["bcl2fastq"]["SampleSheet"]))

localrules: all, create_fastq_list

rule all:
    input:
        getOutput() if projectNum else ""

rule bcl2fastq2:
    input:
        input = config["bcl2fastq"]["SampleSheet"]
    params:
        barcode_mismatches = config["bcl2fastq"]["barcode_mismatch"],
        threads = config["bcl2fastq"]["threads"],
        infolder = getParentDir,
        additionalOptions=[" "+config["bcl2fastq"]["options"],""][len(config["bcl2fastq"]["options"])>0],
        # add param for bcl convert location
        bcl_convert_path = config["bcl_convert_path"]
        out = config["bcl2fastq"]["OutputFolder"]
    log:
        config["bcl2fastq"]["OutputFolder"]+"/logs/e_bcl.log"
    output:
        bcl2fastq2Output = config["bcl2fastq"]["OutputFolder"]+"/Stats/Stats.json"
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
        {params.bcl_convert_path} --bcl-input-directory {params.infolder} --sample-sheet {input[0]} {params.additionalOptions} --no-lane-splitting --output-directory {params.out} --bcl-num-decompression-threads {params.threads} --bcl-conversion-threads {params.threads} --bcl-num-compression-threads {params.threads} --bcl-num-parallel-tiles {params.threads}
       # bcl2fastq -R {params.infolder} --sample-sheet {input[0]} {params.additionalOptions}--no-lane-splitting --barcode-mismatches {params.barcode_mismatches} -o {params.out} --interop-dir {params.out}  -r {params.threads} -p {params.threads} 2> {log[0]}
        cp {input[0]} {params.out}
        """

rule create_fastq_list:
    input:
        rules.bcl2fastq2.output
    params:
        out=config["bcl2fastq"]["OutputFolder"]
    output:
        fastq_list_file=config["bcl2fastq"]["OutputFolder"]+"/fastq_infiles_list.tx"
    shell:
        " cd {params.out} && shopt -s globstar; ls -d ** | grep '.fastq.gz' | grep -v 'Undetermined' >{output}"



onsuccess:
    print("Workflow finished without errors")

onerror:
    print("An error occurred during runtime")

onstart:
    print("Setting up and running pipeline")
