# imports
from snakemake.utils import min_version, validate
import os
import sys
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum, adaptersToStringParams
import os
import pandas as pd


# Snakemake configs and setup
min_version("6.4.0")

configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
samplesheet = config["bcl2fastq"]["SampleSheet"]
outputfolder = config["bcl2fastq"]["OutputFolder"]

# Get the fastq.gz samples
if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
    fastqs = list(samples_dataframe.iloc[:, 0].values)
    sample_names = [fastq[:-9] for fastq in fastqs]


# Validate input from config
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

    if config["gocryptfs"]["gocryptfs_active"]:
        all.extend(expand("{out}/encrypted/gocryptfs.conf",out=outputfolder))
    else:
        all.extend(expand("{out}/{name}.fastq.gz",out=outputfolder,name=sample_names))
        all.extend(expand("{out}/multiqc_report_{prj}.html",out=outputfolder,prj=projectNum))
        all.extend(expand("{out}/allFastQs_from_{prj}.checksums.sha256",out=outputfolder,prj=projectNum))    

    if not validRun or not os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
        all = list()
        print("It seems like the bcl2fastq2 run didnt finish properly or the fastq_infiles_list.tx doesnt exist")
    return all


def getSHA256(wildcards):
    sha256s = list()
    sha256s.extend(expand("{out}/sha256/{name}.checksums.sha256",out=outputfolder,name=sample_names))
    return sha256s


def getFastQCs(wildcards):
    fastQCs = list()
    fastQCs.extend(expand("{out}/fastqc/{name}_fastqc.{ext}",out=outputfolder,name=sample_names,ext=["html","zip"]))
    if config["cutadapt"]["cutadapt_active"]:
        fastQCs.extend(expand("{out}/fastqc_untrimmed/{name}_untrimmed_fastqc.{ext}",out=outputfolder,name=sample_names,ext=["html","zip"]))
    return fastQCs


def getPath(wildcards):
    """
    Returns the path in the samplename wildcard in case there is a path.
    Whether there is a path depends on the samplesheet+bcl2fastq2
    """
    file = "{name}".format(name=wildcards["name"])
    if "/" not in file:
        return ""
    file_split = file.split('/')[:-1]
    path = '/'.join(file_split)
    return path[1:] if path.startswith('/') else path

# local rules
localrules: all, mergeSHA256, gocryptfs



rule all:
    input:
        getOutput()


# conditional rule if cutadapt is set to True
if config["cutadapt"]["cutadapt_active"]:
    rule FastQC_untrimmed:
        input:
            samples = outputfolder+"/{name}.fastq.gz"
        params:
            fastQC=config["bcl2fastq"]["FastQC_version"],
            path=getPath,
            html = outputfolder+"/fastqc_untrimmed/{name}_fastqc.html",
            zip = outputfolder+"/fastqc_untrimmed/{name}_fastqc.zip",
            out=outputfolder
        output:
            zip = expand("{out}/fastqc_untrimmed/{{name}}_untrimmed_fastqc.zip", out=outputfolder),
            html = expand("{out}/fastqc_untrimmed/{{name}}_untrimmed_fastqc.html", out=outputfolder)
        threads: 1
        message:
            "Run untrimmed FastQC"
        shell:
            """
            mkdir -p {params.out}/fastqc_untrimmed/{params.path}
            chmod ago+rwx -R {params.out}
            {params.fastQC} -q --threads {threads} {input} -o {params.out}/fastqc_untrimmed/{params.path}
            mv {params.zip} {output.zip} 
            mv {params.html} {output.html}
            """
           

    rule cutadapt:
        input:
            samples = outputfolder+"/{name}.fastq.gz"
        params:
            adapters=adaptersToStringParams(config["cutadapt"]["adapters"],config["cutadapt"]["adapter_type"]),
            otherParams=config["cutadapt"]["other_params"],
            path=getPath,
            out=outputfolder
        output:
            expand("{out}/trimmed/{{name}}.fastq.gz", out=outputfolder)
        threads: 1
        message:
            "Run cutadapt"
        envmodules:
            config["cutadapt"]["cutadapt_version"]
        shell:
            """
            mkdir -p {params.out}/trimmed; mkdir -p {params.out}/trimmed/report; mkdir -p {params.out}/trimmed/report/{params.path}
            cutadapt {params.adapters} {params.otherParams} -o {output} {input} > {params.out}/trimmed/report/{wildcards.name}.txt
            """


rule FastQC:
    input:
        samples = outputfolder+"/{name}.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else rules.cutadapt.output
    params:
        fastQC=config["bcl2fastq"]["FastQC_version"],
        path=getPath,
        out=outputfolder
    output:
        expand("{out}/fastqc/{{name}}_fastqc.{ext}", out=outputfolder, ext=["html","zip"])
    threads: 1
    message:
        "Run FastQC"
    shell:
        """
        mkdir -p {params.out}/fastqc/{params.path}
        chmod ago+rwx -R {params.out}
        {params.fastQC} -q --threads {threads} {input} -o {params.out}/fastqc/{params.path}
        """


rule SHA256:
    input:
        samples = outputfolder+"/{name}.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else rules.cutadapt.output
    output:
        expand("{out}/sha256/{{name}}.checksums.{ext}", out=outputfolder, ext=["sha256"])
    threads: 1
    message:
        "Run SHA256"
    shell:
        """
        (sha256sum {input} | cut -d ' ' -f 1| tr -d '\n') > {output}
        echo -n -e " $(basename {input}) \n" >> {output}
        """


rule mergeSHA256:
    input:
        getSHA256
    output:
        "{out}/allFastQs_from_{prj}.checksums.sha256".format(out=outputfolder,prj=projectNum)
    threads: 1
    message:
        "Run mergeSHA256"
    shell:
        """
        cat {input} >> {output}
        """


rule multiqc:
    input:
        getFastQCs
    output:
        "{out}/multiqc_report_{prj}.html".format(out=outputfolder,prj=projectNum)
    params:
        output=outputfolder
    group: "multiqc"
    message:
        "Run MultiQC"
    envmodules:
        config["bcl2fastq"]["MultiQC_version"]
    shell:
        """
        chmod ago+rwx -R {params.output}
        multiqc --filename {output} -q -d --ignore-samples Undetermined* -x Undetermined* -m 'fastqc' {params.output} -m 'cutadapt' {params.output} -o {params.output}
        chmod ago+rwx -R {params.output}
        """


rule gocryptfs:
    input:
        rules.multiqc.output,
        rules.mergeSHA256.output
    output:
        "{out}/encrypted/gocryptfs.conf".format(out=outputfolder)
    params:
        output=outputfolder
    message: "Run Gocryptfs"
    envmodules:
        config["gocryptfs"]["gocryptfs_version"]
    shell:
        """
        mkdir -p {params.output}/encrypted/
        gocryptfs -init {params.output}/encrypted/
        mkdir -p {params.output}/plain/
        gocryptfs {params.output}/encrypted/ {params.output}/plain/
        rsync -av --progress {params.output}/* {params.output}/plain/ --exclude encrypted/ --exclude plain/ --remove-source-files
        fusermount -u {params.output}/plain/
        rm -r {params.output}/plain/
        shopt -s extglob
        rm -rf {params.output}/!(encrypted)/
        """


onsuccess:
    print("Workflow finished without errors")
    if config["cutadapt"]["cutadapt_active"]:
        print("Cutadapt ran with the following adapters and settings: %s" % adaptersToStringParams(config["cutadapt"]["adapters"],config["cutadapt"]["adapter_type"]))

onerror:
    print("An error occurred during runtime")

onstart:
    print("Setting up and running pipeline")
