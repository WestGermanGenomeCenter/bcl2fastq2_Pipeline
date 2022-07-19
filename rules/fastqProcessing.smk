# imports
from snakemake.utils import min_version, validate
import os
import sys
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum, adaptersToStringParams
import pandas as pd
p = os.path.abspath(".")

# Snakemake configs and setup
min_version("6.4.0")

configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
samplesheet = config["bcl2fastq"]["SampleSheet"]
outputfolder = config["bcl2fastq"]["OutputFolder"]

### include rules ###

include: "qcShort.smk"

# Get the fastq.gz samples
sample_names = list()
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


def getFiles():
    files = list()
    for sample in sample_names:
        trySplit = sample.split(os.sep)[-1]
        if trySplit:
            sampleName = trySplit.split("_R")[0]
            # make fastq/669-7_S7_R1_001_deduped.Aligned.sortedByCoord.out.bam into 669-1_S1/deduped.Aligned.sortedByCoord.out.bam
            files.append(sampleName)
    return files





# Set expected pipeline output
def getOutput():
    all = list()

    if config["gocryptfs"]["gocryptfs_active"]:
        all.extend(expand("{out}/encrypted/gocryptfs.conf",out=outputfolder))
    else:
        all.extend(expand("{out}/{name}.fastq.gz",out=outputfolder,name=sample_names))
        all.extend(expand("{out}/multiqc_report_{prj}.html",out=outputfolder,prj=projectNum))
        all.extend(expand("{out}/allFastQs_from_{prj}.checksums.sha256",out=outputfolder,prj=projectNum))

    if config["rseqc"]["rseqc_active"]:
        all.extend(expand("{out}/counts/all.tsv",out=outputfolder))

    if config["umi_tools"]["umi_tools_active"]:
        all.extend(expand("{out}/{name}.umis-extracted.fastq.gz",out=outputfolder,name=sample_names))
        all.extend(expand("{out}/star/{name}/deduped.Aligned.sortedByCoord.out.bam",out=outputfolder,name=sample_names))

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
        for sample in sample_names:
            relativePath = os.path.dirname(sample)+"/" if "/" in sample else ""
            fastQCs.extend(expand("{out}/fastqc_untrimmed/{path}untrimmed_{file}_fastqc.{ext}",out=outputfolder,path=relativePath,file=os.path.split(sample)[-1],ext=["html","zip"]))
    if config["rseqc"]["rseqc_active"]:
        for sample in sample_names:
            fastQCs.extend(expand("{out}/star/{file}/Aligned.sortedByCoord.out.bam",out=outputfolder,file=sample.split("/")[-1].split("_R")[0]))
            fastQCs.extend(expand("{out}/star/{file}/ReadsPerGene.out.tab",out=outputfolder,file=sample.split("/")[-1].split("_R")[0]))
    return fastQCs

def isSingleEnd() -> bool:
    """
    Returns wether the fastqs are single-end=True or paired-end=False
    """
    R1 = list()
    R2 = list()
    for sample in sample_names:
        if sample.split("_R")[1].startswith("1"):
            R1.append(sample)
        else:
            R2.append(sample)
    if len(R1)!=len(R2):
        return True
    return False



def getSamples(wildcards):
    if isSingleEnd():
        for sample in sample_names:
            if str(wildcards) in sample:
                if config["umi_tools"]["umi_tools_active"]:
                    return expand("{out}/{file}.umis-extracted.fastq.gz",out=outputfolder,file=sample)
                else:
                    return expand("{out}/{file}.fastq.gz",out=outputfolder,file=sample)
    sampleName = str(wildcards).split("_R")[0]
    R1and2 = list()
    for sample in sample_names:
        trySplit = sample.split(os.sep)[-1]
        if trySplit:
            if trySplit.startswith(sampleName):
                R1and2.append(sample)
    if config["umi_tools"]["umi_tools_active"]:
        return expand("{out}/{file}.umis-extracted.fastq.gz",out=outputfolder,file=R1and2)
    else:
        return expand("{out}/{file}.fastq.gz",out=outputfolder,file=R1and2)



def getPath(wildcards):
    """
    Returns the path in the samplename wildcard in case there is a path.
    Whether there is a path depends on the samplesheet+bcl2fastq2
    """
    if len(set(wildcards))>1:
        file = "{path}".format(path=wildcards["path"])
    else:
        file = "{name}".format(name=wildcards["name"])
    if "/" not in file:
        return ""
    file_split = file.split('/')[:-1]
    path = '/'.join(file_split)
    return path[1:] if path.startswith('/') else path

def getFastqs(wildcards):
    samples = getSamples(wildcards)
    if len(samples)>1:
        return " ".join(samples)
    return samples

# local rules
localrules: all, mergeSHA256, gocryptfs


rule all:
    input:
        getOutput() if projectNum else ""

# conditional rule if cutadapt is set to True
if config["cutadapt"]["cutadapt_active"]:
    rule FastQC_untrimmed:
        input:
            samples = outputfolder+"/{path}{file}.fastq.gz"
        params:
            path=getPath,
            html = outputfolder+"/fastqc_untrimmed/{path}{file}_fastqc.html",
            zip = outputfolder+"/fastqc_untrimmed/{path}{file}_fastqc.zip",
            out=outputfolder
        output:
            zip = expand("{out}/fastqc_untrimmed/{{path}}untrimmed_{{file}}_fastqc.zip", out=outputfolder),
            html = expand("{out}/fastqc_untrimmed/{{path}}untrimmed_{{file}}_fastqc.html", out=outputfolder)
        threads: config["others"]["fastQC_threads"]
        log:
            outputfolder+"/logs/fastqc/untrimmedFastQC_{path}{file}.log"
        conda:
            p+"/envs/fastqc.yaml"
        message:
            "Run untrimmed FastQC"
        shell:
            """
            mkdir -p {params.out}/fastqc_untrimmed/{params.path}
            chmod ago+rwx -R {params.out}/fastqc_untrimmed/ || :
            unset command_not_found_handle
            fastqc -q --threads {threads} {input} -o {params.out}/fastqc_untrimmed/{params.path} >> {log} 2>&1
            mv {params.zip} {output.zip}
            mv {params.html} {output.html}
            """


    rule cutadapt:
        input:
            outputfolder+"/{name}.fastq.gz" if not config["umi_tools"]["umi_tools_active"] else outputfolder+"/{name}.umis-extracted.fastq.gz" 
        params:
            adapters=adaptersToStringParams(config["cutadapt"]["adapters"],config["cutadapt"]["adapter_type"]),
            otherParams=config["cutadapt"]["other_params"],
            path=getPath,
            out=outputfolder
        output:
            expand("{out}/trimmed/{{name}}.fastq.gz", out=outputfolder)
        threads: config["cutadapt"]["cutadapt_threads"]
        log:
            outputfolder+"/logs/cutadapt/trimmed_{name}.log"
        message:
            "Run cutadapt"
        conda:
            p+"/envs/cutadapt.yaml"
        shell:
            """
            mkdir -p {params.out}/trimmed; mkdir -p {params.out}/trimmed/report; mkdir -p {params.out}/trimmed/report/{params.path}
            cutadapt {params.adapters} {params.otherParams} -o {output} {input} 2> {log}
            """

if config["rseqc"]["rseqc_active"]:
    rule align:
        input:# need to change this soon to be umi-tools aware
            fastqs=getSamples
        output:
            # see STAR manual for additional output files
            bams=expand("{out}/star/{{file}}/Aligned.sortedByCoord.out.bam",out=outputfolder), #this needs to be deduped 
            tabs=expand("{out}/star/{{file}}/ReadsPerGene.out.tab",out=outputfolder)
        wildcard_constraints:
            file="(.*(?:_).*)"
        log:
            outputfolder+"/logs/star/{file}.log"
        params:
            # optional parameters
            extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {} {}".format(
                  config["rseqc"]["gtf_file"], config["rseqc"]["other_params"]),
            outp = lambda w: os.path.dirname(outputfolder+"/star/"+w.file+"/")+"/",
            genDir = config["rseqc"]["genomeDir"],
            fastqs=lambda w: getFastqs(w.file)
        threads: config["rseqc"]["STAR_threads"]
        conda:
            p+"/envs/STAR.yaml"
        shell:
            "STAR {params.extra} --genomeDir {params.genDir} --runThreadN {threads} --readFilesIn {params.fastqs} --readFilesCommand zcat --outFileNamePrefix {params.outp} --outStd Log {log}"


rule FastQC:
    input:
        samples = outputfolder+"/{name}.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else rules.cutadapt.output
    params:
        path=getPath,
        out=outputfolder
    output:
        expand("{out}/fastqc/{{name}}_fastqc.{ext}", out=outputfolder, ext=["html","zip"])
    threads: config["others"]["fastQC_threads"]
    log:
        outputfolder+"/logs/fastqc/FastQC_{name}.log"
    conda:
        p+"/envs/fastqc.yaml"
    message:
        "Run FastQC"
    shell:
        """
        mkdir -p {params.out}/fastqc/{params.path}
        chmod ago+rwx -R {params.out}/fastqc/ || :
        fastqc -q --threads {threads} {input} -o {params.out}/fastqc/{params.path} >> {log} 2>&1
        """


# the umi extract is wanted before the fastqc happens, we dont need umi qc scores

if config["umi_tools"]["umi_tools_active"]:
    rule umi_extract:
        input:
            samples = outputfolder+"/{name}.fastq.gz" 

        params:
            path=getPath,
            umi_ptrn=config["umi_tools"]["pattern"],
            unzipped_files=outputfolder+"/{name}.fastq",
            extended_name=config["umi_tools"]["extended_name"],
            outputf=outputfolder+"/umi_extract",
            out=outputfolder
        log:
            outputfolder+"/logs/{name}.umi_extract.log"
        conda:
            p+"/envs/umi_tools.yaml"
        output:
            expand("{out}/{{name}}.umis-extracted.fastq.gz",out=outputfolder)
        shell:
            """
            mkdir -p {params.outputf}
            umi_tools extract --stdin={input} --bc-pattern={params.umi_ptrn} --log={log} --stdout={output}
            """ 
    
    rule umi_dedup:
        input:
            outputfolder+"/star/{name}/Aligned.sortedByCoord.out.bam"
        params:
            path=getPath,
            out=outputfolder,
            sorted_bam=outputfolder+"/star/{name}/sorted.Aligned.sortedByCoord.out.bam"
        log:
            outputfolder+"/logs/{name}.umi_dedup.log"
        conda:
            p+"/envs/umi_tools.yaml"
        threads:
            config["umi_tools"]["threads"]
        output:
            outputfolder+"/star/{name}/deduped.Aligned.sortedByCoord.out.bam"
        shell:
            """
            samtools sort --threads {threads} {input} -o {params.sorted_bam}
            samtools index {params.sorted_bam}
            umi_tools dedup -I {params.sorted_bam} --output-stats={log} -S {output}
            """


    rule featurecounts:
        input:
            expand("{out}/star/{name}/deduped.Aligned.sortedByCoord.out.bam",out=outputfolder,name=sample_names)
        params:
            gtf=config["rseqc"]["gtf_file"],
            stranded=config["umi_tools"]["stranded"],
            outdir=outputfolder+"/counts",

        threads:
            config["umi_tools"]["threads"]
        conda:
            p+"/envs/umi_tools.yaml"
        output:
            counts_file=outputfolder+"/counts/all.tsv",
        shell:
            """
            mkdir -p {params.outdir}
            featureCounts -a {params.gtf} -o {output.counts_file} {input} -T {threads} -g gene_id -s {params.stranded}
            """




rule SHA256:
    input:
        samples = outputfolder+"/{name}.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else rules.cutadapt.output
    output:
        expand("{out}/sha256/{{name}}.checksums.{ext}", out=outputfolder, ext=["sha256"])
    threads: config["others"]["SHA256_threads"]
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
    log:
        outputfolder+"/logs/MultiQC.log"
    message:
        "Run MultiQC"
    conda:
        p+"/envs/multiqc.yaml"
    shell:
        """
        chmod ago+rwx -R {params.output} || :
        multiqc --filename {output} --fn_as_s_name -q --ignore-samples Undetermined* -x Undetermined* -m "star" {params.output} -m 'rseqc' {params.output} -m 'fastqc' {params.output} -m 'cutadapt' {params.output} -o {params.output} >> {log} 2>&1
        chmod ago+rwx -R {params.output} || :
        """


rule gocryptfs:
    input:
        rules.multiqc.output,
        rules.mergeSHA256.output
    output:
        "{out}/encrypted/gocryptfs.conf".format(out=outputfolder)
    params:
        output=outputfolder
    log:
        outputfolder+"/logs/gocryptfs.log"
    message: "Run Gocryptfs"
    envmodules:
        config["gocryptfs"]["gocryptfs_version"]
    shell:
        """
        mkdir -p {params.output}/encrypted/
        gocryptfs -init {params.output}/encrypted/
        mkdir -p {params.output}/plain/
        echo 'Password:'
        gocryptfs {params.output}/encrypted/ {params.output}/plain/ >> {log} 2>&1
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

