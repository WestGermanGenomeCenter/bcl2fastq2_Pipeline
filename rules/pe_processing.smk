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

#include: "qc.smk"
#include: "qc_pe.smk"
#include: "common.smk"
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



def isSingleEnd() -> bool:
    """
    Returns wether the fastqs are single-end=True or paired-end=False
    """
    R1 = list()
    R2 = list()
    for sample in sample_names:
        if sample.split("_R")[1].startswith("1"):
            R1.append(sample)
            #print("issingleend: attaching to R1:")
            #print(sample)
            only_sample=sample.replace('_R1', '')
                      # outputfolder+"/trimmed/{file}_trimmed.fastq.gz"
            #read1=expand("{out}/trimmed/{file}_{read}_trimmed.fastq.gz",out=outputfolder,
            #print("returning pe out read1:")
            #print(read1)
        else:
            R2.append(sample)
            #print("issingleend: attaching to R2:")
            #print(sample)
    if len(R1)!=len(R2):
        #print("R1 and R2 are not the same length, returning isSingleEnd=True")
        return True
#    print ("R1 and R2 are same length, returning isSingleEnd=False")
    else:
        return False

if isSingleEnd() == False:
    rule cutadapt_pe:
        input:
            in_1=outputfolder+"/untrimmed_fastq/{short}_R1.fastq.gz",
            in_2=outputfolder+"/untrimmed_fastq/{short}_R2.fastq.gz"
        params:
            #fastq_input=get_pe_trimming,
            #input_fastqs=getcutadapt_pe_string,
            adapters=adaptersToStringParams(config["cutadapt"]["adapters"],config["cutadapt"]["adapter_type"]),
            otherParams=config["cutadapt"]["other_params"],
            out=outputfolder,
            sha_sum_file=outputfolder+"/trimmed/sha256sums_fastqfiles_trimmed.sha256"
        output:
            fq1=outputfolder+"/trimmed/{short}_R1_trimmed.fastq.gz",
            fq2=outputfolder+"/trimmed/{short}_R2_trimmed.fastq.gz"


        threads: config["cutadapt"]["cutadapt_threads"]
        #log:
         #   outputfolder+"/logs/cutadapt/{short}_trimmed.log"
        message:
            "Run cutadapt"
        conda:
            p+"/envs/cutadapt.yaml"
        shell:
            """
            mkdir -p {params.out}/trimmed
            cutadapt {params.adapters} {params.otherParams} -o {output.fq1} -p {output.fq2} {input.in_1} {input.in_2} 
            sha256sum {output} >>{params.sha_sum_file}
            """



    rule align_pe: # aligning pe can only input trimmed  or umi_extracted data, no untrimmed data allowed!
        input:
            pe1=outputfolder + "/umi_extract/{short}_R1.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R1_trimmed.fastq.gz",
            pe2=outputfolder + "/umi_extract/{short}_R2.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R2_trimmed.fastq.gz"
        output:
            bams=outputfolder + "/star/{short}_pe_Aligned.sortedByCoord.out.bam", #this needs to be deduped 
            tabs=outputfolder + "/star/{short}_pe_ReadsPerGene.out.tab"
#        wildcard_constraints:
#            file="(.*(?:_).*)"
        log:
            outputfolder+"/logs/star/{short}.log"
        params:
            # optional parameters
            extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {} {}".format(
                config["rseqc"]["gtf_file"], config["rseqc"]["other_params"]),
            #outp = lambda w: os.path.dirname(outputfolder+"/star/"+w.file+"/"+w.file)+"/",
            prefix = outputfolder + "/star/{short}_pe_",
            genDir = config["rseqc"]["genomeDir"],
            qc_dir = outputfolder + "/qc",
            rseqc_dir = outputfolder+"/qc/rseqc",
            #fastqs=lambda w: getFastqs(w.file),

        threads: config["rseqc"]["STAR_threads"]
        conda:
            p+"/envs/STAR.yaml"
        shell:
            """
            mkdir -p {params.qc_dir}
            mkdir -p {params.rseqc_dir}
            STAR {params.extra} --genomeDir {params.genDir} --runThreadN {threads} --readFilesIn {input.pe1} {input.pe2} --readFilesCommand zcat --outFileNamePrefix {params.prefix} --outStd Log {log}
            chmod ago+rwx -R {output} >> {log} 2>&1
            """



    rule FastQC_pe:
        input:  # need to test this if this works, also the input fun of kraken rule 
        #             outputfolder+"/trimmed/{short}_R1_trimmed.fastq.gz",
            outputfolder+"/trimmed/{short}_R1_trimmed.fastq.gz",
            outputfolder+"/trimmed/{short}_R2_trimmed.fastq.gz"
           # outputfolder + "/umi_extract/{file}.umis-extracted.fastq.gz"  if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{file}_trimmed.fastq.gz
        params:
            #path=getPath,
            out=outputfolder,
            fastp_report=outputfolder+"/fastqc/{short}_after_filtering_fastp.html",
            fastp_json= outputfolder+"/fastqc/{short}_after_filtering_fastp.json"

        output:
            outputfolder+"/fastqc/{short}_R1_trimmed_fastqc.html",
            outputfolder+"/fastqc/{short}_R1_trimmed_fastqc.zip" ,
            outputfolder+"/fastqc/{short}_R2_trimmed_fastqc.html",
            outputfolder+"/fastqc/{short}_R2_trimmed_fastqc.zip" ,

        threads: config["others"]["fastQC_threads"]
        log:
            outputfolder+"/logs/fastqc/FastQC_{short}.log"
        conda:
            p+"/envs/fastqc.yaml"
        message:
            "Run FastQC"
        shell:
            """
            mkdir -p {params.out}/fastqc/ >> {log} 2>&1
            fastqc -q --threads {threads} {input} -o {params.out}/fastqc/ >> {log} 2>&1
            fastp -i {input} -h {params.fastp_report} -j {params.fastp_json} >> {log} 2>&1
            chmod ago+rwx -R {params.out}/fastqc/ >> {log} 2>&1
            """



    if config["umi_tools"]["umi_tools_active"]:# umi can be used without cutadapt, but thats not the default
        rule umi_extract_pe:
            input:
                in_1=outputfolder+"/untrimmed_fastq/{short}_R1.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else outputfolder+"/trimmed/{short}_R1_trimmed.fastq.gz",
                in_2=outputfolder+"/untrimmed_fastq/{short}_R2.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else outputfolder+"/trimmed/{short}_R2_trimmed.fastq.gz"                
            params:
                #path=getPath,
                umi_ptrn=lambda wc:config["umi_tools"]["pattern"], # need the lambda pseudo-fun for correct curly bracket in pattern recognition
                #unzipped_files=outputfolder+"/{file}.fastq",
                ptrn_2=lambda wc:config["umi_tools"]["pattern2"],
                #extended_name=config["umi_tools"]["extended_name"],
                outputf=outputfolder+"/umi_extract",
                logfolder=outputfolder+"/logs/umi_tools/",
                checksum_file=expand(outputfolder+"/umi_extract/checksums_umis_extracted_{prj}.sha256", prj=projectNum)
            log:
                outputfolder+"/logs/umi_tools/{short}.umi_extract.log"
            conda:
                p+"/envs/umi_tools.yaml"
            output:
                out1=outputfolder + "/umi_extract/{short}_R1.umis-extracted.fastq.gz",
                out2=outputfolder + "/umi_extract/{short}_R2.umis-extracted.fastq.gz"
            shell:
                """
                mkdir -p {params.outputf}
                mkdir -p {params.logfolder}
                umi_tools extract --stdin={input.in_1} --extract-method=regex --bc-pattern={params.umi_ptrn} --log={log} --stdout={output.out1} --read2-in={input.in_2} --read2-out={output.out2} --bc-pattern2={params.ptrn_2}
                chmod ago+rwx -R {params.outputf}
                sha256sum {output} >>{params.checksum_file}
                """ 

# umi_tools extract --extract-method=string
#--bc-pattern=[PATTERN] --bc-pattern2=[PATTERN]
#--read2-in=[FASTQIN] --read2-out=[FASTQOUT] -L extract.log [OPTIONS]

# add umi-tools if all works out
# test,test,test - even with se data











































































































