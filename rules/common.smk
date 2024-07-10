import os
import pandas as pd
p = os.path.abspath(".")
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum


configfile: "config.yaml"
outputfolder = config["bcl2fastq"]["OutputFolder"]

# Check for projectNum in samplesheet name
projectNum=validateProjectNum(samplesheet)
if not projectNum:
    print("There was no project number to be found in the samplesheet and/or the name of the samplesheet")

sample_names = list()
if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
    fastqs = list(samples_dataframe.iloc[:, 0].values)
    sample_names = [fastq[:-9] for fastq in fastqs]

def getFiles():
    files = list()
    for sample in sample_names:
        trySplit = sample.split(os.sep)[-1]
        if trySplit:
            sampleName = trySplit.split("_R")[0]
            files.append(sampleName)
    return files


# Get the fastq.gz samples
if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
    fastqs = list(samples_dataframe.iloc[:, 0].values)
    sample_names = [fastq[:-9] for fastq in fastqs]




def isSingleEnd() -> bool:
    """
    Returns wether the fastqs are single-end=True or paired-end=False
    """
    R1 = list()
    R2 = list()
    for sample in sample_names:
        if sample.split("_R")[1].startswith("1"):
            R1.append(sample)
            only_sample=sample.replace('_R1', '')
        else:
            R2.append(sample)
    if len(R1)!=len(R2):
        return True
    else:
        return False


def getSample_names_post_mapping():# maybe wildcards ne
	if isSingleEnd () == True:
		return sample_names
	else:
		pe_samplenames = list()
		for sample in sample_names:
			if sample.split("_R")[1].startswith("1"):
				only_sample=sample.replace('_R1','_pe')
				pe_samplenames.append(only_sample)
		return pe_samplenames



def getSamples(wildcards):
    
    if isSingleEnd() == True :
        #print("Single end input, thus aligning files:")
        for sample in sample_names:
            if str(wildcards) in sample:
                if config["umi_tools"]["umi_tools_active"]:
                    return expand("{out}/umi_extract/{file}.umis-extracted.fastq.gz",out=outputfolder,file=sample)
                elif config["cutadapt"]["cutadapt_active"]:# continue here - we have 3 cases: raw, cutadapt, umi_tools
                    return expand("{out}/trimmed/{file}_trimmed.fastq.gz",out=outputfolder,file=sample)
                else:
                    return expand("{out}/untrimmed_fastq/{file}.fastq.gz",out=outputfolder,file=sample)
                
    else:   
        #print("found paired end:")
        sampleName = str(wildcards).split("_R")[0]
        R1and2 = list()
        for sample in sample_names:
            if sample.split("_R")[1].startswith("1"):
                R1and2.append(sample)
            else:
                R1and2.append(sample)


        if config["umi_tools"]["umi_tools_active"]:
       	    return expand("{out}/umi_extract/{file}.umis-extracted.fastq.gz",out=outputfolder,file=R1and2)
        elif config["cutadapt"]["cutadapt_active"]:# continue here - we have 3 cases: raw, cutadapt, umi_tools
            returning_files=expand("{out}/trimmed/{file}_trimmed.fastq.gz",out=outputfolder,file=R1and2)
            #print ("pe, after cutadapt, returning files:")
            #print(returning_files)
            return expand("{out}/trimmed/{file}_trimmed.fastq.gz",out=outputfolder,file=R1and2)
        else:
       	    return expand("{out}/untrimmed_fastq/{file}.fastq.gz",out=outputfolder,file=R1and2)


def get_pe_pairingsheet():
    # get samplenamespostmapping, get r1 , get r2 and save all onto a file
    # use the old structure
    R1 = list()
    R2 = list() 
    pe_samplenames = list ()
    for sample in sample_names:
        if sample.split("_R")[1].startswith("1"):
            R1.append(sample)
            only_sample=sample.replace('_R1','_pe')
            pe_samplenames.append(only_sample)
        else:
            R2.append(sample)

    df = pd.DataFrame(list(zip(R1, R2, pe_samplenames)),columns=['read1','read2','pe_samplename'])

    df.to_csv(outputfolder+"/pe_samples.tsv", sep="\t", index=False)


only_sample=list() # if paired end, this spits out a  list of samplenames without the _R1, and only half of the fastq files 
for sample_full in sample_names:
    if sample_full.split("_R")[1].startswith("1"):
        only_sample_single=sample_full.replace('_R1', '')
        only_sample.append(only_sample_single)

def getFastqs(wildcards):
    samples = getSamples(wildcards)
    if len(samples)>1:
        return " ".join(samples)
    return samples
