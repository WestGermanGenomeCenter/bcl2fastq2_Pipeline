# RSEQC
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

localrules: count_matrix

if not config["umi_tools"]["umi_tools_active"]:# since if umis are used, we need to use featurecounts
    rule count_matrix:
        input:
            expand(outputfolder+"/star/{file}/ReadsPerGene.out.tab", file=getFiles())
        output:
            outputfolder+"/counts/all.tsv"
        params:
            samples=getFiles(),
            strand=config["rseqc"]["strandedness"]
        conda:
            p+"/envs/pandas.yaml"
        script:
            "../scripts/count-matrix.py"

