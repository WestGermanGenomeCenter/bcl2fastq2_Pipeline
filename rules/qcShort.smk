# mapping
#import os
#import pandas as pd
#p = os.path.abspath(".")
#sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
#from helper import validateSamplesheet, validateOutput, validateProjectNum
#
#
#configfile: "config.yaml"
#outputfolder = config["demux"]["OutputFolder"]
#
#
#sample_names = list()
#if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
#    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
#    fastqs = list(samples_dataframe.iloc[:, 0].values)
#    sample_names = [fastq[:-9] for fastq in fastqs if "_I1_" not in fastq and "_I2_" not in fastq] # now excluding index reads in the sample names list
#
#def isSingleEnd() -> bool:
#    """
#    Returns wether the fastqs are single-end=True or paired-end=
#    """
#    R1 = list()
#    R2 = list()
#    for sample in sample_names:
#        if "_R" in sample:
#            if sample.split("_R")[1].startswith("1"):
#                R1.append(sample)
#            else:
#                R2.append(sample)
#        elif "_I" in sample:
#                print (" detected index reads, skipping these...")
#    if len(R1)!=len(R2):
#        return True
#    else:
#        return False
#
#
#def getSample_names_post_mapping():# maybe wildcards ne
#	if isSingleEnd () == True:
#		return sample_names
#	else:
#		pe_samplenames = list()
#		for sample in sample_names:
#			if sample.split("_R")[1].startswith("1"):
#				only_sample=sample.replace('_R1','_pe')
#				pe_samplenames.append(only_sample)
#		return pe_samplenames
#
##
## Check for projectNum in samplesheet name
#projectNum=validateProjectNum(samplesheet)
#if not projectNum:
#    print("There was no project number to be found in the samplesheet and/or the name of the samplesheet")
#
##sample_names = list()
##if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
##    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
##    fastqs = list(samples_dataframe.iloc[:, 0].values)
##    sample_names = [fastq[:-9] for fastq in fastqs]
#
#def getFiles():
#    files = list()
#    for sample in sample_names:
#        trySplit = sample.split(os.sep)[-1]
#        if trySplit:
#            sampleName = trySplit.split("_R")[0]
#            files.append(sampleName)
#    return files
#
#
#
#
#
#
#localrules: count_matrix
# now included in fastqProcessing, might need to salvage parts from here?
#if not config["umi_tools"]["umi_tools_active"]:# since if umis are used, we need to use featurecounts
#    rule count_matrix:
#        input:
#            expand(outputfolder+"/star/{file}/ReadsPerGene.out.tab", file=getSample_names_post_mapping())
#        output:
#            outputfolder+"/counts/all.tsv"
#        resources:
#            threads=lambda wildcards, attempt: attempt * 1,
#            time_hrs=lambda wildcards, attempt: attempt * 1,
#            mem_gb=lambda wildcards, attempt: attempt * 4
#        params:
#            samples=getFiles(),
#            strand="1"
#        conda:
#            p+"/envs/pandas.yaml"
#        script:
#            "../scripts/count-matrix.py"
#
