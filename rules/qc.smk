## RSEQC
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
    
rule rseqc_gtf2bed:
    input:
        config["rseqc"]["gtf_file"],
    output:
        bed=outputfolder+"/qc/rseqc/annotation.bed",
        db=temp(outputfolder+"/qc/rseqc/annotation.db"),
    log:
        outputfolder+"/logs/rseqc_gtf2bed.log",
    conda:
        p+"/envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam=outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam",
        bed=outputfolder+"/qc/rseqc/annotation.bed",
    output:
        outputfolder+"/qc/rseqc/{file}.junctionanno.junction.bed"
    priority: 1
    log:
        outputfolder+"/logs/rseqc/rseqc_junction_annotation/{file}.log"
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        p+"/envs/rseqc.yaml"
    shell:
        """
        junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} >> {log} 2>&1
        [ ! -f {output} ] && touch {output} || :
        """


rule rseqc_junction_saturation:
    input:
        bam=outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam",
        bed=outputfolder+"/qc/rseqc/annotation.bed",
    output:
        outputfolder+"/qc/rseqc/{file}.junctionsat.junctionSaturation_plot.pdf"
    priority: 1
    log:
        outputfolder+"/logs/rseqc/rseqc_junction_saturation/{file}.log"
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        p+"/envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        ">> {log} 2>&1"


rule rseqc_stat:
    input:
        outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam"
    output:
        outputfolder+"/qc/rseqc/{file}.stats.txt"
    priority: 1
    log:
        outputfolder+"/logs/rseqc/rseqc_stat/{file}.log"
    conda:
        p+"/envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam=outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam",
        bed=outputfolder+"/qc/rseqc/annotation.bed",
    output:
        outputfolder+"/qc/rseqc/{file}.infer_experiment.txt"
    priority: 1
    log:
        outputfolder+"/logs/rseqc/rseqc_infer/{file}.log"
    conda:
        p+"/envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam=outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam",
        bed=outputfolder+"/qc/rseqc/annotation.bed",
    output:
        outputfolder+"/qc/rseqc/{file}.inner_distance_freq.inner_distance.txt"
    priority: 1
    log:
        outputfolder+"/logs/rseqc/rseqc_innerdis/{file}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        p+"/envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} >> {log} 2>&1"


rule rseqc_readdis:
    input:
        bam=outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam",
        bed=outputfolder+"/qc/rseqc/annotation.bed",
    output:
        outputfolder+"/qc/rseqc/{file}.readdistribution.txt"
    priority: 1
    log:
        outputfolder+"/logs/rseqc/rseqc_readdis/{file}.log"
    conda:
        p+"/envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam"
    output:
        outputfolder+"/qc/rseqc/{file}.readdup.DupRate_plot.pdf"
    priority: 1
    log:
        outputfolder+"/logs/rseqc/rseqc_readdup/{file}.log"
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        p+"/envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} >> {log} 2>&1"


rule rseqc_readgc:
    input:
        outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam"
    output:
        outputfolder+"/qc/rseqc/{file}.readgc.GC_plot.pdf"
    priority: 1
    log:
        outputfolder+"/logs/rseqc/rseqc_readgc/{file}.log"
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        p+"/envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} >> {log} 2>&1"


rule rseqc_done:
    input:
        expand(
            outputfolder+"/star/{file}/Aligned.sortedByCoord.out.bam",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/qc/rseqc/{file}.junctionanno.junction.bed",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/qc/rseqc/{file}.junctionsat.junctionSaturation_plot.pdf",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/qc/rseqc/{file}.infer_experiment.txt",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/qc/rseqc/{file}.stats.txt",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/qc/rseqc/{file}.inner_distance_freq.inner_distance.txt",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/qc/rseqc/{file}.readdistribution.txt",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/qc/rseqc/{file}.readdup.DupRate_plot.pdf",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/qc/rseqc/{file}.readgc.GC_plot.pdf",
            file=getFiles(),
        ),
        expand(
            outputfolder+"/logs/rseqc/rseqc_junction_annotation/{file}.log",
            file=getFiles(),
        ),
        expand(outputfolder+"/counts/all.tsv"),
    output:
        outputfolder+"/qc/rseqc/done.flag",
    shell:
        """
        touch {output}
        """

