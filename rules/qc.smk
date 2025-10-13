## RSEQC
import os
import pandas as pd
p = os.path.abspath(".")
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum


configfile: "config.yaml"
outputfolder = config["demux"]["OutputFolder"]


projectNum=validateProjectNum(samplesheet)

#def getSample_names_post_mapping():
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
#
#sample_names = list()
#if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
#    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
#    fastqs = list(samples_dataframe.iloc[:, 0].values)
#    sample_names = [fastq[:-9] for fastq in fastqs]
#


sample_names = list()
if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
    fastqs = list(samples_dataframe.iloc[:, 0].values)
    sample_names = [fastq[:-9] for fastq in fastqs if "_I1_" not in fastq and "_I2_" not in fastq] # now excluding index reads in the sample names list


rule rseqc_gtf2bed:
    input:
        config["mapping"]["gtf_file"],
    output:
        bed=outputfolder+"/rseqc/annotation.bed",
        db=temp(outputfolder+"/rseqc/annotation.db")
    params:
    	qc_dir=outputfolder+"",
    	rseqc_dir=outputfolder+"/rseqc",
    log:
        outputfolder+"/logs/rseqc_gtf2bed.log",
    conda:
        p+"/envs/gffutils.yaml"
    message: "converting gtf to rseqc-bed..."

    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: attempt * 4    
    script:
        "../scripts/gtf2bed.py"



rule rseqc_main:
    input:
        bam=outputfolder + "/star/{file}_Aligned.sortedByCoord.out.bam" if not config["umi_tools"]["umi_tools_active"] else outputfolder+"/star/{file}_deduped.Aligned.sortedByCoord.out.bam",
        bed=outputfolder+"/rseqc/annotation.bed",
    output:
        junc_anno=outputfolder+"/rseqc/{file}.junctionanno.junction.bed",
        junc_sat=outputfolder+"/rseqc/{file}.junctionsat.junctionSaturation_plot.r",
        inner_dis=outputfolder+"/rseqc/{file}.inner_distance_freq.inner_distance_freq.txt"
    log:
        outputfolder+"/logs/rseqc/{file}.log"
    params:
        read_dis=outputfolder+"/rseqc/{file}.readdistribution.txt",
        infer_ex=outputfolder+"/rseqc/{file}.infer_experiment.txt",
        stats=outputfolder+"/rseqc/{file}.stats.txt",
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix_juncanno=outputfolder+"/rseqc/{file}.junctionanno",
        prefix_juncsat=outputfolder+"/rseqc/{file}.junctionsat",
        prefix_innerdis=outputfolder+"/rseqc/{file}.inner_distance_freq",
        nvc_outpre=outputfolder+"/rseqc/{file}_nucleotide_composition_bias",
        tin_prefix=outputfolder+"/rseqc/{file}_tin_.bed",
        genebody_pre=outputfolder+"/rseqc/{file}_genebody_coverage",
        fpkm_out=outputfolder+"/rseqc/{file}_rseqc_fpkm_count.tsv",
        rseqc_dir=outputfolder+"/rseqc",
        log_dir=outputfolder+"/logs/rseqc/"

    conda:
        p+"/envs/rseqc.yaml"
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 3,
        mem_gb=lambda wildcards, attempt: attempt * 12
    message: "RSeQC  Part 1/2 ..."

    shell:
        """
        mkdir -p {params.rseqc_dir} && mkdir -p {params.log_dir} >> {log} 2>&1
        junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix_juncanno} >> {log} 2>&1
        bam_stat.py -i {input.bam} > {params.stats} 2> {log}
        junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix_juncsat}  >> {log} 2>&1
        infer_experiment.py -r {input.bed} -i {input.bam} > {params.infer_ex} 2> {log}
        read_NVC.py -i {input.bam} -o {params.nvc_outpre} >> {log} 2>&1
        FPKM_count.py -r {input.bed} -i {input.bam} -o {params.fpkm_out} >> {log} 2>&1
        inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix_innerdis} >> {log} 2>&1
        """




rule rseqc_minor:
    input:
        bam=outputfolder + "/star/{file}_Aligned.sortedByCoord.out.bam" if not config["umi_tools"]["umi_tools_active"] else outputfolder+"/star/{file}_deduped.Aligned.sortedByCoord.out.bam",
        bed=outputfolder+"/rseqc/annotation.bed",
    output:
        dups=outputfolder+"/rseqc/{file}.readdup.DupRate_plot.r",
        gc=outputfolder+"/rseqc/{file}.readgc.GC_plot.r"

    log:
        outputfolder+"/logs/rseqc/{file}.log"
    params:
        read_dis=outputfolder+"/rseqc/{file}.readdistribution.txt",
        infer_ex=outputfolder+"/rseqc/{file}.infer_experiment.txt",
        stats=outputfolder+"/rseqc/{file}.stats.txt",
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        tin_prefix=outputfolder+"/rseqc/{file}_tin_.bed",
        genebody_pre=outputfolder+"/rseqc/{file}_genebody_coverage",
        fpkm_out=outputfolder+"/rseqc/{file}_rseqc_fpkm_count.tsv",
        prefix_readdup=outputfolder+"/rseqc/{file}.readdup",
        prefix_readgc=outputfolder+"/rseqc/{file}.readgc",
        prefix_bodycov=outputfolder+"/rseqc/{file}.genebody",
        rseqc_dir=outputfolder+"/rseqc",
        log_dir=outputfolder+"/logs/rseqc/"

    conda:
        p+"/envs/rseqc.yaml"
    message:
        "RSeqQC: Part 2..."
    resources:
        threads=lambda wildcards, attempt: attempt * 2,
        time_hrs=lambda wildcards, attempt: attempt * 3,
        mem_gb=lambda wildcards, attempt: attempt * 12
    shell:
        """
        mkdir -p {params.rseqc_dir} && mkdir -p {params.log_dir} >> {log} 2>&1
        read_distribution.py -r {input.bed} -i {input.bam} > {params.read_dis} 2> {log}
        read_duplication.py -i {input.bam} -o {params.prefix_readdup} >> {log} 2>&1
        read_GC.py -i {input.bam} -o {params.prefix_readgc} >> {log} 2>&1
        #geneBody_coverage.py -r {input.bed} -i {input.bam}  -o {params.prefix_bodycov}  >> {log} 2>&1
        """








rule rseqc_done:
    input:
        expand(
            outputfolder + "/star/{file}_Aligned.sortedByCoord.out.bam" if not config["umi_tools"]["umi_tools_active"] else outputfolder+"/star/{file}_deduped.Aligned.sortedByCoord.out.bam",
            file=getSample_names_post_mapping(),
        ),
        expand(
            outputfolder+"/rseqc/{file}.junctionanno.junction.bed",
            file=getSample_names_post_mapping(),
        ),
        expand(
            outputfolder+"/rseqc/{file}.junctionsat.junctionSaturation_plot.r",
            file=getSample_names_post_mapping(),
        ),
        expand(
            outputfolder+"/rseqc/{file}.inner_distance_freq.inner_distance_freq.txt",
            file=getSample_names_post_mapping(),
        ),
        expand(
            outputfolder+"/rseqc/{file}.readdup.DupRate_plot.r",
            file=getSample_names_post_mapping(),
        ),
        expand(
            outputfolder+"/rseqc/{file}.readgc.GC_plot.r",
            file=getSample_names_post_mapping(),
        ),

    output:
        outputfolder+"/rseqc/done.flag",
    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: attempt * 1
    shell:
        """
        touch {output}
        """


