# imports
from snakemake.utils import min_version, validate
import os
import sys
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum, adaptersToStringParams
import pandas as pd
p = os.path.abspath(".")

# Snakemake configs and setup
min_version("8.28.0") # use-envmodules are not working in versions of 5.10.0 and below for clusters

configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
samplesheet = config["demux"]["SampleSheet"]
outputfolder = config["demux"]["OutputFolder"]

### include rules ###

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
    else:
        return False

if isSingleEnd() == False:


# add the second read, make mapping dependent on it
# need to make sortmerna put out pe aswell before this works


    def get_mapping_input_pe1(wildcards):
        if config["sortmerna"]["sortmerna_active"]:
            mapping_se_fastq1=outputfolder+"/sortmerna/{short}_R1_001_non-ribosomal_rna.fq.gz" # maybe this works
        elif config["umi_tools"]["umi_tools_active"]:
            mapping_se_fastq1=outputfolder + "/umi_extract/{short}_R1_001.fastq.gz"
        elif config ["cutadapt"]["cutadapt_active"]:
            mapping_se_fastq1=outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz"
        else:
            mapping_se_fastq1=outputfolder+"/untrimmed_fastq/{short}_R1_001.fastq.gz",
        return mapping_se_fastq1

    def get_mapping_input_pe2(wildcards):
        if config["sortmerna"]["sortmerna_active"]:
            mapping_se_fastq2=outputfolder+"/sortmerna/{short}_R2_001_non-ribosomal_rna.fq.gz" # maybe this works
        elif config["umi_tools"]["umi_tools_active"]:
            mapping_se_fastq2=outputfolder + "/umi_extract/{short}_R2_001.fastq.gz"
        elif config ["cutadapt"]["cutadapt_active"]:
            mapping_se_fastq2=outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"
        else:
            mapping_se_fastq2=outputfolder+"/untrimmed_fastq/{short}_R2_001.fastq.gz",
        return mapping_se_fastq2


    rule cutadapt_pe:
        input:
            in_1=outputfolder+"/untrimmed_fastq/{short}_R1_001.fastq.gz",
            in_2=outputfolder+"/untrimmed_fastq/{short}_R2_001.fastq.gz"
        params:
            adapters=adaptersToStringParams(config["cutadapt"]["adapters"],config["cutadapt"]["adapter_type"]),
            otherParams=config["cutadapt"]["other_params"],
            out=outputfolder,
            log_folder=outputfolder+"/logs/cutadapt/",
            sha_sum_file=outputfolder+"/trimmed/sha256sums_fastqfiles_trimmed.sha256"
        output:
            fq1=outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
            fq2=outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"

        resources:
            threads=lambda wildcards, attempt: attempt * 1,
            time_hrs=lambda wildcards, attempt: attempt * 1,
            mem_gb=lambda wildcards, attempt: attempt * 12
        log:
            outputfolder+"/logs/cutadapt/{short}_trimmed.log"
        message:
            "Run cutadapt"
        conda:
            p+"/envs/cutadapt.yaml"
        shell:
            """
            mkdir -p {params.out}/trimmed
            mkdir -p {params.log_folder}
            cutadapt {params.adapters} {params.otherParams} -o {output.fq1} -p {output.fq2} {input.in_1} {input.in_2} >> {log} 2>&1
            sha256sum {output} >>{params.sha_sum_file}
            """






    rule align_pe: # aligning pe can only input trimmed  or umi_extracted data, no untrimmed data allowed!
        input:
            pe1 = get_mapping_input_pe1,
            pe2 = get_mapping_input_pe2
            #pe1=outputfolder + "/umi_extract/{short}_R1.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R1_trimmed.fastq.gz",
            #pe2=outputfolder + "/umi_extract/{short}_R2.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R2_trimmed.fastq.gz"
        output:
            bams=outputfolder + "/star/{short}_pe_001_Aligned.sortedByCoord.out.bam", #this needs to be deduped 
            tabs=outputfolder + "/star/{short}_pe_001_ReadsPerGene.out.tab"
#        wildcard_constraints:
#            file="(.*(?:_).*)"
        log:
            outputfolder+"/logs/star/{short}.log"
        params:
            extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {} {}".format(
                config["mapping"]["gtf_file"], config["mapping"]["other_params"]),
            prefix = outputfolder + "/star/{short}_pe_001_",
            genDir = config["mapping"]["genomeDir"],
            mapping_dir = outputfolder+"/star",

        resources:
            threads=lambda wildcards, attempt: attempt * 12,
            time_hrs=lambda wildcards, attempt: attempt * 2,
            mem_gb=lambda wildcards, attempt: attempt * 40
        conda:
            p+"/envs/STAR.yaml"
        shell:
            """
            mkdir -p {params.mapping_dir}
            STAR {params.extra} --genomeDir {params.genDir} --runThreadN {resources.threads} --readFilesIn {input.pe1} {input.pe2} --readFilesCommand zcat --outFileNamePrefix {params.prefix} --outStd Log {log}
            chmod ago+rwx -R {params.mapping_dir} >> {log} 2>&1
            """

    
    rule umi_dedup_pe:
        input:
            outputfolder+"/star/{short}_pe_Aligned.sortedByCoord.out.bam"
        params:
            out=outputfolder,
            sorted_bam=outputfolder+"/star/{short}_pe_sorted.Aligned.sortedByCoord.out.bam",
            logfolder=outputfolder+"/logs/umi_tools/"
        log:
            outputfolder+"/logs/umi_tools/{short}_pe.umi_dedup.log"
        conda:
            p+"/envs/umi_tools.yaml"
        message: "umi_tools dedup: removing duplicate reads from mapping output..."

        resources:
            threads=lambda wildcards, attempt: attempt * 12,
            time_hrs=lambda wildcards, attempt: attempt * 4,
            mem_gb=lambda wildcards, attempt: attempt * 8

        output:
            outputfolder+"/star/{short}_deduped.Aligned.sortedByCoord.out.bam"
        shell:
            """
            mkdir -p {params.logfolder}
            samtools sort --threads {resources.threads} {input} -o {params.sorted_bam} 2>{log}
            samtools index {params.sorted_bam} 2>{log}
            umi_tools dedup -I {params.sorted_bam} --output-stats={log} -S {output} 2>{log}
            #chmod ago+rwx {params.out} 2>{log}
            """

    rule FastQC_pe:
        input:  # need to test this if this works, also the input fun of kraken rule 
            outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
            outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"
        params:
            out=outputfolder,
            fastp_report=outputfolder+"/fastqc/{short}_after_filtering_fastp.html",
            fastp_json= outputfolder+"/fastqc/{short}_after_filtering_fastp.json"

        output:
            outputfolder+"/fastqc/{short}_R1_001_trimmed_fastqc.html",
            outputfolder+"/fastqc/{short}_R1_001_trimmed_fastqc.zip" ,
            outputfolder+"/fastqc/{short}_R2_001_trimmed_fastqc.html",
            outputfolder+"/fastqc/{short}_R2_001_trimmed_fastqc.zip" ,

        resources:
            threads=lambda wildcards, attempt: attempt * 1,
            time_hrs=lambda wildcards, attempt: attempt * 1,
            mem_gb=lambda wildcards, attempt: attempt * 4
        log:
            outputfolder+"/logs/fastqc/FastQC_{short}.log"
        conda:
            p+"/envs/fastqc.yaml"
        message:
            "Run FastQC"
        shell:
            """
            mkdir -p {params.out}/fastqc/ >> {log} 2>&1
            fastqc -q --threads {resources.threads} {input} -o {params.out}/fastqc/ >> {log} 2>&1
            fastp -i {input} -h {params.fastp_report} -j {params.fastp_json} >> {log} 2>&1
            chmod ago+rwx -R {params.out}/fastqc/ >> {log} 2>&1
            """

#FastQC untrimmed missing


    if not config["umi_tools"]["umi_tools_active"]:# since if umis are used, we need to use featurecounts
        rule count_matrix_pe:
            input:
                expand(outputfolder+"/star/{file}/ReadsPerGene.out.tab", file=getSample_names_post_mapping())
            output:
                outputfolder+"/counts/all.tsv"
            resources:
                threads=lambda wildcards, attempt: attempt * 1,
                time_hrs=lambda wildcards, attempt: attempt * 1,
                mem_gb=lambda wildcards, attempt: attempt * 4
            params:
                samples=getSample_names_post_mapping(),
                strand="True"
            conda:
                p+"/envs/pandas.yaml"
            script:
                "../scripts/count-matrix.py"



    rule FastQC_untrimmed_pe:
        input:
            outputfolder+"/untrimmed_fastq/{short}_R1_001.fastq.gz",
            outputfolder+"/untrimmed_fastq/{short}_R2_001.fastq.gz"
        params:
        #    path=getPath,
            fastp_report = outputfolder+"/fastqc_untrimmed/{short}_untrimmed_fastp.html",
            fastp_json = outputfolder+"/fastqc_untrimmed/{short}_untrimmed_fastp.json",
            html = outputfolder+"/fastqc_untrimmed/{short}_fastqc.html",
            zip = temp(outputfolder+"/fastqc_untrimmed/{short}_fastqc.zip"),
            out = outputfolder
        output:
            outputfolder+"/fastqc_untrimmed/{short}_R1_001_fastqc.html",
            outputfolder+"/fastqc_untrimmed/{short}_R1_001_fastqc.zip" ,
            outputfolder+"/fastqc_untrimmed/{short}_R2_001_fastqc.html",
            outputfolder+"/fastqc_untrimmed/{short}_R2_001_fastqc.zip" ,
        resources:
            threads=lambda wildcards, attempt: attempt * 1,
            time_hrs=lambda wildcards, attempt: attempt * 1,
            mem_gb=lambda wildcards, attempt: attempt * 4
        log:
            outputfolder+"/logs/fastqc/untrimmedFastQC{short}.log"
        conda:
            p+"/envs/fastqc.yaml"
        message:
            "Run untrimmed FastQC"
        shell:
            """
            mkdir -p {params.out}/fastqc_untrimmed/
            chmod ago+rwx -R {params.out}/fastqc_untrimmed/ || :
            unset command_not_found_handle
            fastqc -q --threads {resources.threads} {input} -o {params.out}/fastqc_untrimmed/ >> {log} 2>&1
            fastp -i {input} -h {params.fastp_report} -j {params.fastp_json}>> {log} 2>&1
            mv {params.zip} {output.zip}
            mv {params.html} {output.html}
            """


    rule crypt4gh_pe:
        input:
            pe1=outputfolder + "/umi_extract/{short}_R1.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
            pe2=outputfolder + "/umi_extract/{short}_R2.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"
        params:
            priv_key=config["crypt4gh"]["crypt4gh_own_private_key"],
            pub_key=config["crypt4gh"]["crypt4gh_client_public_key_dir"],# here the collaborators an our many public keys, all that can decrypt later
            c4gh_folder=outputfolder+"/encrypted_fastqs",
            c4gh_logfolder=outputfolder+"/logs/crypt4gh",
            checksum_file=outputfolder+"/encrypted_fastqs/checksums_encrypted_fastqs.sha256",
        log:
            outputfolder+"/logs/crypt4gh/crypt4gh_encryption_{short}.log"
        message: "Running crypt4gh encryption on trimmed (and umis extracted if enabled) .fastq.gz files..."
        conda:
            p+"/envs/crypt4gh.yaml"
        output:
            c4gh_out1=outputfolder+"/encrypted_fastqs/{short}_R1_001_processed.fastq.gz.c4gh",
            c4gh_out2=outputfolder+"/encrypted_fastqs/{short}_R2_001_processed.fastq.gz.c4gh"
        resources:
            threads=lambda wildcards, attempt: attempt * 1,
            time_hrs=lambda wildcards, attempt: attempt * 2,
            mem_gb=lambda wildcards, attempt: attempt * 6 

        shell:
            """
            mkdir -p {params.c4gh_folder}
            mkdir -p {params.c4gh_logfolder}
            crypt4gh encrypt --sk {params.priv_key} `bash scripts/get_c4gh_string.sh {params.pub_key}` < {input.pe1} > {output.c4gh_out1} 2>{log}
            crypt4gh encrypt --sk {params.priv_key} `bash scripts/get_c4gh_string.sh {params.pub_key}` < {input.pe2} > {output.c4gh_out2} 2>{log}
            sha256sum {output} >>{params.checksum_file}
            """


    rule diamond_pe:
       input:
           in_1=outputfolder + "/umi_extract/{short}_R1.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
           in_2=outputfolder + "/umi_extract/{short}_R2.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"

           #outputfolder+"/umi_extract/{file}.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{file}_trimmed.fastq.gz"   # only possible if cutadapt and/or umi_tools are active
       output:
           log_file1=outputfolder+"/diamond/{short}_R1_001/diamond.log",
           log_file2=outputfolder+"/diamond/{short}_R2_001/diamond.log"

       params:
           out_folder=outputfolder+"/diamond/",
           out_sample_folder=outputfolder+"/diamond/{short}", # run diamond here
           diamond_ref_file=config["diamond"]["diamond_ref_file"],
           textfile1=outputfolder+"/diamond/{short}/{short}_r1_diamond_output.txt",
           textfile2=outputfolder+"/diamond/{short}/{short}_r2_diamond_output.txt",
           diamond_summary_file=outputfolder+"/diamond/{short}_diamond_output_summary.txt",
       log:
           outputfolder+"/diamond/{short}.diamond.log"
       conda:
           p+"/envs/diamond.yaml"
       message: "diamond: classifying Sequence data with peptide similarity..."
       resources:
           threads=lambda wildcards, attempt: attempt * 12,
           time_hrs=lambda wildcards, attempt: attempt * 2,
           mem_gb=lambda wildcards, attempt: attempt * 24

       shell:
           """
           mkdir -p {params.out_folder} 2>{log}
           mkdir -p {params.out_sample_folder}
           cd {params.out_sample_folder}
           diamond blastx --threads {resources.threads} --db {params.diamond_ref_file} --out {params.textfile1} --query {input.in_1} --outfmt 6 sscinames staxids sskingdoms skingdoms sphylums --taxon-k 1 --max-target-seqs 1 --log >> {log} 2>&1
           diamond blastx --threads {resources.threads} --db {params.diamond_ref_file} --out {params.textfile2} --query {input.in_2} --outfmt 6 sscinames staxids sskingdoms skingdoms sphylums --taxon-k 1 --max-target-seqs 1 --log >> {log} 2>&1
           cat {params.textfile1} |sort |uniq -c | sort -n | tac >>{params.diamond_summary_file}
           cat {params.textfile2} |sort |uniq -c | sort -n | tac >>{params.diamond_summary_file}

           """

    rule sortmerna_pe:
        input:
            sort_1=outputfolder + "/umi_extract/{short}_R1.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
            sort_2=outputfolder + "/umi_extract/{short}_R2.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"
            #outputfolder+"/umi_extract/{file}.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{file}_trimmed.fastq.gz"
        params:
            ref_string=lambda wc:config["sortmerna"]["sortmerna_reference_list"],
            fq_rrna_string=outputfolder+"/sortmerna/{short}_ribosomal_rna",
            fq_rrna_free_string=outputfolder+"/sortmerna/{short}_non-ribosomal_rna",
            log_folder=outputfolder+"/logs/sortmerna/",
            folder_sort=outputfolder+"/sortmerna/",
            workdir=outputfolder+"/sortmerna/{short}_sortmerna",
            fq_rrna1=outputfolder+"/sortmerna/{short}_ribosomal_rna_fwd.fq.gz",
            fq_rrna_free1=outputfolder+"/sortmerna/{short}_non-ribosomal_rna_fwd.fq.gz",
            fq_rrna2=outputfolder+"/sortmerna/{short}_ribosomal_rna_rev.fq.gz",
            fq_rrna_free2=outputfolder+"/sortmerna/{short}_non-ribosomal_rna_rev.fq.gz",

# what the names should be
            out_rrna1=outputfolder+"/sortmerna/{short}_R1_001_ribosomal_rna.fq.gz",
            out_rrna_free1=outputfolder+"/sortmerna/{short}_R1_001_non-ribosomal_rna.fq.gz",
            out_rrna2=outputfolder+"/sortmerna/{short}_R2_001_ribosomal_rna.fq.gz",
            out_rrna_free2=outputfolder+"/sortmerna/{short}_R2_001_non-ribosomal_rna.fq.gz",

        output:
        # mapping_se_fastq2=outputfolder+"/sortmerna/{short}_R2_001_non-ribosomal_rna.fq.gz"
            out_fq_rrna1=outputfolder+"/sortmerna/{short}_R1_001_ribosomal_rna.fq.gz",
            out_fq_rrna_free1=outputfolder+"/sortmerna/{short}_R1_001_non-ribosomal_rna.fq.gz",
            out_fq_rrna2=outputfolder+"/sortmerna/{short}_R2_001_ribosomal_rna.fq.gz",
            out_fq_rrna_free2=outputfolder+"/sortmerna/{short}_R2_001_non-ribosomal_rna.fq.gz",
        log:
            outputfolder+"/logs/sortmerna/rrna_removal_{short}.log"
        resources:
            threads=lambda wildcards, attempt: attempt * 24,
            time_hrs=lambda wildcards, attempt: attempt * 2,
            mem_gb=lambda wildcards, attempt: attempt * 24
        conda:
            p+"/envs/sortmerna.yaml"
        message:"sortmerna: removing rRNA reads from trimmed/ umi-moved .fastq files"
        shell:
            """
            rm -rf {params.workdir}
            mkdir -p {params.log_folder} 
            mkdir -p {params.folder_sort} 2>{log}
            sortmerna {params.ref_string} --reads {input.sort_1} --reads {input.sort_2} --threads {resources.threads} --workdir {params.workdir} --aligned {params.fq_rrna_string} --fastx --out2 --other {params.fq_rrna_free_string} --paired_in >> {log} 2>&1
            mv {params.fq_rrna1} {params.out_rrna1}
            mv {params.fq_rrna2} {params.out_rrna2}
            mv {params.fq_rrna_free1} {params.out_rrna_free1}
            mv {params.fq_rrna_free2} {params.out_rrna_free2}
            """


    if config["umi_tools"]["umi_tools_active"]:# umi can be used without cutadapt, but thats not the default
        rule umi_extract_pe:
            input:
                in_1=outputfolder+"/untrimmed_fastq/{short}_R1.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
                in_2=outputfolder+"/untrimmed_fastq/{short}_R2.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"                
            params:
                umi_ptrn=lambda wc:config["umi_tools"]["pattern"], # need the lambda pseudo-fun for correct curly bracket in pattern recognition
                ptrn_2=lambda wc:config["umi_tools"]["pattern2"],
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
                umi_tools extract --stdin={input.in_1} --extract-method=regex --bc-pattern={params.umi_ptrn} --log={log} --stdout={output.out1} --read2-in={input.in_2} --read2-out={output.out2} --bc-pattern2={params.ptrn_2} >> {log} 2>&1
                chmod ago+rwx -R {params.outputf}
                sha256sum {output} >>{params.checksum_file}
                """ 
