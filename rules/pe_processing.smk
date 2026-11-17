# imports
from snakemake.utils import min_version
import os
import sys
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum, adaptersToStringParams
import pandas as pd
p = os.path.abspath(".")

# Snakemake configs and setup
min_version("8.28.0") # use-envmodules are not working in versions of 5.10.0 and below for clusters

configfile: "config.yaml"
samplesheet = config["demux"]["SampleSheet"]
outputfolder = config["demux"]["OutputFolder"]

# Validate input from config
def validateBefore(outputfolder):
    success = validateSamplesheet(samplesheet)
    outputfolder = validateOutput(outputfolder)
    return True if success and outputfolder is not None else False

validRun = validateBefore(outputfolder)

def get_pe_pairingsheet():
    parts = [(s, splitReadTag(s)) for s in sample_names]
    r1 = {t[0]: s for s, t in parts if t[2] == 1}
    r2 = {t[0]: s for s, t in parts if t[2] == 2}
    rows = [(r1[b], r2[b], b + "_pe_001") for b in sorted(r1) if b in r2]
    df = pd.DataFrame(rows, columns=['read1', 'read2', 'pe_samplename'])
    df.to_csv(outputfolder+"/pe_samples.tsv", sep="\t", index=False)

if isSingleEnd() == False:

    def _pe_mapping_input(wildcards, read):
        short = wildcards.short
        if config["sortmerna"]["sortmerna_active"]:
            return f"{outputfolder}/sortmerna/{short}_R{read}_001_non-ribosomal_rna.fq.gz"
        if config["umi_tools"]["umi_tools_active"]:
            return f"{outputfolder}/umi_extract/{short}_R{read}_001.umis-extracted.fastq.gz"
        if config["cutadapt"]["cutadapt_active"]:
            return f"{outputfolder}/trimmed/{short}_R{read}_001_trimmed.fastq.gz"
        return f"{outputfolder}/untrimmed_fastq/{short}_R{read}_001.fastq.gz"

    def get_mapping_input_pe1(wildcards):
        return _pe_mapping_input(wildcards, 1)

    def get_mapping_input_pe2(wildcards):
        return _pe_mapping_input(wildcards, 2)


    rule cutadapt_pe:
        input:
            in_1=outputfolder+"/untrimmed_fastq/{short}_R1_001.fastq.gz",
            in_2=outputfolder+"/untrimmed_fastq/{short}_R2_001.fastq.gz",
            nonempty_fastqs=outputfolder+"/inputfiles_checkpoint_completed.out"
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
            rm -f {output.fq1} {output.fq2}
            cutadapt {params.adapters} {params.otherParams} -o {output.fq1} -p {output.fq2} {input.in_1} {input.in_2} >> {log} 2>&1
            sha256sum {output} >>{params.sha_sum_file}
            """


    rule align_pe: # aligning pe can only input trimmed  or umi_extracted data, no untrimmed data allowed!
        input:
            pe1 = get_mapping_input_pe1,
            pe2 = get_mapping_input_pe2
        output:
            bams=outputfolder + "/star/{short}_pe_001_Aligned.sortedByCoord.out.bam", #this needs to be deduped 
            tabs=outputfolder + "/star/{short}_pe_001_ReadsPerGene.out.tab"
        log:
            outputfolder+"/logs/star/{short}.log"
        params:
            extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {} {}".format(
                config["mapping"]["gtf_file"], config["mapping"]["other_params"]),
            prefix = outputfolder + "/star/{short}_pe_001_",
            genDir = config["mapping"]["genomeDir"],
            mapping_dir = outputfolder+"/star",
            work_dir = outputfolder + "/star/{short}_pe_001__STARgenome",


        resources:
            threads=lambda wildcards, attempt: attempt * 12,
            time_hrs=lambda wildcards, attempt: attempt * 2,
            mem_gb=lambda wildcards, attempt: attempt * 40
        conda:
            p+"/envs/STAR.yaml"
        shell:
            """
            mkdir -p {params.mapping_dir}
            STAR {params.extra} --genomeDir {params.genDir} --runThreadN {resources.threads} --readFilesIn {input.pe1} {input.pe2} --readFilesCommand zcat --outFileNamePrefix {params.prefix} --outStd Log >> {log} 2>&1
            chmod ago+rwx -R {params.mapping_dir} >> {log} 2>&1
            rm -rf {params.work_dir} >> {log} 2>&1
            """


    if config["umi_tools"]["umi_tools_active"]:
        # the single-end umi_dedup rule is not gated by isSingleEnd and would match the same output files
        ruleorder: umi_dedup_pe > umi_dedup

        rule umi_dedup_pe:
            input:
                outputfolder+"/star/{short}_pe_001_Aligned.sortedByCoord.out.bam"
            params:
                sorted_bam=outputfolder+"/star/{short}_pe_001_sorted.Aligned.sortedByCoord.out.bam",
                logfolder=outputfolder+"/logs/umi_tools/",
                output_stats=outputfolder+"/logs/umi_tools/{short}_pe_dedup_stats.out"
            log:
                outputfolder+"/logs/umi_tools/{short}_pe.umi_dedup.log"
            conda:
                p+"/envs/umi_tools.yaml"
            message: "umi_tools dedup (paired): removing duplicate read pairs from mapping output..."
            resources:
                threads=lambda wildcards, attempt: attempt * 12,
                time_hrs=lambda wildcards, attempt: attempt * 4,
                mem_gb=lambda wildcards, attempt: attempt * 8
            output:
                outputfolder+"/star/{short}_pe_001_deduped.Aligned.sortedByCoord.out.bam"
            shell:
                """
                rm -f {output}
                mkdir -p {params.logfolder}
                samtools sort --threads {resources.threads} {input} -o {params.sorted_bam} 2>{log}
                samtools index {params.sorted_bam} 2>>{log}
                umi_tools dedup -I {params.sorted_bam} --paired --output-stats={params.output_stats} -S {output} --log={log}.dedup 2>>{log}
                """


    rule FastQC_pe:
        input:  # need to test this if this works, also the input fun of kraken rule 
            outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
            outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"
        params:
            out=outputfolder,
            fastp_report=outputfolder+"/fastqc/{short}_after_filtering_fastp.html",
            fastp_json= outputfolder+"/fastqc/{short}_after_filtering_fastp.json",
            fastqc_out_folder= outputfolder+"/fastqc",
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
            fastp -i {input[0]} -I {input[1]} -h {params.fastp_report} -j {params.fastp_json} >> {log} 2>&1
            chmod ago+rwx -R {params.out}/fastqc/ >> {log} 2>&1
            """


    rule FastQC_untrimmed_pe:
        input:
            r1 = outputfolder+"/untrimmed_fastq/{short}_R1_001.fastq.gz",
            r2 = outputfolder+"/untrimmed_fastq/{short}_R2_001.fastq.gz"
        params:
            fastp_report = outputfolder+"/fastqc_untrimmed/{short}_untrimmed_fastp.html",
            fastp_json = outputfolder+"/fastqc_untrimmed/{short}_untrimmed_fastp.json",
            html1 = outputfolder+"/fastqc_untrimmed/{short}_R1_001_fastqc.html",
            zip1 = temp(outputfolder+"/fastqc_untrimmed/{short}_R1_001_fastqc.zip"),
            html2 = outputfolder+"/fastqc_untrimmed/{short}_R2_001_fastqc.html",
            zip2 = temp(outputfolder+"/fastqc_untrimmed/{short}_R2_001_fastqc.zip"),           
            out = outputfolder
        output:
            r1_h = outputfolder+"/fastqc_untrimmed/untrimmed_{short}_R1_001_fastqc.html",
            r1_z = outputfolder+"/fastqc_untrimmed/untrimmed_{short}_R1_001_fastqc.zip" ,
            r2_h = outputfolder+"/fastqc_untrimmed/untrimmed_{short}_R2_001_fastqc.html",
            r2_z = outputfolder+"/fastqc_untrimmed/untrimmed_{short}_R2_001_fastqc.zip" ,
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
            fastqc -q --threads {resources.threads} {input.r1} -o {params.out}/fastqc_untrimmed/ >> {log} 2>&1
            fastqc -q --threads {resources.threads} {input.r2} -o {params.out}/fastqc_untrimmed/ >> {log} 2>&1
            fastp -i {input.r1} -I {input.r2} -h {params.fastp_report} -j {params.fastp_json} >> {log} 2>&1
            mv {params.zip1} {output.r1_z}
            mv {params.html1} {output.r1_h}
            mv {params.zip2} {output.r2_z}
            mv {params.html2} {output.r2_h}
            """


    rule crypt4gh_pe:
        input:
            pe1=outputfolder + "/umi_extract/{short}_R1_001.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
            pe2=outputfolder + "/umi_extract/{short}_R2_001.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"
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
            crypt4gh encrypt --sk {params.priv_key} `bash scripts/get_c4gh_string.sh {params.pub_key}` < {input.pe2} > {output.c4gh_out2} 2>>{log}
            sha256sum {output} >>{params.checksum_file}
            """


    rule sortmerna_pe:
        input:
            sort_1=outputfolder + "/umi_extract/{short}_R1_001.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
            sort_2=outputfolder + "/umi_extract/{short}_R2_001.umis-extracted.fastq.gz" if config["umi_tools"]["umi_tools_active"] else outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"
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
            out_rrna1=outputfolder+"/sortmerna/{short}_R1_001_ribosomal_rna.fq.gz",
            out_rrna_free1=outputfolder+"/sortmerna/{short}_R1_001_non-ribosomal_rna.fq.gz",
            out_rrna2=outputfolder+"/sortmerna/{short}_R2_001_ribosomal_rna.fq.gz",
            out_rrna_free2=outputfolder+"/sortmerna/{short}_R2_001_non-ribosomal_rna.fq.gz",

        output:
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
        # the single-end umi_extract rule is not gated by isSingleEnd and would match the same output files
        ruleorder: umi_extract_pe > umi_extract

        rule umi_extract_pe:
            input:
                in_1=outputfolder+"/untrimmed_fastq/{short}_R1_001.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else outputfolder+"/trimmed/{short}_R1_001_trimmed.fastq.gz",
                in_2=outputfolder+"/untrimmed_fastq/{short}_R2_001.fastq.gz" if not config["cutadapt"]["cutadapt_active"] else outputfolder+"/trimmed/{short}_R2_001_trimmed.fastq.gz"                
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
                out1=outputfolder + "/umi_extract/{short}_R1_001.umis-extracted.fastq.gz",
                out2=outputfolder + "/umi_extract/{short}_R2_001.umis-extracted.fastq.gz"
            shell:
                """
                mkdir -p {params.outputf}
                mkdir -p {params.logfolder}
                umi_tools extract --stdin={input.in_1} --extract-method=regex --bc-pattern={params.umi_ptrn} --log={log} --stdout={output.out1} --read2-in={input.in_2} --read2-out={output.out2} --bc-pattern2={params.ptrn_2} >> {log} 2>&1
                chmod ago+rwx -R {params.outputf}
                sha256sum {output} >>{params.checksum_file}
                """