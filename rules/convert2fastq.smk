from snakemake.utils import min_version
import os
import sys
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum
import pandas as pd
import re

min_version("8.28.0") # use-envmodules are not working in versions of 5.10.0 and below for clusters
p = os.path.abspath(".")
#Snakemake configs and setup
configfile: "config.yaml"
samplesheet = config["demux"]["SampleSheet"]
outputfolder = config["demux"]["OutputFolder"]

include: "common.smk"


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




sample_names = list()
if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
    fastqs = list(samples_dataframe.iloc[:, 0].values)
    sample_names = [fastq[:-9] for fastq in fastqs if "_I1_" not in fastq and "_I2_" not in fastq] # now excluding index reads in the sample names list








# Check for projectNum in samplesheet name
projectNum=validateProjectNum(samplesheet)
if not projectNum:
    print("There was no project number to be found in the samplesheet and/or the name of the samplesheet")

# Set expected pipeline output
def getOutput():
    all = list()
    all.extend(expand("{out}/fastq_infiles_list.tx",out=outputfolder))
    all.extend(expand("{out}/untrimmed_fastq/{num}_sha256sums_fastqfiles.sha256",out=outputfolder,num=projectNum))# projectNum+"_sha256sums_fastqfiles.txt"
    all.extend(expand("{out}/{num}_software_environment_mqc_versions.yml",out=outputfolder,num=projectNum))
    all.extend(expand("{out}/{num}_config.yaml",out=outputfolder,num=projectNum))
    if config["interop_plots"]["interop_plots_active"]:
        all.extend(expand("{out}/interop_plots/interop_plots_done.flag",out=outputfolder))
    if config["prevent_revcomp"]["prevent_revcomp_active"]:
        all.extend(expand("{out}/revcomp_prevented.flag",out=outputfolder))
        all.extend(expand("{out}/reversed_revcomp_prevented.flag",out=outputfolder))
    if config["skip_demux"]["skip_demux_active"]:
        all.extend(expand("{out}/skipped_demuxing.flag",out=outputfolder))
    if not config["skip_demux"]["skip_demux_active"]:
        if not config["demux"]["use_bcl2fastq"]:
            all.extend(expand("{out}/SampleSheet_check.out",out=outputfolder))# this uses convert to check,
    if config["demux"]["use_bcl2fastq"]:
        all.extend(expand("{out}/used_bcl2fastq2.flag",out=outputfolder))

    return all

def getParentDir(wildcards):
    return os.path.dirname(str(config["demux"]["SampleSheet"]))

localrules: all, create_fastq_list

rule all:
    input:
        getOutput()

# this rule has no parameter/flag/file that it needs to be executed before demultiplexing
# this works just because the ressources for this job are tiny, such that it will 
# all the time start before the demultiplexing starts
rule prevent_revcomp:
    input:
        sample_sheet=config["demux"]["SampleSheet"]
    params:
        folder_of_runxml=getParentDir,
        script_to_fix="scripts/prevent_revcomp_nexts2.sh",
        out = config["demux"]["OutputFolder"]
    output:
        config["demux"]["OutputFolder"]+"/revcomp_prevented.flag" 
    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: attempt * 2
    conda:
        p+"/envs/demux.yaml"
    shell:
        """
        mkdir -p {params.out}
        bash {params.script_to_fix} {params.folder_of_runxml}
        touch {output}
        """



rule reverse_prevent_revcomp: # this rule only exists to reverse the edit of the runinfo.xml if the read3 revcomp option was activated
    input:
        demux_done=config["demux"]["OutputFolder"]+"/Reports/Adapter_Metrics.csv",
        sample_sheet=config["demux"]["SampleSheet"],
        step1_complete=config["demux"]["OutputFolder"]+"/revcomp_prevented.flag" 

    params:
        folder_of_runxml=getParentDir,
        script_to_fix="scripts/reverse_revcomp_nexts2.sh",
        out = config["demux"]["OutputFolder"]
    output:
        config["demux"]["OutputFolder"]+"/reversed_revcomp_prevented.flag" 
    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: attempt * 2
    conda:
        p+"/envs/demux.yaml"
    shell:
        """
        mkdir -p {params.out}
        bash {params.script_to_fix} {params.folder_of_runxml}
        touch {output}
        """



if not config["skip_demux"]["skip_demux_active"]:

    rule check_inputs:
        input:
            samplesheet = config["demux"]["SampleSheet"]
        params:
            infolder = getParentDir,
            bcl_convert_path = config["demux"]["bcl_convert_path"],
            output_dir = config["demux"]["OutputFolder"],
            Logs_path= config["demux"]["OutputFolder"] + "/Logs",
            script_for_checking= p+"/scripts/bclconvert_samplesheet_check.py",
            script_out=config["demux"]["OutputFolder"]+"/SampleSheet_check_script.out",
            samshee_out=config["demux"]["OutputFolder"]+"/SampleSheet_check_samshee.out",
        output:
            samplesheet_errors = config["demux"]["OutputFolder"]+"/SampleSheet_check.out"
        message:
            "Checking the samplesheet with Bcl Convert and the python script"
        resources:
            threads=lambda wildcards, attempt: attempt * 1,
            time_hrs=lambda wildcards, attempt: attempt * 1,
            mem_gb=lambda wildcards, attempt: attempt * 4
        conda:
            p+"/envs/demux.yaml"
        shell:
            """
            mkdir -p {params.output_dir}
            {params.bcl_convert_path} --bcl-input-directory {params.infolder} --sample-sheet {input.samplesheet} --output-directory {params.output_dir} --bcl-validate-sample-sheet-only true  --force 2>&1>>{output}
            cat {params.Logs_path}/*.log >>{output} 
            rm -rf {params.output_dir}/Reports {params.output_dir}/Logs # just empty?
            python {params.script_for_checking} {input.samplesheet} >>{params.script_out} 2>/dev/null # should not be the reason for this rule to fail
            """


    rule demux:
        input:
            samplesheet = config["demux"]["SampleSheet"],
            samplesheet_errors = config["demux"]["OutputFolder"]+"/SampleSheet_check.out"
        params:
            infolder = getParentDir,
            additionalOptions=config["demux"]["options"],
            bcl_convert_path = config["demux"]["bcl_convert_path"],
            outfastqs = config["demux"]["OutputFolder"] + "/*.fastq.gz", # we want to exclude the undetermined maybe here, maybe after the multiqc
            output_dir = config["demux"]["OutputFolder"],
            undetermined = config["demux"]["OutputFolder"]  + "/untrimmed_fastq/" + "Undetermined*.fastq.gz",
            out_fastqs_dir = config["demux"]["OutputFolder"] + "/untrimmed_fastq/",
            Logs_path= config["demux"]["OutputFolder"] + "/Logs",
            bcl_log_path=config["demux"]["OutputFolder"] + "/logs/bcl_convert_logs"
        log:
            config["demux"]["OutputFolder"]+"/logs/e_bcl.log"
        output:
            demux2Output = config["demux"]["OutputFolder"]+"/Reports/Adapter_Metrics.csv"
        message:
            "Running bcl convert"
        resources:
            threads=lambda wildcards, attempt:  12 + (attempt * 4) ,
            time_hrs=lambda wildcards, attempt: attempt * 4,
            mem_gb=lambda wildcards, attempt: round(config["demux"]["demux_try1_gb"] * 0.7 ) +(config["demux"]["demux_try1_gb"] * ( 0.3 * attempt) )
        conda:
            p+"/envs/demux.yaml"
        shell:
            """
            mkdir -p {params.output_dir} >> {log} 2>&1
            chmod ago+rwx -R {params.output_dir} >> {log} 2>&1
            rm -rf {params.out_fastqs_dir}  >> {log} 2>&1 # to prevent 2 times demuxing into the same folder
            rm -f {params.outfastqs}  >> {log} 2>&1 
            {params.bcl_convert_path} --bcl-input-directory {params.infolder} --sample-sheet {input.samplesheet} {params.additionalOptions}  --output-directory {params.output_dir} --force  --bcl-num-conversion-threads {resources.threads} >> {log} 2>&1
            cp -f {input[0]} {params.output_dir} >> {log} 2>&1
            mkdir -p {params.out_fastqs_dir} && mv {params.outfastqs} {params.out_fastqs_dir} >> {log} 2>&1
            mv -f {params.undetermined} {params.output_dir} 2>/dev/null
            mv -f {params.Logs_path} {params.bcl_log_path} >> {log} 2>&1
            """

if config["demux"]["use_bcl2fastq"]: # allow miseq also into the mix
    rule bcl2fastq2:
        input:
            input = config["demux"]["SampleSheet"]
        params:
            barcode_mismatches = config["demux"]["bcl2fastq_mismatches"],
            infolder = getParentDir,
            additionalOptions=config["demux"]["bcl2fastq_options"],
            out = config["demux"]["OutputFolder"],
            fastq_destination=config["demux"]["OutputFolder"]+"/untrimmed_fastq/",
            exec_path = config["demux"]["bcl2fastq2_path"]
        log:
            config["demux"]["OutputFolder"]+"/logs/e_bcl.log"
        output:
            bcl2fastq2Output = config["demux"]["OutputFolder"]+"/used_bcl2fastq2.flag"
        message:
            "Running bcl2fastq2 - this is for older samplesheets"  
        resources:
            threads=lambda wildcards, attempt: attempt * 24,
            time_hrs=lambda wildcards, attempt: attempt * 1,
            mem_gb=lambda wildcards, attempt: attempt * 32
        conda:
            p+"/envs/bcl2fastq2.yaml"

        shell:
            """
            mkdir -p {params.out} >> {log} 2>&1
            chmod ago+rwx -R {params.out} >> {log} 2>&1
            rm -rf {params.fastq_destination} >> {log} 2>&1
            {params.exec_path} -R {params.infolder} --sample-sheet {input[0]} {params.additionalOptions} --barcode-mismatches {params.barcode_mismatches} -o {params.out} --interop-dir {params.out}  -r {resources.threads} -p {resources.threads} 2> {log[0]}
            cp -f {input[0]} {params.out} >> {log} 2>&1
            mkdir -p {params.fastq_destination} >> {log} 2>&1
            mv -f {params.out}/*/*.fastq.gz {params.fastq_destination} >> {log} 2>&1
            touch {output}
           """

if config["skip_demux"]["skip_demux_active"]:

    rule skip_demux:
        input:
        params:
            dir_w_fastq=config["skip_demux"]["fastq_folder"],
            out_folder=outputfolder,
            reports_dir=outputfolder +"/Reports",
            reports_file=outputfolder +"/Reports/Adapter_Metrics.csv",
            untrimmed_fastq_folder=outputfolder +"/untrimmed_fastq/",
        output:
            flagfile=outputfolder+"/skipped_demuxing.flag",
        resources:
            threads=lambda wildcards, attempt: attempt * 1,
            time_hrs=lambda wildcards, attempt: attempt * 1,
            mem_gb=lambda wildcards, attempt: attempt * 2
        conda:
            p+"/envs/demux.yaml"
        shell:
            """
            mkdir -p {params.out_folder}
            mkdir -p {params.reports_dir}
            mkdir -p {params.untrimmed_fastq_folder}
            touch {params.reports_file}
            mv -f {params.dir_w_fastq}/*.gz {params.untrimmed_fastq_folder}
            touch {output.flagfile}
            """


def getlist_input():# 
    if config["skip_demux"]["skip_demux_active"]:
        return rules.skip_demux.output
    else:
        if config["demux"]["use_bcl2fastq"]:
            return rules.bcl2fastq2.output
        else:
            return rules.demux.output




rule create_fastq_list:
    input:
        fastqs_there= getlist_input()
    params:
        out_fastqs_dir = config["demux"]["OutputFolder"] + "/untrimmed_fastq/"
    output:
        fastq_list_file=config["demux"]["OutputFolder"]+"/fastq_infiles_list.tx"
    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: attempt * 2
    conda:
        p+"/envs/demux.yaml"
    shell:
        " cd {params.out_fastqs_dir} && shopt -s globstar; ls -d ** | grep '.fastq.gz' | grep -v 'Undetermined' >{output}"

rule create_checksums:
    input:
        fastq_list_file=config["demux"]["OutputFolder"]+"/fastq_infiles_list.tx"
    params:
        out=config["demux"]["OutputFolder"]+"/untrimmed_fastq/"
    output:
        checksum_file=config["demux"]["OutputFolder"]+"/untrimmed_fastq/"+projectNum+"_sha256sums_fastqfiles.sha256"
    resources:
        threads=lambda wildcards, attempt: attempt * 16,
        time_hrs=lambda wildcards, attempt: attempt * 2,
        mem_gb=lambda wildcards, attempt: attempt * 8
    conda:
        p+"/envs/demux.yaml"
    shell:
        "cd {params.out} && shopt -s globstar; ls -d ** | grep '.fastq.gz' | grep -v 'Undetermined' | xargs -n 1 -P {resources.threads} sha256sum >>{output}"


rule copy_software_env:
    input:
        versions_file=config["software_versions_file"]
    params:
        out=config["demux"]["OutputFolder"] # not used
    output:
        env_file=config["demux"]["OutputFolder"]+"/"+projectNum+"_software_environment_mqc_versions.yml",
        config_file=config["demux"]["OutputFolder"]+"/"+projectNum+"_config.yaml"
    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: attempt * 1
    conda:
        p+"/envs/demux.yaml"
    shell:
        "cp {input.versions_file} {output.env_file} && cp config.yaml {output.config_file}"

# add interop_more.py : python interop_more.py ../../Nextseq2000-Run324/260113_VL00286_324_AAHLCG3M5 .
rule interop_plots:
    input:
        fastqs_there= getlist_input()
    output:
        done_flag=outputfolder+"/interop_plots/interop_plots_done.flag"
    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: attempt * 22
    params:
        tilt_script=p+"/scripts/tilt_report.py",
        tilt_pdf=config["demux"]["OutputFolder"]+"/interop_plots/TiltReport.pdf",
        script=p+"/scripts/interop_minimal.py",
        # interop_dir=config["interop_plots"]["interop_binaries_dir"],
        log_dir=config["demux"]["OutputFolder"]+"/logs/interop/",
        output_dir=config["demux"]["OutputFolder"]+"/interop_plots",
        script_interop=p+"/scripts/interop_plots.sh",
        cluster_pdf_script=p+"/scripts/cluster_pdf.py",
        cluster_pdf=config["demux"]["OutputFolder"]+"/interop_plots/RunReport.pdf",
        cluster_file=outputfolder+"/interop_plots/Sequencing_report.pdf",
        demux_py=p+"/scripts/demux_report.py",
        demux_pdf=outputfolder+"/interop_plots/Demux_report.pdf",
        reports_dir=outputfolder+"/Reports",
        raw_files_place=getParentDir        

    conda:
        p+"/envs/interop.yaml"
    log:
        outputfolder+"/logs/interop/"+projectNum+"_interop_plots.log"
    shell:
        """ 
        mkdir -p {params.log_dir}
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        bash {params.script_interop} {params.raw_files_place} {params.output_dir}  >> /dev/null 2>&1 
        # add the get_cluster_numbers.py script, redirect errors to /dev/null
        # alternatively to the script, one can also interop_plot_q_score path/to/runfolder | gnuplot. the interop package is already in anaconda: illumina-interop
        # interop via pip wont install on the hpc, so conda install must be enough for now. supressing all errors because of this - should work on any modern machine
        python {params.script} -i {params.raw_files_place} -o {params.cluster_file}  >> /dev/null 2>&1
        python {params.cluster_pdf_script} {params.raw_files_place} {params.cluster_pdf} >> {log} 2>&1
        python {params.demux_py} -i {params.reports_dir} -o {params.demux_pdf} >> {log} 2>&1
        python {params.tilt_script} {params.raw_files_place} {params.tilt_pdf} >> {log} 2>&1
        touch {output}
        """

onsuccess:
    print("Workflow finished without errors")

onerror:
    print("An error occurred during runtime")

onstart:
    print("Setting up and running pipeline")
