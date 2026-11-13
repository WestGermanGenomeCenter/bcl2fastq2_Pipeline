from snakemake.utils import min_version
from snakemake.exceptions import WorkflowError
import os
import pandas as pd

min_version("8.28.0") # use-envmodules are not working in versions of 5.10.0 and below for clusters
p = os.path.abspath(".")
#Snakemake configs and setup
configfile: "config.yaml"
samplesheet = config["demux"]["SampleSheet"]   # with skip_demux only the filename (project number) is used, the file itself need not exist
outputfolder = config["demux"]["OutputFolder"]

SKIP_DEMUX = config["skip_demux"]["skip_demux_active"]
# idemux works on top of bclconvert output OR on externally demultiplexed fastqs (skip_demux)
USE_IDEMUX = config["demux"].get("use_idemux", False)


def validate_idemux_sheet(path):
    """Fail at parse time (seconds) instead of after hours of bcl convert."""
    if not path or not os.path.isfile(str(path)):
        raise WorkflowError(f"use_idemux is on, but idemux_samplesheet is empty or not a file: {path!r}")
    df = pd.read_csv(path, dtype=str).fillna("")
    missing = {"pool", "sample_name", "i1"} - set(df.columns)
    if missing:
        raise WorkflowError(f"idemux sheet {path}: missing columns {sorted(missing)}")
    errors = []
    dup_names = df.loc[df["sample_name"].duplicated(), "sample_name"].tolist()
    if dup_names:
        errors.append(f"duplicated sample_name (must be unique across all pools): {dup_names}")
    bad_names = df.loc[~df["sample_name"].str.fullmatch(r"[A-Za-z0-9_.-]+"), "sample_name"].tolist()
    if bad_names:
        errors.append(f"sample_name with characters unsafe for file names: {bad_names}")
    bad_pools = df.loc[~df["pool"].str.fullmatch(r"[A-Za-z0-9_.-]+"), "pool"].tolist()
    if bad_pools:
        errors.append(f"empty or invalid pool names: {bad_pools}")
    bad_i1 = df.loc[~df["i1"].str.fullmatch(r"[ACGTN]+"), ["pool", "sample_name", "i1"]]
    if len(bad_i1):
        errors.append("empty or invalid i1 sequence:\n" + bad_i1.to_string(index=False))
    for pool, g in df.groupby("pool"):
        if g["i1"].duplicated().any():
            errors.append(f"pool {pool}: duplicated i1 barcode")
        if g["i1"].str.len().nunique() > 1:
            print(f"WARNING: pool {pool} has i1 barcodes of different lengths")
    if errors:
        raise WorkflowError("idemux sheet problems:\n" + "\n".join(errors))

if USE_IDEMUX:
    validate_idemux_sheet(config["demux"].get("idemux_samplesheet"))

# provides projectNum, sample_names, isSingleEnd(), getSample_names_post_mapping(), getSamples(), ...
include: "common.smk"


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
    if SKIP_DEMUX:
        all.extend(expand("{out}/skipped_demuxing.flag",out=outputfolder))
    else:
        all.extend(expand("{out}/SampleSheet_check.out",out=outputfolder))

    return all

def getParentDir(wildcards):
    return os.path.dirname(str(config["demux"]["SampleSheet"]))

localrules: all, create_fastq_list

rule all:
    input:
        getOutput()

# demux depends on this flag when prevent_revcomp_active is set (see rule demux)
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


# ----------------------------------------------------------------------------
# primary demultiplexing: bcl convert  (rules: check_inputs, demux)
# ----------------------------------------------------------------------------
if not SKIP_DEMUX:

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
            {params.bcl_convert_path} --bcl-input-directory {params.infolder} --sample-sheet {input.samplesheet} --output-directory {params.output_dir} --bcl-validate-sample-sheet-only true  --force >>{output} 2>&1
            cat {params.Logs_path}/*.log >>{output} 
            rm -rf {params.output_dir}/Reports {params.output_dir}/Logs # just empty?
            python {params.script_for_checking} {input.samplesheet} >>{params.script_out} 2>/dev/null # should not be the reason for this rule to fail
            """

    rule demux:
        input:
            samplesheet = config["demux"]["SampleSheet"],
            samplesheet_errors = config["demux"]["OutputFolder"]+"/SampleSheet_check.out",
            # makes the revcomp edit a real dependency instead of relying on tiny resources
            revcomp = [config["demux"]["OutputFolder"]+"/revcomp_prevented.flag"] if config["prevent_revcomp"]["prevent_revcomp_active"] else []
        params:
            infolder = getParentDir,
            additionalOptions = config["demux"]["options"],
            bcl_convert_path = config["demux"]["bcl_convert_path"],
            output_dir = config["demux"]["OutputFolder"],
            outfastqs = config["demux"]["OutputFolder"] + "/*.fastq.gz",
            out_fastqs_dir = config["demux"]["OutputFolder"] + "/untrimmed_fastq/"
        log:
            config["demux"]["OutputFolder"]+"/logs/e_bcl.log"
        output:
            demux2Output = config["demux"]["OutputFolder"]+"/Reports/Adapter_Metrics.csv"
        message:
            "Running bcl convert"
        resources:
            threads=lambda wildcards, attempt:  12 + (attempt * 4),
            time_hrs=lambda wildcards, attempt: attempt * 4,
            mem_gb=lambda wildcards, attempt: round(config["demux"]["demux_try1_gb"] * 0.7 ) + (config["demux"]["demux_try1_gb"] * ( 0.3 * attempt) )
        conda:
            p+"/envs/demux.yaml"
        shell:
            """
            exec >> {log} 2>&1
            mkdir -p {params.output_dir}
            chmod ago+rwx -R {params.output_dir}
            # clear leftovers so nothing is demultiplexed into old output
            rm -rf {params.out_fastqs_dir}
            rm -f {params.outfastqs}
            {params.bcl_convert_path} --bcl-input-directory {params.infolder} --sample-sheet {input.samplesheet} {params.additionalOptions} --output-directory {params.output_dir} --force --bcl-num-conversion-threads {resources.threads}
            """


# ----------------------------------------------------------------------------
# demultiplexing was done elsewhere (e.g. on-device): hard-link the fastqs in place
# no samplesheet file is needed in this mode
# ----------------------------------------------------------------------------
if SKIP_DEMUX:

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
            set -euo pipefail
            shopt -s nullglob
            mkdir -p {params.out_folder} {params.reports_dir}
            touch {params.reports_file}

            # reset: only removes links, the originals stay untouched
            rm -rf {params.untrimmed_fastq_folder}

            files=({params.dir_w_fastq}/*.fastq.gz)
            if [ "${{#files[@]}}" -eq 0 ]; then
                echo "ERROR: no *.fastq.gz in {params.dir_w_fastq}"; exit 1
            fi
            for f in "${{files[@]}}"; do
                d={params.out_folder}/$(basename "$f")
                if [ ! "$f" -ef "$d" ]; then
                    ln -f "$f" "$d" || {{ echo "ERROR: hard link failed for $f (different filesystem, or not owner/writable: fs.protected_hardlinks)"; exit 1; }}
                fi
            done
            touch {output.flagfile}
            """


# upstream output of the primary demultiplexing, whichever way it happened
def upstream_demux_output():
    if SKIP_DEMUX:
        return rules.skip_demux.output
    return rules.demux.output


# ----------------------------------------------------------------------------
# organize: shared by bcl convert and skip_demux
#  - sample fastqs -> untrimmed_fastq/, Undetermined stays in the output root
#  - fails with a clear message if there are no sample fastqs
#  - the samplesheet is only copied when there is one (not with skip_demux)
# ----------------------------------------------------------------------------
rule organize_fastqs:
    input:
        upstream = upstream_demux_output(),
        samplesheet = [] if SKIP_DEMUX else config["demux"]["SampleSheet"]
    params:
        out = config["demux"]["OutputFolder"],
        fq_dir = config["demux"]["OutputFolder"] + "/untrimmed_fastq"
    log:
        config["demux"]["OutputFolder"]+"/logs/organize_fastqs.log"
    output:
        flag = config["demux"]["OutputFolder"]+"/demux_organized.flag"
    message:
        "Sorting demultiplexed fastqs into place"
    resources:
        threads=lambda wildcards, attempt: attempt * 1,
        time_hrs=lambda wildcards, attempt: attempt * 1,
        mem_gb=lambda wildcards, attempt: attempt * 2
    conda:
        p+"/envs/demux.yaml"
    shell:
        """
        set -euo pipefail
        exec >> {log} 2>&1
        mkdir -p {params.fq_dir}

        # sample fastqs -> untrimmed_fastq (Undetermined stays in the output root)
        find {params.out} -maxdepth 1 -name '*.fastq.gz' ! -name 'Undetermined*' -exec mv -f -t {params.fq_dir} {{}} +

        if [ -z "$(find {params.fq_dir} -name '*.fastq.gz' ! -name 'Undetermined*' -print -quit)" ]; then
            echo "ERROR: no sample fastq.gz in {params.fq_dir} (everything undetermined? check samplesheet / index settings)"
            exit 1
        fi

        if [ -n "{input.samplesheet}" ]; then cp -f {input.samplesheet} {params.out}/; fi
        if [ -d {params.out}/Logs ]; then
            rm -rf {params.out}/logs/bcl_convert_logs
            mkdir -p {params.out}/logs
            mv -f {params.out}/Logs {params.out}/logs/bcl_convert_logs
        fi
        touch {output.flag}
        """


# ----------------------------------------------------------------------------
# idemux: i1 inline barcode demultiplexing of the pool fastqs from the primary demux
#  - pool fastqs are parked in <OutputFolder>/<primary_demux_dir>/ (default primary_demux_out)
#  - per-sample fastqs land in untrimmed_fastq/ named like bcl convert output
#  - undetermined reads go to <OutputFolder>/Undetermined_<pool>_R*.fastq.gz
#  - pools not listed in the idemux sheet are passed through unchanged
# ----------------------------------------------------------------------------
if USE_IDEMUX:

    rule idemux:
        input:
            primary = rules.organize_fastqs.output,
            idemux_sheet = config["demux"]["idemux_samplesheet"]
        params:
            out_dir     = config["demux"]["OutputFolder"],
            fastq_dir   = config["demux"]["OutputFolder"] + "/untrimmed_fastq",
            idemux_dir  = config["demux"]["OutputFolder"] + "/idemux",
            primary_dir = config["demux"]["OutputFolder"] + "/" + config["demux"].get("primary_demux_dir", "primary_demux_out"),
            i1_start    = config["demux"].get("idemux_i1_start", 11),
            options     = config["demux"].get("idemux_options", ""),
            exec_path   = config["demux"].get("idemux_path", "idemuxCPP")
        log:
            config["demux"]["OutputFolder"] + "/logs/idemux.log"
        output:
            flag = config["demux"]["OutputFolder"] + "/used_idemux.flag"
        message:
            "Running idemuxCPP (i1 inline barcodes) on the pool fastqs"
        resources:
            threads = lambda wildcards, attempt: 8 + attempt * 4,
            time_hrs = lambda wildcards, attempt: attempt * 2,
            mem_gb = lambda wildcards, attempt: attempt * 16
        conda:
            p + "/envs/idemux.yaml"   # needs idemuxcpp from bioconda
        shell:
            """
            set -euo pipefail
            exec >> {log} 2>&1
            shopt -s nullglob

            FQ=$(realpath {params.fastq_dir})
            OUT=$(realpath -m {params.idemux_dir})
            PRIM=$(realpath -m {params.primary_dir})
            ROOT=$(realpath {params.out_dir})
            SHEET=$(realpath {input.idemux_sheet})
            MANIFEST=$OUT/moved_samples.txt
            mkdir -p "$OUT" "$PRIM"

            # --- rerun handling ---
            # drop sample files of a previous idemux run
            if [ -s "$MANIFEST" ]; then xargs -a "$MANIFEST" -I{{}} rm -f "$FQ/{{}}"; fi
            : > "$MANIFEST"
            # restore parked pool fastqs, unless the primary demux was rerun (fresh file exists -> parked one is stale)
            for f in "$PRIM"/*.fastq.gz; do
                b=$(basename "$f")
                if [ -e "$FQ/$b" ]; then rm -f "$f"; else mv "$f" "$FQ/"; fi
            done

            # --- every pool in the idemux sheet must have fastqs ---
            for pool in $(tail -n +2 "$SHEET" | tr -d '\\r' | cut -d, -f1 | grep -v '^$' | sort -u); do
                files=("$FQ/${{pool}}"_S*_R1_001.fastq.gz)
                if [ "${{#files[@]}}" -eq 0 ]; then
                    echo "ERROR: pool '$pool' is in the idemux sheet but has no R1 fastq in $FQ"
                    exit 1
                fi
            done

            cd "$FQ"
            for r1 in $(ls | grep -E '_R1_001\\.fastq\\.gz$' | grep -v Undetermined); do
                key=$(echo "$r1" | sed -E 's/_R1_001\\.fastq\\.gz$//')              # Pool1_S1_L001 or Pool1_S1
                pool=$(echo "$key" | sed -E 's/_S[0-9]+(_L[0-9]+)?$//')             # Pool1
                lane=$(echo "$key" | sed -nE 's/.*_S[0-9]+(_L[0-9]+)$/\\1/p')        # _L001 or empty
                r2=$(echo "$r1" | sed 's/_R1_001/_R2_001/')
                psheet="$OUT/$pool.sheet.csv"

                # per-pool sheet in idemux format, i7/i5 empty (CRLF-safe)
                tr -d '\\r' < "$SHEET" | awk -F, -v p="$pool" 'NR==1{{print "sample_name,i7,i5,i1"; next}} $1==p{{print $2",,,"$3}}' > "$psheet"
                if [ "$(wc -l < "$psheet")" -le 1 ]; then
                    echo "pool $pool not in idemux sheet -> passed through unchanged"
                    rm -f "$psheet"; continue
                fi

                mkdir -p "$OUT/$key"
                if [ -f "$r2" ]; then
                    {params.exec_path} --r1 "$r1" --r2 "$r2" --sample-sheet "$psheet" --out "$OUT/$key" \\
                        --i1-start {params.i1_start} -r 2 -p {resources.threads} -w {resources.threads} {params.options}
                else
                    {params.exec_path} --r1 "$r1" --sample-sheet "$psheet" --out "$OUT/$key" \\
                        --i1-read 1 --i1-start {params.i1_start} -r 1 -p {resources.threads} -w {resources.threads} {params.options}
                fi

                # rename idemux output to bcl convert style, lane kept so lanes cannot overwrite each other
                for f in "$OUT/$key"/*.fastq.gz; do
                    b=$(basename "$f")
                    case "$b" in
                        undetermined_*) mv -f "$f" "$ROOT/Undetermined_${{key}}_${{b#undetermined_}}" ;;
                        *_R1.fastq.gz)  nb="${{b%_R1.fastq.gz}}_S0${{lane}}_R1_001.fastq.gz"; mv -f "$f" "$FQ/$nb"; echo "$nb" >> "$MANIFEST" ;;
                        *_R2.fastq.gz)  nb="${{b%_R2.fastq.gz}}_S0${{lane}}_R2_001.fastq.gz"; mv -f "$f" "$FQ/$nb"; echo "$nb" >> "$MANIFEST" ;;
                        *)              echo "WARNING: unexpected idemux output name $b, moved as is"; mv -f "$f" "$FQ/"; echo "$b" >> "$MANIFEST" ;;
                    esac
                done

                # park the pool fastqs so they don't end up in the list/checksums
                mv -f "$r1" "$PRIM/"
                if [ -f "$r2" ]; then mv -f "$r2" "$PRIM/"; fi
            done
            touch {output.flag}
            """


def getlist_input():
    if USE_IDEMUX:
        return rules.idemux.output
    return rules.organize_fastqs.output


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
        "cd {params.out_fastqs_dir} && find . -type f -name '*.fastq.gz' ! -name 'Undetermined*' -printf '%P\\n' | sort > {output}"

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
        """
        LIST=$(realpath {input.fastq_list_file})
        cd {params.out}
        xargs -a "$LIST" -n 1 -P {resources.threads} sha256sum | sort -k2 > {output}
        """


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
        script=p+"/scripts/test.py",
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
        python {params.script} {params.raw_files_place} {params.cluster_file}  >> /dev/null 2>&1
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