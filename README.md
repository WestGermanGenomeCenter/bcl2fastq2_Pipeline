# Demultiplexing & RNA-Seq Pipeline for Illumina Platforms

This Snakemake 8.28+ pipeline automates Illumina BCL demultiplexing (bcl-convert or bcl2fastq2) and provides downstream quality control, contamination detection, and RNA-Seq analysis on HPC environments (PBS Pro / HILBERT).

## Table of Contents

1. Overview & Key Features
2. Prerequisites & System Requirements
3. Deployment & Quick Start
4. Running the Pipeline & Execution Options
5. Pipeline Stages & Architecture
6. The check_fastq_files Rule & check_input_fastqs.sh
7. Paired-End (PE) vs. Single-End (SE) Layout Rules
8. Automation & Autostart Cronjob Setup
9. Complete Tool Reference & Setup Instructions
10. Rule Registry (Snakefiles & .smk includes)
11. Scripts Directory Guide (scripts/)
12. Checksums, Cryptography & MultiQC Reports
13. Configuration Reference (config.yaml)
14. rMATS Group File Setup
15. Troubleshooting & Common Pitfalls
16. Known Code Issues & Verification Notes

---

## 1. Overview & Key Features

### Core Capabilities

- **Dual Demultiplexer Support**: bcl-convert (recommended for NovaSeq X / NextSeq 2000) or bcl2fastq2 (legacy / MiSeq fallback). Controlled by `demux.use_bcl2fastq`.
- **Flexible Data Input**: raw BCL run directories or pre-demultiplexed FASTQs via `skip_demux.skip_demux_active`.
- **Read Layout Agnostic**: PE and SE supported. Layout is inferred at parse time from FASTQ filenames via `isSingleEnd()`, using the list cached in `fastq_infiles_list.tx`.
- **UMI Handling**: single/dual UMI regex patterns (`umi_tools.pattern`, `umi_tools.pattern2`); UMIs moved to read headers (`umi_extract/`); BAM deduplication post-mapping (`umi_dedup` -> `*_deduped.Aligned.sortedByCoord.out.bam`). The PE + UMI path is not fully wired, see section 16.
- **Contamination Screening**: Kraken2, BioBloom Tools, BLAST+, Kaiju, Sourmash. (Diamond exists only as a commented-out `diamond_pe` rule; there is no active Diamond rule.)
- **End-to-End RNA-Seq**: STAR (`--quantMode GeneCounts`) -> Qualimap / RSeQC / SAMtools -> FeatureCounts (`counts/all.tsv`) -> StringTie (3 stages) + gffcompare -> rMATS (`rmats/summary.txt`).

### Advanced Features

- **Data Security & Encryption**: per-file crypt4gh with one aggregate `encrypted_fastqs/checksums_encrypted_fastqs.sha256`; or directory-wide gocryptfs (which reduces the output folder to only `encrypted/`).
- **Automated Data Transfer**: `rsync -rzvP` via rule `transfer_files` (config `transfer.transfer_to`, `transfer.transfer_as`), not scp.
- **Reproducibility & Auditing**: `<projectNum>_software_environment_mqc_versions.yml`, `<projectNum>_config.yaml`, per-run `filelist_project_<timestamp>.filelist` and `checksums_<timestamp>.sha256`, `<projectNum>_sha256sums_fastqfiles.sha256`, per-stage rulegraph PDFs and Snakemake HTML reports (wrapper-level).
- **Smart Scaling**: rules use attempt-scaled `threads`, `time_hrs`, `mem_gb` resources. Demux threads = `12 + attempt * 4`; demux memory is `round(demux_try1_gb * 0.7) + demux_try1_gb * 0.3 * attempt` (approx. 96 GB attempt 1, approx. 125 GB attempt 2 at default `demux_try1_gb: 96`). Kraken2 (`178 + attempt * 12`) and Kaiju (`88 + attempt * 40`) use custom memory formulas; sourmash is fixed at 1 thread.

---

## 2. Prerequisites & System Requirements

### Compute & Storage

- HPC Environment: PBS Pro (tested on HILBERT HPC) or local multi-core.
- Memory (attempt 1; scales on retry):
  - Demux: `demux_try1_gb` default 96 GB; set 126+ for NovaSeq X.
  - Kraken2: `178 + attempt * 12` GB -> approx. 190 GB attempt 1.
  - Kaiju: `88 + attempt * 40` GB.
  - STAR align: requests `attempt * 40` GB. (STAR genome indexing is not part of this repo; a pre-built index is expected, typically ~30-35 GB RAM to build, organism-dependent.)
  - Qualimap: requests `attempt * 12` GB from the scheduler but starts Java with a hardcoded `--java-mem-size=16G` (see section 16).
- Storage: raw BCLs, intermediate FASTQs, aligned BAMs, and (if enabled) `encrypted/`. Several rules leave large intermediates behind (section 16).

### Required Software

1. Conda / Mamba. The repo ships `smk8_environment.yml` and per-tool `envs/*.yaml`.
2. Snakemake >= 8.28.0 (`min_version("8.28.0")` in both entry Snakefiles).
3. bcl-convert and/or bcl2fastq2 binaries from Illumina (credentials required). Paths: `demux.bcl_convert_path`, `demux.bcl2fastq2_path`.
4. Qualimap executable (manual install or `--no-deps` on HPC; see `envs/install_Qualimap_manually.sh`). Path: `mapping.qualimap_dir`. The rule calls `{qualimap_dir}/qualimap` directly.
5. A PBS Snakemake profile named `pbs` (see section 3; the wrapper calls `--profile pbs`).
6. gocryptfs (`gocryptfs.gocryptfs_version`, e.g. `gocryptfs/1.7.1`) as a cluster module if `gocryptfs.gocryptfs_active: True` (loaded via `envmodules:`; needs cluster module support).
7. graphviz (`dot`), needed by `runPipeline.sh` for the rulegraph PDFs.

---

## 3. Deployment & Quick Start

### HILBERT HPC Users

```bash
git clone https://github.com/WestGermanGenomeCenter/bcl2fastq2_Pipeline.git
cd bcl2fastq2_Pipeline
bash setup_pipeline.sh
```

`setup_pipeline.sh` creates the `smk8` conda env and copies a ready-made PBS Snakemake profile from a HILBERT-internal path (`/gpfs/project/projects/bmfz_gtl/software/snakemake_profile/pbs`) to `~/.config/snakemake/pbs`. This only works for HILBERT users.

### Non-HILBERT / Custom Environments

```bash
git clone https://github.com/WestGermanGenomeCenter/bcl2fastq2_Pipeline.git
cd bcl2fastq2_Pipeline
conda env create -n smk8 -f smk8_environment.yml
conda activate smk8
# create a PBS snakemake profile named "pbs" (the snakemake-utils/ directory
# ships PBS profile components: pbs-submit.py, pbs-status.py, pbs-jobscript.sh)
# edit config.yaml
bash runPipeline.sh
```

Note: `runPipeline.sh` hardcodes `source ~/conda/etc/profile.d/conda.sh` and `conda activate smk8`. Adjust both if your conda lives elsewhere or the env has a different name.

There is no top-level `Snakefile`. The pipeline has two entry Snakefiles, both driven through the wrappers:

- `rules/convert2fastq.smk`: demultiplexing stage (includes `common.smk`)
- `rules/fastqProcessing.smk`: processing/analysis stage (includes `pe_processing.smk`, `common.smk`, `qc.smk`; a leftover `# include: "qcShort.smk"` is commented out and harmless)

### Repository Layout (observed)

- `rules/`: the two entry Snakefiles plus `common.smk`, `pe_processing.smk`, `qc.smk`.
- `envs/`: one YAML per tool: `Qualimap.yaml`, `STAR.yaml`, `bcl2fastq2.yaml`, `biobloom.yaml`, `blast.yaml`, `counts.yaml`, `crypt4gh.yaml`, `cutadapt.yaml`, `demux.yaml`, `fastqc.yaml`, `gffcompare.yaml`, `gffutils.yaml`, `interop.yaml`, `jellyfish.yaml`, `kaiju.yaml`, `kraken2.yaml`, `multiqc.yaml`, `pandas.yaml`, `pca_plot.yaml`, `preseq.yaml`, `rmats.yaml`, `rseqc.yaml`, `rustqc.yaml`, `samtools.yaml`, `sortmerna.yaml`, `sourmash.yaml`, `stringtie.yaml`, `transfer.yaml`, `umi_tools.yaml` (plus `tools.yaml` and helper scripts for the software-version list).
- `scripts/`: see section 11.
- `smk8_environment.yml`: Snakemake driver env.
- `autostart/`: cron-based autonomous execution (see section 8).
- `schemas/config.schema.yaml`: legacy schema; not referenced by any Snakefile and out of date with the current config.yaml (do not rely on it for validation).
- `examples/`: example config, rMATS group file, and rendered reports/graphs.
- `snakemake-utils/`: PBS profile components (submit/status scripts, jobscript, test workflow).
- `runPipeline.sh`, `runPipeline_local.sh`, `runPipeline_barebone.sh`, `setup_pipeline.sh`, `setup_pipeline_local.sh`.

---

## 4. Running the Pipeline & Execution Options

### Standard HPC Execution (runPipeline.sh)

```bash
screen -dRR pipeline_run
qsub -I -l select=1:ncpus=2:mem=8gb,walltime=48:00:00 OR ssh snakemake-node (there are no interactive jobs needed for snakemake itself)
conda activate smk8
cd folder/of/runPipeline.sh
bash runPipeline.sh
```

Detach: `Ctrl + A`, then `D`. Reattach: `screen -r pipeline_run`.

### Execution Flags & Variants

- Skip wrapper-level checksum generation (faster):
  ```bash
  bash runPipeline.sh --no-checksums
  ```
  Caveat: this only skips the wrapper-level `checksums_<timestamp>.sha256` and `filelist_project_<timestamp>.filelist`. The demux stage itself still runs rule `create_checksums` (per-project `untrimmed_fastq/<projectNum>_sha256sums_fastqfiles.sha256`); that rule is unconditional in `convert2fastq.smk`.
- Local workstation:
  ```bash
  bash runPipeline_local.sh
  ```
  (runs `snakemake --cores all --sdm conda --jobs 7/4` instead of the PBS profile, and uses `nice`.)
- Barebones (direct Snakemake calls, no reports/checksums):
  ```bash
  bash runPipeline_barebone.sh
  ```

### What runPipeline.sh Actually Does

1. Reads the output folder from `config.yaml` with `grep OutputFolder config.yaml -A1 | tail -n 1 | awk '{print $1}'`, i.e. it takes the line **after** the `OutputFolder` match. The path must therefore sit on the line following the key; verify against your `config.yaml` (a wrong value makes `mkdir -p $out` and all wrapper outputs land in the wrong place).
2. Prints which `config.yaml` lines contain `True` / `False` (a plain `grep`, so any line with those words is listed).
3. If `$out/*config.yaml` exists (e.g. `<projectNum>_config.yaml` from an earlier run), moves old `*_snakemake_report*.html`, `*_rulegraph*.pdf`, `*.sha256` and `*.filelist` from the output root into `<out>/logs/previous_executions/`. Config files are the trigger, they are **not** moved.
4. Writes per-stage rulegraph PDFs: `convert2fastq_rulegraph.<timestamp>.pdf` and `fastqProcessing_rulegraph.<timestamp>.pdf` (via `snakemake --forceall --rulegraph | dot -Tpdf`).
5. Runs `convert2fastq.smk` and then `fastqProcessing.smk` with `--profile pbs`, each followed by `--report` -> `convert2fastq_snakemake_report.<timestamp>.html` / `fastqProcessing_snakemake_report.<timestamp>.html`.
6. Unless `--no-checksums`: `filelist_project_<timestamp>.filelist` and `checksums_<timestamp>.sha256` over the whole output folder.

The timestamp format is `%Y_%m_%d_%I_%M_%p`; there are no plain `filelist.txt` / `checksums.sha256` / `dag.pdf` files.

**The wrapper does not stop on errors** (no `set -e`, no exit-code checks). If the demux stage fails, the processing stage, the reports and the checksums still run, and the script ends with "Completed the run." That message is not a success indicator; check the Snakemake logs.

Note: `copy_software_env` (rule, demux stage) copies the versions file to `<projectNum>_software_environment_mqc_versions.yml` and `config.yaml` to `<projectNum>_config.yaml` in the output folder.

---

## 5. Pipeline Stages & Architecture

### Stage 1: Demultiplexing (`rules/convert2fastq.smk`)

```text
Illumina RunFolder (BCLs) + SampleSheet.csv
                  |
                  v
 check_inputs (bcl-convert --bcl-validate-sample-sheet-only true
               + scripts/bclconvert_samplesheet_check.py)
               [only requested when bcl-convert is used and skip_demux is off]
                  |
                  v
 bcl-convert / bcl2fastq2  ->  *.fastq.gz -> moved to untrimmed_fastq/
                               (bcl-convert: Undetermined*.fastq.gz moved to output root;
                                the SampleSheet is also copied into the output folder)
                  |
                  v
 create_fastq_list  ->  fastq_infiles_list.tx (excludes Undetermined)
                  |
                  v
 create_checksums   ->  untrimmed_fastq/<projectNum>_sha256sums_fastqfiles.sha256
                  |
                  v
 copy_software_env  ->  <projectNum>_software_environment_mqc_versions.yml
                        <projectNum>_config.yaml
                  |
                  v
 interop_plots (if interop_plots_active)  ->  interop_plots/ PDFs + done flag
```

Optional around demux: `prevent_revcomp` / `reverse_prevent_revcomp` (RunInfo.xml patch, off by default; `prevent_revcomp` has no input ordering it before demux, it only wins the race because its resources are tiny), `skip_demux` (moves pre-demultiplexed `*.gz` files from `skip_demux.fastq_folder` into `untrimmed_fastq/`).

### Stage 2: Processing, Contamination & RNA-Seq (`rules/fastqProcessing.smk`)

```text
untrimmed_fastq/
      |
      v
 FastQC + fastp (raw)  ->  fastqc_untrimmed/  ->  multiqc_raw
      |
      v
 cutadapt           ->  trimmed/{file}_trimmed.fastq.gz
      |                 (the only rule that consumes the check_fastq_files flag)
      v
 umi_tools extract  ->  umi_extract/{file}.umis-extracted.fastq.gz     (if UMIs)
      |
      v
 sortmerna          ->  sortmerna/{file}_non-ribosomal_rna.fq.gz      (if active)
      |                 (+ {file}_ribosomal_rna.fq.gz)
      v
 STAR align         ->  star/{file}_Aligned.sortedByCoord.out.bam
      |                 star/{file}_ReadsPerGene.out.tab
      |
      +-> umi_dedup (if UMIs) -> star/{file}_deduped.Aligned.sortedByCoord.out.bam
      +-> samtools / preseq / Qualimap or RSeQC / rustqc
      +-> StringTie (3 stages) + gffcompare, rMATS
      +-> FeatureCounts -> counts/all.tsv -> pca_plot (PCA/heatmap/reports)
      +-> star_result_table
      |
      v
 multiqc_report_complete_<projectNum>.html

Read-level contamination / k-mer rules run on FASTQs, independent of alignment:
 kraken2, kaiju, sourmash : trimmed FASTQ if cutadapt_active, else untrimmed FASTQ
 blast, jellyfish, biobloom: umi_extract FASTQ if UMIs on, else trimmed FASTQ
 (none of them reads the sortmerna output; preseq is BAM-based, not FASTQ-based)
```

`blast`, `jellyfish`, `biobloom`, `sortmerna`, `FastQC` (trimmed) and `crypt4gh` always request `trimmed/...` (or `umi_extract/...`). The code comments say they require cutadapt and/or umi_tools to be active; with both off they still request `trimmed/`, which makes Snakemake run the `cutadapt` rule anyway.

**BAM-input switch**: when `umi_tools.umi_tools_active: True`, the following rules consume `star/{file}_deduped.Aligned.sortedByCoord.out.bam` instead of `star/{file}_Aligned.sortedByCoord.out.bam`: `samtools`, `preseq`, `stringtie`, `stringtie_estimate`, `Qualimap`, `featurecounts`, `star_result_table`, `rmats_create_groupfiles_from_csv`, `rustqc`, and the RSeQC rules in `qc.smk` (`rseqc_main`, `rseqc_minor`, `rseqc_done`). `align` itself takes FASTQs, and `pca_plot` takes `counts/all.tsv`, so neither switches inputs. Enabling UMIs therefore changes the BAM input of 12 rules (plus `umi_dedup` producing the deduped BAMs).

**BAM lifetime**: in the SE `align` rule the BAM is declared `temp()`, so Snakemake deletes `star/{file}_Aligned.sortedByCoord.out.bam` once all consuming rules have finished. With UMIs off no BAM remains in `star/` afterwards; with UMIs on, only the (non-temp) deduped BAM stays. `align_pe` does **not** declare its BAM `temp()`. `umi_dedup` additionally leaves `star/{file}_sorted.Aligned.sortedByCoord.out.bam` behind.

---

## 6. The check_fastq_files Rule & check_input_fastqs.sh

```python
checkpoint check_fastq_files:
    input:
        checksum_file=config["demux"]["OutputFolder"] + "/untrimmed_fastq/"
                      + projectNum + "_sha256sums_fastqfiles.sha256"
    output:
        outputfolder + "/inputfiles_checkpoint_completed.out"
    params:
        dir_to_check = outputfolder + "/untrimmed_fastq/"
    shell:
        "bash scripts/check_input_fastqs.sh {params.dir_to_check} >>{output}"
```

- Input: per-project checksum file (from `create_checksums`), not the FASTQ list directly.
- Output: `inputfiles_checkpoint_completed.out`.
- Action: runs `scripts/check_input_fastqs.sh` on `untrimmed_fastq/`.
- It is declared with `checkpoint`, but nothing in the workflow calls `checkpoints.check_fastq_files.get(...)`, so it behaves like an ordinary rule (no DAG re-evaluation).
- **Who depends on it**: only `cutadapt` and `cutadapt_pe` (input `nonempty_fastqs=`). `getOutput()` no longer lists the flag (the line is commented out). The sanity check therefore only runs if cutadapt runs; with cutadapt not needed, it is skipped.

### check_input_fastqs.sh validation steps

1. Directory argument provided and exists.
2. Directory non-empty.
3. At most ONE zero-byte `*.fastq.gz` is tolerated (emits a warning); more than one fails the check.
4. Every `.fastq.gz` must be >= 80 bytes (filters truncated/corrupted output).
5. On success the message is appended to `inputfiles_checkpoint_completed.out`.

---

## 7. Paired-End (PE) vs. Single-End (SE) Layout Rules

- Layout detection: `isSingleEnd()`, evaluated at parse time; PE/SE blocks gate `cutadapt`, `align`, `crypt4gh`, `FastQC`, `FastQC_untrimmed`, `sortmerna`, and `get_mapping_input` in `fastqProcessing.smk`, and the whole PE rule set in `pe_processing.smk`. Note: if even one R2 file is missing, the entire set is treated as single-end.
- The detector classifies a name by the text between the **first** `_R` and the next `_R` (`name.split("_R")[1]`). Consequences:
  - Sample names must not contain `_R` before the read tag (`S1_Rep1_R1_001`, `Lib_RNA_R1_001`): all such files are classified as R2, and the whole run is silently treated as single-end.
  - A FASTQ name without any `_R` raises an `IndexError` at parse time in `pe_processing.smk` (its `isSingleEnd()` has no guard).
  - With an empty `fastq_infiles_list.tx` the counts are 0 = 0, i.e. "PE".
- Sample name derivation:
  - `sample_names` = FASTQ names from `fastq_infiles_list.tx` minus `.fastq.gz`, excluding `_I1_`/`_I2_` index reads.
  - `getSample_names_post_mapping()` returns `sample_names` for SE, and for PE replaces `_R1` with `_pe` (so `X_R1_001` becomes `X_pe_001`). PE mapping outputs are `star/<sample>_pe_001_...`, not `<sample>_R1_...`.
- Naming: PE inputs must be `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`. Index reads (I1/I2) are ignored after demux.
- cutadapt: SE trims one file per job; PE (`cutadapt_pe`) trims R1+R2 simultaneously, applying all adapter pairs from `config["cutadapt"]["adapters"]` with types from `adapter_type` (default `b,b,b,b`), plus `other_params` (default `--minimum-length 25 --quality-cutoff 5 -n 2`). Output: `trimmed/{file}_trimmed.fastq.gz` (SE) / `trimmed/{short}_R1_001_trimmed.fastq.gz` + R2 (PE).
- umi_tools: `pattern` for R1. The PE rule (`umi_extract_pe`) always passes `--read2-in`/`--read2-out` **and** `--bc-pattern2={pattern2}`; there is no branch for "R2 without UMI", so `pattern2` has to be set for PE. Outputs: `umi_extract/{file}.umis-extracted.fastq.gz` (SE) / `{short}_R1.umis-extracted.fastq.gz` + R2 (PE). See section 16 for the inconsistencies in the PE + UMI path.
- STAR: `--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile <gtf> <mapping.other_params> --readFilesCommand zcat`; outputs `star/{file}_Aligned.sortedByCoord.out.bam` and `star/{file}_ReadsPerGene.out.tab` (SE, BAM is `temp()`) or `star/{short}_pe_001_Aligned...` (PE).
- FeatureCounts: `-g gene_id` + `config["mapping"]["stranded"]` (e.g. `-s 0` unstranded, `-s 1` forward, `-s 2` reverse). That string is the **only** place extra featureCounts options can go: the rule does NOT add `-p` automatically for PE, and `mapping.other_params` is passed to STAR, not to featureCounts. For PE put the flags into `mapping.stranded`, e.g. `stranded: "-p --countReadPairs -s 2"` (`--countReadPairs` is required with Subread >= 2.0.2 to count fragments; check the version in `envs/counts.yaml`).
- rMATS: `mapping.rmats_paired` is passed as `-t`. Example value: `"single --libType fr-secondstrand"`. Snakemake does not quote it, so it ends up as separate shell words; keep it as one YAML string and verify with a dry run.

---

## 8. Automation & Autostart Cronjob Setup

`/autostart` enables autonomous execution for core-facility sequencers via `execute_pipeline.sh` plus one of three cron wrappers:

- `autostart_cronjob_bclconvert.sh`: hardcoded for the NextSeq 2000 watch directory.
- `autostart_novaseq.sh`: same logic, hardcoded for a NovaSeq X watch directory.
- `autostart_portable.sh`: parameterized version (`--watch-dir`, `--pipeline`, `--snakemake-dir`, `--execution-script`, `--env-script`, `--session-prefix`); requires an environment script that sets PATH/modules/conda.

1. **Manual config required**: autostart does not construct config.yaml. An operator must place a populated `config.yaml` in each run folder (the wrapper copies it to the pipeline root before launching).
2. **Directory monitoring**: the cron wrapper scans the instrument storage for run folders that:
   - contain `CopyComplete.txt` (no RTAComplete.txt check is performed),
   - contain exactly `config.yaml` (not `*_config.yaml`),
   - do not contain the lockfile `.pipeline_launched`,
   - while no job named with prefix `demx_pipe_auto` (PBS job name from `execute_pipeline.sh`) is in the `qstat` queue.
3. **Launch**: creates `.pipeline_launched` (with launch date + session), copies the run's `config.yaml` to the pipeline root, waits ~30 s, and submits `execute_pipeline.sh` via `qsub -q default`.

Caveats:

- `execute_pipeline.sh` activates conda env `smk9` (not `smk8`) and hardcodes the pipeline directory; adjust both for your site.
- `autostart_cronjob_bclconvert.sh` and `autostart_novaseq.sh` use `wait 30`, which in bash only waits for background jobs of the current shell; it does NOT sleep. `autostart_portable.sh` correctly uses `sleep 30`.
- The wrappers also hardcode the invoking user's PATH to make `qstat` work from a minimal cron environment.

---

## 9. Complete Tool Reference & Setup Instructions

| Tool | Config Key | Purpose | Docs | Source |
| --- | --- | --- | --- | --- |
| bcl-convert | demux.bcl_convert_path | Primary demux for NovaSeq X / NextSeq 2000 | [Illumina](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) | Proprietary |
| bcl2fastq2 | demux.bcl2fastq2_path | Legacy demux for MiSeq / older flowcells | [Validation guide](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2_val_guide_15051736_02.pdf) | Proprietary |
| FastQC | (untrimmed/trimmed QC) | Per-base quality metrics | [Babraham](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | [GitHub](https://github.com/s-andrews/FastQC) |
| Fastp | - | Fast quality filtering, HTML/JSON reports (run alongside FastQC) | [GitHub](https://github.com/OpenGene/fastp) | [GitHub](https://github.com/OpenGene/fastp) |
| Cutadapt | cutadapt | Adapter trimming, length/quality filters | [Docs](https://cutadapt.readthedocs.io/) | [GitHub](https://github.com/marcelm/cutadapt) |
| STAR | mapping | Splice-aware RNA-Seq alignment (--quantMode GeneCounts) | [Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) | [GitHub](https://github.com/alexdobin/STAR) |
| Qualimap | mapping.qualimap_on, mapping.qualimap_dir | Post-alignment BAM QC | [Qualimap](http://qualimap.bioinfo.cipf.es/) | [Bitbucket](https://bitbucket.org/kokonech/qualimap) |
| RSeQC | mapping.rseqc_active | Alternative post-alignment QC | [Docs](https://rseqc.sourceforge.net/) | [GitHub (Monash mirror)](https://github.com/MonashBioinformaticsPlatform/RSeQC) |
| FeatureCounts (Subread) | mapping | Assign BAM reads to GTF features (counts/all.tsv) | [User guide](https://subread.sourceforge.net/SubreadUsersGuide.pdf) | SourceForge |
| StringTie | mapping.stringtie_on | Assemble/quantify transcripts (3 stages) | [Docs](https://ccb.jhu.edu/software/stringtie/) | [GitHub](https://github.com/gpertea/stringtie) |
| gffcompare | - | Compare StringTie GTF to reference | [Docs](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) | [GitHub](https://github.com/gpertea/gffcompare) |
| rMATS | mapping.rmats_on | Differential alternative splicing | [GitHub](https://github.com/XingLab/rmats-turbo) | [GitHub](https://github.com/XingLab/rmats-turbo) |
| Preseq | preseq | Library complexity / saturation | [GitHub](https://github.com/smithlabcode/preseq) | [GitHub](https://github.com/smithlabcode/preseq) |
| bedtools | - | bamtobed for Preseq | [Docs](https://bedtools.readthedocs.io/) | [GitHub](https://github.com/arq5x/bedtools2) |
| Kraken2 | kraken2 | k-mer taxonomic classification | [Docs](https://ccb.jhu.edu/software/kraken2/) | [GitHub](https://github.com/DerrickWood/kraken2) |
| Sourmash | sourmash | MinHash signature comparison | [Docs](https://sourmash.readthedocs.io/) | [GitHub](https://github.com/sourmash-bio/sourmash) |
| BioBloom Tools | biobloom | Bloom-filter organism filtering | [GitHub](https://github.com/bcgsc/biobloom) | [GitHub](https://github.com/bcgsc/biobloom) |
| BLAST+ | blast | blastn on subsampled reads | [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK279690/) | NCBI |
| Diamond | - | Protein-level alignment; only a commented-out `diamond_pe` rule exists, no active rule and no env YAML | [GitHub](https://github.com/bbuchfink/diamond) | [GitHub](https://github.com/bbuchfink/diamond) |
| Kaiju | kaiju | Protein-level metagenomic classification | [GitHub](https://github.com/bioinformatics-centre/kaiju) | [GitHub](https://github.com/bioinformatics-centre/kaiju) |
| Jellyfish | jellyfish | 21-mer frequency distributions | [Docs](https://www.cbcb.umd.edu/software/jellyfish/) | [GitHub](https://github.com/gmarcais/Jellyfish) |
| SortMeRNA | sortmerna | rRNA removal | [Docs](https://sortmerna.readthedocs.io/) | [GitHub](https://github.com/biocore/sortmerna) |
| UMI-tools | umi_tools | UMI extraction + BAM dedup | [Docs](https://umi-tools.readthedocs.io/) | [GitHub](https://github.com/CGATOxford/UMI-tools) |
| RustQC | rustqc | Alternative RNA QC (needs sambamba markdup) | [Docs](https://seqeralabs.github.io/RustQC/) | [GitHub](https://github.com/seqeralabs/RustQC) |
| sambamba | - | markdup pre-step for rustqc | [Docs](https://lomereiter.github.io/sambamba/) | [GitHub](https://github.com/biod/sambamba) |
| Crypt4gh | crypt4gh | GA4GH per-file public-key encryption | [Docs](https://crypt4gh.readthedocs.io/) | [GitHub](https://github.com/ga4gh/crypt4gh) |
| Gocryptfs | gocryptfs | Encrypted overlay directory | [Docs](https://nuetzlich.net/gocryptfs/) | [GitHub](https://github.com/rfjakob/gocryptfs) |

---

## 10. Rule Registry (Snakefiles & .smk includes)

### rules/convert2fastq.smk (demux stage)

`localrules: all, create_fastq_list`

| Rule | Purpose / Output |
| --- | --- |
| all | collects `getOutput()`: fastq list, checksums, software env copy + config copy, interop flag (if active), revcomp flags (if active), `skipped_demuxing.flag` (if skip_demux), `SampleSheet_check.out` (bcl-convert, no skip_demux), `used_bcl2fastq2.flag` (if use_bcl2fastq) |
| prevent_revcomp | patches RunInfo.xml via `scripts/prevent_revcomp_nexts2.sh` -> `revcomp_prevented.flag` |
| reverse_prevent_revcomp | reverts the patch -> `reversed_revcomp_prevented.flag` |
| check_inputs | bcl-convert samplesheet validation + `scripts/bclconvert_samplesheet_check.py` -> `SampleSheet_check.out` (defined unless skip_demux; requested by `all` only when bcl-convert is used) |
| demux | bcl-convert run -> `Reports/Adapter_Metrics.csv`, FASTQs moved to `untrimmed_fastq/`, Undetermined moved to output root, bcl logs to `logs/bcl_convert_logs`, SampleSheet copied into the output folder |
| bcl2fastq2 | legacy demux (only defined if `use_bcl2fastq`) -> `used_bcl2fastq2.flag`, FASTQs to `untrimmed_fastq/` |
| skip_demux | (only defined if `skip_demux_active`) moves all `*.gz` from `skip_demux.fastq_folder` -> `skipped_demuxing.flag` |
| create_fastq_list | `fastq_infiles_list.tx` (recursive glob, excludes Undetermined) |
| create_checksums | `untrimmed_fastq/<projectNum>_sha256sums_fastqfiles.sha256` (parallel `xargs -P {resources.threads}`) |
| copy_software_env | versions TSV -> `<projectNum>_software_environment_mqc_versions.yml`, config -> `<projectNum>_config.yaml` |
| interop_plots | InterOp plots + reports (`Sequencing_report.pdf`, `RunReport.pdf`, `Demux_report.pdf`, `TiltReport.pdf`) -> `interop_plots/interop_plots_done.flag`. Output of `interop_plots.sh` and `test.py` is sent to `/dev/null`, so their failure does not fail the rule |

### rules/fastqProcessing.smk (processing stage)

`localrules: all, gocryptfs`

| Rule | Gate | Purpose / Output |
| --- | --- | --- |
| all | - | collects `getOutput()`. With gocryptfs active it requests only `encrypted/gocryptfs.conf`; otherwise the untrimmed FASTQs and `multiqc_report_<projectNum>.html`. `multiqc_report_complete_<projectNum>.html` is always appended. If `fastq_infiles_list.tx` is missing the list is reset to that single target (section 15) |
| check_fastq_files | checkpoint (used as a plain rule) | `inputfiles_checkpoint_completed.out`; consumed only by cutadapt |
| cutadapt | SE | `trimmed/{file}_trimmed.fastq.gz` |
| align | SE | STAR -> `star/{file}_Aligned.sortedByCoord.out.bam` (`temp()`) + `star/{file}_ReadsPerGene.out.tab` |
| crypt4gh | SE | `encrypted_fastqs/{file}_processed.fastq.gz.c4gh` |
| FastQC_untrimmed | SE | FastQC + fastp on raw FASTQs -> `fastqc_untrimmed/` |
| FastQC | SE | FastQC + fastp on trimmed/UMI-extracted FASTQs -> `fastqc/` |
| sortmerna | SE | rRNA split -> `sortmerna/{file}_[non-]ribosomal_rna.fq.gz` |
| umi_extract | umi_tools on | `umi_extract/{file}.umis-extracted.fastq.gz` (not layout-gated) |
| umi_dedup | umi_tools on | `star/{file}_deduped.Aligned.sortedByCoord.out.bam` (not layout-gated; leaves `star/{file}_sorted.Aligned...bam` behind) |
| samtools | - | index/stats/coverage/plot-bamstats -> `samtools/{file}_samtools_stats.stats` |
| preseq | - | c_curve on BAM-derived BED -> `preseq/{file}_c_curve_out.txt` |
| Qualimap | - | `qualimap/{file}/qualimapReport.html` (hardcoded `-p strand-specific-forward`, `--java-mem-size=16G`) |
| featurecounts | - | `counts/all.tsv` |
| pca_plot | - | pca.R + heatmap.py + counts_report.py -> declared output `counts/pca.png`; also writes `counts/<projectNum>_CountsReport.pdf` and `counts/<projectNum>_Counts_Report.html` (not declared as outputs) |
| star_result_table | - | summarize_star_out.py -> `star/star_logs_summarized.csv` |
| stringtie | - | `stringtie/transcripts_{file}.gtf` (+ `abundances_{file}.tab`) |
| stringtie_merge | - | `stringtie/transcripts_all_merged.gtf` |
| stringtie_estimate | - | `stringtie/transcript_abundances_{file}.gtf` |
| gffcompare | - | `gffcompare/stringtie_to_ref_comp_{file}.stats` |
| rmats_create_groupfiles_from_csv | rmats on | generate_rmats_groups.sh -> `rmats/group[1,2]_bams.txt` |
| run_rmats | rmats on | rmats.py -> `rmats/summary.txt` |
| kraken2 | only defined if kraken2 on | `kraken2/{file}.report` |
| kaiju | - | kaiju + kaiju2table -> `kaiju/kaiju_species_table_{file}.txt` |
| blast | - | subsample + blastn + filter_blast.sh -> `blast/blast_report_{file}.tsv` (+ `blast_summary_{file}.tsv`) |
| jellyfish | - | 21-mer histogram -> `jellyfish/{file}_trimmed_jf.hist` |
| biobloom | - | biobloomcategorizer -> `biobloom/biobloom_results_{file}_summary.tsv` |
| sourmash | - | sketch + gather -> `sourmash/{file}_sourmash.csv` (`sourmash gather` errors are ignored, the command ends in `\|\| true`) |
| rustqc | - | sambamba markdup + rustqc rna -> `rustqc/featurecounts/{file}.biotype_counts.tsv` |
| transfer_files | transfer on | rsync -> `transfer/transfer_log.log` (destination = `transfer_as` + `transfer_to`) |
| multiqc_raw | - | MultiQC on `fastqc_untrimmed/` -> `multiqc_report_<projectNum>.html` |
| multiqc | - | MultiQC on all results -> `multiqc_report_complete_<projectNum>.html` |
| gocryptfs | gocryptfs on | interactive gocryptfs -> `encrypted/gocryptfs.conf` |

### rules/pe_processing.smk (PE only, `isSingleEnd() == False`)

`cutadapt_pe`, `align_pe`, `umi_dedup_pe`, `FastQC_pe`, `FastQC_untrimmed_pe`, `crypt4gh_pe`, `sortmerna_pe`, `umi_extract_pe` (only if UMIs on), and `count_matrix_pe` (only if UMIs off). `count_matrix_pe` is a ReadsPerGene-based count-matrix builder (`scripts/count-matrix.py`) targeting the same `counts/all.tsv` as `featurecounts`; its input path does not exist in the flat `star/` layout, so `featurecounts` is what actually produces the matrix (section 16). `diamond_pe` exists only commented out.

### rules/qc.smk (RSeQC, `mapping.rseqc_active`)

`rseqc_gtf2bed` (scripts/gtf2bed.py -> `rseqc/annotation.bed`), `rseqc_main` (junction_annotation, bam_stat, junction_saturation, infer_experiment, read_NVC, FPKM_count, inner_distance), `rseqc_minor` (read_distribution, read_duplication, read_GC; geneBody_coverage is commented out), `rseqc_done` (`rseqc/done.flag`).

### rules/common.smk (helpers, no rules)

`getFiles()`, `isSingleEnd()`, `getSample_names_post_mapping()`, `getSamples()`, `getFastqs()`; reads `fastq_infiles_list.tx`, defines `outputfolder` and `projectNum`. `isSingleEnd()` / `getSample_names_post_mapping()` / the `sample_names` block are also defined in `pe_processing.smk` (and the latter two again in `convert2fastq.smk`), so there are several copies to keep in sync (section 16). `getSamples()` / `getFastqs()` / `getFiles()` are not used by any active rule.

---

## 11. Scripts Directory Guide (scripts/)

Scripts actually present and referenced by the current rules:

- `check_input_fastqs.sh`: validates `untrimmed_fastq/` (non-empty, max one zero-byte file, all >= 80 bytes). Called by `check_fastq_files`.
- `bclconvert_samplesheet_check.py`: parses/validates SampleSheet.csv before bcl-convert; called by `check_inputs`.
- `prevent_revcomp_nexts2.sh` / `reverse_revcomp_nexts2.sh`: patch/revert RunInfo.xml to prevent unexpected index-2 reverse complementation. Called by `prevent_revcomp` / `reverse_prevent_revcomp`.
- `interop_plots.sh`: InterOp plot wrapper; called by `interop_plots`.
- `test.py`: produces `interop_plots/Sequencing_report.pdf` (generically named; verify intent).
- `cluster_pdf.py`: produces `interop_plots/RunReport.pdf`.
- `demux_report.py`: produces `interop_plots/Demux_report.pdf`.
- `tilt_report.py`: produces `interop_plots/TiltReport.pdf`.
- `gtf2bed.py`: GTF to BED for RSeQC; called by `rseqc_gtf2bed`.
- `generate_rmats_groups.sh`: builds `rmats/group1_bams.txt` and `rmats/group2_bams.txt` from `mapping.rmats_group_file` (called directly, so it needs the executable bit).
- `pca.R`: PCA on `counts/all.tsv`; writes `counts/pca.png` (called by `pca_plot`).
- `heatmap.py`: writes `counts/<projectNum>_CountsReport.pdf` (called by `pca_plot`).
- `counts_report.py`: writes `counts/<projectNum>_Counts_Report.html` (called by `pca_plot`).
- `summarize_star_out.py`: writes `star/star_logs_summarized.csv` (called by `star_result_table`).
- `fastq_to_fasta_subsampling.sh` / `filter_blast.sh`: subsample FASTQ -> FASTA and summarize blastn output for the `blast` rule. Note: the rule invokes them from the configured `blast.blast_script_dir`, NOT from the repo's `scripts/` folder; keep copies in sync.
- `get_c4gh_string.sh`: builds the `crypt4gh encrypt --recipient_pk ...` string from `crypt4gh_client_public_key_dir`; called by `crypt4gh`/`crypt4gh_pe`.
- `helper.py`: `validateSamplesheet`, `validateOutput`, `validateProjectNum`, `adaptersToStringParams`.
- `count-matrix.py`: builds a count matrix from STAR ReadsPerGene tabs; referenced only by the PE-only `count_matrix_pe` rule (and a commented-out SE `count_matrix`).

Present in the repo but not referenced by any current rule: `interop_minimal.py`, `interop_more.py`, `interop_more_old.py`, `clusters.py`, `get_cluster_numbers.py`, `fastq_to_fasta.sh` (only named in an unused `blast` param), `jobscript.sh`, `tubemap.py`, `versions.py`. Older scripts (`interop_plots.py`, `clusters_pdf.py`, `prevent_revcomp.py`) do not exist in the repo.

---

## 12. Checksums, Cryptography & MultiQC Reports

### Checksums & File Lists

1. Raw demux checksums: `create_checksums` writes `untrimmed_fastq/<projectNum>_sha256sums_fastqfiles.sha256` (sha256sum over all non-Undetermined FASTQs, parallelized with `xargs -P`).
2. Trimmed checksums: `cutadapt`/`cutadapt_pe` appends to `trimmed/sha256sums_fastqfiles_trimmed.sha256`.
3. UMI-extracted checksums: `umi_extract`/`umi_extract_pe` appends to `umi_extract/checksums_umis_extracted_<projectNum>.sha256`.
4. Crypt4gh checksums: `crypt4gh`/`crypt4gh_pe` appends to `encrypted_fastqs/checksums_encrypted_fastqs.sha256`: one aggregate file, not per-file `.c4gh.sha256`.
5. Wrapper-level: `filelist_project_<timestamp>.filelist` and `checksums_<timestamp>.sha256` from `runPipeline.sh` (unless `--no-checksums`).

### MultiQC Reports

1. `multiqc_report_<projectNum>.html`: rule `multiqc_raw`; runs inside `fastqc_untrimmed/` with `-m fastqc`, `--fullnames`, `--ignore-samples Undetermined*`, `-x Undetermined*`, `--no-data-dir`; also writes `logs/fastq_file_list.tsv`.
2. `multiqc_report_complete_<projectNum>.html`: rule `multiqc`; aggregates `getFastQCs()` (raw+trimmed FastQC/Fastp, STAR, samtools, FeatureCounts, contamination, Qualimap, StringTie, gffcompare, rMATS, rustqc, etc.) with `--fullnames` and `--no-data-dir`, then deletes empty `.log` files (list in `logs/deleted_empty_files.out`).

---

## 13. Configuration Reference (config.yaml)

Only two things are mandatory per run: a samplesheet with a project number in the filename, and an output folder. Trailing slashes: required only on `biobloom_filter_dir`, `transfer.transfer_to` and `crypt4gh_client_public_key_dir` (the last one is consumed by `scripts/get_c4gh_string.sh`); all other paths must NOT end with `/`.

| Section | Key | Notes |
| --- | --- | --- |
| demux | SampleSheet | full path; must be inside the sequencing run folder; project number must be part of the samplesheet filename (e.g. 1399_Samplesheet.csv) |
| demux | OutputFolder | complete dir path, no trailing /; `runPipeline.sh` reads it from the line after the key (section 4) |
| demux | options | extra bcl-convert CLI parameters (e.g. `--no-lane-splitting true`) |
| demux | bcl_convert_path | path to the bcl-convert executable, no trailing / |
| demux | demux_try1_gb | RAM for first demux attempt (96; 126+ for NovaSeq X) |
| demux | use_bcl2fastq | True for MiSeq / if bcl-convert errors |
| demux | bcl2fastq2_path, bcl2fastq_mismatches, bcl2fastq_options | bcl2fastq2 settings |
| skip_demux | skip_demux_active, fastq_folder | use pre-demultiplexed FASTQs (all `*.gz` in the folder are **moved**); no trailing / |
| prevent_revcomp | prevent_revcomp_active | patch RunInfo.xml against read-3 revcomp; keep inactive unless you know better |
| interop_plots | interop_plots_active | adds SAV-like plots + shareable report PDFs (needs the Illumina run folder) |
| cutadapt | cutadapt_active, adapters, adapter_type, other_params | one adapter_type per adapter (`a`/`g`/`b`) |
| mapping | mapping_active, gtf_file, genomeDir, other_params, stranded | STAR settings; genomeDir no trailing /; `other_params` goes to STAR only; `stranded` is the string appended to the featureCounts call (`-s` value, and `-p ...` for PE) |
| mapping | stringtie_on, rseqc_active, qualimap_on, qualimap_parameters, qualimap_dir | post-mapping QC; enable rseqc OR qualimap, not both |
| mapping | rmats_on, rmats_read_length, rmats_paired, rmats_group_file | rMATS settings (read length = sequenced length minus UMI/spacer/adapter) |
| rustqc | rustqc_active | experimental section; needs sambamba |
| preseq | preseq_active | sequencing saturation estimate |
| umi_tools | umi_tools_active, pattern, pattern2 | regex patterns need single AND double quotes in YAML; `pattern2` is always used for PE |
| blast | blast_active, blast_db, blast_script_dir | quick species proxy; script dir must contain the subsampling/filter scripts |
| kraken2 | kraken2_active, kraken2_database | shared DB; jobs are slow to start (~190 GB RAM) |
| kaiju | kaiju_active, kaiju_db, kaiju_nodes_dmp, kaiju_names_dmp | metagenomic fallback |
| jellyfish | jellyfish_active | 21-mers; for whole-genome assays |
| biobloom | biobloom_active, biobloom_filter_dir, biobloom_filters | one filter per species; filter dir needs trailing / |
| sourmash | sourmash_active, sourmash_eukaryote_ref, sourmash_microbial_ref | k-mer classification |
| sortmerna | sortmerna_active, sortmerna_reference_list | rRNA counting/removal; needs cutadapt and/or umi_tools |
| crypt4gh | crypt4gh_active, crypt4gh_own_private_key, crypt4gh_client_public_key_dir | human data; public-key dir needs trailing / |
| transfer | transfer_active, transfer_to, transfer_as | login-less connection required; transfer_to needs trailing /; target is `transfer_as` + `transfer_to` |
| gocryptfs | gocryptfs_active, gocryptfs_version | only when absolutely needed, in a screen session (password typed manually) |
| others | (empty placeholder section) | |
| software_versions_file | | TSV of tool versions; update manually (see envs/update_software_env_list.sh) |

---

## 14. rMATS Group File Setup

When `mapping.rmats_on: True`, set `mapping.rmats_group_file` to a group file. The rule `rmats_create_groupfiles_from_csv` calls `scripts/generate_rmats_groups.sh` to produce `rmats/group1_bams.txt` (group A) and `rmats/group2_bams.txt` (group B) from the STAR BAM directory. `run_rmats` then consumes those two lists. The rule also waits for the samtools stats of all samples, so no samtools temp files are picked up.

Example (`examples/example_rmats_group_file.txt`) for 7 samples in two comparison groups (A and B):

```text
*S1*.fastq.gz,*S2*.fastq.gz,*S3*.fastq.gz:*S4*.fastq.gz,*S5*.fastq.gz,*S6*.fastq.gz,*S7*.fastq.gz
```

Patterns must match the demultiplexed (or skip_demux) sample filenames.

---

## 15. Troubleshooting & Common Pitfalls

### Trailing Slash Errors

- Issue: rule failure or missing file errors.
- Fix: remove trailing slashes (`/`) from `OutputFolder`, `fastq_folder`, `genomeDir`, `qualimap_dir`, `bcl_convert_path`, `bcl2fastq2_path`. Only `biobloom_filter_dir`, `transfer_to` and `crypt4gh_client_public_key_dir` need trailing slashes.

### SampleSheet Naming

- Issue: pre-flight check fails or demux target cannot be resolved.
- Fix: ensure SampleSheet.csv resides inside the sequencing run directory and contains a numeric Project ID in its filename (e.g. `1399_Samplesheet.csv`). `validateProjectNum()` is called; `validateSamplesheet`/`validateOutput` are imported but only called by `pe_processing.smk` (`validateBefore`), so the demux stage relies on `validateProjectNum` alone.

### Near-silent "success" with a missing fastq list

- Issue: if `fastq_infiles_list.tx` does not exist, `getOutput()` resets the target list and prints "It seems like the demuxing run didnt finish properly...". `rule all` then reduces to little more than the `multiqc_report_complete_<projectNum>.html` target, so Snakemake can "succeed" while producing almost nothing.
- Fix: verify `untrimmed_fastq/` and `fastq_infiles_list.tx` exist before trusting a green run; check the log for the demuxing message above.

### "Completed the run." does not mean success

- Issue: `runPipeline.sh` ignores exit codes of the individual `snakemake` calls and always continues to the next stage, the reports and the checksums.
- Fix: check the Snakemake log of each stage (or add `set -e` / exit-code checks to the wrapper).

### Wrong layout detected (PE treated as SE)

- Issue: all reads are processed as single-end, or parsing fails with `IndexError: list index out of range`.
- Cause: layout detection splits names on `_R` (section 7). Names such as `S1_Rep1_R1_001` are classified as all-R2 -> single-end; names without `_R` raise an `IndexError` in `pe_processing.smk`.
- Fix: rename samples so that the first `_R` in the file name is the read tag (`_R1_001` / `_R2_001`), i.e. no `_R...` earlier in the name: `S1_Rep1_R1_001` -> `S1_rep1_R1_001`.

### Locked Snakemake Directory

- Issue: `Error: Directory cannot be locked.`
- Fix: run `snakemake -s rules/convert2fastq.smk --unlock` and/or `snakemake -s rules/fastqProcessing.smk --unlock` from the pipeline root (the wrappers do not pass `--unlock` through).

### Memory Starvation on Large Flowcells

- Issue: bcl-convert or kraken2 jobs killed by PBS OOM killer.
- Fix: raise `demux.demux_try1_gb` to 126+ for NovaSeq X; kraken2 needs ~178-190 GB on attempt 1.

### Qualimap strandedness is overridden

- Issue: `mapping.qualimap_parameters` is inserted into the command line BEFORE a hardcoded `-p strand-specific-forward`, so a user-supplied `-p` is overridden. Note the shipped config.yaml default is `-p strand-specific-reverse`, which is silently ignored.
- Fix: this is a code-level issue in `rule Qualimap`. Either accept the hardcoded forward strand, or patch the rule to not hardcode `-p`. Do not rely on `qualimap_parameters` to set strandedness.

### FeatureCounts -p for PE

- Issue: the `featurecounts` rule does not add `-p` automatically for PE.
- Fix: add the PE flags to `mapping.stranded` (the only string appended to the featureCounts call), e.g. `stranded: "-p --countReadPairs -s 2"`. Do **not** put them in `mapping.other_params`; that string goes to STAR and would break the alignment.

### rMATS -t argument

- Issue: `mapping.rmats_paired` is passed as one string to `-t`, e.g. `-t single --libType fr-secondstrand`; quoting matters.
- Fix: keep it as one YAML string as in the example; verify with a dry run.

### crypt4gh recipient string

- Issue: crypt4gh needs one `--recipient_pk` per public key.
- Fix: `scripts/get_c4gh_string.sh` builds the recipient string from `crypt4gh_client_public_key_dir`; ensure that directory contains all collaborator public keys (trailing `/` required).

### gocryptfs output is destructive

- Issue: when `gocryptfs.gocryptfs_active: True`, the rule copies the whole output into an encrypted mount with `--remove-source-files` and then deletes everything outside `encrypted/` (`rm -rf {output}/!(encrypted)/`). It is a `localrule`, prompts for a password, and runs after `multiqc`.
- Fix: only enable after verifying you have a copy; run it in a `screen` session (password must be typed manually); requires the `gocryptfs` env module.

---

## 16. Known Code Issues & Verification Notes

1. **PE path inconsistencies in `pe_processing.smk`** (verify with a dry run before relying on PE, especially with UMIs):
   - `umi_dedup_pe` consumes `star/{short}_pe_Aligned.sortedByCoord.out.bam`, but `align_pe` produces `star/{short}_pe_001_Aligned.sortedByCoord.out.bam` (missing `_001` in the input). The SE `umi_dedup` rule is not layout-gated and its pattern also matches the PE BAM names, so it is most likely what actually runs for PE.
   - `umi_extract_pe` consumes `untrimmed_fastq/{short}_R1.fastq.gz` / `{short}_R2.fastq.gz` (when cutadapt is off), but demux naming is `*_R1_001.fastq.gz`.
   - `umi_extract_pe` writes `umi_extract/{short}_R1.umis-extracted.fastq.gz`, while `get_mapping_input_pe1/2` (used by `align_pe` when UMIs are on and sortmerna is off) request `umi_extract/{short}_R1_001.fastq.gz`; no rule produces that name.
   - `getOutput()` requests `umi_extract/{X_R1_001,X_R2_001}.umis-extracted.fastq.gz` per file; these names match the single-file SE `umi_extract` rule (R1 `pattern` only), not `umi_extract_pe`.
   - For PE with UMIs, `getFastQCs()` requests `fastqc/{file}.umis-extracted_fastqc.*`, but `FastQC_pe` only produces `*_trimmed_fastqc.*`.
   - `count_matrix_pe` reads `star/{file}/ReadsPerGene.out.tab`: a nested `{file}` path that does not match the flat `star/` layout (and it competes with `featurecounts` for `counts/all.tsv`).
2. **Layout detection is fragile**: `isSingleEnd()` uses `name.split("_R")[1]`; names containing an earlier `_R` are misclassified as single-end, and names without `_R` raise `IndexError` in `pe_processing.smk` (section 7). The helper functions exist in several copies (`common.smk`, `pe_processing.smk`; `convert2fastq.smk` redefines `getSample_names_post_mapping()` without the `"_R" in sample` guard); `common.smk` also builds `sample_names` twice.
3. **`temp()` inside `params:` has no effect** (Snakemake only honors it on outputs). These intermediates therefore stay on disk: `kraken2/{file}.kraken2`, `kaiju/{file}_kaiju_result.out`, `biobloom/{file}.fastq` (uncompressed FASTQ), `blast/subsampled_{file}.fa`, `rustqc/{file}.bam` (markdup BAM), `sourmash/{file}_sketch.sig.gz`, `preseq/{file}_bam2bed.bed` and `star/sorted_{file}_Aligned.sortedByCoord.out.bam`. (`jellyfish/{file}.fastq` is a real `temp()` output and is removed.)
4. **Qualimap memory**: the rule requests `attempt * 12` GB from the scheduler but starts Java with a hardcoded `--java-mem-size=16G`, so the job can exceed its PBS memory request.
5. **`wait 30` in autostart**: `autostart_cronjob_bclconvert.sh` and `autostart_novaseq.sh` use `wait 30`, which does not sleep; `autostart_portable.sh` correctly uses `sleep 30`.
6. **Env mismatch in autostart**: `autostart/execute_pipeline.sh` activates `smk9`, while the wrappers/docs use `smk8`.
7. **`schemas/config.schema.yaml` is stale** (references old `bcl2fastq:`/`threads` keys) and is not wired into any Snakefile; do not use it for validation.
8. **Diamond**: only a commented-out `diamond_pe` rule exists in `pe_processing.smk`; there is no active Diamond rule and no `envs/diamond.yaml`.
9. **bash login-shell warning**: the demux rules run `chmod ago+rwx ...` (a bash-specific `chmod` alias seen on HILBERT login shells); on other systems this may fail; the original author's environment accepts it. `FastQC_untrimmed` unsets `command_not_found_handle` for the same reason.
10. **rMATS alternative**: the code notes Spladder ([https://spladder.readthedocs.io](https://spladder.readthedocs.io)) as a drop-in replacement if rMATS causes problems.