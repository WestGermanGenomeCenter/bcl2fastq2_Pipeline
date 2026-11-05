# Demultiplexing &amp; RNA-Seq Pipeline for Illumina Platforms

This Snakemake 8.28+ pipeline automates Illumina BCL demultiplexing (bcl-convert or bcl2fastq2) and provides downstream quality control, contamination detection, and RNA-Seq analysis on HPC environments (PBS Pro / HILBERT).

## Table of Contents

1. Overview &amp; Key Features
2. Prerequisites &amp; System Requirements
3. Deployment &amp; Quick Start
4. Running the Pipeline &amp; Execution Options
5. Pipeline Stages &amp; Architecture
6. The check\_fastq\_files Checkpoint &amp; check\_input\_fastqs.sh
7. Paired-End (PE) vs. Single-End (SE) Layout Rules
8. Automation &amp; Autostart Cronjob Setup
9. Complete Tool Reference &amp; Setup Instructions
10. Rule Registry (Snakefiles &amp; .smk includes)
11. Scripts Directory Guide (scripts/)
12. Checksums, Cryptography &amp; MultiQC Reports
13. Configuration Reference (config.yaml)
14. rMATS Group File Setup
15. Troubleshooting &amp; Common Pitfalls
16. Known Code Issues &amp; Verification Notes

---

## 1. Overview &amp; Key Features

### Core Capabilities

- **Dual Demultiplexer Support**: bcl-convert (recommended for NovaSeq X / NextSeq 2000) or bcl2fastq2 (legacy / MiSeq fallback). Controlled by `demux.use_bcl2fastq`.
- **Flexible Data Input**: raw BCL run directories or pre-demultiplexed FASTQs via `skip_demux.skip_demux_active`.
- **Read Layout Agnostic**: PE and SE supported. Layout inferred at parse time from FASTQ filenames via `isSingleEnd()` and cached in `fastq_infiles_list.tx`.
- **UMI Handling**: single/dual UMI regex patterns (`umi_tools.pattern`, `umi_tools.pattern2`); UMIs moved to read headers (`umi_extract/`); BAM deduplication post-mapping (`umi_dedup` -&gt; `*_deduped.Aligned.sortedByCoord.out.bam`).
- **Contamination Screening**: Kraken2, BioBloom Tools, BLAST+, Kaiju, Sourmash. (Diamond exists only as a commented-out `diamond_pe` rule; there is no active Diamond rule.)
- **End-to-End RNA-Seq**: STAR (`--quantMode GeneCounts`) -&gt; Qualimap / RSeQC / SAMtools -&gt; FeatureCounts (`counts/all.tsv`) -&gt; StringTie (3 stages) + gffcompare -&gt; rMATS (`rmats/summary.txt`).

### Advanced Features

- **Data Security &amp; Encryption**: per-file crypt4gh with one aggregate `encrypted_fastqs/checksums_encrypted_fastqs.sha256`; or directory-wide gocryptfs (which reduces the output folder to only `encrypted/`).
- **Automated Data Transfer**: `rsync -rzvP` via rule `transfer_files` (config `transfer.transfer_to`, `transfer.transfer_as`) — not scp.
- **Reprocibility &amp; Auditing**: `<projectNum>_software_environment_mqc_versions.yml`, `<projectNum>_config.yaml`, per-run `filelist_project_<timestamp>.filelist` and `checksums_<timestamp>.sha256`, `<projectNum>_sha256sums_fastqfiles.sha256`, per-stage rulegraph PDFs and Snakemake HTML reports (wrapper-level).
- **Smart Scaling**: rules use attempt-scaled `threads`, `time_hrs`, `mem_gb` resources. demux threads = `12 + attempt * 4`; demux memory is `round(demux_try1_gb * 0.7) + demux_try1_gb * 0.3 * attempt` (approx. 96 GB attempt 1, approx. 125 GB attempt 2 at default `demux_try1_gb: 96`). Kraken2 (`178 + attempt * 12`) and Kaiju (`88 + attempt * 40`) use custom memory formulas; sourmash is fixed at 1 thread.

---

## 2. Prerequisites &amp; System Requirements

### Compute &amp; Storage

- HPC Environment: PBS Pro (tested on HILBERT HPC) or local multi-core.
- Memory (attempt 1; scales on retry):
  - Demux: `demux_try1_gb` default 96 GB; set 126+ for NovaSeq X.
  - Kraken2: `178 + attempt * 12` GB -&gt; approx. 190 GB attempt 1.
  - Kaiju: `88 + attempt * 40` GB.
  - STAR align: requests `attempt * 40` GB. (STAR genome indexing is not part of this repo; a pre-built index is expected — typically \~30-35 GB RAM to build, organism-dependent.)
- Storage: raw BCLs, intermediate FASTQs, aligned BAMs, and (if enabled) `encrypted/`.

### Required Software

1. Conda / Mamba. The repo ships `smk8_environment.yml` and per-tool `envs/*.yaml`.
2. Snakemake &gt;= 8.28.0 (`min_version("8.28.0")` in both entry Snakefiles).
3. bcl-convert and/or bcl2fastq2 binaries from Illumina (credentials required). Paths: `demux.bcl_convert_path`, `demux.bcl2fastq2_path`.
4. Qualimap executable (manual install or `--no-deps` on HPC; see `envs/install_Qualimap_manually.sh`). Path: `mapping.qualimap_dir` — the rule calls `{qualimap_dir}/qualimap` directly.
5. A PBS Snakemake profile named `pbs` (see section 3 — the wrapper calls `--profile pbs`).
6. gocryptfs (`gocryptfs.gocryptfs_version`, e.g. `gocryptfs/1.7.1`) as a cluster module if `gocryptfs.gocryptfs_active: True` (loaded via `envmodules:`; needs cluster module support).

---

## 3. Deployment &amp; Quick Start

### HILBERT HPC Users

```bash
git clone https://github.com/WestGermanGenomeCenter/bcl2fastq2_Pipeline.git
cd bcl2fastq2_Pipeline
bash setup_pipeline.sh
```

`setup_pipeline.sh` creates the `smk8` conda env and copies a ready-made PBS Snakemake profile from a HILBERT-internal path (`/gpfs/project/projects/bmfz_gtl/software/snakemake_profile/pbs`) to `~/.config/snakemake/pbs` — this only works for HILBERT users.

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

There is no top-level `Snakefile`. The pipeline has two entry Snakefiles, both driven through the wrappers:

- `rules/convert2fastq.smk` — demultiplexing stage (includes `common.smk`)
- {-`rules/fastqProcessing.smk` — processing/analysis stage (includes `pe_processing.smk`, `common.smk`, `qc.smk`)-}
- {+`rules/fastqProcessing.smk` — processing/analysis stage (includes `pe_processing.smk`, `common.smk`, `qc.smk`; also contains a stale `include: "qcShort.smk"` — see section 16)+}

### Repository Layout (observed)

- `rules/` — the two entry Snakefiles plus `common.smk`, `pe_processing.smk`, `qc.smk`.
- `envs/` — one YAML per tool: `Qualimap.yaml`, `STAR.yaml`, `bcl2fastq2.yaml`, `biobloom.yaml`, `blast.yaml`, `counts.yaml`, `crypt4gh.yaml`, `cutadapt.yaml`, `demux.yaml`, `fastqc.yaml`, `gffcompare.yaml`, `gffutils.yaml`, `interop.yaml`, `jellyfish.yaml`, `kaiju.yaml`, `kraken2.yaml`, `multiqc.yaml`, `pandas.yaml`, `pca_plot.yaml`, `preseq.yaml`, `rmats.yaml`, `rseqc.yaml`, `rustqc.yaml`, `samtools.yaml`, `sortmerna.yaml`, `sourmash.yaml`, `stringtie.yaml`, `transfer.yaml`, `umi_tools.yaml` (plus `tools.yaml` and helper scripts for the software-version list).
- `scripts/` — see section 11.
- `smk8_environment.yml` — Snakemake driver env.
- `autostart/` — cron-based autonomous execution (see section 8).
- `schemas/config.schema.yaml` — legacy schema; not referenced by any Snakefile and out of date with the current config.yaml (do not rely on it for validation).
- `examples/` — example config, rMATS group file, and rendered reports/graphs.
- `snakemake-utils/` — PBS profile components (submit/status scripts, jobscript, test workflow).
- `runPipeline.sh`, `runPipeline_local.sh`, `runPipeline_barebone.sh`, `setup_pipeline.sh`, `setup_pipeline_local.sh`.

---

## 4. Running the Pipeline &amp; Execution Options

### Standard HPC Execution (runPipeline.sh)

```bash
screen -dRR pipeline_run
qsub -I -l select=1:ncpus=2:mem=8gb,walltime=48:00:00
conda activate smk8
bash runPipeline.sh
```

Detach: `Ctrl + A`, then `D`. Reattach: `screen -r pipeline_run`.

### Execution Flags &amp; Variants

- Skip wrapper-level checksum generation (faster):
  ```bash
  bash runPipeline.sh --no-checksums
  ```
  Caveat: this only skips the wrapper-level \`checksums\_&lt;timestamp&gt;.sha256\` and \`filelist\_project\_&lt;timestamp&gt;.filelist\`. The demux stage itself still runs rule \`create\_checksums\` (per-project \`untrimmed\_fastq/&lt;projectNum&gt;\_sha256sums\_fastqfiles.sha256\`) — that rule is unconditional in \`convert2fastq.smk\`.
- Local workstation:
  ```bash
  bash runPipeline_local.sh
  ```
  (runs \`snakemake --cores all --sdm conda --jobs 7/4\` instead of the PBS profile, and uses \`nice\`.)
- Barebones (direct Snakemake calls, no reports/checksums):
  ```bash
  bash runPipeline_barebone.sh
  ```

### What runPipeline.sh Actually Adds

1. Per-stage rulegraph PDFs: `convert2fastq_rulegraph.<timestamp>.pdf` and `fastqProcessing_rulegraph.<timestamp>.pdf` (via `snakemake --forceall --rulegraph | dot -Tpdf`; `dot` from graphviz is required).
2. Per-stage Snakemake HTML reports: `convert2fastq_snakemake_report.<timestamp>.html` and `fastqProcessing_snakemake_report.<timestamp>.html`.
3. `filelist_project_<timestamp>.filelist` and `checksums_<timestamp>.sha256` over the whole output folder (unless `--no-checksums`).
4. Moves files from a previous execution of the same output folder into `<out>/logs/previous_executions/` (reports, rulegraphs, checksums, filelists, configs).
5. Prints which config.yaml options are enabled/disabled before execution.

Note: `copy_software_env` (rule, demux stage) copies the versions file to `<projectNum>_software_environment_mqc_versions.yml` and `config.yaml` to `<projectNum>_config.yaml` in the output folder. Wrapper file names carry a `%Y_%m_%d_%I_%M_%p` timestamp — there are no plain `filelist.txt`/`checksums.sha256`/`dag.pdf` files.

---

## 5. Pipeline Stages &amp; Architecture

### Stage 1: Demultiplexing (`rules/convert2fastq.smk`)

```text
Illumina RunFolder (BCLs) + SampleSheet.csv
                  |
                  v
 check_inputs (bcl-convert --bcl-validate-sample-sheet-only true
               + scripts/bclconvert_samplesheet_check.py)
                  |
                  v
 bcl-convert / bcl2fastq2  ->  *.fastq.gz -> moved to untrimmed_fastq/
                               (Undetermined*.fastq.gz moved to output root)
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

Optional around demux: `prevent_revcomp` / `reverse_prevent_revcomp` (RunInfo.xml patch, off by default), `skip_demux` (move pre-demultiplexed FASTQs into `untrimmed_fastq/`).

### Stage 2: Processing, Contamination &amp; RNA-Seq (`rules/fastqProcessing.smk`)

```text
untrimmed_fastq/
      |
      v
 CHECKPOINT check_fastq_files (scripts/check_input_fastqs.sh)
      |
      v
 FastQC + Fastp untrimmed (fastqc_untrimmed/) -> multiqc_raw
      |
      v
 cutadapt           ->  trimmed/{file}_trimmed.fastq.gz
      |
      v
 umi_tools extract  ->  umi_extract/{file}.umis-extracted.fastq.gz
      |
      v
 sortmerna          ->  sortmerna/{file}_non-ribosomal_rna.fq.gz
      |                (+ {file}_ribosomal_rna.fq.gz)
      |
  +---+---------------+---------------------+
  |                   |                     |
  v                   v                     v
Contamination       Preseq / Jellyfish   STAR align
(Kraken2/BioBloom/  (preseq/, jellyfish/) -> star/{file}_Aligned.sortedByCoord.out.bam
 BLAST/Kaiju/                               star/{file}_ReadsPerGene.out.tab
 Sourmash)                                   |
                                             v
                                  Qualimap / RSeQC / SAMtools
                                             |
                                             v
                                  umi_tools dedup (if UMIs) ->
                                  star/{file}_deduped.Aligned.sortedByCoord.out.bam
                                             |
                                             v
                      FeatureCounts (counts/all.tsv) + PCA/heatmap/reports
                      + StringTie (3 stages) + gffcompare
                                             |
                                             v
                      multiqc_report_complete_<projectNum>.html
```

**BAM-input switch** — when `umi_tools.umi_tools_active: True`, the following rules consume `star/{file}_deduped.Aligned.sortedByCoord.out.bam` instead of `star/{file}_Aligned.sortedByCoord.out.bam`: `samtools`, `preseq`, `stringtie`, `stringtie_estimate`, `Qualimap`, `featurecounts`, `star_result_table`, `rmats_create_groupfiles_from_csv`, `rustqc`, and the RSeQC rules in `qc.smk` (`rseqc_main`, `rseqc_minor`, `rseqc_done`). `align` itself takes FASTQs, and `pca_plot` takes `counts/all.tsv`, so neither switches inputs. Enabling `umi_tools.umi_tools_active: True` therefore changes the BAM input of \~10 rules (plus `umi_dedup` producing the deduped BAMs).

---

## 6. The check\_fastq\_files Checkpoint &amp; check\_input\_fastqs.sh

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

### check\_input\_fastqs.sh validation steps

1. Directory argument provided and exists.
2. Directory non-empty.
3. At most ONE zero-byte `*.fastq.gz` is tolerated (emits a warning); more than one fails the checkpoint.
4. Every `.fastq.gz` must be &gt;= 80 bytes (filters truncated/corrupted output).
5. On success the message is appended to `inputfiles_checkpoint_completed.out`.

Downstream rules depend on the flag via `nonempty_fastqs=` (e.g. `cutadapt`, `cutadapt_pe`) and via `getOutput()` in rule `all`.

---

## 7. Paired-End (PE) vs. Single-End (SE) Layout Rules

- Layout detection: `isSingleEnd()`, evaluated at parse time; PE/SE blocks gate `cutadapt`, `align`, `crypt4gh`, `FastQC`, `FastQC_untrimmed`, `sortmerna`, and `get_mapping_input` in `fastqProcessing.smk`, and the whole PE rule set in `pe_processing.smk`. Note: if even one R2 file is missing, the entire set is treated as single-end.
- Sample name derivation:
  - `sample_names` = FASTQ basenames minus `.fastq.gz`, excluding `_I1_`/`_I2_` index reads.
  - `getSample_names_post_mapping()` returns `sample_names` for SE, and for PE replaces `_R1` with `_pe`. Mapping outputs use `<sample>_pe_...`, not `<sample>_R1_...`.
- Naming: PE inputs must be `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`. Index reads (I1/I2) are ignored after demux.
- cutadapt: SE trims one file per job; PE (`cutadapt_pe`) trims R1+R2 simultaneously, applying all adapter pairs from `config["cutadapt"]["adapters"]` with types from `adapter_type` (default `b,b,b,b`), plus `other_params` (default `--minimum-length 25 --quality-cutoff 5 -n 2`). Output: `trimmed/{file}_trimmed.fastq.gz` (SE) / `trimmed/{short}_R1_001_trimmed.fastq.gz` + R2 (PE).
- umi\_tools: `pattern` for R1, optional `pattern2` for R2 in dual-UMI setups (PE uses `--read2-in`/`--read2-out` + `--bc-pattern2`); outputs `umi_extract/{file}.umis-extracted.fastq.gz` (SE) / `{short}_R1.umis-extracted.fastq.gz` + R2 (PE).
- STAR: `--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts`; outputs `star/{file}_Aligned.sortedByCoord.out.bam` and `star/{file}_ReadsPerGene.out.tab`.
- FeatureCounts: `-g gene_id` + `config["mapping"]["stranded"]` (`-s 0` unstranded, `-s 1` forward, `-s 2` reverse); the rule does NOT add `-p` automatically for PE — add `-p` to `mapping.other_params` (or `stranded`) yourself when running PE.
- rMATS: `mapping.rmats_paired` is passed as `-t`. Example value: `"single --libType fr-secondstrand"` — Snakemake splits it into separate shell words, so keep it as one string and verify with a dry run.

---

## 8. Automation &amp; Autostart Cronjob Setup

`/autostart` enables autonomous execution for core-facility sequencers via `execute_pipeline.sh` plus one of three cron wrappers:

- `autostart_cronjob_bclconvert.sh` — hardcoded for the NextSeq 2000 watch directory.
- `autostart_novaseq.sh` — same logic, hardcoded for a NovaSeq X watch directory.
- `autostart_portable.sh` — parameterized version (`--watch-dir`, `--pipeline`, `--snakemake-dir`, `--execution-script`, `--env-script`, `--session-prefix`); requires an environment script that sets PATH/modules/conda.

1. **Manual config required**: autostart does not construct config.yaml. An operator must place a populated `config.yaml` in each run folder (the wrapper copies it to the pipeline root before launching).
2. **Directory monitoring**: the cron wrapper scans the instrument storage for run folders that:
  - contain `CopyComplete.txt` (no RTAComplete.txt check is performed),
  - contain exactly `config.yaml` (not `*_config.yaml`),
  - do not contain the lockfile `.pipeline_launched`,
  - while no job named with prefix `demx_pipe_auto` (PBS job name from `execute_pipeline.sh`) is in the `qstat` queue.
3. **Launch**: creates `.pipeline_launched` (with launch date + session), copies the run's `config.yaml` to the pipeline root, waits \~30 s, and submits `execute_pipeline.sh` via `qsub -q default`.

Caveats:

- `execute_pipeline.sh` activates conda env `smk9` (not `smk8`) and hardcodes the pipeline directory — adjust both for your site.
- `autostart_cronjob_bclconvert.sh` and `autostart_novaseq.sh` use `wait 30`, which in bash only waits for background jobs of the current shell — it does NOT sleep. `autostart_portable.sh` correctly uses `sleep 30`.
- The wrappers also hardcode the invoking user's PATH to make `qstat` work from a minimal cron environment.

---

## 9. Complete Tool Reference &amp; Setup Instructions


| Tool                    | Config Key                                  | Purpose                                                                                                 | Docs                                                                                                                                                                                                                                                                                                                    | Source                                                                                           |
| ----------------------- | ------------------------------------------- | ------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| bcl-convert             | demux.bcl\_convert\_path                    | Primary demux for NovaSeq X / NextSeq 2000                                                              | [https://support.illumina.com/sequencing/sequencing\_software/bcl-convert.html](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html)                                                                                                                                                           | Proprietary                                                                                      |
| bcl2fastq2              | demux.bcl2fastq2\_path                      | Legacy demux for MiSeq / older flowcells                                                                | [https://support.illumina.com/content/dam/illumina-support/documents/documentation/software\_documentation/bcl2fastq/bcl2fastq2\_val\_guide\_15051736\_02.pdf](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2_val_guide_15051736_02.pdf) | Proprietary                                                                                      |
| FastQC                  | (untrimmed/trimmed QC)                      | Per-base quality metrics                                                                                | [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)                                                                                                                                                                                                | [https://github.com/s-andrews/FastQC](https://github.com/s-andrews/FastQC)                       |
| Fastp                   | -                                           | Fast quality filtering, HTML/JSON reports (run alongside FastQC)                                        | [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)                                                                                                                                                                                                                                                  | [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)                           |
| Cutadapt                | cutadapt                                    | Adapter trimming, length/quality filters                                                                | [https://cutadapt.readthedocs.io/](https://cutadapt.readthedocs.io/)                                                                                                                                                                                                                                                    | [https://github.com/marcelm/cutadapt](https://github.com/marcelm/cutadapt)                       |
| STAR                    | mapping                                     | Splice-aware RNA-Seq alignment (--quantMode GeneCounts)                                                 | [https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)                                                                                                                                                                                    | [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)                           |
| Qualimap                | mapping.qualimap\_on, mapping.qualimap\_dir | Post-alignment BAM QC                                                                                   | [http://qualimap.bioinfo.cipf.es/](http://qualimap.bioinfo.cipf.es/)                                                                                                                                                                                                                                                    | [https://bitbucket.org/kokonech/qualimap](https://bitbucket.org/kokonech/qualimap)               |
| RSeQC                   | mapping.rseqc\_active                       | Alternative post-alignment QC                                                                           | [https://rseqc.sourceforge.net/](https://rseqc.sourceforge.net/)                                                                                                                                                                                                                                                        | [https://github.com/davidemms/RSeQC](https://github.com/davidemms/RSeQC)                         |
| FeatureCounts (Subread) | mapping                                     | Assign BAM reads to GTF features (counts/all.tsv)                                                       | [https://subread.sourceforge.net/SubreadUsersGuide.pdf](https://subread.sourceforge.net/SubreadUsersGuide.pdf)                                                                                                                                                                                                          | SourceForge                                                                                      |
| StringTie               | mapping.stringtie\_on                       | Assemble/quantify transcripts (3 stages)                                                                | [https://ccb.jhu.edu/software/stringtie/](https://ccb.jhu.edu/software/stringtie/)                                                                                                                                                                                                                                      | [https://github.com/gpertea/stringtie](https://github.com/gpertea/stringtie)                     |
| gffcompare              | -                                           | Compare StringTie GTF to reference                                                                      | [https://ccb.jhu.edu/software/stringtie/gffcompare.shtml](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)                                                                                                                                                                                                      | [https://github.com/gpertea/gffcompare](https://github.com/gpertea/gffcompare)                   |
| rMATS                   | mapping.rmats\_on                           | Differential alternative splicing                                                                       | [https://github.com/XingLab/rmats-turbo](https://github.com/XingLab/rmats-turbo)                                                                                                                                                                                                                                        | [https://github.com/XingLab/rmats-turbo](https://github.com/XingLab/rmats-turbo)                 |
| Preseq                  | preseq                                      | Library complexity / saturation                                                                         | [https://github.com/smithlabcode/preseq](https://github.com/smithlabcode/preseq)                                                                                                                                                                                                                                        | [https://github.com/smithlabcode/preseq](https://github.com/smithlabcode/preseq)                 |
| bedtools                | -                                           | bamtobed for Preseq                                                                                     | [https://bedtools.readthedocs.io/](https://bedtools.readthedocs.io/)                                                                                                                                                                                                                                                    | [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2)                         |
| Kraken2                 | kraken2                                     | k-mer taxonomic classification                                                                          | [https://ccb.jhu.edu/software/kraken2/](https://ccb.jhu.edu/software/kraken2/)                                                                                                                                                                                                                                          | [https://github.com/DerrickWood/kraken2](https://github.com/DerrickWood/kraken2)                 |
| Sourmash                | sourmash                                    | MinHash signature comparison                                                                            | [https://sourmash.readthedocs.io/](https://sourmash.readthedocs.io/)                                                                                                                                                                                                                                                    | [https://github.com/sourmash-bio/sourmash](https://github.com/sourmash-bio/sourmash)             |
| BioBloom Tools          | biobloom                                    | Bloom-filter organism filtering                                                                         | [https://github.com/bcgsc/biobloom](https://github.com/bcgsc/biobloom)                                                                                                                                                                                                                                                  | [https://github.com/bcgsc/biobloom](https://github.com/bcgsc/biobloom)                           |
| BLAST+                  | blast                                       | blastn on subsampled reads                                                                              | [https://www.ncbi.nlm.nih.gov/books/NBK279690/](https://www.ncbi.nlm.nih.gov/books/NBK279690/)                                                                                                                                                                                                                          | NCBI                                                                                             |
| Diamond                 | -                                           | Protein-level alignment — only a commented-out `diamond_pe` rule exists; no active rule and no env YAML | [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)                                                                                                                                                                                                                                            | [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)                     |
| Kaiju                   | kaiju                                       | Protein-level metagenomic classification                                                                | [https://github.com/bioinformatics-centre/kaiju](https://github.com/bioinformatics-centre/kaiju)                                                                                                                                                                                                                        | [https://github.com/bioinformatics-centre/kaiju](https://github.com/bioinformatics-centre/kaiju) |
| Jellyfish               | jellyfish                                   | 21-mer frequency distributions                                                                          | [https://www.cbcb.umd.edu/software/jellyfish/](https://www.cbcb.umd.edu/software/jellyfish/)                                                                                                                                                                                                                            | [https://github.com/gmarcais/Jellyfish](https://github.com/gmarcais/Jellyfish)                   |
| SortMeRNA               | sortmerna                                   | rRNA removal                                                                                            | [https://sortmerna.readthedocs.io/](https://sortmerna.readthedocs.io/)                                                                                                                                                                                                                                                  | [https://github.com/biocore/sortmerna](https://github.com/biocore/sortmerna)                     |
| UMI-tools               | umi\_tools                                  | UMI extraction + BAM dedup                                                                              | [https://umi-tools.readthedocs.io/](https://umi-tools.readthedocs.io/)                                                                                                                                                                                                                                                  | [https://github.com/CGATOxford/UMI-tools](https://github.com/CGATOxford/UMI-tools)               |
| rustqc                  | rustqc                                      | Alternative RNA QC (needs sambamba markdup)                                                             | [https://github.com/rustqc/rustqc](https://github.com/rustqc/rustqc)                                                                                                                                                                                                                                                    | [https://github.com/rustqc/rustqc](https://github.com/rustqc/rustqc)                             |
| sambamba                | -                                           | markdup pre-step for rustqc                                                                             | [https://lomereiter.github.io/sambamba/](https://lomereiter.github.io/sambamba/)                                                                                                                                                                                                                                        | [https://github.com/biod/sambamba](https://github.com/biod/sambamba)                             |
| Crypt4gh                | crypt4gh                                    | GA4GH per-file public-key encryption                                                                    | [https://crypt4gh.readthedocs.io/](https://crypt4gh.readthedocs.io/)                                                                                                                                                                                                                                                    | [https://github.com/ga4gh/crypt4gh](https://github.com/ga4gh/crypt4gh)                           |
| Gocryptfs               | gocryptfs                                   | Encrypted overlay directory                                                                             | [https://nuetzlich.net/gocryptfs/](https://nuetzlich.net/gocryptfs/)                                                                                                                                                                                                                                                    | [https://github.com/rfjakob/gocryptfs](https://github.com/rfjakob/gocryptfs)                     |


---

## 10. Rule Registry (Snakefiles &amp; .smk includes)

### rules/convert2fastq.smk (demux stage)

`localrules: all, create_fastq_list`


| Rule                      | Purpose / Output                                                                                                                                                |
| ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| all                       | collects `getOutput()` (fastq list, checksums, software env copy, interop flag, revcomp flags, bcl2fastq2 flag, SampleSheet\_check.out)                         |
| prevent\_revcomp          | patches RunInfo.xml via `scripts/prevent_revcomp_nexts2.sh` -&gt; `revcomp_prevented.flag`                                                                      |
| reverse\_prevent\_revcomp | reverts the patch -&gt; `reversed_revcomp_prevented.flag`                                                                                                       |
| check\_inputs             | bcl-convert samplesheet validation + `scripts/bclconvert_samplesheet_check.py` -&gt; `SampleSheet_check.out`                                                    |
| demux                     | bcl-convert run -&gt; `Reports/Adapter_Metrics.csv`, FASTQs moved to `untrimmed_fastq/`, Undetermined moved to output root, bcl logs to `logs/bcl_convert_logs` |
| bcl2fastq2                | legacy demux -&gt; `used_bcl2fastq2.flag`, FASTQs to `untrimmed_fastq/`                                                                                         |
| skip\_demux               | moves pre-demultiplexed FASTQs -&gt; `skipped_demuxing.flag`                                                                                                    |
| create\_fastq\_list       | `fastq_infiles_list.tx` (recursive glob, excludes Undetermined)                                                                                                 |
| create\_checksums         | `untrimmed_fastq/<projectNum>_sha256sums_fastqfiles.sha256` (parallel `xargs -P {threads}`)                                                                     |
| copy\_software\_env       | versions TSV -&gt; `<projectNum>_software_environment_mqc_versions.yml`, config -&gt; `<projectNum>_config.yaml`                                                |
| interop\_plots            | InterOp plots + reports (Sequencing\_report.pdf, RunReport.pdf, Demux\_report.pdf, TiltReport.pdf) -&gt; `interop_plots/interop_plots_done.flag`                |


### rules/fastqProcessing.smk (processing stage)

`localrules: all, gocryptfs`


| Rule                                 | Gate          | Purpose / Output                                                                                                                  |
| ------------------------------------ | ------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| all                                  | -             | collects `getOutput()`                                                                                                            |
| check\_fastq\_files                  | checkpoint    | `inputfiles_checkpoint_completed.out`                                                                                             |
| cutadapt                             | SE            | `trimmed/{file}_trimmed.fastq.gz`                                                                                                 |
| align                                | SE            | STAR -&gt; `star/{file}_Aligned.sortedByCoord.out.bam` + `ReadsPerGene.out.tab`                                                   |
| crypt4gh                             | SE            | `encrypted_fastqs/{file}_processed.fastq.gz.c4gh`                                                                                 |
| FastQC\_untrimmed                    | SE            | FastQC + fastp on raw FASTQs -&gt; `fastqc_untrimmed/`                                                                            |
| FastQC                               | SE            | FastQC + fastp on trimmed/UMI-extracted FASTQs -&gt; `fastqc/`                                                                    |
| sortmerna                            | SE            | rRNA split -&gt; `sortmerna/{file}_[non-]ribosomal_rna.fq.gz`                                                                     |
| umi\_extract                         | umi\_tools on | `umi_extract/{file}.umis-extracted.fastq.gz`                                                                                      |
| umi\_dedup                           | umi\_tools on | `star/{file}_deduped.Aligned.sortedByCoord.out.bam`                                                                               |
| samtools                             | -             | index/stats/coverage/plot-bamstats -&gt; `samtools/{file}_samtools_stats.stats`                                                   |
| preseq                               | -             | c\_curve on BAM-derived BED -&gt; `preseq/{file}_c_curve_out.txt`                                                                 |
| Qualimap                             | -             | `qualimap/{file}/qualimapReport.html` (hardcoded `-p strand-specific-forward`)                                                    |
| featurecounts                        | -             | `counts/all.tsv`                                                                                                                  |
| pca\_plot                            | -             | pca.R + heatmap.py + counts\_report.py -&gt; `counts/pca.png`, `<projectNum>_CountsReport.pdf`, `<projectNum>_Counts_Report.html` |
| star\_result\_table                  | -             | summarize\_star\_out.py -&gt; `star/star_logs_summarized.csv`                                                                     |
| stringtie                            | -             | `stringtie/transcripts_{file}.gtf` (+ abundances .tab)                                                                            |
| stringtie\_merge                     | -             | `stringtie/transcripts_all_merged.gtf`                                                                                            |
| stringtie\_estimate                  | -             | `stringtie/transcript_abundances_{file}.gtf`                                                                                      |
| gffcompare                           | -             | `gffcompare/stringtie_to_ref_comp_{file}.stats`                                                                                   |
| rmats\_create\_groupfiles\_from\_csv | rmats on      | generate\_rmats\_groups.sh -&gt; `rmats/group[1,2]_bams.txt`                                                                      |
| run\_rmats                           | rmats on      | rmats.py -&gt; `rmats/summary.txt`                                                                                                |
| kraken2                              | kraken2 on    | `kraken2/{file}.report`                                                                                                           |
| kaiju                                | -             | kaiju + kaiju2table -&gt; `kaiju/kaiju_species_table_{file}.txt`                                                                  |
| blast                                | -             | subsample + blastn + filter\_blast.sh -&gt; `blast/blast_report_{file}.tsv`, `blast_summary_{file}.tsv`                           |
| jellyfish                            | -             | 21-mer histogram -&gt; `jellyfish/{file}_trimmed_jf.hist`                                                                         |
| biobloom                             | -             | biobloomcategorizer -&gt; `biobloom/biobloom_results_{file}_summary.tsv`                                                          |
| sourmash                             | -             | sketch + gather -&gt; `sourmash/{file}_sourmash.csv`                                                                              |
| rustqc                               | -             | sambamba markdup + rustqc rna -&gt; `rustqc/featurecounts/{file}.biotype_counts.tsv`                                              |
| transfer\_files                      | transfer on   | rsync -&gt; `transfer/transfer_log.log`                                                                                           |
| multiqc\_raw                         | -             | MultiQC on `fastqc_untrimmed/` -&gt; `multiqc_report_<projectNum>.html`                                                           |
| multiqc                              | -             | MultiQC on all results -&gt; `multiqc_report_complete_<projectNum>.html`                                                          |
| gocryptfs                            | gocryptfs on  | interactive gocryptfs -&gt; `encrypted/gocryptfs.conf`                                                                            |


### rules/pe\_processing.smk (PE only, `isSingleEnd() == False`)

`cutadapt_pe`, `align_pe`, `umi_dedup_pe`, `FastQC_pe`, `FastQC_untrimmed_pe`, `crypt4gh_pe`, `sortmerna_pe`, `umi_extract_pe`, `count_matrix_pe` (featureCounts alternative when UMIs are off; `diamond_pe` exists only commented out).

### rules/qc.smk (RSeQC, `mapping.rseqc_active`)

`rseqc_gtf2bed` (scripts/gtf2bed.py -&gt; `rseqc/annotation.bed`), `rseqc_main` (junction\_annotation, bam\_stat, junction\_saturation, infer\_experiment, read\_NVC, FPKM\_count, inner\_distance), `rseqc_minor` (read\_distribution, read\_duplication, read\_GC), `rseqc_done` (`rseqc/done.flag`).

### rules/common.smk (helpers, no rules)

`getFiles()`, `isSingleEnd()`, `getSample_names_post_mapping()`, `getSamples()`, `getFastqs()`; reads `fastq_infiles_list.tx` and defines `outputfolder`.

---

## 11. Scripts Directory Guide (scripts/)

Scripts actually present and referenced by the current rules:

- `check_input_fastqs.sh` — validates `untrimmed_fastq/` (non-empty, max one zero-byte file, all &gt;= 80 bytes). Called by `check_fastq_files`.
- `bclconvert_samplesheet_check.py` — parses/validates SampleSheet.csv before bcl-convert; called by `check_inputs`.
- `prevent_revcomp_nexts2.sh` / `reverse_revcomp_nexts2.sh` — patch/revert RunInfo.xml to prevent unexpected index-2 reverse complementation. Called by `prevent_revcomp` / `reverse_prevent_revcomp`.
- `interop_plots.sh` — InterOp plot wrapper; called by `interop_plots`.
- `test.py` — produces `interop_plots/Sequencing_report.pdf` (generically named; verify intent).
- `cluster_pdf.py` — produces `interop_plots/RunReport.pdf`.
- `demux_report.py` — produces `interop_plots/Demux_report.pdf`.
- `tilt_report.py` — produces `interop_plots/TiltReport.pdf`.
- `gtf2bed.py` — GTF to BED for RSeQC; called by `rseqc_gtf2bed`.
- `generate_rmats_groups.sh` — builds `rmats/group1_bams.txt` and `rmats/group2_bams.txt` from `mapping.rmats_group_file`.
- `pca.R` — PCA on `counts/all.tsv`; writes `counts/pca.png` (called by `pca_plot`).
- `heatmap.py` — writes `<projectNum>_CountsReport.pdf` (called by `pca_plot`).
- `counts_report.py` — writes `<projectNum>_Counts_Report.html` (called by `pca_plot`).
- `summarize_star_out.py` — writes `star/star_logs_summarized.csv` (called by `star_result_table`).
- `fastq_to_fasta_subsampling.sh` / `filter_blast.sh` — subsample FASTQ -&gt; FASTA and summarize blastn output for the `blast` rule. Note: the rule invokes them from the configured `blast.blast_script_dir`, NOT from the repo's `scripts/` folder — keep copies in sync.
- `get_c4gh_string.sh` — builds the `crypt4gh encrypt --recipient_pk ...` string from `crypt4gh_client_public_key_dir`; called by `crypt4gh`/`crypt4gh_pe`.
- `helper.py` — `validateSamplesheet`, `validateOutput`, `validateProjectNum`, `adaptersToStringParams`.
- `count-matrix.py` — builds `counts/all.tsv` from STAR ReadstPerGene tabs; referenced only by the (commented-out / PE-only) `count_matrix` rules.

Present in the repo but not referenced by any current rule: `interop_minimal.py`, `interop_more.py`, `interop_more_old.py`, `clusters.py`, `get_cluster_numbers.py`, `fastq_to_fasta.sh`, `jobscript.sh`, `tubemap.py`, `versions.py`. Older scripts (`interop_plots.py`, `clusters_pdf.py`, `prevent_revcomp.py`) do not exist in the repo.

---

## 12. Checksums, Cryptography &amp; MultiQC Reports

### Checksums &amp; File Lists

1. Raw demux checksums: `create_checksums` writes `untrimmed_fastq/<projectNum>_sha256sums_fastqfiles.sha256` (sha256sum over all non-Undetermined FASTQs, parallelized with `xargs -P {threads}`).
2. Trimmed checksums: `cutadapt`/`cutadapt_pe` appends to `trimmed/sha256sums_fastqfiles_trimmed.sha256`.
3. UMI-extracted checksums: `umi_extract`/`umi_extract_pe` appends to `umi_extract/checksums_umis_extracted_<projectNum>.sha256`.
4. Crypt4gh checksums: `crypt4gh`/`crypt4gh_pe` appends to `encrypted_fastqs/checksums_encrypted_fastqs.sha256` — one aggregate file, not per-file `.c4gh.sha256`.
5. Wrapper-level: `filelist_project_<timestamp>.filelist` and `checksums_<timestamp>.sha256` from `runPipeline.sh` (unless `--no-checksums`).

### MultiQC Reports

1. `multiqc_report_<projectNum>.html` — rule `multiqc_raw`; runs inside `fastqc_untrimmed/` with `-m fastqc`, `--fullnames`, `--ignore-samples Undetermined*`, `-x Undetermined*`, `--no-data-dir`; also writes `logs/fastq_file_list.tsv`.
2. `multiqc_report_complete_<projectNum>.html` — rule `multiqc`; aggregates `getFastQCs()` (raw+trimmed FastQC/Fastp, STAR, samtools, FeatureCounts, contamination, Qualimap, StringTie, gffcompare, rMATS, rustqc, etc.) with `--fullnames` and `--no-data-dir`, then deletes empty `.log` files (list in `logs/deleted_empty_files.out`).

---

## 13. Configuration Reference (config.yaml)

Only two things are mandatory per run: a samplesheet with a project number in the filename, and an output folder. Trailing slashes: allowed/required ONLY on `biobloom_filter_dir` and `transfer.transfer_to`; all other paths must NOT end with `/`.


| Section                  | Key                                                                               | Notes                                                                                                                                     |
| ------------------------ | --------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| demux                    | SampleSheet                                                                       | full path; must be inside the sequencing run folder; project number must be part of the samplesheet filename (e.g. 1399\_Samplesheet.csv) |
| demux                    | OutputFolder                                                                      | complete dir path, no trailing /                                                                                                          |
| demux                    | options                                                                           | extra bcl-convert CLI parameters (e.g. `--no-lane-splitting true`)                                                                        |
| demux                    | bcl\_convert\_path                                                                | path to the bcl-convert executable, no trailing /                                                                                         |
| demux                    | demux\_try1\_gb                                                                   | RAM for first demux attempt (96; 126+ for NovaSeq X)                                                                                      |
| demux                    | use\_bcl2fastq                                                                    | True for MiSeq / if bcl-convert errors                                                                                                    |
| demux                    | bcl2fastq2\_path, bcl2fastq\_mismatches, bcl2fastq\_options                       | bcl2fastq2 settings                                                                                                                       |
| skip\_demux              | skip\_demux\_active, fastq\_folder                                                | use pre-demultiplexed FASTQs; folder must contain the `.fastq.gz` files (no trailing /)                                                   |
| prevent\_revcomp         | prevent\_revcomp\_active                                                          | patch RunInfo.xml against read-3 revcomp; keep inactive unless you know better                                                            |
| interop\_plots           | interop\_plots\_active                                                            | adds SAV-like plots + shareable report PDFs (needs the Illumina run folder)                                                               |
| cutadapt                 | cutadapt\_active, adapters, adapter\_type, other\_params                          | one adapter\_type per adapter (`a`/`g`/`b`)                                                                                               |
| mapping                  | mapping\_active, gtf\_file, genomeDir, other\_params, stranded                    | STAR settings; genomeDir no trailing /; stranded is the featureCounts `-s` value                                                          |
| mapping                  | stringtie\_on, rseqc\_active, qualimap\_on, qualimap\_parameters, qualimap\_dir   | post-mapping QC; enable rseqc OR qualimap, not both                                                                                       |
| mapping                  | rmats\_on, rmats\_read\_length, rmats\_paired, rmats\_group\_file                 | rMATS settings (read length = sequenced length minus UMI/spacer/adapter)                                                                  |
| rustqc                   | rustqc\_active                                                                    | experimental section; needs sambamba                                                                                                      |
| preseq                   | preseq\_active                                                                    | sequencing saturation estimate                                                                                                            |
| umi\_tools               | umi\_tools\_active, pattern, pattern2                                             | regex patterns need single AND double quotes in YAML                                                                                      |
| blast                    | blast\_active, blast\_db, blast\_script\_dir                                      | quick species proxy; script dir must contain the subsampling/filter scripts                                                               |
| kraken2                  | kraken2\_active, kraken2\_database                                                | shared DB; jobs are slow to start (\~190 GB RAM)                                                                                          |
| kaiju                    | kaiju\_active, kaiju\_db, kaiju\_nodes\_dmp, kaiju\_names\_dmp                    | metagenomic fallback                                                                                                                      |
| jellyfish                | jellyfish\_active                                                                 | 21-mers; for whole-genome assays                                                                                                          |
| biobloom                 | biobloom\_active, biobloom\_filter\_dir, biobloom\_filters                        | one filter per species; filter dir needs trailing /                                                                                       |
| sourmash                 | sourmash\_active, sourmash\_eukaryote\_ref, sourmash\_microbial\_ref              | k-mer classification                                                                                                                      |
| sortmerna                | sortmerna\_active, sortmerna\_reference\_list                                     | rRNA counting/removal; needs cutadapt and/or umi\_tools                                                                                   |
| crypt4gh                 | crypt4gh\_active, crypt4gh\_own\_private\_key, crypt4gh\_client\_public\_key\_dir | human data; public-key dir needs trailing /                                                                                               |
| transfer                 | transfer\_active, transfer\_to, transfer\_as                                      | login-less connection required; transfer\_to needs trailing /                                                                             |
| gocryptfs                | gocryptfs\_active, gocryptfs\_version                                             | only when absolutely needed, in a screen session (password typed manually)                                                                |
| others                   | (empty placeholder section)                                                       |                                                                                                                                           |
| software\_versions\_file |                                                                                   | TSV of tool versions; update manually (see envs/update\_software\_env\_list.sh)                                                           |


---

## 14. rMATS Group File Setup

When `mapping.rmats_on: True`, set `mapping.rmats_group_file` to a group file. The rule `rmats_create_groupfiles_from_csv` calls `scripts/generate_rmats_groups.sh` to produce `rmats/group1_bams.txt` (group A) and `rmats/group2_bams.txt` (group B) from the STAR BAM directory. `run_rmats` then consumes those two lists.

Example (`examples/example_rmats_group_file.txt`) for 7 samples in two comparison groups (A and B):

```text
*S1*.fastq.gz,*S2*.fastq.gz,*S3*.fastq.gz:*S4*.fastq.gz,*S5*.fastq.gz,*S6*.fastq.gz,*S7*.fastq.gz
```

Patterns must match the demultiplexed (or skip\_demux) sample filenames.

---

## 15. Troubleshooting &amp; Common Pitfalls

### Trailing Slash Errors

- Issue: rule failure or missing file errors.
- Fix: remove trailing slashes (`/`) from `OutputFolder`, `fastq_folder`, `genomeDir`, `qualimap_dir`, `bcl_convert_path`, `bcl2fastq2_path`. Only `biobloom_filter_dir` and `transfer_to` require trailing slashes.

### SampleSheet Naming

- Issue: pre-flight check fails or demux target cannot be resolved.
- Fix: ensure SampleSheet.csv resides inside the sequencing run directory and contains a numeric Project ID in its filename (e.g. `1399_Samplesheet.csv`). `validateProjectNum()` is called; `validateSamplesheet`/`validateOutput` are imported but only called by `pe_processing.smk` (`validateBefore`) — the demux stage relies on `validateProjectNum` alone.

### Near-silent "success" with a missing fastq list

- Issue: if `fastq_infiles_list.tx` does not exist, `getOutput()` resets the target list and prints "It seems like the demuxing run didnt finish properly...". `rule all` then reduces to little more than the `multiqc_report_complete_<projectNum>.html` target, so Snakemake can "succeed" while producing almost nothing.
- Fix: verify `untrimmed_fastq/` and `fastq_infiles_list.tx` exist before trusting a green run; check the log for the demuxing message above.

### Locked Snakemake Directory

- Issue: `Error: Directory cannot be locked.`
- Fix: run `snakemake -s rules/convert2fastq.smk --unlock` and/or `snakemake -s rules/fastqProcessing.smk --unlock` from the pipeline root (the wrappers do not pass `--unlock` through).

### Memory Starvation on Large Flowcells

- Issue: bcl-convert or kraken2 jobs killed by PBS OOM killer.
- Fix: raise `demux.demux_try1_gb` to 126+ for NovaSeq X; kraken2 needs \~178-190 GB on attempt 1.

### Qualimap strandedness is overridden

- Issue: `mapping.qualimap_parameters` is inserted into the command line BEFORE a hardcoded `-p strand-specific-forward`, so a user-supplied `-p` is overridden. Note the shipped config.yaml default is `-p strand-specific-reverse` — which is silently ignored.
- Fix: this is a code-level issue in `rule Qualimap`. Either accept the hardcoded forward strand, or patch the rule to not hardcode `-p`. Do not rely on `qualimap_parameters` to set strandedness.

### FeatureCounts -p for PE

- Issue: the `featurecounts` rule does not add `-p` automatically for PE.
- Fix: add `-p` to `mapping.other_params` (or the stranded string) when PE is active, together with the strandedness flag in `mapping.stranded`.

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

## 16. Known Code Issues &amp; Verification Notes

1. **{-PE path inconsistencies in `pe_processing.smk`** (verify/patch before relying on PE with UMIs):-}
2. {+**Stale include**: `rules/fastqProcessing.smk` still contains `include: "qcShort.smk"`, but `qcShort.smk` has been removed from the project. A fresh clone fails at parse time until this line is deleted/commented out.+}
3. **{+PE path inconsistencies in `pe_processing.smk`** (verify/patch before relying on PE with UMIs):+}
  - `umi_dedup_pe` consumes `star/{short}_pe_Aligned.sortedByCoord.out.bam`, but `align_pe` produces `star/{short}_pe_001_Aligned.sortedByCoord.out.bam` (missing `_001` in the input).
  - `umi_extract_pe` consumes `untrimmed_fastq/{short}_R1.fastq.gz` / `{short}_R2.fastq.gz` (without cutadapt), but demux naming is `*_R1_001.fastq.gz`.
  - `count_matrix_pe` reads `star/{file}/ReadsPerGene.out.tab` — with a nested `{file}` path that does not match the flat `star/` layout.
4. **{-`wait 30` in autostart**: `autostart_cronjob_bclconvert.sh` and `autostart_novaseq.sh` use `wait 30`, which does not sleep; `autostart_portable.sh` correctly uses `sleep 30`.-}
5. {-**Env mismatch in autostart**: `autostart/execute_pipeline.sh` activates `smk9`, while the wrappers/docs use `smk8`.-}
6. **{-`schemas/config.schema.yaml` is stale** (references old `bcl2fastq:`/`threads` keys) and is not wired into any Snakefile — do not use it for validation.-}
7. {-**Diamond**: only a commented-out `diamond_pe` rule exists in `pe_processing.smk`; there is no active Diamond rule and no `envs/diamond.yaml`.-}
8. {-**bash login-shell warning**: the demux rules run `chmod ago+rwx ...` (a bash-specific `chmod` alias seen on HILBERT login shells); on other systems this may fail — the original author's environment accepts it. `FastQC_untrimmed` unsets `command_not_found_handle` for the same reason.-}
9. {-**rMATS alternative**: the code notes Spladder ([https://spladder.readthedocs.io](https://spladder.readthedocs.io)) as a drop-in replacement if rMATS causes problems.-}3. **{+`wait 30` in autostart**: `autostart_cronjob_bclconvert.sh` and `autostart_novaseq.sh` use `wait 30`, which does not sleep; `autostart_portable.sh` correctly uses `sleep 30`.+}
10. {+**Env mismatch in autostart**: `autostart/execute_pipeline.sh` activates `smk9`, while the wrappers/docs use `smk8`.+}
11. **{+`schemas/config.schema.yaml` is stale** (references old `bcl2fastq:`/`threads` keys) and is not wired into any Snakefile — do not use it for validation.+}
12. {+**Diamond**: only a commented-out `diamond_pe` rule exists in `pe_processing.smk`; there is no active Diamond rule and no `envs/diamond.yaml`.+}
13. {+**bash login-shell warning**: the demux rules run `chmod ago+rwx ...` (a bash-specific `chmod` alias seen on HILBERT login shells); on other systems this may fail — the original author's environment accepts it. `FastQC_untrimmed` unsets `command_not_found_handle` for the same reason.+}
14. {+**rMATS alternative**: the code notes Spladder ([https://spladder.readthedocs.io](https://spladder.readthedocs.io)) as a drop-in replacement if rMATS causes problems.+}