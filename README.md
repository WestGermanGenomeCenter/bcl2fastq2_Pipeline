# Demultiplexing Pipeline for Illumina Machines with extensive QC


This Pipeline uses Snakemake 8 on a PBS Pro HPC Cluster to Demultiplex Illumina Sequencing Runs, either with bcl2fastq2 or the newer bclconvert.


#  Deployment
if you are HILBERT HPC user:

 `bash setup_pipeline.sh`

if not: 
 - create your snakemake profile if needed
 - create a conda env with `conda env create -n smk8 -f smk8_environment.yml`
 - `conda activate smk8 ` <- this env should then include snakemake 8.x.x


then:
 - edit the config.yaml file. This is where all the settings are for this Pipeline. Only use tools that you need!

then:

`bash runPipeline.sh`

the runPipeline scripts creates 2 files after the pipeline is done:
- a list of all files inside the output folder (ls -alth)
- a list of checksums of all files inside the output folder (sha256sum)
if you do not want / need this, you can disable this by running:

`bash runPipeline.sh --no-checksums`


when not using the HPC, use:

`bash runPipeline_local.sh`

There is also a minimal option to execute the pipeline, without the nice-to-haves of the runPipeline scripts.
This is for testing, or if you do not like the options and things added outside of the bare snakemake pipeline:

`bash runPipeline_barebone.sh`

# Examples
The Pipeline is split into two main Parts:
- the demultiplexing (rules/convert2fastq.smk)
- all other things (fastqProcessing.smk + all other .smk files in rules/)

For example DAGs and Reports, check the examples/ folder inside this repo

The demultiplexing part creates its own DAG, in a minimal setting this looks like:
![alt text](convert_dag.png)


Possible configuration DAG for the second part of the Pipeline:
![alt text](dag_fastqProcessing.png)


If enabling less and working with pe data, the DAG can also look like:
![alt text](small_dag.png)

Only enabling trimming with cutadapt given pe data yields this dag:
![alt text](mini_dag.png)


# General remarks: The Pipeline 
- is setup to always use the config.yaml in the root folder
- expects bcl2fastq2/bclconvert to be somewhere as executable (needs to be set in config.yaml)
- always does some basic steps like sha256sums / fastqc + multiqc reports
- is best to use one pipeline run per project and one project per sequencing run. That means if you have 2 or more projects in one flowcell/ sequencing run, you can eiter 
    -- split the samplesheet by project, run the pipeline once for each project but get a big undetermined.fastq.gz file or
    -- run once, but then mapping and all other after-mapping qc steps can only be made with one reference genome, and not one reference genome per project
    -- in essence, split the samplesheet if you want to do mapping or even more and if you only want to do basic qc on all the included projects one samplesheet is also a solution
    -- you can also demultiplex once, move/ copy the .fastq files and then use the skip_demux option to run the pipeline with one set of parameters for each project. This is the most ressource-efficient way for now.
- can deal with .fastq files as input (if demultiplexing already happened), for this enable the skip_demux part in the config.yaml
- can deal with paired and and single read sequencing setups. With pe setups the naming scheme is strict (see config.yaml). Index read files are ignored after initial demultiplexing.
- if qc steps are enabled that require trimming of adapters, but trimming is disabled the pipeline will trimm regardless. This is to ensure the proper input quality of data for each step.
- The major focus here is demultiplexing of MiSeq/Nextseq2000/NovaseqX and do RNASeq analysis until transcript count tables are created
- can deal with UMIs, even in paired end setups (check config.yaml)
- collects qc of raw sequencing output in a multiqc file, and output of all analysis steps in a multiqc_report_complete file.
- general deletion of data does not happen. You need to cleanup on your own.
- expects a RunFolder as input that includes a SampleSheet. The SampleSheet needs to have a number (ProjectID) in the Filename. The RunFolder is where the demultiplexing then starts. The Samplesheet cannot be outside the root folder of the RunFolder created by Illuminas machines. See config.yaml for an example.
- SampleSheet and OutputFolder are the only 2 required inputs in the config.yaml, all other steps can be disabled.
- currently uses RSeQC as post-mapping QC for RNASeq runs, but only the for us useful modules of RSeQC. Alternatively Qualimap can be used. Here the installation is not correctly handled by conda (needs to be installed with the --no-deps parameter on the HPC to avoid a failing install). Qualimap needs to be available as a executable for now. This might change in the future (check the config.yaml)
- includes many QC steps to identify contamination, from most reliable to least: Kraken , BioBloom, BLAST, Diamond. Please be aware that either results are dependent on Quality of the input and the Database used.
- automatically scales up the ressources for each job after a failed attempt. after 3 attempts it gives up - you can configure this in your snakemake profile.
- many results of analyses are shown in the multiqc report complete file, but a lot of tools also put out more data then multiqc catches. Check the output folder after each run for more!

# Setup steps needed
- install the correct star version and create a reference genome, include a annotation .gtf file that you can use later in the pipeline during this step. Reference genomes can be downloaded from ncbi, ensembl or anywhere else. Download .fa(sta) and corresponding .gtl file, then do `star --runMode GenomeGenerate ` . the resulting Folder and gtf file can then be used in the config.yaml
- download binaries for bcl2fastq2 and/or bclconvert. This is sadly only available on the illumina website with a login. The complete path to the executables can then be filled in to the config.yaml
- create / download kraken2 / biobloomtools / blast / diamond / sortmerna databases. Most tools have downloadable databases available. Best is to start there and then adjust to your needs. The file paths can then be filled into the config.yaml. For biobloomtools, each organism needs one .bf file that need to be chained one after the other. All biobloomfilter .bf files need to be in the same, biobloom_filter_dir that also needs to be filled into the config.yaml file. Sortmerna needs also a special way: the reference list needs to be split by --ref for each specific file like so:
config.yaml: sortmerna, sortmerna_reference_list: '--ref /path/to/file1 --ref /path/to/file2 --ref /path/to/file3' ... and so on 
- all software should be handed by snakemake / conda except: 
    - bcl2fastq2
    - bclconvert
    - Qualimap (since the conda install is broken on our HPC)
- blast_script_dir, interop_script_path are the folder scripts/ in this repo, so filling in only "scripts/" should also work for both
- start with only minimal options enabled, then enable more as you get familiar with how to fill in / use the config.yaml file. If you need inspiration, check the examples/ folder included in this dir.
- kraken2 is okay with UMIs, biobloom, blast and diamond not. This is why kraken2 is started before the UMIs are potentially moved to the .fastq header, giving you a sooner glimpse of potential contamination

# Tools included / steps available
- demultiplex with bcl2fastq2 or bclconvert (needs to be available as executable)
- FastQC, Fastp, sha256sums and MultiQC of raw demultiplexing output 
- ability to automatically reverse-complement the samplesheet info for second index read
- can automatically skip demultiplexing if .fastq.gz files are already somewhere
## Second part
- Illumina Interop and demux stats in second MultiQC report (..report_complete.html)
- cutadapt
- umi_tools
- sortmerna
- FastQC, MultiQC, sha256sums of trimmed/Umis moved to header / rRNA removed .fastq files
- contamination detection with:
    - Kraken2
    - BioBloom
    - BLAST
    - Diamond
    - Kaiju
- Jellyfish for k-mer distribution checks
- encryption of .fastq files with crypt4gh, encryption of all output with gocryptfs also possible
- automatic transport of all output data to other server/ scp adress
## RNA-Seq Workflow in the second part
- Mapping with STAR
- RseQC and samtools-plots for mapping output QC
- preseq for sequencing saturation estimation
- stringtie
- rmats (very simple groups file needs to be made)
- featurecounts for final transcript count table

- can work with pe/se data, can also work with dual UMI setup
- ends in a final MultiQC report that includes a overview of all analysis done and a list of used software
- copies samplesheet / config.yaml / software versions into output folder
- includes creating snakemake rulegraphs and reports for each run
- includes complete filelist and checksums of all created files as part of runPipeline(_local).sh. This can be switched off.


# The runPipeline scripts
The snakemake pipeline itself is executed with the runPipeline(local).sh script.
This adds useful information like DAGs, snakemake reports and checksums.
If needed, this can be ignored and only the snakemake itself can be executed.
requirements for the runPipeline.sh script to work properly:
  + is executed in an interactive job that itself runs inside a screen
  + a conda environment "smk8" needs to be available to the executing user with snakemake 8.x installed in that environment
  + Execution happens on the local HPC Hilbert
  + the pipeline is also tested with snakemake9+, results do not change


### The rmats group file
There is an example rmats group file in the examples/  folder: example_rmats_group_file.txt - this fits for 7 samples in 2 groups (A and B), for now single read only:

`*S1*.fastq.gz,
*S2*.fastq.gz,
*S3*.fastq.gz,
*S4*.fastq.gz,
*S5*.fastq.gz,
*S6*.fastq.gz,
*S7*.fastq.gz`

need to be files that come from demultiplexing (or skip-demux step).



_______________________________________________________________________________________________________
_______________________________________________________________________________________________________
_______________________________________________________________________________________________________
_______________________________________________________________________________________________________
_______________________________________________________________________________________________________



# A more extensive documentation, now with chapters

## Table of Contents
- [Overview](#overview)
- [Key Features](#key-features)
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites](#prerequisites)
- [Installation & Setup](#installation--setup)
- [Configuration Guide](#configuration-guide)
- [Running the Pipeline](#running-the-pipeline)
- [Pipeline Stages](#pipeline-stages)
- [Example Workflows](#example-workflows)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

---

## Overview

This Snakemake-based pipeline provides automated demultiplexing of Illumina sequencing data (BCL to FASTQ conversion) with extensive quality control 
optional downstream RNA-Seq analysis. Built for high-throughput environments, it supports both bcl2fastq2 and the newer bclconvert tools.

**Target Users:** Bioinformaticians, core facility operators, and researchers working with Illumina sequencing platforms (MiSeq, NextSeq 2000, NovaSe

**Current Version:** Snakemake 8.x compatible (also tested with Snakemake 9+)

---

## Key Features

### Core Capabilities
- **Dual demultiplexer support**: bcl2fastq2 or bclconvert
- **Automatic quality control**: FastQC, MultiQC, and comprehensive QC metrics
- **Flexible input**: Works with BCL files or pre-demultiplexed FASTQ files
- **Read type agnostic**: Handles single-end (SE) and paired-end (PE) data
- **UMI handling**: Supports unique molecular identifiers, including dual-UMI setups
- **Contamination detection**: Kraken2, BioBloom, BLAST, Diamond, Kaiju
- **RNA-Seq workflow**: Complete pipeline from FASTQ to count tables

### Advanced Features
- **Data security**: File encryption with crypt4gh or gocryptfs
- **Comprehensive reporting**: Automated MultiQC reports with detailed metrics
- **Automatic retry**: Smart resource scaling after failed job attempts
- **Data transfer**: Automated output delivery via SCP
- **Reproducibility**: Version tracking, checksums, and complete audit trails

### Automation
- **cronjob**: setup a cronjob for auto-execution of a pipeline
---

## Pipeline Architecture

The pipeline consists of two main stages:

### Stage 1: Demultiplexing
Converts BCL files to FASTQ format with quality control.

![Demultiplexing DAG](convert_dag.png)
*Figure 1: Minimal demultiplexing workflow showing BCL conversion and initial QC steps*

### Stage 2: Analysis & QC
Performs trimming, contamination checks, mapping, and transcript quantification.

![Full Pipeline DAG](dag_fastqProcessing.png)
*Figure 2: Complete pipeline with all optional modules enabled (paired-end data)*

#### Simplified Workflows

**Paired-end with fewer modules enabled:**

![Simplified PE DAG](small_dag.png)
*Figure 3: Streamlined paired-end workflow*

**Minimal trimming-only workflow:**

![Minimal DAG](mini_dag.png)
*Figure 4: Cutadapt-only configuration for paired-end data*

---

## Prerequisites

### Required Software
1. **Conda/Mamba** (for environment management)
2. **Snakemake 8.x** or higher
3. **bcl2fastq2** OR **bclconvert** (download from Illumina - requires login)

### Optional Software (depending on enabled features)
- **STAR** (for RNA-Seq mapping)
- **Qualimap** (alternative post-mapping QC - note: conda install requires `--no-deps` flag)
- **Kraken2** (contamination detection)
- **BioBloom Tools** (contamination detection)
- **BLAST+** (contamination detection)
- **Diamond** (contamination detection)

### System Requirements
- **Storage**: Large disk space for BCL files, FASTQ outputs, and intermediate files
- **Memory**: Varies by sample size; STAR indexing can require 30+ GB RAM
- **Compute**: Multi-core system or HPC cluster recommended


---

## Installation & Setup

### Quick Start (HILBERT HPC Users)

```bash
git clone https://github.com/WestGermanGenomeCenter/bcl2fastq2_Pipeline.git
cd bcl2fastq2_Pipeline
bash setup_pipeline.sh
```

### Manual Setup (Other Systems)

#### 1. Clone the Repository
```bash
git clone https://github.com/WestGermanGenomeCenter/bcl2fastq2_Pipeline.git
cd bcl2fastq2_Pipeline
```

#### 2. Create Snakemake Environment
```bash
conda env create -n smk8 -f smk8_environment.yml
conda activate smk8
```

Verify installation:
```bash
snakemake --version  # Should show 8.x.x or higher
```

#### 3. Install bcl2fastq2/bclconvert
Download binaries from [Illumina Support](https://support.illumina.com/) (requires account).

```bash
# Example installation (adjust paths as needed)
chmod +x bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
sudo rpm -ivh bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm

# Verify installation
bcl2fastq --version
```

Update `config.yaml` with the full path to the executable.

#### 4. Set Up Reference Genomes (for RNA-Seq)

**Download reference files:**
```bash
# Example: Human genome from ENSEMBL
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

gunzip *.gz
```

**Generate STAR index:**
```bash
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.110.gtf \
     --runThreadN 16 \
     --sjdbOverhang 99  # adjust based on read length
```

#### 5. Set Up Contamination Databases (Optional)

**Kraken2:**
```bash
# Download pre-built standard database
kraken2-build --standard --db kraken2_standard_db

# Or build custom database
kraken2-build --download-taxonomy --db custom_db
kraken2-build --download-library bacteria --db custom_db
kraken2-build --build --db custom_db
```

**BioBloom:**
```bash
# Download or create .bf files for each organism
# Place all .bf files in a single directory
# In config.yaml, chain them: filter1.bf filter2.bf filter3.bf
```

**BLAST:**
```bash
# Download NT database
update_blastdb.pl --decompress nt

# Or build custom database
makeblastdb -in sequences.fasta -dbtype nucl -out custom_blast_db
```

**Diamond:**
```bash
# Download NCBI nr database
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz -d nr
```

**SortMeRNA:**
```bash
# Download rRNA databases
wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
tar -xzf database.tar.gz
```

---

## Configuration Guide

All pipeline behavior is controlled through `config.yaml`. This file must be edited before running the pipeline.

### Required Settings

#### Minimal Configuration
```yaml
##########################
# REQUIRED SETTINGS
##########################

# Path to sample sheet (must be in RunFolder or subdirectory)
SampleSheet: "/path/to/RunFolder/SampleSheet_12345.csv"

# Output directory for all results
OutputFolder: "/path/to/output"

# Demultiplexer choice
use_bcl2fastq: true    # true = bcl2fastq2, false = bclconvert
bcl2fastq_path: "/usr/local/bin/bcl2fastq"
bclconvert_path: "/usr/local/bin/bcl-convert"


```

### Sample Sheet Requirements

**Important Notes:**
- Sample sheet **must** contain a project ID number in the filename: `SampleSheet_12345.csv`
- Sample sheet **must** be inside the RunFolder or a subdirectory
- For paired-end data, file naming must follow pattern: `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`
- Index reads (I1/I2) are automatically ignored after demultiplexing

**Example Sample Sheet (IEM format):**
```csv
[Header]
IEMFileVersion,4
Investigator Name,John Doe
Experiment Name,RNASeq_Experiment
Date,2024-02-04
Workflow,GenerateFASTQ
Application,RNA-Seq
Assay,TruSeq
Description,Paired-end RNA sequencing
Chemistry,Default

[Reads]
151
151

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,Plate1,A01,D701,ATTACTCG,D501,TATAGCCT,12345,Control_rep1
Sample2,Sample2,Plate1,A02,D702,TCCGGAGA,D501,TATAGCCT,12345,Control_rep2
Sample3,Sample3,Plate1,A03,D703,CGCTCATT,D502,ATAGAGGC,12345,Treatment_rep1
```

### Optional Features

#### Skip Demultiplexing (Use Pre-existing FASTQ Files)
```yaml
skip_demux:
  skip_demux_active: False
  fastq_folder: "/gpfs/project/projects/bmfz_gtl/projects/1162_Konkat/Konkat+Trimmed/fastq" #<- inside that folder are the .fastq.gz files expected then (untrimmed) # / not allowed at the end
```

#### Adapter Trimming
```yaml

cutadapt: 
  cutadapt_active: True
  adapters: 
    - AGATCGGAAGAGCACACGTCTGAACTCCAGTCA    
    - CTGTCTCTTATACACATCT
    - AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    - ACACTCTTTCCCTACACGACGCTCTTCCGATCT

  adapter_type: 
    - b # b for can be at both ends, a and g are for adapters expected on specific ends. b makes it more general, filtering a little better
    - b
    - b
    - b
  other_params: "--minimum-length 25 --quality-cutoff 5 -n 2" # very conservative, should work with any assay.

```

#### UMI Handling
```yaml
# umi_tools is used in two steps to first move the umis from a defined position of the trimmed read to the read header inside the .fastq.gz files
umi_tools: # for this to work, cutadapt needs to be active 
  umi_tools_active: True
 # pattern: '".+(?P<umi_1>.{8})(?P<discard_1>G{3,5})"' # TaKaRa: 8nt umi at 5' + 3-5Gs at the end  of the read {s<=2}
  pattern: '"(?P<umi_1>.{6})(?P<discard_1>TATA){s<=1}.+"' # Lexogen: 6nt umi + tata = 10nt UMI , allowing one mismatch in the spacer 
 # pattern: '"(?P<umi_1>.{6})(?P<discard_1>TATA).+"' # Lexogen: 6nt umi + tata = 10nt UMI , always with single and double quotes! 
  pattern2: '""' # umi pattern for the second read if needed, only possible in paired-end sequencing setups

```

#### rRNA Removal (SortMeRNA)
```yaml


# sortmerna counts rRNA from input files, splitting the .fastq in the process.  # needs cutadapt and/ or umi_tools to be active
sortmerna:
  sortmerna_active: True
  sortmerna_reference_list: '--ref /db/silva-bac-16s.fasta --ref /db/silva-arc-16s.fasta --ref /db/silva-euk-18s.fasta'
```

#### Contamination Detection
```yaml
# Kraken2 (recommended - works with UMIs)
kraken2:
  kraken2_active: True   #
  kraken2_database: "/gpfs/project/databases/Kraken2-2022-09-28/kraken_db_plus" # also possible: /gpfs/project/databases/Kraken2_refseq+hfmv # also possible: /gpfs/project/projects/bmfz_gtl/software/Kraken

# BioBloom (fast, no UMI support)
biobloom:
  biobloom_active: True
  biobloom_filter_dir: "/path/to/bloom_filters"
  biobloom_filters: "ecoli.bf human.bf mouse.bf"

# BLAST (slow, thorough)
blast: # for this to work, cutadapt needs to be active
  blast_active: False    #
  blast_db: "/path/to/nt"
  blast_script_dir: "scripts/"

# Diamond (fast protein-level search)
diamond:
  diamond_active: False  #
  diamond_ref_file: "/gpfs/project/projects/bmfz_gtl/databases/diamond/swissprot_protein_db.dmnd"

# Kaiju (fast k-mer based)
# kaiju is for metagenomic data only, if kraken/biobloom/blast/diamond give not a clear result, try this here
kaiju:
  kaiju_active: False
  kaiju_db: "/gpfs/project/projects/bmfz_gtl/software/kaiju/kaiju_db_refseq_nr.fmi" 
  kaiju_nodes_dmp: "/gpfs/project/projects/bmfz_gtl/software/kaiju/nodes.dmp"
  kaiju_names_dmp: "/gpfs/project/projects/bmfz_gtl/software/kaiju/names.dmp"

```

#### RNA-Seq Mapping & Quantification
```yaml

mapping: # transcriptome mapping happens here with STAR
  mapping_active: True
#  gtf_file: "/gpfs/project/projects/bmfz_gtl/software/reference_genomes/hg38/star_2710b_hg38_ensembl/Homo_sapiens.GRCh38.111.gtf" 
#  genomeDir: "/gpfs/project/projects/bmfz_gtl/software/reference_genomes/hg38/star_2710b_hg38_ensembl" # / not allowed at the end
  gtf_file: "/gpfs/project/projects/bmfz_gtl/software/reference_genomes/mouse_ensembl/Mus_musculus.GRCm39.111.gtf"
  genomeDir: "/gpfs/project/projects/bmfz_gtl/software/reference_genomes/mouse_ensembl"
  other_params: ""
  stranded: "-s 2 " # featurecounts parameter -s , adapt this to paired ("-p") if paired- end sequencing, else "-s 1" for forward stranded, "-s 0" for unstranded and "-s 2" for reverse strandedness"
  stringtie_on: True # if mapping is disabled, but one of these options is enabled - mapping will happen anyway
  rseqc_active: False #  qc of mapping, somewhat redundant to qualimap:	rseqc and qualimap do the same,	just in	different ways:	for now	please just enable one of the two, i would prefer qualimap
  qualimap_on: True # for this to work you need to manually create Qualimap env with qualimap=2.3 [--no-deps] installed, or use the local installment:
  qualimap_parameters: "-p strand-specific-reverse" # for qualimap to work properly, you need to set the strandedness, otherwise it assumes unstrandedness: -p strand-specific-reverse strand-specific-forward, strand-specific-reverse or non-strand-specific
  qualimap_dir: "/gpfs/project/projects/bmfz_gtl/software/qualimap/qualimap_2_3/qualimap_v2.3" # an executable of qualimap needs to be in that folder, / not allowed at the end
  rmats_on: True
  rmats_read_length: 100 # length of input reads. sequenced length minus umi, minus spacer, minus adapter.
  rmats_paired: "single --libType fr-secondstrand" # Type of read used in the analysis: either "paired" for paired-end data or "single" for single-end data. Default: paired +   --libType {fr-unstranded,fr-firststrand,fr-secondstrand}                         Library type. Use fr-firststrand or fr-secondstrand for strand-specific data. Only relevant to the prep step, not the post step. Default: fr-unstranded
  rmats_group_file: "/gpfs/project/projects/bmfz_gtl/projects/1399/star/groups.txt"

```

#### Encryption
```yaml
# Per-file encryption (crypt4gh)
# crypt4gh encryption of umi_extract/cutadapt output files. needs both files to be present on the same filesystem as files # use only for human data, only transfer the files in the /encrypted_fastqs/ folder 
crypt4gh:
  crypt4gh_active: False
  crypt4gh_own_private_key: "/home/daric102/daric102_nop.sec"
  crypt4gh_client_public_key_dir: "/gpfs/project/projects/bmfz_gtl/software/c4gh_public_keys/" # recipient and all gtl public keys that can decrypt later / at the end needed 


# use this only when absolutely needed, in a screen session! pw needs to be manually typed in after all analysis are done.
gocryptfs:
  gocryptfs_active: False
  gocryptfs_version: gocryptfs/1.7.1
```

#### Data Transfer
```yaml
transfer:
  transfer_active: False
  transfer_to: "136.69.264.68:/mnt/data/pipelines/" # / at the end needed 
  transfer_as: "daniel@"

```

### Full Configuration Examples

Check the `examples/` directory for complete configuration files:
- `config_minimal.yaml` - Demultiplexing + basic QC only
- `config_rnaseq_full.yaml` - Complete RNA-Seq pipeline
- `config_contamination_check.yaml` - Focus on contamination detection

---

## Running the Pipeline

### Standard Execution (HPC with PBS Pro)

```bash
# Activate conda environment
conda activate smk8

# Run with all bells and whistles (recommended)
bash runPipeline.sh

# Run without checksum generation (faster)
bash runPipeline.sh --no-checksums

# Local execution (no cluster)
bash runPipeline_local.sh

# Bare-bones execution (testing only)
bash runPipeline_barebone.sh
```

### What runPipeline.sh Does

The wrapper script adds useful features beyond basic Snakemake execution:

1. **Generates workflow visualizations:**
   - DAG graph of all pipeline steps
   - Rule graph showing dependencies
   
2. **Creates comprehensive reports:**
   - Snakemake HTML report with timing, resource usage
   - Complete file listing (`ls -alth`)
   - SHA256 checksums of all outputs

3. **Environment documentation:**
   - Copies `config.yaml` to output directory
   - Exports software versions

4. **Safety checks:**
   - Validates config file syntax
   - Checks for required executables

**Requirements for runPipeline.sh:**
- Must be executed in an **interactive job** inside a **screen** session
- Conda environment `smk8` must exist with Snakemake 8.x

**Example HPC submission:**
```bash


# Inside the HPC shell, start screen
screen -dRR pipeline_run

# Inside screen, start interactive job
qsub -I -l select=1:ncpus=2:mem=8gb,walltime=48:00:00


# run
bash runPipeline.sh

# Detach from screen: Ctrl+A, then D
# Reattach later: screen -r pipeline_run
```

### Direct Snakemake Execution

For advanced users who want full control:

```bash
snakemake \
  --snakefile Snakefile \
  --configfile config.yaml \
  --cores 32 \
  --use-conda \
  --cluster "qsub -l select=1:ncpus={threads}:mem={resources.mem_mb}mb,walltime={resources.walltime}" \
  --jobs 50 \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds
```

### Dry Run (Test Without Execution)

```bash
snakemake --snakefile Snakefile --configfile config.yaml --dry-run
```

### Unlock After Crash

If pipeline crashes and locks the working directory:

```bash
snakemake --snakefile Snakefile --unlock
```

---

## Automation of start
The pipeline can be made to autostart, if:
  - there is access to a crontab
  - the pc /vm  where the crontab is active has direct access to the illumina outputfolder

To automate the pipeline:
  - edit *both* files in the dir /autostart : paths, usernames and job parameters need to fit to your environment
  - add a cronjob executing the autostart_cronjob_bclconvert.sh

#### How the script works:
It scans a folder of folders (The folder where all illumina runfolders are) for runs (folders)that:
  - have a CopyComplete.txt in a folder. (The folder can also contain multiple SampleSheets, only the one in the config.yaml gets used).
  - have a config.yaml in **the same folder**. This needs to be a copy of the config.yaml from the root of this pipeline. aldkf234_config.yaml also works, as does config_oiquzreoizewr123213.yaml. The config.yaml gets copied from the runfolder to the root of this pipeline. In this config.yaml the SampleSheet and the outputfolder need to be configured. This will be the config.yaml that gets used during pipeline execution, so enable/ disable / set all things in the config.yaml according to your needs
  - have NOT a flag file  (.pipeline_launched) in **the same folder** that gets put there on the first automatic execution of the pipeline. Do a 'ls -alth' to find the flag
  - as indirect requirement: the pipeline needs the samplesheet in the illumina runfolder that it should demultiplex, of course only if skip_demux is not active
The script also only starts the pipeline if a job of the same pipeline autostart (its called "demx_pipe_auto") is not already / still running. it checks with qstat for the name of the pipeline job that it would execute.

You can also run the starter script manually to check if your setup is ok: 'bash autostart_cronjob_bclconvert.sh'

Then the scripts submits the jobscript in the same folder: autostart/execute_pipeline.sh

To summarize:
- edit the bash scripts according to your setup (folders, username)
- add to your crontab

To make a pipeline auto-run:
- a runfolder needs to have: CopyComplete.txt, config.yaml, no flagfile (that gets put there on the first automatic execution)
- the config.yaml should have your desired pipeline parameters, as these are the settings the pipeline will be executed with

Pitfalls:
- the pipeline will be executed if you analyzed the data already with the pipeline manually before - as long as the previously stated files are as expected - because the manual pipeline execution does not create the flagfile. 
- if you use multiple machines, you need to setup multiple instances of the cronjob. For this the fastest path is to copy the autostart_cronjob_bclconvert.sh and edit the WATCH_DIR. Then add another cronjob with the copied file that watches the other dir. 


## Pipeline Stages

### Stage 1: Demultiplexing

**Input:** Illumina RunFolder with BCL files  
**Output:** FASTQ files per sample

**Steps:**
1. **BCL Conversion** - bcl2fastq2 or bclconvert
2. **InterOp Parsing** - Extract run metrics
3. **Initial QC** - FastQC on raw FASTQ files
4. **Checksum Generation** - SHA256 for data integrity
5. **MultiQC Report** - Aggregated QC metrics

**Key Outputs:**
- `*.fastq.gz` - Demultiplexed sequence files
- `Undetermined_*.fastq.gz` - Reads that didn't match any index
- `multiqc_report.html` - Quality control summary
- `demux_stats/` - Demultiplexing statistics

### Stage 2: Processing & Analysis

#### 2A. Pre-Processing

**UMI Extraction (if enabled):**
- Moves UMI sequences to read headers
- Compatible with paired-end and single-end
- Dual UMI support for R1 + R2

**Adapter Trimming (cutadapt):**
- Removes sequencing adapters
- Quality trimming
- Length filtering
- Generates trimming statistics

**rRNA Removal (SortMeRNA):**
- Filters ribosomal RNA reads
- Saves non-rRNA reads for downstream analysis
- Quantifies rRNA contamination

**Post-trim QC:**
- FastQC on processed files
- MultiQC aggregation

#### 2B. Contamination Detection

**Kraken2** (Recommended - fast, UMI-compatible):
- Taxonomic classification
- Species-level identification
- Contamination percentage

**BioBloom** (Very fast, k-mer based):
- Screens against known organisms
- No UMI support - run before UMI extraction
- Low memory footprint

**BLAST** (Thorough but slow):
- Homology search against NT database
- Detailed alignment information
- Best for unknown contaminants

**Diamond** (Fast protein-level):
- Translated search against NR database
- Good for viral contamination
- Much faster than BLASTX

**Kaiju** (K-mer based, protein-level):
- Fast taxonomic classification
- Good for metagenomics
- Low memory usage

#### 2C. RNA-Seq Workflow

**Mapping (STAR):**
- Splice-aware alignment
- Generates BAM files sorted by coordinate

**Post-Mapping QC:**

*RSeQC Suite:*
- `bam_stat.py` - Basic BAM statistics
- `read_distribution.py` - Reads per genomic feature
- `read_duplication.py` - PCR duplication rates
- `junction_annotation.py` - Splice junction analysis
- `infer_experiment.py` - Strand specificity detection

*Qualimap (Alternative):*
- Comprehensive BAM QC
- Gene coverage profiles
- Insert size distributions

*SAMtools:*
- Alignment statistics
- Coverage depth
- Flagstat metrics

**Sequencing Saturation (Preseq):**
- Predicts library complexity
- Estimates future sequencing yield
- Guides additional sequencing decisions

**UMI Deduplication:**
- `umi_tools dedup` - Removes PCR duplicates using UMI information
- Separate handling for single vs. dual UMI

**Transcript Quantification:**

*StringTie:*
- Transcript assembly
- Novel transcript discovery
- Abundance estimation (TPM, FPKM)

*featureCounts:*
- Gene-level read counting
- Output compatible with DESeq2/edgeR
- Handles multi-mapping reads

**Alternative Splicing (rMATS):**
- Detects differential splicing events
- Compares experimental groups
- Requires group file (see `examples/example_rmats_group_file.txt`)

**Final MultiQC:**
- Aggregates all analysis steps
- `multiqc_report_complete.html`

---

## Example Workflows

### Example 1: Basic Demultiplexing Only

**Use Case:** Just need FASTQ files, no downstream analysis

**config.yaml:**
```yaml
SampleSheet: "/data/RunFolder_20240204/SampleSheet_001.csv"
OutputFolder: "/results/Project_001"
use_bcl2fastq: true
bcl2fastq_path: "/usr/local/bin/bcl2fastq"

# Disable all optional features

```

**Expected Runtime:** 2-6 hours (depends on flowcell size)

**Outputs:**
- FASTQ files
- FastQC reports
- MultiQC summary
- Checksums



---

## Troubleshooting

### Common Issues

#### Issue: "bcl2fastq not found"

**Cause:** bcl2fastq/bclconvert not in PATH or incorrect path in config

**Solution:**
```bash
# Find the executable
which bcl2fastq
# or
find /usr -name bcl2fastq 2>/dev/null

# Update config.yaml with full path
bcl2fastq_path: "/usr/local/bin/bcl2fastq"
```

---

#### Issue: "Sample sheet validation failed"

**Cause:** Sample sheet format errors or missing project ID

**Common fixes:**
- Ensure filename contains a number: `SampleSheet_12345.csv`
- Check for Windows line endings: `dos2unix SampleSheet.csv`
- Verify all required columns are present
- Check for trailing commas or extra columns

**Validation command:**
```bash
# Check line endings
file SampleSheet.csv

# View in hex to see hidden characters
xxd SampleSheet.csv | less
```

---


---

#### Issue: Kraken2 database errors

**Cause:** Database not properly built or path incorrect

**Solution:**
```bash
# Verify database integrity
kraken2-inspect --db /path/to/kraken_db | head

# Rebuild if necessary
kraken2-build --build --db /path/to/kraken_db
```

---

#### Issue: Pipeline stalls on cluster

**Cause:** Job submission limits or cluster queue issues

**Solution:**
```bash
# Check cluster job status
qstat -u $USER

# Check Snakemake job submission
tail -f .snakemake/log/*.log

# Reduce concurrent jobs
snakemake --jobs 10  # instead of 50
```

---

#### Issue: UMI extraction fails

**Cause:** Incorrect UMI pattern or read structure

**Solution:**
```bash
# Verify read structure with first read
zcat fastq/Sample1_R1_001.fastq.gz | head -2

# Output example:
# @M00123:1:000000000-ABCDE:1:1101:10000:1000 1:N:0:ATCACG
# NNNNNNNNATCGATCGATCG...
#  ^^^^^^^^ UMI

# Adjust UMI pattern in config
umi_pattern: "NNNNNNNN"  # 8 N's for 8bp UMI
```

---

### Performance Optimization







#### Reduce Disk I/O
```yaml


# Disable unnecessary QC
run_blast: false  # BLAST is very slow
run_diamond: false
```

#### Parallelize Across Samples
```bash
# Increase concurrent jobs (if cluster allows)
snakemake --jobs 100
```

#### Use Faster Storage
- Place RunFolder on SSD if possible
- Use local scratch for intermediate files
- Output to high-performance filesystem

---

### Debugging Tips

#### Dry Run Before Execution
```bash
snakemake --dry-run --printshellcmds
```

#### Check Specific Rule
```bash
snakemake --until fastqc_raw --dry-run
```

#### Examine Log Files
```bash
# Snakemake logs
tail -f .snakemake/log/*.log

# Individual job logs (if using cluster)
tail -f logs/demultiplex/*.log
```

#### Visualize Workflow
```bash
# Generate DAG
snakemake --dag | dot -Tpng > dag.png

# Generate rule graph
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

#### Force Re-run Specific Rule
```bash
snakemake --forcerun fastqc_raw
```

#### Clean Failed Jobs
```bash
# Remove incomplete outputs
snakemake --delete-all-output

# Or just unlock
snakemake --unlock
```

---



## FAQ

**Q: Can I use this pipeline without a cluster?**  
A: Yes, use `runPipeline_local.sh` for single-machine execution. Adjust `--cores` based on your system.

**Q: How do I handle index hopping?**  
A: Check `Undetermined_*.fastq.gz` files. High undetermined reads suggest index hopping or incorrect sample sheet. Consider increasing 
`--minimum-trimmed-read-length` in bcl2fastq or using `--filter-barcode` in bclconvert.

**Q: Can I demultiplex NovaSeq runs with this?**  
A: Yes, both NovaSeq 6000 and NovaSeq X are supported. Use `bclconvert` for NovaSeq X (bcl2fastq2 doesn't support it).

**Q: What if my samples have different read lengths?**  
A: bcl2fastq2/bclconvert handle variable lengths automatically. Make sure your SampleSheet `[Reads]` section specifies the maximum length.

**Q: How much disk space do I need?**  
A: Rule of thumb: 3-4x the size of BCL files (for BCL + FASTQ + intermediate files). Use `du -sh RunFolder/` to estimate.

**Q: Can I run multiple pipelines simultaneously?**  
A: Yes, copy/paste the complete pipeline to a different second folder just to be sure. Also ensure each uses a different `OutputFolder` and avoid cluster job limits.

**Q: How do I cite this pipeline?**  
A: Include the GitHub repository URL and the specific commit hash.

**Q: Is there a Docker/Singularity container?**  
A: No, but conda environments make installation straightforward on most systems.

**Q: What about 10x Genomics data?**  
A: This pipeline is for bulk RNA-Seq. For 10x single-cell data, use Cell Ranger instead.

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Commit your changes: `git commit -m 'Add amazing feature'`
4. Push to branch: `git push origin feature/amazing-feature`
5. Open a Pull Request

**Guidelines:**
- Follow existing code style
- Update documentation
- Use descriptive commit messages

**Bug Reports:**
- Use GitHub Issues
- Include config.yaml (redact sensitive paths)
- Attach relevant log files
- Describe expected vs. actual behavior

---

## Citation

If you use this pipeline in your research, please cite:

```
West German Genome Center bcl2fastq2_Pipeline (2026)
https://github.com/WestGermanGenomeCenter/bcl2fastq2_Pipeline
```

Additionally, please cite the tools used in your specific workflow. References are provided in the final MultiQC report.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

**West German Genome Center**  
Website: [https://www.wggc.de](https://www.wggc.de)  
GitHub: [https://github.com/WestGermanGenomeCenter](https://github.com/WestGermanGenomeCenter)

For support:
- Open a GitHub Issue
- Check existing Issues for similar problems
- Consult the Troubleshooting section above

---

## Acknowledgments

This pipeline builds upon excellent open-source tools including:
- Snakemake workflow management
- bcl2fastq2 and bclconvert (Illumina)
- FastQC, MultiQC, cutadapt
- STAR aligner, RSeQC, StringTie
- Kraken2, BLAST, Diamond for contamination detection
- And many more listed in environment files

Thanks to all contributors and the bioinformatics community!

---

**Last Updated:** February 4, 2026  
**Pipeline Version:** See `git log` for commit history