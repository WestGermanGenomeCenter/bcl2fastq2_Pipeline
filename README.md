# bcl2fastq2Pipe

a snakemake based workflow for demultiplexing Illumina runs with bclconvert


automatically does our most-used application, RNASeq analysis with a lot of QC before and after mapping.


For best results, use Ensembl reference genome and annotation!


## features
- bclconvert/bcl2fastq2
- interop_plots
- checksum creation
- fastqc, fastp
- cutadapt
- umi_tools
- kraken2
- blastn
- qualimap
- rseqc
- samtools stats, coverage, plot
- stringtie
- prevent weird auto-reverse complement in read 3 if samplesheet created on the wrong Machine
- STAR as main mapper
- 2 multiqc reports: simple + all analysis inclusive
- auto-transfer of folder to other servers
- gocryptfs
- preseq

  # Structure:
  the workflow is actually 3 parts:
  1. demultiplexing, interop plot, checksums of untrimmed data
  2. fastqProcessing: all other parts, excluding
  3. RseQC: only active if rseqc_modules_active= True in config.yaml

# Usage:
git pull

edit config.yaml

bash runPipeline.sh
