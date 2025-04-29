#!/bin/bash
input_fastq_gz=$1
output_fasta=$2
zcat $1 | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' |head -n 80000 |tail -n 40000 >$2

