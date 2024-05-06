#!/bin/bash
projects_dir=/gpfs/project/projects/bmfz_gtl/projects/
head_projects=`ls -lth $projects_dir | head -n 25`

echo "type in only your project folder name. "
echo "folder MUST include .fastq.gz files"
echo "below the newest projects:"
echo "$head_projects"


read -r -p "project dir: " dir|| exit 100
mkdir -p $projects_dir/$dir/untrimmed_fastq

cd $projects_dir/$dir && shopt -s globstar; ls -f1 *.gz  | grep '.fastq.gz' | grep -v 'Undetermined' >$projects_dir/$dir/fastq_infiles_list.tx

mv *.gz untrimmed_fastq/
mv untrimmed_fastq/*Undetermined* $projects_dir/$dir/ 
echo "created fastq_infiles_list.tx"
echo "moved .fastq.gz files to untrimmed_fastq folder"
mkdir $projects_dir/$dir/Reports/ -p
touch $projects_dir/$dir/Reports/Demultiplex_Stats.csv
echo "this file was created by the skip_bcl.sh script that enables pipeline usage with already demultiplexed data.">$projects_dir/$dir/Reports/Demultiplex_Stats.csv
