#!/bin/bash
input_report=$1
output_file=$2
echo "number_of_best_hits hit" >$2
cat $1 | awk '!x[$1]++' | awk '{print $2}'  | grep -v BLASTN |sort | uniq -c | sort -n | tac >>$2
