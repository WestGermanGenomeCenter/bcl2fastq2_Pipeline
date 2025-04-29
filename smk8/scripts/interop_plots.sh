#!/bin/bash
dir_of_cmds=$1
dir_of_data=$2
output_dir=$3

cd $3;

nice $1/plot_qscore_histogram $2 | gnuplot
nice $1/plot_qscore_heatmap $2 | gnuplot
nice $1/plot_flowcell $2 | gnuplot
nice $1/plot_by_lane $2 | gnuplot
nice $1/plot_by_cycle $2 | gnuplot
nice $1/plot_sample_qc $2 | gnuplot
