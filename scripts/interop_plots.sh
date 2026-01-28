#!/bin/bash
#dir_of_cmds=$1
dir_of_data=$1
output_dir=$2

cd $2;

nice interop_plot_qscore_histogram $1 | gnuplot
nice interop_plot_qscore_heatmap $1 | gnuplot
nice interop_plot_flowcell $1 | gnuplot
nice interop_plot_by_lane $1 | gnuplot
nice interop_plot_by_cycle $1 | gnuplot
nice interop_plot_sample_qc $1 | gnuplot
nice interop_imaging_table $1 >imaging_table.csv
nice interop_index-summary $1 --csv=1 >illumina-interop_index_summary.csv
nice interop_dumptext $1 >interop_textdump.csv
nice interop_summary $1 --csv=1 >interop_summary.csv