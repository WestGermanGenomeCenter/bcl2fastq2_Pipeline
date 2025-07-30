#!/bin/bash

# Check if a directory argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/your/folder"
    exit 1
fi

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"


# Check for empty files
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

if [ -z "$empty_files" ]; then
    echo "No empty files found in $DIRECTORY."
else
    echo "Empty files found in $DIRECTORY:$empty_files"
    exit 1
fi















#set -e 
## Check if exactly two arguments are provided (output file and input directory)
#if [ "$#" -ne 2 ]; then
#    echo "Usage: $0 output_file input_directory"
#    exit 1
#fi
#
## First argument is the output file
#output_file="$1"
## Second argument is the input directory
#input_directory="$2"
#
## Clear the output file if it exists
#> "$output_file"
#
## Check if the input directory exists
#if [ ! -d "$input_directory" ]; then
#    echo "Error: $input_directory is not a valid directory."
#    exit 1
#fi
#
## Loop through each .fastq.gz file in the input directory
#for file in "$input_directory"/*.fastq.gz; do
#    if [ -f "$file" ]; then  # Check if it is a file
#        line_count=$(wc -l < "$file")  # Get the number of lines
#        if [ "$line_count" -eq 0 ]; then
#            trap echo "$file: Error - File is empty." >> "$output_file" ERR
#            exit 1
#        else
#            echo "$file: $line_count lines." >> "$output_file"
#        fi
#    else
#        echo "No .fastq.gz files found in $input_directory." >> "$output_file"
#        break  # Exit the loop if no files are found
#    fi
#done
#
#echo "Results written to $output_file."