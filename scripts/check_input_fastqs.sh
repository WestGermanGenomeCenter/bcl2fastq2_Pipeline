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
