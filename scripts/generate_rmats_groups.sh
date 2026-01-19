#!/bin/bash
# Usage: ./generate_rmats_groups.sh <csv_file> <bam_dir> <group_name> <output_file>

CSV=$1
BAM_DIR=$2
TARGET_GROUP=$3
OUTPUT_FILE=$4

# 1. Extract sample IDs for the specific group (assumes CSV format: sample,group)
# We use awk to skip the header and match the second column
SAMPLES=$(awk -F',' -v grp="$TARGET_GROUP" '$2 == grp {print $1}' "$CSV")

# 2. Find the full paths for each sample
PATHS=()
for ID in $SAMPLES; do
    # Search for the file containing the sample ID
    FILE=$(find "$BAM_DIR" -name "*${ID}*.bam" | grep -v ".tmp." | head -n 1)
    
    if [ -z "$FILE" ]; then
        echo "Error: BAM file for sample $ID not found in $BAM_DIR" >&2
        exit 1
    fi
    
    # Get the absolute path
    PATHS+=($(readlink -f "$FILE"))
done

# 3. Join array with commas and write to the output file
IFS=',' 
echo "${PATHS[*]}" > "$OUTPUT_FILE"

echo "Successfully created $OUTPUT_FILE for group $TARGET_GROUP"