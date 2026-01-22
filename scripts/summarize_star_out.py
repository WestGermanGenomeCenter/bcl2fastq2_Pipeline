#!/usr/bin/env python3
"""
Parse STAR aligner Log.final.out files and generate a CSV summary.
Usage: python star_log_parser.py <folder_path> <output.csv>
"""

import os
import re
import csv
import sys
from pathlib import Path

def parse_star_log(filepath):
    """Parse a single STAR Log.final.out file and extract all metrics."""
    data = {}
    data['sample_name'] = Path(filepath).name.replace('_Log.final.out', '')
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if '|' in line:
                parts = line.split('|')
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()
                    
                    # Remove % signs and convert to float where applicable
                    value_clean = value.replace('%', '').strip()
                    
                    # Store the value
                    data[key] = value_clean
    
    return data

def main():
    if len(sys.argv) != 3:
        print("Usage: python star_log_parser.py <folder_path> <output.csv>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    output_csv = sys.argv[2]
    
    # Find all Log.final.out files
    log_files = list(Path(folder_path).glob('*Log.final.out'))
    
    if not log_files:
        print(f"No *Log.final.out files found in {folder_path}")
        sys.exit(1)
    
    print(f"Found {len(log_files)} STAR log files")
    
    # Parse all files
    all_data = []
    for log_file in sorted(log_files):
        print(f"Parsing {log_file.name}...")
        data = parse_star_log(log_file)
        all_data.append(data)
    
    # Get all unique keys (column headers)
    all_keys = set()
    for data in all_data:
        all_keys.update(data.keys())
    
    # Define column order (sample_name first, then chronological order from the log)
    ordered_columns = [
        'sample_name',
        'Started job on',
        'Started mapping on',
        'Finished on',
        'Mapping speed, Million of reads per hour',
        'Number of input reads',
        'Average input read length',
        'Uniquely mapped reads number',
        'Uniquely mapped reads %',
        'Average mapped length',
        'Number of splices: Total',
        'Number of splices: Annotated (sjdb)',
        'Number of splices: GT/AG',
        'Number of splices: GC/AG',
        'Number of splices: AT/AC',
        'Number of splices: Non-canonical',
        'Mismatch rate per base, %',
        'Deletion rate per base',
        'Deletion average length',
        'Insertion rate per base',
        'Insertion average length',
        'Number of reads mapped to multiple loci',
        '% of reads mapped to multiple loci',
        'Number of reads mapped to too many loci',
        '% of reads mapped to too many loci',
        'Number of reads unmapped: too many mismatches',
        '% of reads unmapped: too many mismatches',
        'Number of reads unmapped: too short',
        '% of reads unmapped: too short',
        'Number of reads unmapped: other',
        '% of reads unmapped: other',
        'Number of chimeric reads',
        '% of chimeric reads'
    ]
    
    # Add any additional columns that might exist
    for key in sorted(all_keys):
        if key not in ordered_columns:
            ordered_columns.append(key)
    
    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=ordered_columns)
        writer.writeheader()
        for data in all_data:
            writer.writerow(data)
    
    print(f"\nSuccessfully created {output_csv} with {len(all_data)} samples")

if __name__ == '__main__':
    main()