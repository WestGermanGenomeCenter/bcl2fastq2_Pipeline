#!/usr/bin/env python3
"""
Minimal script to output total clusters and PF clusters from Illumina run folder.

Usage:
    python get_clusters_simple.py /path/to/RunFolder

Output:
    Line 1: Total clusters
    Line 2: Clusters passing filter
"""

import sys
from pathlib import Path

def main():
    if len(sys.argv) != 2:
        print("Usage: python get_clusters_simple.py /path/to/RunFolder", file=sys.stderr)
        sys.exit(1)
    
    run_folder = sys.argv[1]
    
    try:
        from interop import py_interop_run_metrics
    except ImportError:
        print("ERROR: illumina-interop not installed. Install with: pip install interop", file=sys.stderr)
        sys.exit(1)
    
    # Validate run folder
    run_folder_path = Path(run_folder).resolve()
    if not run_folder_path.exists() or not run_folder_path.is_dir():
        print(f"ERROR: Invalid run folder: {run_folder}", file=sys.stderr)
        sys.exit(1)
    
    # Read metrics
    run_metrics = py_interop_run_metrics.run_metrics()
    
    try:
        run_metrics.read(str(run_folder_path))
    except Exception as e:
        print(f"ERROR: Failed to read run metrics: {e}", file=sys.stderr)
        sys.exit(1)
    
    tile_metrics = run_metrics.tile_metric_set()
    
    # Sum clusters across all tiles
    total_clusters = sum(tile_metrics.at(i).cluster_count() for i in range(tile_metrics.size()))
    pf_clusters = sum(tile_metrics.at(i).cluster_count_pf() for i in range(tile_metrics.size()))
    
    # Output just the two numbers
    print("Possible Clusters (# Tiles):" total_clusters)
    print("Clusters passing Filters:" pf_clusters)

if __name__ == "__main__":
    main()