#!/usr/bin/env python3
"""
Minimal script to output cluster counts from Illumina run folder.

Usage:
    python get_clusters_simple.py /path/to/RunFolder

Output:
    Line 1: Total clusters
    Line 2: Clusters passing filter (PF)
    Line 3: Non-PF clusters (failed filter)
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
    total_clusters = 0
    pf_clusters = 0
    
    for i in range(tile_metrics.size()):
        tile = tile_metrics.at(i)
        total_clusters += tile.cluster_count()
        pf_clusters += tile.cluster_count_pf()
    
    # Calculate non-PF clusters (failed filter)
    non_pf_clusters = total_clusters - pf_clusters
    
    # Output the three numbers
    print(int(total_clusters),"Total possible Clusters (# Wells)")
    print(int(pf_clusters),"Clusters passing Filters")
    print(int(non_pf_clusters),"Clusters not passed")

if __name__ == "__main__":
    main()