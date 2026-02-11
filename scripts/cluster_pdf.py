#!/usr/bin/env python3
"""
Comprehensive Illumina InterOp Analysis and QC Report Generator

Combines cluster counting, quality metrics, and detailed QC plots
into a single comprehensive PDF report.

Usage:
    python interop_comprehensive.py <run_folder> <output_pdf>
    
Example:
    python interop_comprehensive.py /data/NovaSeq_Run_001 comprehensive_report.pdf

Requirements:
    pip install interop matplotlib seaborn pandas numpy
"""

import sys
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.gridspec import GridSpec
    from datetime import datetime
except ImportError as e:
    print(f"ERROR: Missing required package: {e}", file=sys.stderr)
    print("\nInstall with:", file=sys.stderr)
    print("  pip install matplotlib seaborn pandas numpy", file=sys.stderr)
    sys.exit(1)

# Try to import interop - works with both pip and conda installations
try:
    # Try pip installation first (interop package)
    from interop import py_interop_run_metrics, py_interop_run
    import interop.core as ic
    print("Using interop from pip installation", file=sys.stderr)
except ImportError:
    try:
        # Try conda installation (illumina-interop package)
        from interop import py_interop_run_metrics, py_interop_run
        from interop import core as ic
        print("Using interop from conda installation", file=sys.stderr)
    except ImportError:
        print("ERROR: illumina-interop not installed", file=sys.stderr)
        print("\nInstall with ONE of:", file=sys.stderr)
        print("  pip install interop", file=sys.stderr)
        print("  conda install -c bioconda illumina-interop", file=sys.stderr)
        sys.exit(1)

# Styling
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150

# Color palettes
LANE_COLORS = {
    1: '#1f77b4', 2: '#ff7f0e', 3: '#2ca02c', 4: '#d62728',
    5: '#9467bd', 6: '#8c564b', 7: '#e377c2', 8: '#7f7f7f',
}

CHANNEL_COLORS = {
    'A': '#2ca02c', 'C': '#1f77b4', 'G': '#000000', 'T': '#d62728',
}

PF_COLOR = '#2ecc71'
NON_PF_COLOR = '#e74c3c'
WARN_COLOR = '#f39c12'


def get_cluster_counts(run_metrics):
    """Extract cluster counts (total, PF, non-PF) per lane and overall."""
    tile_metrics = run_metrics.tile_metric_set()
    
    total_clusters = 0
    pf_clusters = 0
    lane_data = {}
    
    for tile_idx in range(tile_metrics.size()):
        tile = tile_metrics.at(tile_idx)
        lane = tile.lane()
        tile_total = tile.cluster_count()
        tile_pf = tile.cluster_count_pf()
        
        total_clusters += tile_total
        pf_clusters += tile_pf
        
        if lane not in lane_data:
            lane_data[lane] = {'total': 0, 'pf': 0}
        lane_data[lane]['total'] += tile_total
        lane_data[lane]['pf'] += tile_pf
    
    non_pf_clusters = total_clusters - pf_clusters
    
    for lane in lane_data:
        lane_data[lane]['non_pf'] = lane_data[lane]['total'] - lane_data[lane]['pf']
        lane_data[lane]['percent_pf'] = (
            (lane_data[lane]['pf'] / lane_data[lane]['total'] * 100) 
            if lane_data[lane]['total'] > 0 else 0
        )
    
    return {
        'total': int(total_clusters),
        'pf': int(pf_clusters),
        'non_pf': int(non_pf_clusters),
        'percent_pf': (pf_clusters / total_clusters * 100) if total_clusters > 0 else 0,
        'lanes': lane_data
    }


def load_imaging_table(run_folder):
    """Load imaging table with error handling."""
    try:
        ar = ic.imaging(str(run_folder))
        df = pd.DataFrame(ar)
        
        # Print available columns for debugging
        print(f"  → Imaging columns available: {', '.join(df.columns.tolist())}")
        
        int_cols = ['Lane', 'Tile', 'Tile Number', 'Read', 'Cycle', 
                    'Cycle Within Read', 'Swath', 'Surface']
        for col in int_cols:
            if col in df.columns:
                df[col] = df[col].astype(int)
        
        return df
    except Exception as e:
        print(f"  WARNING: Could not load imaging table: {e}", file=sys.stderr)
        return None


def has_valid_data(data):
    """Check if data contains valid non-zero values."""
    if data is None or len(data) == 0:
        return False
    
    if isinstance(data, pd.Series):
        valid_data = data.dropna()
        if len(valid_data) == 0 or (valid_data == 0).all():
            return False
        if valid_data.std() < 1e-10:
            return False
    elif isinstance(data, pd.DataFrame):
        valid_data = data.dropna(how='all')
        if len(valid_data) == 0:
            return False
        numeric_cols = valid_data.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0 and (valid_data[numeric_cols] == 0).all().all():
            return False
    
    return True


# ============================================================================
# PDF PLOTTING FUNCTIONS
# ============================================================================

def plot_summary_page(pdf, run_folder, cluster_counts, run_metrics):
    """Page 1: Summary statistics with run info."""
    fig = plt.figure(figsize=(8.5, 11))
    fig.suptitle(f'Comprehensive Sequencing Run Report\n{Path(run_folder).name}', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    ax = fig.add_subplot(111)
    ax.axis('off')
    
    total = cluster_counts['total']
    pf = cluster_counts['pf']
    non_pf = cluster_counts['non_pf']
    pct_pf = cluster_counts['percent_pf']
    
    # Build summary text
    summary = f"""
RUN INFORMATION:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Run Folder:    {run_folder}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

"""
    
    # Add run info if available
    try:
        run_info = run_metrics.run_info()
        summary += f"""Flowcell ID:   {run_info.flowcell().name()}
Total Cycles:  {run_info.total_cycles()}
Number Reads:  {run_info.reads().size()}

"""
    except:
        pass
    
    summary += f"""
OVERALL CLUSTER STATISTICS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Total Clusters:                {total:,}
Clusters Passing Filter (PF):  {pf:,}
Clusters Failing Filter:       {non_pf:,}

Percent PF:                    {pct_pf:.2f}%
Percent Non-PF:                {100 - pct_pf:.2f}%

QUALITY ASSESSMENT:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
    
    if pct_pf >= 80:
        summary += "✓ EXCELLENT - >80% clusters passing filter\n"
    elif pct_pf >= 70:
        summary += "✓ GOOD - 70-80% clusters passing filter\n"
    else:
        summary += "⚠ WARNING - <70% clusters passing filter\n"
    
    if (100 - pct_pf) > 30:
        summary += "⚠ High failure rate detected\n"
    
    summary += "\n\nPER-LANE BREAKDOWN:\n"
    summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
    
    for lane in sorted(cluster_counts['lanes'].keys()):
        ld = cluster_counts['lanes'][lane]
        summary += f"Lane {lane}:\n"
        summary += f"  Total:    {ld['total']:>15,}\n"
        summary += f"  PF:       {ld['pf']:>15,}  ({ld['percent_pf']:.2f}%)\n"
        summary += f"  Non-PF:   {ld['non_pf']:>15,}\n\n"
    
    ax.text(0.05, 0.95, summary, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_cluster_distribution(pdf, cluster_counts):
    """Page 2: Cluster distribution - stacked bars and percent PF."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5))
    fig.suptitle('Cluster Distribution Across Lanes', fontsize=14, fontweight='bold')
    
    lanes = sorted(cluster_counts['lanes'].keys())
    lane_labels = [f'L{l}' for l in lanes]
    
    total_counts = [cluster_counts['lanes'][l]['total'] for l in lanes]
    pf_counts = [cluster_counts['lanes'][l]['pf'] for l in lanes]
    non_pf_counts = [cluster_counts['lanes'][l]['non_pf'] for l in lanes]
    
    x = np.arange(len(lanes))
    width = 0.35
    
    # Stacked bar chart
    ax1.bar(x, pf_counts, width, label='PF Clusters', color=PF_COLOR, alpha=0.8)
    ax1.bar(x, non_pf_counts, width, bottom=pf_counts, label='Non-PF', color=NON_PF_COLOR, alpha=0.8)
    ax1.set_xlabel('Lane', fontsize=11)
    ax1.set_ylabel('Cluster Count', fontsize=11)
    ax1.set_title('Stacked PF vs Non-PF\nShows quality distribution per lane', fontsize=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels(lane_labels)
    ax1.legend(fontsize=9)
    ax1.grid(axis='y', alpha=0.3)
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x/1e6):.0f}M'))
    
    # Percent PF by lane
    percent_pf = [cluster_counts['lanes'][l]['percent_pf'] for l in lanes]
    colors = [PF_COLOR if p >= 80 else WARN_COLOR if p >= 70 else NON_PF_COLOR for p in percent_pf]
    
    ax2.bar(x, percent_pf, width*2, color=colors, alpha=0.8)
    ax2.axhline(y=80, color='green', linestyle='--', linewidth=1, alpha=0.5, label='Excellent (80%)')
    ax2.axhline(y=70, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='Good (70%)')
    ax2.set_xlabel('Lane', fontsize=11)
    ax2.set_ylabel('Percent Passing Filter (%)', fontsize=11)
    ax2.set_title('Percent PF by Lane\nHigher is better', fontsize=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels(lane_labels)
    ax2.set_ylim([0, 100])
    ax2.legend(fontsize=8)
    ax2.grid(axis='y', alpha=0.3)
    
    for i, v in enumerate(percent_pf):
        ax2.text(i, v + 2, f'{v:.1f}%', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_percent_occupied_vs_pf(pdf, imaging_df):
    """Page 3: % Occupied vs % Pass Filter scatter plot - QC diagnostic."""
    if imaging_df is None:
        print("    (No imaging data available)", file=sys.stderr)
        return
    
    if '% Occupied' not in imaging_df.columns or '% Pass Filter' not in imaging_df.columns:
        print(f"    (Missing required columns for % Occupied vs % PF plot)", file=sys.stderr)
        return
    
    if not has_valid_data(imaging_df['% Occupied']) or not has_valid_data(imaging_df['% Pass Filter']):
        print(f"    (Data is invalid/empty)", file=sys.stderr)
        return
    
    try:
        plot_data = imaging_df.copy()
        
        # Get data ranges for auto-scaling (don't force start at 0)
        occupied_vals = plot_data['% Occupied'].dropna()
        pf_vals = plot_data['% Pass Filter'].dropna()
        
        x_min = max(0, occupied_vals.min() - 2)
        x_max = min(100, occupied_vals.max() + 2)
        y_min = max(0, pf_vals.min() - 2)
        y_max = min(100, pf_vals.max() + 2)
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        fig.suptitle('% Occupied vs % Pass Filter (by Lane)', fontsize=14, fontweight='bold')
        
        lanes = sorted(plot_data['Lane'].unique())
        
        # Plot each lane with consistent color
        for lane in lanes:
            lane_data = plot_data[plot_data['Lane'] == lane]
            if len(lane_data) > 0:
                color = LANE_COLORS.get(lane, f'C{lane-1}')
                ax.scatter(lane_data['% Occupied'], lane_data['% Pass Filter'],
                          label=f'Lane {lane}', alpha=0.6, s=25, color=color, edgecolors='none')
        
        # Reference lines
        ax.axhline(y=80, color='orange', linestyle='--', linewidth=1.5, alpha=0.5, label='80% PF threshold')
        ax.axvline(x=95, color='green', linestyle='--', linewidth=1.5, alpha=0.5, label='95% Occupied threshold')
        
        ax.set_xlabel('% Occupied', fontsize=12)
        ax.set_ylabel('% Pass Filter', fontsize=12)
        ax.set_title('Each point = one tile/cycle; good runs cluster in upper-right', fontsize=10)
        ax.legend(fontsize=9, loc='lower right', framealpha=0.9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])
        ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating % Occupied vs % PF plot: {e})", file=sys.stderr)
        return


def plot_percent_occupied_vs_pf_by_surface(pdf, imaging_df):
    """Page 4: % Occupied vs % Pass Filter scatter plot colored by Surface."""
    if imaging_df is None:
        print("    (No imaging data available)", file=sys.stderr)
        return
    
    if '% Occupied' not in imaging_df.columns or '% Pass Filter' not in imaging_df.columns:
        print(f"    (Missing required columns for % Occupied vs % PF plot)", file=sys.stderr)
        return
    
    if 'Surface' not in imaging_df.columns:
        print(f"    (No 'Surface' column found)", file=sys.stderr)
        return
    
    if not has_valid_data(imaging_df['% Occupied']) or not has_valid_data(imaging_df['% Pass Filter']):
        print(f"    (Data is invalid/empty)", file=sys.stderr)
        return
    
    try:
        plot_data = imaging_df.copy()
        
        # Get data ranges for auto-scaling
        occupied_vals = plot_data['% Occupied'].dropna()
        pf_vals = plot_data['% Pass Filter'].dropna()
        
        x_min = max(0, occupied_vals.min() - 2)
        x_max = min(100, occupied_vals.max() + 2)
        y_min = max(0, pf_vals.min() - 2)
        y_max = min(100, pf_vals.max() + 2)
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        fig.suptitle('% Occupied vs % Pass Filter (by Surface)', fontsize=14, fontweight='bold')
        
        surfaces = sorted(plot_data['Surface'].unique())
        
        # Define colors for surfaces (typically 1 and 2 for flow cells)
        surface_colors = {
            1: '#3498db',  # Blue
            2: '#e74c3c',  # Red
        }
        
        # Plot each surface with different color
        for surface in surfaces:
            surface_data = plot_data[plot_data['Surface'] == surface]
            if len(surface_data) > 0:
                color = surface_colors.get(surface, f'C{surface}')
                ax.scatter(surface_data['% Occupied'], surface_data['% Pass Filter'],
                          label=f'Surface {surface}', alpha=0.6, s=25, color=color, edgecolors='none')
        
        # Reference lines
        ax.axhline(y=80, color='orange', linestyle='--', linewidth=1.5, alpha=0.5, label='80% PF threshold')
        ax.axvline(x=95, color='green', linestyle='--', linewidth=1.5, alpha=0.5, label='95% Occupied threshold')
        
        ax.set_xlabel('% Occupied', fontsize=12)
        ax.set_ylabel('% Pass Filter', fontsize=12)
        ax.set_title('Each point = one tile/cycle; helps identify surface-specific issues', fontsize=10)
        ax.legend(fontsize=9, loc='lower right', framealpha=0.9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])
        ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating % Occupied vs % PF by Surface plot: {e})", file=sys.stderr)
        return


def plot_qscore_heatmap(pdf, imaging_df):
    """Page 4: Q30 heatmap across cycles and lanes."""
    if imaging_df is None or '%>= Q30' not in imaging_df.columns:
        print("    (No Q30 data available)", file=sys.stderr)
        return
    
    if not has_valid_data(imaging_df['%>= Q30']):
        print("    (Q30 data is invalid/empty)", file=sys.stderr)
        return
    
    # Pivot data
    try:
        pivot = imaging_df.pivot_table(values='%>= Q30', index='Lane', columns='Cycle', aggfunc='mean')
        
        if pivot.empty or not has_valid_data(pivot):
            print("    (Pivot table is empty)", file=sys.stderr)
            return
        
        fig, ax = plt.subplots(1, 1, figsize=(11, 6))
        fig.suptitle('Q30 Percentage Heatmap', fontsize=14, fontweight='bold')
        
        sns.heatmap(pivot, cmap='RdYlGn', vmin=70, vmax=95, cbar_kws={'label': '% >= Q30'},
                    ax=ax, linewidths=0.5, linecolor='white')
        ax.set_xlabel('Cycle', fontsize=12)
        ax.set_ylabel('Lane', fontsize=12)
        ax.set_title('Higher values (green) indicate better quality', fontsize=10)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating heatmap: {e})", file=sys.stderr)
        return


def plot_qscore_by_cycle(pdf, imaging_df):
    """Page 5: Q-score trends by cycle for each lane."""
    if imaging_df is None or '%>= Q30' not in imaging_df.columns:
        print("    (No Q30 data available)", file=sys.stderr)
        return
    
    if not has_valid_data(imaging_df['%>= Q30']):
        print("    (Q30 data is invalid/empty)", file=sys.stderr)
        return
    
    try:
        fig, ax = plt.subplots(1, 1, figsize=(11, 6))
        fig.suptitle('Q30 Percentage by Cycle', fontsize=14, fontweight='bold')
        
        lanes = sorted(imaging_df['Lane'].unique())
        plotted = False
        
        for lane in lanes:
            lane_data = imaging_df[imaging_df['Lane'] == lane]
            cycle_mean = lane_data.groupby('Cycle')['%>= Q30'].mean()
            
            if has_valid_data(cycle_mean):
                color = LANE_COLORS.get(lane, f'C{lane-1}')
                ax.plot(cycle_mean.index, cycle_mean.values, 
                       label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                plotted = True
        
        if not plotted:
            print("    (No valid lane data to plot)", file=sys.stderr)
            plt.close(fig)
            return
        
        ax.axhline(y=80, color='green', linestyle='--', linewidth=1, alpha=0.5, label='80% threshold')
        ax.set_xlabel('Cycle', fontsize=12)
        ax.set_ylabel('% Bases >= Q30', fontsize=12)
        ax.set_title('Q30 represents 99.9% base call accuracy; higher is better', fontsize=10)
        ax.legend(fontsize=9, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 100])
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating Q30 plot: {e})", file=sys.stderr)
        return


def plot_error_rate(pdf, imaging_df):
    """Page 6: Quality metrics - No Calls percentage by cycle."""
    if imaging_df is None:
        print("    (No imaging data available)", file=sys.stderr)
        return
    
    # Use '% No Calls' as a quality metric (lower is better)
    if '% No Calls' not in imaging_df.columns:
        print(f"    (No '% No Calls' column found)", file=sys.stderr)
        return
    
    if not has_valid_data(imaging_df['% No Calls']):
        print(f"    ('% No Calls' data is invalid/empty)", file=sys.stderr)
        return
    
    try:
        #  Filter to Read 1 only to avoid duplicate cycles
        if 'Read' in imaging_df.columns:
            plot_data = imaging_df[imaging_df['Read'] == 1]
            if len(plot_data) == 0:
                plot_data = imaging_df
        else:
            plot_data = imaging_df
        
        fig, ax = plt.subplots(1, 1, figsize=(11, 6))
        fig.suptitle('No Calls Percentage by Cycle', fontsize=14, fontweight='bold')
        
        lanes = sorted(plot_data['Lane'].unique())
        plotted = False
        
        for lane in lanes:
            lane_data = plot_data[plot_data['Lane'] == lane]
            cycle_mean = lane_data.groupby('Cycle')['% No Calls'].mean()
            
            if has_valid_data(cycle_mean):
                color = LANE_COLORS.get(lane, f'C{lane-1}')
                ax.plot(cycle_mean.index, cycle_mean.values, 
                       label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                plotted = True
        
        if not plotted:
            print("    (No valid lane data)", file=sys.stderr)
            plt.close(fig)
            return
        
        ax.set_xlabel('Cycle', fontsize=12)
        ax.set_ylabel('% No Calls', fontsize=12)
        ax.set_title('Lower values indicate better base calling accuracy', fontsize=10)
        ax.legend(fontsize=9, loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating No Calls plot: {e})", file=sys.stderr)
        return




def plot_phasing_offset_slope(pdf, imaging_df):
    """Phasing Offset and Slope by lane."""
    if imaging_df is None:
        print("    (No imaging data available)", file=sys.stderr)
        return
    
    # Check for Phasing Offset and Phasing Slope columns
    has_offset = 'Phasing Offset' in imaging_df.columns
    has_slope = 'Phasing Slope' in imaging_df.columns
    
    if not (has_offset or has_slope):
        print(f"    (No Phasing Offset or Slope columns found)", file=sys.stderr)
        return
    
    try:
        # Filter to Read 1 only
        if 'Read' in imaging_df.columns:
            plot_data = imaging_df[imaging_df['Read'] == 1].copy()
            if len(plot_data) == 0:
                plot_data = imaging_df
        else:
            plot_data = imaging_df
        
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
        fig.suptitle('Phasing Offset and Slope by Cycle', fontsize=14, fontweight='bold')
        
        lanes = sorted(plot_data['Lane'].unique())
        
        # Plot Phasing Offset
        if has_offset and has_valid_data(plot_data['Phasing Offset']):
            ax = axes[0]
            plotted = False
            for lane in lanes:
                lane_data = plot_data[plot_data['Lane'] == lane]
                cycle_mean = lane_data.groupby('Cycle')['Phasing Offset'].mean()
                
                if has_valid_data(cycle_mean):
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                    plotted = True
            
            if plotted:
                ax.set_xlabel('Cycle', fontsize=11)
                ax.set_ylabel('Phasing Offset', fontsize=11)
                ax.set_title('Phasing Offset\nLower is better', fontsize=10)
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
            else:
                ax.axis('off')
        else:
            axes[0].axis('off')
        
        # Plot Phasing Slope
        if has_slope and has_valid_data(plot_data['Phasing Slope']):
            ax = axes[1]
            plotted = False
            for lane in lanes:
                lane_data = plot_data[plot_data['Lane'] == lane]
                cycle_mean = lane_data.groupby('Cycle')['Phasing Slope'].mean()
                
                if has_valid_data(cycle_mean):
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                    plotted = True
            
            if plotted:
                ax.set_xlabel('Cycle', fontsize=11)
                ax.set_ylabel('Phasing Slope', fontsize=11)
                ax.set_title('Phasing Slope\nShould be stable', fontsize=10)
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
            else:
                ax.axis('off')
        else:
            axes[1].axis('off')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating phasing offset/slope plot: {e})", file=sys.stderr)
        return




def plot_intensity_metrics(pdf, imaging_df):
    """Page 8: Signal intensity by channel using Corrected values."""
    if imaging_df is None:
        print("    (No imaging data available)", file=sys.stderr)
        return
    
    # Check for Corrected/A, Corrected/C, Corrected/G, Corrected/T
    channels = ['A', 'C', 'G', 'T']
    corrected_cols = [f'Corrected/{ch}' for ch in channels]
    
    # Check if all corrected columns exist
    if not all(col in imaging_df.columns for col in corrected_cols):
        print(f"    (Missing Corrected intensity columns)", file=sys.stderr)
        return
    
    try:
        # Filter to Read 1 only
        if 'Read' in imaging_df.columns:
            plot_data = imaging_df[imaging_df['Read'] == 1].copy()
            if len(plot_data) == 0:
                plot_data = imaging_df
        else:
            plot_data = imaging_df
        
        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        fig.suptitle('Corrected Intensity by Channel', fontsize=14, fontweight='bold')
        axes = axes.flatten()
        
        plotted_any = False
        
        for idx, channel in enumerate(channels):
            ax = axes[idx]
            col_name = f'Corrected/{channel}'
            
            if has_valid_data(plot_data[col_name]):
                lanes = sorted(plot_data['Lane'].unique())
                
                for lane in lanes:
                    lane_data = plot_data[plot_data['Lane'] == lane]
                    cycle_mean = lane_data.groupby('Cycle')[col_name].mean()
                    
                    if has_valid_data(cycle_mean):
                        color = LANE_COLORS.get(lane, f'C{lane-1}')
                        ax.plot(cycle_mean.index, cycle_mean.values, 
                               label=f'L{lane}', linewidth=1.5, alpha=0.7, color=color)
                        plotted_any = True
                
                ax.set_xlabel('Cycle', fontsize=10)
                ax.set_ylabel('Corrected Intensity', fontsize=10)
                ax.set_title(f'Channel {channel}', fontsize=11, fontweight='bold', 
                            color=CHANNEL_COLORS.get(channel, 'black'))
                ax.legend(fontsize=8, loc='best')
                ax.grid(True, alpha=0.3)
            else:
                ax.axis('off')
        
        if not plotted_any:
            print("    (No valid intensity data to plot)", file=sys.stderr)
            plt.close(fig)
            return
        
        fig.text(0.5, 0.02, 'Consistent intensity across cycles indicates stable sequencing chemistry', 
                 ha='center', fontsize=9, style='italic')
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating intensity plot: {e})", file=sys.stderr)
        return


def plot_fwhm_metrics(pdf, imaging_df):
    """Page 9: FWHM by color channel (green/blue for 2-channel sequencing)."""
    if imaging_df is None:
        print("    (No imaging data available)", file=sys.stderr)
        return
    
    # Check for Fwhm/green and Fwhm/blue (2-channel chemistry)
    has_fwhm = 'Fwhm/green' in imaging_df.columns or 'Fwhm/blue' in imaging_df.columns
    
    if not has_fwhm:
        print(f"    (No FWHM columns found)", file=sys.stderr)
        return
    
    try:
        # Filter to Read 1 only
        if 'Read' in imaging_df.columns:
            plot_data = imaging_df[imaging_df['Read'] == 1].copy()
            if len(plot_data) == 0:
                plot_data = imaging_df
        else:
            plot_data = imaging_df
        
        fig, ax = plt.subplots(1, 1, figsize=(11, 6))
        fig.suptitle('FWHM (Full Width Half Maximum) by Color Channel', fontsize=14, fontweight='bold')
        
        plotted = False
        
        # Plot green channel in green color
        if 'Fwhm/green' in plot_data.columns and has_valid_data(plot_data['Fwhm/green']):
            lanes = sorted(plot_data['Lane'].unique())
            for lane in lanes:
                lane_data = plot_data[plot_data['Lane'] == lane]
                cycle_mean = lane_data.groupby('Cycle')['Fwhm/green'].mean()
                
                if has_valid_data(cycle_mean):
                    # Use green color for green channel
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'L{lane} Green', linewidth=2, alpha=0.7, 
                           color='green', linestyle='-')
                    plotted = True
        
        # Plot blue channel in blue color
        if 'Fwhm/blue' in plot_data.columns and has_valid_data(plot_data['Fwhm/blue']):
            lanes = sorted(plot_data['Lane'].unique())
            for lane in lanes:
                lane_data = plot_data[plot_data['Lane'] == lane]
                cycle_mean = lane_data.groupby('Cycle')['Fwhm/blue'].mean()
                
                if has_valid_data(cycle_mean):
                    # Use blue color for blue channel
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'L{lane} Blue', linewidth=2, alpha=0.7, 
                           color='blue', linestyle='--')
                    plotted = True
        
        if not plotted:
            print("    (No valid FWHM data to plot)", file=sys.stderr)
            plt.close(fig)
            return
        
        ax.set_xlabel('Cycle', fontsize=12)
        ax.set_ylabel('FWHM', fontsize=12)
        ax.set_title('Lower FWHM values indicate better focus and cluster definition', fontsize=10)
        ax.legend(fontsize=9, loc='best', ncol=2)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating FWHM plot: {e})", file=sys.stderr)
        return


def plot_phasing_prephasing(pdf, imaging_df):
    """Page 10: Phasing and prephasing weight metrics."""
    if imaging_df is None:
        print("    (No imaging data available)", file=sys.stderr)
        return
    
    # Use 'Phasing Weight' and 'Prephasing Weight'
    has_phasing = 'Phasing Weight' in imaging_df.columns and has_valid_data(imaging_df['Phasing Weight'])
    has_prephasing = 'Prephasing Weight' in imaging_df.columns and has_valid_data(imaging_df['Prephasing Weight'])
    
    if not (has_phasing or has_prephasing):
        print(f"    (No Phasing/Prephasing Weight columns found)", file=sys.stderr)
        return
    
    try:
        # Filter to Read 1 only
        if 'Read' in imaging_df.columns:
            plot_data = imaging_df[imaging_df['Read'] == 1].copy()
            if len(plot_data) == 0:
                plot_data = imaging_df
        else:
            plot_data = imaging_df
        
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
        fig.suptitle('Phasing and Prephasing Weight Metrics', fontsize=14, fontweight='bold')
        
        lanes = sorted(plot_data['Lane'].unique())
        
        # Phasing Weight
        if has_phasing:
            ax = axes[0]
            plotted = False
            for lane in lanes:
                lane_data = plot_data[plot_data['Lane'] == lane]
                cycle_mean = lane_data.groupby('Cycle')['Phasing Weight'].mean()
                
                if has_valid_data(cycle_mean):
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                    plotted = True
            
            if plotted:
                ax.set_xlabel('Cycle', fontsize=11)
                ax.set_ylabel('Phasing Weight', fontsize=11)
                ax.set_title('Phasing Weight\nLower is better', fontsize=10)
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
            else:
                ax.axis('off')
        else:
            axes[0].axis('off')
        
        # Prephasing Weight
        if has_prephasing:
            ax = axes[1]
            plotted = False
            for lane in lanes:
                lane_data = plot_data[plot_data['Lane'] == lane]
                cycle_mean = lane_data.groupby('Cycle')['Prephasing Weight'].mean()
                
                if has_valid_data(cycle_mean):
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                    plotted = True
            
            if plotted:
                ax.set_xlabel('Cycle', fontsize=11)
                ax.set_ylabel('Prephasing Weight', fontsize=11)
                ax.set_title('Prephasing Weight\nLower is better', fontsize=10)
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
            else:
                ax.axis('off')
        else:
            axes[1].axis('off')
        
        # Save if at least one plot was created
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"    (Error creating phasing plot: {e})", file=sys.stderr)
        return


def generate_comprehensive_report(run_folder, output_pdf):
    """Generate comprehensive PDF report with all metrics."""
    run_folder_path = Path(run_folder).resolve()
    
    if not run_folder_path.exists():
        print(f"ERROR: Run folder does not exist: {run_folder}", file=sys.stderr)
        sys.exit(1)
    
    if not (run_folder_path / "InterOp").exists():
        print(f"ERROR: InterOp directory not found: {run_folder_path}/InterOp", file=sys.stderr)
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print("COMPREHENSIVE ILLUMINA RUN ANALYSIS")
    print(f"{'='*80}\n")
    print(f"Run folder: {run_folder}")
    print(f"Output PDF: {output_pdf}")
    
    # Load metrics
    print("\n[1/3] Loading InterOp metrics...")
    run_metrics = py_interop_run_metrics.run_metrics()
    
    try:
        run_metrics.read(str(run_folder_path))
        print("  ✓ Successfully loaded run metrics")
    except Exception as e:
        print(f"  ✗ ERROR: Failed to load metrics: {e}", file=sys.stderr)
        sys.exit(1)
    
    print("\n[2/3] Extracting data...")
    cluster_counts = get_cluster_counts(run_metrics)
    print(f"  ✓ Cluster counts extracted")
    print(f"    - Total: {cluster_counts['total']:,}")
    print(f"    - PF: {cluster_counts['pf']:,} ({cluster_counts['percent_pf']:.2f}%)")
    print(f"    - Lanes: {len(cluster_counts['lanes'])}")
    
    imaging_df = load_imaging_table(run_folder_path)
    if imaging_df is not None:
        print(f"  ✓ Imaging table loaded ({len(imaging_df)} rows)")
    else:
        print("  ⚠ Imaging table not available (some plots will be skipped)")
    
    # Generate PDF
    print("\n[3/3] Generating PDF report...")
    page_count = 0
    
    with PdfPages(output_pdf) as pdf:
        # Page 1: Summary (always included)
        try:
            print("  → Page 1: Summary statistics")
            plot_summary_page(pdf, run_folder, cluster_counts, run_metrics)
            page_count += 1
        except Exception as e:
            print(f"  ✗ Failed to generate summary page: {e}", file=sys.stderr)
        
        # Page 2: Cluster distribution (always included)
        try:
            print("  → Page 2: Cluster distribution")
            plot_cluster_distribution(pdf, cluster_counts)
            page_count += 1
        except Exception as e:
            print(f"  ✗ Failed to generate cluster distribution: {e}", file=sys.stderr)
        
        # Remaining pages - only if imaging data available
        if imaging_df is not None:
            # Page 3: % Occupied vs % Pass Filter (by Lane)
            try:
                print("  → Page 3: % Occupied vs % Pass Filter (by Lane)")
                plot_percent_occupied_vs_pf(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped % Occupied vs PF (by Lane): {e}", file=sys.stderr)
            
            # Page 4: % Occupied vs % Pass Filter (by Surface)
            try:
                print("  → Page 4: % Occupied vs % Pass Filter (by Surface)")
                plot_percent_occupied_vs_pf_by_surface(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped % Occupied vs PF (by Surface): {e}", file=sys.stderr)
            
            # Page 5: Q30 heatmap
            try:
                print("  → Page 5: Q30 heatmap")
                plot_qscore_heatmap(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped Q30 heatmap: {e}", file=sys.stderr)
            
            # Page 6: Q30 by cycle
            try:
                print("  → Page 6: Q30 by cycle")
                plot_qscore_by_cycle(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped Q30 by cycle: {e}", file=sys.stderr)
            
            # Page 7: Error rate (% No Calls)
            try:
                print("  → Page 7: Error rate by cycle")
                plot_error_rate(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped error rate: {e}", file=sys.stderr)
            
            # Page 8: Signal intensity
            try:
                print("  → Page 8: Signal intensity")
                plot_intensity_metrics(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped signal intensity: {e}", file=sys.stderr)
            
            # Page 9: FWHM metrics
            try:
                print("  → Page 9: FWHM metrics")
                plot_fwhm_metrics(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped FWHM metrics: {e}", file=sys.stderr)
            
            # Page 10: Phasing Offset/Slope
            try:
                print("  → Page 10: Phasing Offset/Slope")
                plot_phasing_offset_slope(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped phasing offset/slope: {e}", file=sys.stderr)
            
            # Page 11: Phasing/Prephasing Weights
            try:
                print("  → Page 11: Phasing/Prephasing Weights")
                plot_phasing_prephasing(pdf, imaging_df)
                page_count += 1
            except Exception as e:
                print(f"  ⚠ Skipped phasing/prephasing: {e}", file=sys.stderr)
        
        # Set PDF metadata
        d = pdf.infodict()
        d['Title'] = f'Comprehensive Run Report - {Path(run_folder).name}'
        d['Author'] = 'interop_comprehensive.py'
        d['Subject'] = 'Illumina Sequencing QC Report'
        d['Keywords'] = 'Illumina, InterOp, QC, Sequencing, Clusters'
        d['CreationDate'] = None
    
    print(f"\n{'='*80}")
    print("REPORT GENERATION COMPLETE")
    print(f"{'='*80}")
    print(f"\nOutput: {output_pdf}")
    print(f"Pages generated: {page_count}")
    print(f"\nSummary:")
    print(f"  Total clusters:     {cluster_counts['total']:,}")
    print(f"  PF clusters:        {cluster_counts['pf']:,} ({cluster_counts['percent_pf']:.2f}%)")
    print(f"  Non-PF clusters:    {cluster_counts['non_pf']:,}")
    print(f"  Number of lanes:    {len(cluster_counts['lanes'])}")
    print()


def main():
    """Main entry point."""
    if len(sys.argv) != 3:
        print("Usage: python interop_comprehensive.py <run_folder> <output_pdf>")
        print("\nExample:")
        print("  python interop_comprehensive.py /data/NovaSeq_Run_001 report.pdf")
        sys.exit(1)
    
    run_folder = sys.argv[1]
    output_pdf = sys.argv[2]
    
    try:
        generate_comprehensive_report(run_folder, output_pdf)
    except Exception as e:
        print(f"\nERROR: Report generation failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()