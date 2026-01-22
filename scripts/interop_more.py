#!/usr/bin/env python3
"""
Illumina InterOp QC Plot Generator using Python
Compatible with NovaSeqX+, NextSeq, and MiSeq
Generates comprehensive quality control plots for sequencing runs

Usage:
    python interop_qc_plots.py <run_folder> <output_dir>

Requirements:
    conda install -c bioconda illumina-interop
    conda install pandas matplotlib seaborn numpy reportlab
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
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    from interop import py_interop_run_metrics, py_interop_run, py_interop_plot
    import interop.core as ic
    from reportlab.lib.pagesizes import letter, A4
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    from datetime import datetime
except ImportError as e:
    print(f"Error: Missing required package: {e}")
    print("\nPlease install required packages:")
    print("  conda install -c bioconda illumina-interop")
    print("  conda install pandas matplotlib seaborn numpy reportlab")
    sys.exit(1)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# Define consistent color palette for lanes
LANE_COLORS = {
    1: '#1f77b4',  # Blue
    2: '#ff7f0e',  # Orange
    3: '#2ca02c',  # Green
    4: '#d62728',  # Red
    5: '#9467bd',  # Purple
    6: '#8c564b',  # Brown
    7: '#e377c2',  # Pink
    8: '#7f7f7f',  # Gray
}

# Define channel colors for FWHM
CHANNEL_COLORS = {
    'A': '#2ca02c',  # Green
    'C': '#1f77b4',  # Blue
    'G': '#000000',  # Black
    'T': '#d62728',  # Red
}

class InterOpQCPlotter:
    """Generate comprehensive QC plots from Illumina InterOp data"""
    
    def __init__(self, run_folder, output_dir):
        self.run_folder = Path(run_folder)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.plot_files = []  # Track generated plot files for PDF
        
        # Load run metrics
        print(f"[INFO] Loading run metrics from: {self.run_folder}")
        self.run_metrics = py_interop_run_metrics.run_metrics()
        
        # Check if run folder has required files
        if not (self.run_folder / "InterOp").exists():
            raise FileNotFoundError(f"InterOp folder not found in {self.run_folder}")
        
        try:
            self.run_metrics.read(str(self.run_folder))
            print(f"[INFO] Successfully loaded run metrics")
        except Exception as e:
            print(f"[ERROR] Failed to load run metrics: {e}")
            raise
        
        # Load imaging data
        try:
            self.imaging_df = self._load_imaging_table()
            print(f"[INFO] Loaded imaging table with {len(self.imaging_df)} rows")
        except Exception as e:
            print(f"[WARN] Could not load imaging table: {e}")
            self.imaging_df = None
        
        # Load summary data
        try:
            self.summary_df = self._load_summary_table()
            print(f"[INFO] Loaded summary table")
        except Exception as e:
            print(f"[WARN] Could not load summary table: {e}")
            self.summary_df = None
    
    def _has_valid_data(self, data, metric_name):
        """Check if data contains valid non-zero values"""
        if data is None or len(data) == 0:
            return False
        
        # Check if all values are NaN, zero, or invalid
        if isinstance(data, pd.Series):
            valid_data = data.dropna()
            if len(valid_data) == 0:
                return False
            if (valid_data == 0).all():
                return False
            # Check for meaningful variance
            if valid_data.std() < 1e-10:
                return False
        elif isinstance(data, pd.DataFrame):
            valid_data = data.dropna(how='all')
            if len(valid_data) == 0:
                return False
            # Check if all numeric columns are zero
            numeric_cols = valid_data.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) > 0:
                if (valid_data[numeric_cols] == 0).all().all():
                    return False
        
        return True
    
    def _load_imaging_table(self):
        """Load imaging table into pandas DataFrame"""
        ar = ic.imaging(str(self.run_folder))
        df = pd.DataFrame(ar)
        
        # Convert appropriate columns to integers
        int_cols = ['Lane', 'Tile', 'Tile Number', 'Read', 'Cycle', 
                    'Cycle Within Read', 'Swath', 'Surface']
        for col in int_cols:
            if col in df.columns:
                df[col] = df[col].astype(int)
        
        return df
    
    def _load_summary_table(self):
        """Load summary table"""
        try:
            summary_data = ic.summary(str(self.run_folder), level='Lane')
            return pd.DataFrame(summary_data)
        except:
            # Try alternative method
            summary_data = ic.summary(str(self.run_folder))
            return pd.DataFrame(summary_data)
    
    def plot_qscore_by_cycle(self):
        """Plot Q-score metrics by cycle"""
        print("[INFO] Generating Q-score by cycle plots...")
        
        if self.imaging_df is None:
            print("[WARN] No imaging data available")
            return
        
        metrics = ['%>= Q30', 'Error Rate']
        available_metrics = [m for m in metrics if m in self.imaging_df.columns]
        
        if not available_metrics:
            print("[WARN] Q-score metrics not available")
            return
        
        for metric in available_metrics:
            try:
                # Check if data is valid before plotting
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    print(f"[WARN] Skipping {metric} - no valid data")
                    continue
                
                fig, ax = plt.subplots(figsize=(12, 6))
                
                lanes = sorted(self.imaging_df['Lane'].unique())
                for lane in lanes:
                    lane_data = self.imaging_df[self.imaging_df['Lane'] == lane]
                    cycle_mean = lane_data.groupby('Cycle')[metric].mean()
                    
                    if not self._has_valid_data(cycle_mean, f"{metric}_lane{lane}"):
                        continue
                    
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                
                # Check if any lines were plotted
                if len(ax.lines) == 0:
                    plt.close()
                    print(f"[WARN] Skipping {metric} - no valid data to plot")
                    continue
                
                ax.set_xlabel('Cycle', fontsize=12)
                ax.set_ylabel(metric, fontsize=12)
                ax.set_title(f'{metric} by Cycle', fontsize=14, fontweight='bold')
                ax.legend()
                ax.grid(True, alpha=0.3)
                
                safe_name = metric.replace('>=', 'ge').replace('%', 'pct').replace(' ', '_')
                plot_path = self.output_dir / f'01_{safe_name}_by_cycle.png'
                plt.tight_layout()
                plt.savefig(plot_path, bbox_inches='tight')
                plt.close()
                self.plot_files.append((f'{metric} by Cycle', str(plot_path)))
                print(f"  ✓ Created {safe_name}_by_cycle.png")
            except Exception as e:
                print(f"[WARN] Could not create {metric} plot: {e}")
    
    def plot_density_metrics(self):
        """Plot cluster density metrics"""
        print("[INFO] Generating density plots...")
        
        if self.imaging_df is None:
            return
        
        density_cols = [col for col in self.imaging_df.columns 
                       if 'Density' in col or 'Cluster Count' in col]
        
        for metric in density_cols:
            try:
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    print(f"[WARN] Skipping {metric} - no valid data")
                    continue
                
                fig, ax = plt.subplots(figsize=(10, 6))
                
                lane_means = self.imaging_df.groupby('Lane')[metric].mean()
                lanes = lane_means.index.tolist()
                colors = [LANE_COLORS.get(lane, f'C{lane-1}') for lane in lanes]
                
                ax.bar(lanes, lane_means.values, alpha=0.7, color=colors)
                
                ax.set_xlabel('Lane', fontsize=12)
                ax.set_ylabel(metric, fontsize=12)
                ax.set_title(f'{metric} by Lane', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3, axis='y')
                
                safe_name = metric.replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '')
                plot_path = self.output_dir / f'02_{safe_name}_by_lane.png'
                plt.tight_layout()
                plt.savefig(plot_path, bbox_inches='tight')
                plt.close()
                self.plot_files.append((f'{metric} by Lane', str(plot_path)))
                print(f"  ✓ Created {safe_name}_by_lane.png")
            except Exception as e:
                print(f"[WARN] Could not create {metric} plot: {e}")
    
    def plot_intensity_metrics(self):
        """Plot intensity metrics by cycle"""
        print("[INFO] Generating intensity plots...")
        
        if self.imaging_df is None:
            return
        
        # Look for intensity-related columns
        intensity_cols = [col for col in self.imaging_df.columns 
                         if any(x in col for x in ['Corrected', 'Called', 'Intensity', 'P90'])]
        
        if not intensity_cols:
            print("[WARN] No intensity metrics found")
            return
        
        # Group by base if available
        base_metrics = {}
        for col in intensity_cols:
            if '/' in col:
                metric_type, base = col.rsplit('/', 1)
                if metric_type not in base_metrics:
                    base_metrics[metric_type] = []
                base_metrics[metric_type].append(col)
        
        for metric_type, cols in base_metrics.items():
            try:
                # Check if any of the columns have valid data
                has_valid = False
                for col in cols:
                    if self._has_valid_data(self.imaging_df[col], col):
                        has_valid = True
                        break
                
                if not has_valid:
                    print(f"[WARN] Skipping {metric_type} - no valid data")
                    continue
                
                fig, ax = plt.subplots(figsize=(12, 6))
                
                cycle_data = self.imaging_df.groupby('Cycle')[cols].mean()
                
                plotted_any = False
                for col in cols:
                    if not self._has_valid_data(cycle_data[col], col):
                        continue
                    
                    base = col.split('/')[-1]
                    color = CHANNEL_COLORS.get(base, None)
                    ax.plot(cycle_data.index, cycle_data[col], 
                           label=base, linewidth=2, alpha=0.8, color=color)
                    plotted_any = True
                
                if not plotted_any:
                    plt.close()
                    print(f"[WARN] Skipping {metric_type} - no valid data to plot")
                    continue
                
                ax.set_xlabel('Cycle', fontsize=12)
                ax.set_ylabel('Intensity', fontsize=12)
                ax.set_title(f'{metric_type} by Cycle', fontsize=14, fontweight='bold')
                ax.legend()
                ax.grid(True, alpha=0.3)
                
                safe_name = metric_type.replace('/', '_').replace(' ', '_')
                plot_path = self.output_dir / f'03_{safe_name}_by_cycle.png'
                plt.tight_layout()
                plt.savefig(plot_path, bbox_inches='tight')
                plt.close()
                self.plot_files.append((f'{metric_type} by Cycle', str(plot_path)))
                print(f"  ✓ Created {safe_name}_by_cycle.png")
            except Exception as e:
                print(f"[WARN] Could not create {metric_type} plot: {e}")
    
    def plot_base_composition(self):
        """Plot base composition by cycle"""
        print("[INFO] Generating base composition plots...")
        
        if self.imaging_df is None:
            return
        
        base_cols = [col for col in self.imaging_df.columns if '% Base' in col]
        
        if not base_cols:
            print("[WARN] Base composition data not available")
            return
        
        try:
            # Check if any columns have valid data
            has_valid = False
            for col in base_cols:
                if self._has_valid_data(self.imaging_df[col], col):
                    has_valid = True
                    break
            
            if not has_valid:
                print("[WARN] Skipping base composition - no valid data")
                return
            
            fig, ax = plt.subplots(figsize=(12, 6))
            
            cycle_data = self.imaging_df.groupby('Cycle')[base_cols].mean()
            
            plotted_any = False
            for col in base_cols:
                if not self._has_valid_data(cycle_data[col], col):
                    continue
                
                base = col.split('/')[-1]
                color = CHANNEL_COLORS.get(base, None)
                ax.plot(cycle_data.index, cycle_data[col], 
                       label=f'Base {base}', linewidth=2, alpha=0.8, color=color)
                plotted_any = True
            
            if not plotted_any:
                plt.close()
                print("[WARN] Skipping base composition - no valid data to plot")
                return
            
            ax.set_xlabel('Cycle', fontsize=12)
            ax.set_ylabel('% Base', fontsize=12)
            ax.set_title('Base Composition by Cycle', fontsize=14, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_ylim([0, 100])
            
            plot_path = self.output_dir / '04_base_composition_by_cycle.png'
            plt.tight_layout()
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()
            self.plot_files.append(('Base Composition by Cycle', str(plot_path)))
            print("  ✓ Created base_composition_by_cycle.png")
        except Exception as e:
            print(f"[WARN] Could not create base composition plot: {e}")
    
    def plot_fwhm_metrics(self):
        """Plot FWHM (focus) metrics"""
        print("[INFO] Generating FWHM plots...")
        
        if self.imaging_df is None:
            return
        
        fwhm_cols = [col for col in self.imaging_df.columns if 'Fwhm' in col or 'FWHM' in col]
        
        if not fwhm_cols:
            print("[WARN] FWHM data not available")
            return
        
        try:
            # Check if any columns have valid data
            has_valid = False
            for col in fwhm_cols:
                if self._has_valid_data(self.imaging_df[col], col):
                    has_valid = True
                    break
            
            if not has_valid:
                print("[WARN] Skipping FWHM - no valid data")
                return
            
            fig, ax = plt.subplots(figsize=(12, 6))
            
            cycle_data = self.imaging_df.groupby('Cycle')[fwhm_cols].mean()
            
            plotted_any = False
            for col in fwhm_cols:
                if not self._has_valid_data(cycle_data[col], col):
                    continue
                
                if '/' in col:
                    channel = col.split('/')[-1].upper()
                    label = f'Channel {channel}'
                    # Map channel to color: A=Green, C=Blue
                    if channel == 'A' or 'green' in col.lower():
                        color = '#2ca02c'  # Green
                    elif channel == 'C' or 'blue' in col.lower():
                        color = '#1f77b4'  # Blue
                    else:
                        color = CHANNEL_COLORS.get(channel, None)
                else:
                    label = col
                    color = None
                
                ax.plot(cycle_data.index, cycle_data[col], 
                       label=label, linewidth=2, alpha=0.8, color=color)
                plotted_any = True
            
            if not plotted_any:
                plt.close()
                print("[WARN] Skipping FWHM - no valid data to plot")
                return
            
            ax.set_xlabel('Cycle', fontsize=12)
            ax.set_ylabel('FWHM', fontsize=12)
            ax.set_title('FWHM (Focus Metric) by Cycle', fontsize=14, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            plot_path = self.output_dir / '05_fwhm_by_cycle.png'
            plt.tight_layout()
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()
            self.plot_files.append(('FWHM by Cycle', str(plot_path)))
            print("  ✓ Created fwhm_by_cycle.png")
        except Exception as e:
            print(f"[WARN] Could not create FWHM plot: {e}")
    
    def plot_occupancy(self):
        """Plot occupancy metrics (for patterned flowcells)"""
        print("[INFO] Generating occupancy plots...")
        
        if self.imaging_df is None:
            return
        
        occ_cols = [col for col in self.imaging_df.columns 
                   if 'Occupancy' in col or 'Occupied' in col]
        
        if not occ_cols:
            print("[WARN] Occupancy data not available (may not be patterned flowcell)")
            return
        
        for metric in occ_cols:
            try:
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    print(f"[WARN] Skipping {metric} - no valid data")
                    continue
                
                fig, ax = plt.subplots(figsize=(12, 6))
                
                lanes = sorted(self.imaging_df['Lane'].unique())
                plotted_any = False
                
                for lane in lanes:
                    lane_data = self.imaging_df[self.imaging_df['Lane'] == lane]
                    cycle_mean = lane_data.groupby('Cycle')[metric].mean()
                    
                    if not self._has_valid_data(cycle_mean, f"{metric}_lane{lane}"):
                        continue
                    
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                    plotted_any = True
                
                if not plotted_any:
                    plt.close()
                    print(f"[WARN] Skipping {metric} - no valid data to plot")
                    continue
                
                ax.set_xlabel('Cycle', fontsize=12)
                ax.set_ylabel(metric, fontsize=12)
                ax.set_title(f'{metric} by Cycle', fontsize=14, fontweight='bold')
                ax.legend()
                ax.grid(True, alpha=0.3)
                
                safe_name = metric.replace('%', 'pct').replace(' ', '_')
                plot_path = self.output_dir / f'06_{safe_name}_by_cycle.png'
                plt.tight_layout()
                plt.savefig(plot_path, bbox_inches='tight')
                plt.close()
                self.plot_files.append((f'{metric} by Cycle', str(plot_path)))
                print(f"  ✓ Created {safe_name}_by_cycle.png")
            except Exception as e:
                print(f"[WARN] Could not create {metric} plot: {e}")
    
    def plot_phasing_metrics(self):
        """Plot phasing and prephasing metrics"""
        print("[INFO] Generating phasing plots...")
        
        if self.imaging_df is None:
            return
        
        phasing_cols = [col for col in self.imaging_df.columns 
                       if 'Phasing' in col or 'Prephasing' in col]
        
        if not phasing_cols:
            print("[WARN] Phasing data not available")
            return
        
        try:
            # Check if any columns have valid data
            has_valid = False
            for col in phasing_cols:
                if self._has_valid_data(self.imaging_df[col], col):
                    has_valid = True
                    break
            
            if not has_valid:
                print("[WARN] Skipping phasing metrics - no valid data")
                return
            
            fig, ax = plt.subplots(figsize=(12, 6))
            
            plotted_any = False
            for metric in phasing_cols:
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    continue
                
                cycle_mean = self.imaging_df.groupby('Cycle')[metric].mean()
                ax.plot(cycle_mean.index, cycle_mean.values, 
                       label=metric, linewidth=2, alpha=0.8, marker='o', markersize=3)
                plotted_any = True
            
            if not plotted_any:
                plt.close()
                print("[WARN] Skipping phasing metrics - no valid data to plot")
                return
            
            ax.set_xlabel('Cycle', fontsize=12)
            ax.set_ylabel('Percentage', fontsize=12)
            ax.set_title('Phasing/Prephasing by Cycle', fontsize=14, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            plot_path = self.output_dir / '07_phasing_by_cycle.png'
            plt.tight_layout()
            plt.savefig(plot_path, bbox_inches='tight')
            plt.close()
            self.plot_files.append(('Phasing by Cycle', str(plot_path)))
            print("  ✓ Created phasing_by_cycle.png")
        except Exception as e:
            print(f"[WARN] Could not create phasing plot: {e}")
    
    def plot_passfilter_distribution(self):
        """Plot pass filter distribution"""
        print("[INFO] Generating pass filter plots...")
        
        if self.imaging_df is None:
            return
        
        pf_cols = [col for col in self.imaging_df.columns if 'Pass Filter' in col]
        
        if not pf_cols:
            print("[WARN] Pass filter data not available")
            return
        
        for metric in pf_cols:
            try:
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    print(f"[WARN] Skipping {metric} - no valid data")
                    continue
                
                fig, axes = plt.subplots(1, 2, figsize=(15, 6))
                
                # By cycle
                lanes = sorted(self.imaging_df['Lane'].unique())
                plotted_cycle = False
                
                for lane in lanes:
                    lane_data = self.imaging_df[self.imaging_df['Lane'] == lane]
                    cycle_mean = lane_data.groupby('Cycle')[metric].mean()
                    
                    if not self._has_valid_data(cycle_mean, f"{metric}_lane{lane}"):
                        continue
                    
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    axes[0].plot(cycle_mean.index, cycle_mean.values, 
                               label=f'Lane {lane}', linewidth=2, alpha=0.8, color=color)
                    plotted_cycle = True
                
                if not plotted_cycle:
                    plt.close()
                    print(f"[WARN] Skipping {metric} - no valid data to plot")
                    continue
                
                axes[0].set_xlabel('Cycle', fontsize=12)
                axes[0].set_ylabel(metric, fontsize=12)
                axes[0].set_title(f'{metric} by Cycle', fontsize=12, fontweight='bold')
                axes[0].legend()
                axes[0].grid(True, alpha=0.3)
                
                # By lane (box plot)
                lane_data_list = [self.imaging_df[self.imaging_df['Lane'] == lane][metric].values 
                                 for lane in lanes]
                bp = axes[1].boxplot(lane_data_list, labels=lanes, patch_artist=True)
                
                # Color boxes by lane
                for patch, lane in zip(bp['boxes'], lanes):
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                axes[1].set_xlabel('Lane', fontsize=12)
                axes[1].set_ylabel(metric, fontsize=12)
                axes[1].set_title(f'{metric} Distribution by Lane', fontsize=12, fontweight='bold')
                axes[1].grid(True, alpha=0.3, axis='y')
                
                plot_path = self.output_dir / '08_passfilter_analysis.png'
                plt.tight_layout()
                plt.savefig(plot_path, bbox_inches='tight')
                plt.close()
                self.plot_files.append((f'{metric} Analysis', str(plot_path)))
                print("  ✓ Created passfilter_analysis.png")
            except Exception as e:
                print(f"[WARN] Could not create pass filter plot: {e}")
    
    def plot_tile_heatmap(self):
        """Generate heatmap of metrics across tiles"""
        print("[INFO] Generating tile heatmaps...")
        
        if self.imaging_df is None or 'Tile' not in self.imaging_df.columns:
            return
        
        # Select key metrics for heatmap
        heatmap_metrics = []
        for metric in ['%>= Q30', '% Pass Filter', 'Cluster Count PF (k)', 'Error Rate']:
            if metric in self.imaging_df.columns:
                heatmap_metrics.append(metric)
        
        if not heatmap_metrics:
            print("[WARN] No suitable metrics for tile heatmap")
            return
        
        for metric in heatmap_metrics:
            try:
                # Get first cycle data for simplicity
                first_cycle = self.imaging_df[self.imaging_df['Cycle'] == 1].copy()
                
                if len(first_cycle) == 0:
                    continue
                
                pivot_data = first_cycle.pivot_table(
                    values=metric,
                    index='Tile',
                    columns='Lane',
                    aggfunc='mean'
                )
                
                fig, ax = plt.subplots(figsize=(12, 10))
                sns.heatmap(pivot_data, annot=False, cmap='RdYlGn', 
                           cbar_kws={'label': metric}, ax=ax)
                ax.set_title(f'{metric} Across Tiles (Cycle 1)', 
                            fontsize=14, fontweight='bold')
                ax.set_xlabel('Lane', fontsize=12)
                ax.set_ylabel('Tile', fontsize=12)
                
                safe_name = metric.replace('>=', 'ge').replace('%', 'pct').replace(' ', '_').replace('(', '').replace(')', '')
                plot_path = self.output_dir / f'09_tile_heatmap_{safe_name}.png'
                plt.tight_layout()
                plt.savefig(plot_path, bbox_inches='tight')
                plt.close()
                self.plot_files.append((f'Tile Heatmap: {metric}', str(plot_path)))
                print(f"  ✓ Created tile_heatmap_{safe_name}.png")
            except Exception as e:
                print(f"[WARN] Could not create tile heatmap for {metric}: {e}")
    
    def generate_summary_report(self):
        """Generate text summary report"""
        print("[INFO] Generating summary report...")
        
        report_path = self.output_dir / 'summary_report.txt'
        
        with open(report_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("ILLUMINA SEQUENCING RUN QC REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Run Folder: {self.run_folder.name}\n")
            f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Run info
            try:
                run_info = self.run_metrics.run_info()
                f.write(f"Flowcell ID: {run_info.flowcell().name()}\n")
                f.write(f"Total Cycles: {run_info.total_cycles()}\n")
                f.write(f"Reads: {run_info.reads().size()}\n")
                
                for i in range(run_info.reads().size()):
                    read = run_info.reads()[i]
                    f.write(f"  Read {i+1}: {read.total_cycles()} cycles")
                    if read.is_index():
                        f.write(" (Index)")
                    f.write("\n")
                f.write("\n")
            except Exception as e:
                f.write(f"Run info unavailable: {e}\n\n")
            
            # Summary statistics
            if self.summary_df is not None:
                f.write("SUMMARY STATISTICS\n")
                f.write("-" * 80 + "\n")
                f.write(self.summary_df.to_string())
                f.write("\n\n")
            
            # Imaging statistics
            if self.imaging_df is not None:
                f.write("IMAGING METRICS SUMMARY\n")
                f.write("-" * 80 + "\n")
                
                numeric_cols = self.imaging_df.select_dtypes(include=[np.number]).columns
                summary_stats = self.imaging_df[numeric_cols].describe()
                f.write(summary_stats.to_string())
                f.write("\n\n")
        
        print(f"  ✓ Created summary_report.txt")
    
    def export_csv_tables(self):
        """Export data tables to CSV"""
        print("[INFO] Exporting CSV tables...")
        
        if self.imaging_df is not None:
            csv_path = self.output_dir / 'imaging_metrics.csv'
            self.imaging_df.to_csv(csv_path, index=False)
            print(f"  ✓ Created imaging_metrics.csv")
        
        if self.summary_df is not None:
            csv_path = self.output_dir / 'summary_metrics.csv'
            self.summary_df.to_csv(csv_path, index=False)
            print(f"  ✓ Created summary_metrics.csv")
    
    def generate_pdf_report(self):
        """Generate comprehensive PDF report with all plots"""
        print("[INFO] Generating PDF report...")
        
        pdf_path = self.output_dir / 'QC_Report.pdf'
        doc = SimpleDocTemplate(str(pdf_path), pagesize=letter,
                               rightMargin=0.5*inch, leftMargin=0.5*inch,
                               topMargin=0.5*inch, bottomMargin=0.5*inch)
        
        story = []
        styles = getSampleStyleSheet()
        
        # Custom styles
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#1f77b4'),
            spaceAfter=30,
            alignment=TA_CENTER
        )
        
        heading_style = ParagraphStyle(
            'CustomHeading',
            parent=styles['Heading2'],
            fontSize=16,
            textColor=colors.HexColor('#2ca02c'),
            spaceAfter=12,
            spaceBefore=12
        )
        
        # Title page
        story.append(Spacer(1, 1*inch))
        story.append(Paragraph("Illumina Sequencing Run", title_style))
        story.append(Paragraph("Quality Control Report", title_style))
        story.append(Spacer(1, 0.5*inch))
        
        # Run information
        run_info_data = [
            ['Run Folder:', self.run_folder.name],
            ['Analysis Date:', datetime.now().strftime('%Y-%m-%d %H:%M:%S')],
        ]
        
        # Add run metrics if available
        try:
            run_info = self.run_metrics.run_info()
            run_info_data.extend([
                ['Flowcell ID:', run_info.flowcell().name()],
                ['Total Cycles:', str(run_info.total_cycles())],
                ['Number of Reads:', str(run_info.reads().size())],
            ])
        except:
            pass
        
        info_table = Table(run_info_data, colWidths=[2*inch, 4*inch])
        info_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
            ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        
        #story.append(info_table)
        story.append(PageBreak())
        
        # Add summary statistics if available
        if self.summary_df is not None and len(self.summary_df) > 0:
            #story.append(Paragraph("Summary Statistics", heading_style))
            story.append(Spacer(1, 0.2*inch))
            
            # Convert summary dataframe to table
            summary_data = [self.summary_df.columns.tolist()] + self.summary_df.values.tolist()
            
            # Limit columns if too many
            if len(summary_data[0]) > 6:
                # Select most important columns
                important_cols = ['Lane', 'Yield', 'Projected Yield', '% >= Q30', 'Error Rate']
                available_cols = [col for col in important_cols if col in self.summary_df.columns]
                if available_cols:
                    summary_df_subset = self.summary_df[available_cols]
                    summary_data = [summary_df_subset.columns.tolist()] + summary_df_subset.values.tolist()
            
            # Format data
            for i in range(len(summary_data)):
                for j in range(len(summary_data[i])):
                    val = summary_data[i][j]
                    if isinstance(val, float):
                        summary_data[i][j] = f'{val:.2f}'
                    else:
                        summary_data[i][j] = str(val)
            
            col_width = 6.5*inch / len(summary_data[0])
            summary_table = Table(summary_data, colWidths=[col_width] * len(summary_data[0]))
            summary_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 8),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)
            ]))
            
            #story.append(summary_table)
            #story.append(PageBreak())
        
        # Add all plots
        if self.plot_files:
            story.append(Paragraph("Quality Control Plots", heading_style))
            story.append(Spacer(1, 0.2*inch))
            
            for plot_title, plot_path in self.plot_files:
                if os.path.exists(plot_path):
                    # Add plot title
                    story.append(Paragraph(plot_title, styles['Heading3']))
                    story.append(Spacer(1, 0.1*inch))
                    
                    # Add image - scale to fit page width
                    img = Image(plot_path, width=6.5*inch, height=4*inch)
                    story.append(img)
                    story.append(Spacer(1, 0.3*inch))
                    
                    # Page break after every 2 plots
                    if self.plot_files.index((plot_title, plot_path)) % 2 == 1:
                        story.append(PageBreak())
        
        # Build PDF
        doc.build(story)
        print(f"  ✓ Created QC_Report.pdf")
        print(f"\n[INFO] PDF report saved to: {pdf_path}")
    
    def generate_all_plots(self):
        """Generate all available QC plots"""
        print("\n" + "=" * 80)
        print("GENERATING QC PLOTS")
        print("=" * 80 + "\n")
        
        self.plot_qscore_by_cycle()
        self.plot_density_metrics()
        self.plot_intensity_metrics()
        self.plot_base_composition()
        self.plot_fwhm_metrics()
        self.plot_occupancy()
        self.plot_phasing_metrics()
        self.plot_passfilter_distribution()
        self.plot_tile_heatmap()
        self.generate_summary_report()
        self.export_csv_tables()
        self.generate_pdf_report()
        
        print("\n" + "=" * 80)
        print("QC PLOT GENERATION COMPLETE")
        print("=" * 80)
        print(f"\nOutput directory: {self.output_dir}")
        print(f"\nGenerated files:")
        for f in sorted(self.output_dir.glob('*')):
            print(f"  - {f.name}")
        print()


def main():
    if len(sys.argv) != 3:
        print("Usage: python interop_qc_plots.py <run_folder> <output_dir>")
        print("\nExample:")
        print("  python interop_qc_plots.py /data/210101_A00000_0001_AHXXXX ./qc_plots")
        sys.exit(1)
    
    run_folder = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Validate run folder
    if not os.path.exists(run_folder):
        print(f"[ERROR] Run folder does not exist: {run_folder}")
        sys.exit(1)
    
    if not os.path.exists(os.path.join(run_folder, "InterOp")):
        print(f"[ERROR] InterOp folder not found in: {run_folder}")
        sys.exit(1)
    
    try:
        plotter = InterOpQCPlotter(run_folder, output_dir)
        plotter.generate_all_plots()
    except Exception as e:
        print(f"\n[ERROR] Failed to generate plots: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()