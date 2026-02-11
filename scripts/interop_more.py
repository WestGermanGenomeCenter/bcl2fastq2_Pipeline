#!/usr/bin/env python3
"""
Illumina InterOp QC Plot Generator using Python
Compatible with NovaSeqX+, NextSeq, NextSeq2000, and MiSeq
Generates comprehensive quality control plots for sequencing runs

Usage:
    python interop_more.py <run_folder> <output_pdf>

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
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from interop import py_interop_run_metrics, py_interop_run
    import interop.core as ic
    from reportlab.lib.pagesizes import letter, landscape
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
    from datetime import datetime
    from io import BytesIO
except ImportError as e:
    print(f"Error: Missing required package: {e}")
    print("\nPlease install required packages:")
    print("  conda install -c bioconda illumina-interop")
    print("  conda install pandas matplotlib seaborn numpy reportlab")
    sys.exit(1)

# Set style similar to Illumina SAV
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9

# Illumina SAV color palette
LANE_COLORS = {
    1: '#0072B2',
    2: '#E69F00',
    3: '#009E73',
    4: '#D55E00',
    5: '#CC79A7',
    6: '#F0E442',
    7: '#56B4E9',
    8: '#999999',
}

CHANNEL_COLORS = {
    'A': '#009E73',
    'C': '#0072B2',
    'G': '#000000',
    'T': '#D55E00',
}

class InterOpQCPlotter:
    """Generate comprehensive QC plots from Illumina InterOp data"""
    
    def __init__(self, run_folder):
        self.run_folder = Path(run_folder)
        self.plot_images = []  # List of (title, description, image_data)
        
        print(f"[INFO] Loading run metrics from: {self.run_folder}")
        self.run_metrics = py_interop_run_metrics.run_metrics()
        
        if not (self.run_folder / "InterOp").exists():
            raise FileNotFoundError(f"InterOp folder not found in {self.run_folder}")
        
        try:
            self.run_metrics.read(str(self.run_folder))
            print(f"[INFO] Successfully loaded run metrics")
        except Exception as e:
            print(f"[ERROR] Failed to load run metrics: {e}")
            raise
        
        try:
            self.imaging_df = self._load_imaging_table()
            print(f"[INFO] Loaded imaging table with {len(self.imaging_df)} rows")
        except Exception as e:
            print(f"[WARN] Could not load imaging table: {e}")
            self.imaging_df = None
        
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
        
        if isinstance(data, pd.Series):
            valid_data = data.dropna()
            if len(valid_data) == 0:
                return False
            if (valid_data == 0).all():
                return False
            if valid_data.std() < 1e-10:
                return False
        elif isinstance(data, pd.DataFrame):
            valid_data = data.dropna(how='all')
            if len(valid_data) == 0:
                return False
            numeric_cols = valid_data.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) > 0:
                if (valid_data[numeric_cols] == 0).all().all():
                    return False
        
        return True
    
    def _load_imaging_table(self):
        """Load imaging table into pandas DataFrame"""
        ar = ic.imaging(str(self.run_folder))
        df = pd.DataFrame(ar)
        
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
            try:
                summary_data = ic.summary(str(self.run_folder))
                return pd.DataFrame(summary_data)
            except:
                return None
    
    def _save_fig_to_image(self, fig):
        """Save figure to PNG image for PDF embedding"""
        buf = BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
        buf.seek(0)
        plt.close(fig)
        return buf.getvalue()
    
    def _get_tile_metrics(self):
        """Extract tile metrics from run_metrics"""
        tile_metric_set = self.run_metrics.tile_metric_set()
        tile_data = []
        
        for i in range(tile_metric_set.size()):
            metric = tile_metric_set.at(i)
            tile_data.append({
                'Lane': metric.lane(),
                'Tile': metric.tile(),
                'Cluster Count': metric.cluster_count(),
                'Cluster Count PF': metric.cluster_count_pf(),
                'Cluster Density': metric.cluster_density(),
                'Cluster Density PF': metric.cluster_density_pf(),
                'Percent PF': metric.percent_pf(),
                'Percent Occupied': metric.percent_occupied() if hasattr(metric, 'percent_occupied') else np.nan,
            })
        
        return pd.DataFrame(tile_data) if tile_data else None
    
    def _get_max_lanes(self):
        """Get maximum number of lanes in flowcell"""
        try:
            return self.run_metrics.run_info().flowcell().lane_count()
        except:
            return 8
    
    def plot_summary_page(self):
        """Create summary page with overall run statistics"""
        print("[INFO] Generating summary page...")
        
        try:
            tile_df = self._get_tile_metrics()
            if tile_df is None:
                return
            
            total_clusters = tile_df['Cluster Count'].sum()
            pf_clusters = tile_df['Cluster Count PF'].sum()
            non_pf_clusters = total_clusters - pf_clusters
            pct_pf = (pf_clusters / total_clusters * 100) if total_clusters > 0 else 0
            
            fig = plt.figure(figsize=(10, 6))
            ax = fig.add_subplot(111)
            ax.axis('off')
            
            # Summary text
            summary_text = f"""
ILLUMINA RUN SUMMARY REPORT
Run Folder: {self.run_folder.name}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

OVERALL CLUSTER STATISTICS
{'─' * 60}

Total Clusters:                    {total_clusters:,}
Clusters Passing Filter (PF):      {pf_clusters:,}
Clusters Failing Filter:           {non_pf_clusters:,}

Percent PF:                        {pct_pf:.2f}%
Percent Non-PF:                    {100 - pct_pf:.2f}%

QUALITY ASSESSMENT
{'─' * 60}
"""
            if pct_pf >= 80:
                summary_text += "✓ EXCELLENT - >80% clusters passing filter\n"
            elif pct_pf >= 70:
                summary_text += "◐ GOOD - 70-80% clusters passing filter\n"
            else:
                summary_text += "✗ WARNING - <70% clusters passing filter\n"
            
            summary_text += f"""
PER-LANE BREAKDOWN
{'─' * 60}
"""
            lane_groups = tile_df.groupby('Lane').agg({
                'Cluster Count': 'sum',
                'Cluster Count PF': 'sum'
            }).reset_index()
            
            for _, row in lane_groups.iterrows():
                lane = int(row['Lane'])
                total = int(row['Cluster Count'])
                pf = int(row['Cluster Count PF'])
                pct = (pf / total * 100) if total > 0 else 0
                summary_text += f"Lane {lane}: Total={total:,}  PF={pf:,}  ({pct:.2f}%)\n"
            
            ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
                   fontsize=10, verticalalignment='top', fontfamily='monospace',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
            
            fig.tight_layout()
            description = "Overall run summary showing cluster statistics and per-lane quality breakdown."
            self.plot_images.append(('Run Summary', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created summary page")
        except Exception as e:
            print(f"[WARN] Could not create summary page: {e}")
    
    def plot_cluster_distribution(self):
        """Plot cluster count distribution across lanes"""
        print("[INFO] Generating cluster distribution plots...")
        
        try:
            tile_df = self._get_tile_metrics()
            if tile_df is None:
                return
            
            lanes = sorted(tile_df['Lane'].unique())
            max_lanes = self._get_max_lanes()
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            lane_groups = tile_df.groupby('Lane').agg({
                'Cluster Count': 'sum',
                'Cluster Count PF': 'sum'
            }).reset_index()
            
            total_counts = lane_groups['Cluster Count'].values
            pf_counts = lane_groups['Cluster Count PF'].values
            non_pf_counts = total_counts - pf_counts
            lane_list = lane_groups['Lane'].values.astype(int)
            
            x = np.arange(len(lane_list))
            width = 0.6
            
            colors_list = [LANE_COLORS.get(lane, f'C{lane-1}') for lane in lane_list]
            
            ax.bar(x, pf_counts, width, label='Passing Filter', color='#2ecc71', alpha=0.8, edgecolor='black', linewidth=1)
            ax.bar(x, non_pf_counts, width, bottom=pf_counts, label='Failing Filter', 
                   color='#e74c3c', alpha=0.8, edgecolor='black', linewidth=1)
            
            ax.set_xlabel('Lane', fontsize=11, fontweight='bold')
            ax.set_ylabel('Cluster Count', fontsize=11, fontweight='bold')
            ax.set_title('Stacked Cluster Distribution: PF vs Non-PF', fontsize=13, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([f'Lane {int(l)}' for l in lane_list])
            ax.legend(loc='upper right', framealpha=0.95)
            ax.grid(True, alpha=0.3, axis='y', linestyle='--')
            ax.set_axisbelow(True)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x/1e6):.0f}M' if x >= 1e6 else f'{int(x/1e3):.0f}K'))
            
            fig.tight_layout()
            description = "Stacked bar chart showing the proportion of clusters passing versus failing quality filters for each lane."
            self.plot_images.append(('Cluster Distribution', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created cluster distribution plot")
        except Exception as e:
            print(f"[WARN] Could not create cluster distribution plot: {e}")
    
    def plot_overall_cluster_pie(self):
        """Plot pie chart showing overall PF vs non-PF distribution"""
        print("[INFO] Generating cluster breakdown pie chart...")
        
        try:
            tile_df = self._get_tile_metrics()
            if tile_df is None:
                return
            
            pf = tile_df['Cluster Count PF'].sum()
            non_pf = tile_df['Cluster Count'].sum() - pf
            total = pf + non_pf
            pct_pf = (pf / total * 100) if total > 0 else 0
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            sizes = [pf, non_pf]
            labels = [f'Passing Filter\n{pf:,}\n({pct_pf:.2f}%)', 
                     f'Failing Filter\n{non_pf:,}\n({100 - pct_pf:.2f}%)']
            colors = ['#2ecc71', '#e74c3c']
            explode = (0.05, 0)
            
            wedges, texts, autotexts = ax.pie(sizes, explode=explode, labels=labels, colors=colors, 
                                               autopct='', shadow=True, startangle=90,
                                               textprops={'fontsize': 11, 'fontweight': 'bold'})
            ax.axis('equal')
            
            fig.tight_layout()
            description = "Overall quality distribution showing the proportion of clusters passing versus failing quality filters across the entire run."
            self.plot_images.append(('Overall Cluster Quality', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created cluster breakdown pie chart")
        except Exception as e:
            print(f"[WARN] Could not create pie chart: {e}")
    
    def plot_lane_balance(self):
        """Plot lane balance showing deviation from mean cluster count"""
        print("[INFO] Generating lane balance plot...")
        
        try:
            tile_df = self._get_tile_metrics()
            if tile_df is None:
                return
            
            lanes = sorted(tile_df['Lane'].unique())
            max_lanes = self._get_max_lanes()
            
            lane_groups = tile_df.groupby('Lane')['Cluster Count'].sum()
            total_counts = [lane_groups.get(lane, 0) for lane in lanes]
            
            mean_count = np.mean(total_counts)
            deviations = [(count - mean_count) / mean_count * 100 for count in total_counts]
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            colors = ['#2ecc71' if abs(d) < 10 else '#f39c12' if abs(d) < 20 else '#e74c3c' 
                     for d in deviations]
            
            x = np.arange(len(lanes))
            ax.bar(x, deviations, color=colors, alpha=0.8, edgecolor='black', linewidth=1)
            ax.axhline(y=0, color='black', linewidth=1)
            ax.axhline(y=10, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, label='±10% threshold')
            ax.axhline(y=-10, color='orange', linestyle='--', linewidth=1.5, alpha=0.7)
            
            ax.set_xlabel('Lane', fontsize=11, fontweight='bold')
            ax.set_ylabel('Deviation from Mean (%)', fontsize=11, fontweight='bold')
            ax.set_title('Lane Balance: Cluster Count Deviation from Mean', fontsize=13, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([f'Lane {int(l)}' for l in lanes])
            ax.legend(loc='upper right', framealpha=0.95)
            ax.grid(True, alpha=0.3, axis='y', linestyle='--')
            ax.set_axisbelow(True)
            
            for i, d in enumerate(deviations):
                ax.text(i, d + (1 if d >= 0 else -1), f'{d:.1f}%', 
                       ha='center', va='bottom' if d >= 0 else 'top', fontsize=9, fontweight='bold')
            
            fig.tight_layout()
            description = "Lane balance analysis showing deviation from mean cluster count; lanes within ±10% indicate well-balanced sequencing."
            self.plot_images.append(('Lane Balance', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created lane balance plot")
        except Exception as e:
            print(f"[WARN] Could not create lane balance plot: {e}")
    
    def plot_percent_pf_by_lane(self):
        """Plot percent PF by lane"""
        print("[INFO] Generating percent PF by lane plot...")
        
        try:
            tile_df = self._get_tile_metrics()
            if tile_df is None:
                return
            
            lanes = sorted(tile_df['Lane'].unique())
            max_lanes = self._get_max_lanes()
            
            lane_groups = tile_df.groupby('Lane').agg({
                'Cluster Count': 'sum',
                'Cluster Count PF': 'sum'
            }).reset_index()
            
            percent_pf = [(row['Cluster Count PF'] / row['Cluster Count'] * 100) 
                         if row['Cluster Count'] > 0 else 0 
                         for _, row in lane_groups.iterrows()]
            lane_list = lane_groups['Lane'].values.astype(int)
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            x = np.arange(len(lane_list))
            colors = ['#2ecc71' if p >= 80 else '#f39c12' if p >= 70 else '#e74c3c' for p in percent_pf]
            
            bars = ax.bar(x, percent_pf, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5, width=0.6)
            ax.axhline(y=80, color='green', linestyle='--', linewidth=1.5, alpha=0.7, label='Excellent (80%)')
            ax.axhline(y=70, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, label='Good (70%)')
            
            ax.set_xlabel('Lane', fontsize=11, fontweight='bold')
            ax.set_ylabel('% Passing Filter', fontsize=11, fontweight='bold')
            ax.set_title('Percentage of Clusters Passing Filter by Lane', fontsize=13, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([f'Lane {int(l)}' for l in lane_list])
            ax.set_ylim([0, 100])
            ax.legend(loc='upper right', framealpha=0.95)
            ax.grid(True, alpha=0.3, axis='y', linestyle='--')
            ax.set_axisbelow(True)
            
            for i, v in enumerate(percent_pf):
                ax.text(i, v + 2, f'{v:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
            
            fig.tight_layout()
            description = "Percent passing filter per lane; >80% indicates excellent quality, 70-80% is good, <70% may indicate issues."
            self.plot_images.append(('% PF by Lane', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created percent PF plot")
        except Exception as e:
            print(f"[WARN] Could not create percent PF plot: {e}")
    
    def plot_cluster_density_heatmap(self):
        """Plot heatmap of cluster density across tiles"""
        print("[INFO] Generating cluster density heatmap...")
        
        try:
            tile_df = self._get_tile_metrics()
            if tile_df is None:
                return
            
            lanes = sorted(tile_df['Lane'].unique())
            num_lanes = len(lanes)
            
            if num_lanes <= 2:
                nrows, ncols = 1, num_lanes
                figsize = (14, 5)
            elif num_lanes <= 4:
                nrows, ncols = 2, 2
                figsize = (14, 10)
            else:
                nrows = (num_lanes + 1) // 2
                ncols = 2
                figsize = (14, 5 * nrows)
            
            fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
            fig.suptitle('Cluster Density (PF Clusters) Distribution by Tile', fontsize=13, fontweight='bold', y=0.995)
            
            for idx, lane in enumerate(lanes):
                row = idx // ncols
                col = idx % ncols
                ax = axes[row, col]
                
                lane_data = tile_df[tile_df['Lane'] == lane].sort_values('Tile')
                tiles = lane_data['Tile'].values.astype(int)
                densities = lane_data['Cluster Count PF'].values
                
                tile_positions = np.arange(len(tiles))
                cmap = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(tiles)))
                
                bars = ax.bar(tile_positions, densities, color=cmap, alpha=0.8, edgecolor='black', linewidth=0.5)
                ax.set_xlabel('Tile Index', fontsize=10, fontweight='bold')
                ax.set_ylabel('PF Cluster Count', fontsize=10, fontweight='bold')
                ax.set_title(f'Lane {int(lane)}', fontsize=11, fontweight='bold')
                ax.grid(True, alpha=0.3, axis='y', linestyle='--')
                ax.set_axisbelow(True)
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x/1e3):.0f}K'))
                
                mean_density = np.mean(densities)
                ax.axhline(y=mean_density, color='blue', linestyle=':', linewidth=1.5, alpha=0.7, label=f'Mean: {mean_density:.0f}')
                ax.legend(fontsize=8, loc='upper right')
            
            for idx in range(num_lanes, nrows * ncols):
                row = idx // ncols
                col = idx % ncols
                axes[row, col].axis('off')
            
            fig.tight_layout()
            description = "Tile-level cluster density heatmap showing uniform distribution indicates good flowcell quality; outliers suggest optical or chemistry issues."
            self.plot_images.append(('Cluster Density Heatmap', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created cluster density heatmap")
        except Exception as e:
            print(f"[WARN] Could not create cluster density heatmap: {e}")
    
    def plot_qscore_by_cycle(self):
        """Plot Q-score metrics by cycle"""
        print("[INFO] Generating Q-score by cycle plots...")
        
        if self.imaging_df is None:
            return
        
        metrics = ['%>= Q30', 'Error Rate']
        available_metrics = [m for m in metrics if m in self.imaging_df.columns]
        
        if not available_metrics:
            return
        
        for metric in available_metrics:
            try:
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    continue
                
                fig, ax = plt.subplots(figsize=(10, 6))
                
                lanes = sorted(self.imaging_df['Lane'].unique())
                plotted_any = False
                
                for lane in lanes:
                    lane_data = self.imaging_df[self.imaging_df['Lane'] == lane]
                    cycle_mean = lane_data.groupby('Cycle')[metric].mean()
                    
                    if not self._has_valid_data(cycle_mean, f"{metric}_lane{lane}"):
                        continue
                    
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'Lane {lane}', linewidth=2.5, alpha=0.8, color=color, marker='o', markersize=4)
                    plotted_any = True
                
                if not plotted_any:
                    plt.close(fig)
                    continue
                
                ax.set_xlabel('Cycle', fontsize=11, fontweight='bold')
                ax.set_ylabel(metric, fontsize=11, fontweight='bold')
                ax.set_title(f'{metric} by Cycle', fontsize=13, fontweight='bold')
                ax.legend(loc='best', framealpha=0.95, ncol=4)
                ax.grid(True, alpha=0.3, linestyle='--')
                ax.set_axisbelow(True)
                fig.tight_layout()
                
                description = f"Per-cycle {metric} shows sequencing quality across all lanes; declining quality indicates run issues."
                self.plot_images.append((f'{metric} by Cycle', description, self._save_fig_to_image(fig)))
                print(f"  ✓ Created {metric} plot")
            except Exception as e:
                print(f"[WARN] Could not create {metric} plot: {e}")
    
    def plot_intensity_distribution(self):
        """Plot intensity metrics by cycle"""
        print("[INFO] Generating intensity plots...")
        
        if self.imaging_df is None:
            return
        
        intensity_cols = [col for col in self.imaging_df.columns 
                         if any(x in col for x in ['Intensity', 'Corrected', 'Called'])]
        
        if not intensity_cols:
            return
        
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            cycle_data = self.imaging_df.groupby('Cycle')[intensity_cols].mean()
            
            plotted_any = False
            for col in intensity_cols:
                if not self._has_valid_data(cycle_data[col], col):
                    continue
                
                label = col
                ax.plot(cycle_data.index, cycle_data[col], 
                       label=label, linewidth=2.5, alpha=0.8, marker='o', markersize=3)
                plotted_any = True
            
            if not plotted_any:
                plt.close(fig)
                return
            
            ax.set_xlabel('Cycle', fontsize=11, fontweight='bold')
            ax.set_ylabel('Intensity (AU)', fontsize=11, fontweight='bold')
            ax.set_title('Intensity Metrics by Cycle', fontsize=13, fontweight='bold')
            ax.legend(loc='best', framealpha=0.95)
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_axisbelow(True)
            fig.tight_layout()
            
            description = "Average signal intensity across all reads and lanes shows overall photometric quality; declining intensity indicates optical degradation."
            self.plot_images.append(('Intensity Metrics', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created intensity plot")
        except Exception as e:
            print(f"[WARN] Could not create intensity plot: {e}")
    
    def plot_base_composition(self):
        """Plot base composition by cycle"""
        print("[INFO] Generating base composition plots...")
        
        if self.imaging_df is None:
            return
        
        base_cols = [col for col in self.imaging_df.columns if '% Base' in col or '%Base' in col]
        
        if not base_cols:
            return
        
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            cycle_data = self.imaging_df.groupby('Cycle')[base_cols].mean()
            
            has_data = False
            for col in base_cols:
                if not self._has_valid_data(cycle_data[col], col):
                    continue
                
                base = col.split('/')[-1] if '/' in col else col
                color = CHANNEL_COLORS.get(base, None)
                ax.plot(cycle_data.index, cycle_data[col], 
                       label=f'Base {base}', linewidth=2.5, alpha=0.8, color=color, marker='o', markersize=4)
                has_data = True
            
            if not has_data:
                plt.close(fig)
                return
            
            ax.set_xlabel('Cycle', fontsize=11, fontweight='bold')
            ax.set_ylabel('Percentage (%)', fontsize=11, fontweight='bold')
            ax.set_title('Base Composition by Cycle', fontsize=13, fontweight='bold')
            ax.legend(loc='best', framealpha=0.95)
            ax.set_ylim([0, 100])
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_axisbelow(True)
            fig.tight_layout()
            
            description = "Per-cycle distribution of the four DNA bases (A, C, G, T) detects nucleotide imbalance or contamination."
            self.plot_images.append(('Base Composition', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created base composition plot")
        except Exception as e:
            print(f"[WARN] Could not create base composition plot: {e}")
    
    def plot_fwhm_metrics(self):
        """Plot FWHM (focus) metrics"""
        print("[INFO] Generating FWHM plots...")
        
        if self.imaging_df is None:
            return
        
        fwhm_cols = [col for col in self.imaging_df.columns if 'Fwhm' in col or 'FWHM' in col]
        
        if not fwhm_cols:
            return
        
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            has_data = False
            for metric in fwhm_cols:
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    continue
                
                cycle_mean = self.imaging_df.groupby('Cycle')[metric].mean()
                ax.plot(cycle_mean.index, cycle_mean.values, 
                       label=metric, linewidth=2.5, alpha=0.8, marker='s', markersize=4)
                has_data = True
            
            if not has_data:
                plt.close(fig)
                return
            
            ax.set_xlabel('Cycle', fontsize=11, fontweight='bold')
            ax.set_ylabel('FWHM', fontsize=11, fontweight='bold')
            ax.set_title('FWHM (Focus Metric) by Cycle', fontsize=13, fontweight='bold')
            ax.legend(loc='best', framealpha=0.95)
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_axisbelow(True)
            fig.tight_layout()
            
            description = "Full Width at Half Maximum measures optical focus quality; higher values indicate worse focus and potential quality issues."
            self.plot_images.append(('FWHM Metrics', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created FWHM plot")
        except Exception as e:
            print(f"[WARN] Could not create FWHM plot: {e}")
    
    def plot_phasing_prephasing(self):
        """Plot phasing and prephasing metrics"""
        print("[INFO] Generating phasing plots...")
        
        if self.imaging_df is None:
            return
        
        phasing_cols = [col for col in self.imaging_df.columns 
                       if 'Phasing' in col or 'Prephasing' in col]
        
        if not phasing_cols:
            return
        
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            has_data = False
            for metric in phasing_cols:
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    continue
                
                cycle_mean = self.imaging_df.groupby('Cycle')[metric].mean()
                ax.plot(cycle_mean.index, cycle_mean.values, 
                       label=metric, linewidth=2.5, alpha=0.8, marker='o', markersize=4)
                has_data = True
            
            if not has_data:
                plt.close(fig)
                return
            
            ax.set_xlabel('Cycle', fontsize=11, fontweight='bold')
            ax.set_ylabel('Percentage (%)', fontsize=11, fontweight='bold')
            ax.set_title('Phasing/Prephasing by Cycle', fontsize=13, fontweight='bold')
            ax.legend(loc='best', framealpha=0.95)
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_axisbelow(True)
            fig.tight_layout()
            
            description = "Phasing and prephasing rates measure template strand synchronization issues that accumulate over sequencing cycles."
            self.plot_images.append(('Phasing Metrics', description, self._save_fig_to_image(fig)))
            print(f"  ✓ Created phasing plot")
        except Exception as e:
            print(f"[WARN] Could not create phasing plot: {e}")
    
    def plot_q40_metrics(self):
        """Plot % clusters with Q >= Q40 if available"""
        print("[INFO] Generating Q40 plots...")
        
        if self.imaging_df is None:
            return
        
        q40_cols = [col for col in self.imaging_df.columns if 'Q40' in col or 'Q>=40' in col]
        
        if not q40_cols:
            print("[WARN] No Q40 metrics available")
            return
        
        for metric in q40_cols:
            try:
                if not self._has_valid_data(self.imaging_df[metric], metric):
                    continue
                
                fig, ax = plt.subplots(figsize=(10, 6))
                
                lanes = sorted(self.imaging_df['Lane'].unique())
                plotted_any = False
                
                for lane in lanes:
                    lane_data = self.imaging_df[self.imaging_df['Lane'] == lane]
                    cycle_mean = lane_data.groupby('Cycle')[metric].mean()
                    
                    if not self._has_valid_data(cycle_mean, f"{metric}_lane{lane}"):
                        continue
                    
                    color = LANE_COLORS.get(lane, f'C{lane-1}')
                    ax.plot(cycle_mean.index, cycle_mean.values, 
                           label=f'Lane {lane}', linewidth=2.5, alpha=0.8, color=color, marker='o', markersize=4)
                    plotted_any = True
                
                if not plotted_any:
                    plt.close(fig)
                    continue
                
                ax.set_xlabel('Cycle', fontsize=11, fontweight='bold')
                ax.set_ylabel(metric, fontsize=11, fontweight='bold')
                ax.set_title(f'{metric} by Cycle', fontsize=13, fontweight='bold')
                ax.legend(loc='best', framealpha=0.95, ncol=4)
                ax.grid(True, alpha=0.3, linestyle='--')
                ax.set_axisbelow(True)
                fig.tight_layout()
                
                description = f"Percentage of clusters achieving Q≥40 quality score by cycle indicates high-confidence base calling performance."
                self.plot_images.append((f'Q40 Metrics', description, self._save_fig_to_image(fig)))
                print(f"  ✓ Created Q40 metrics plot")
            except Exception as e:
                print(f"[WARN] Could not create Q40 plot: {e}")
    
    def generate_pdf_report(self, output_pdf):
        """Generate comprehensive PDF report with all plots"""
        print("[INFO] Generating PDF report...")
        
        doc = SimpleDocTemplate(output_pdf, pagesize=landscape(letter),
                               rightMargin=0.5*inch, leftMargin=0.5*inch,
                               topMargin=0.5*inch, bottomMargin=0.5*inch)
        
        story = []
        styles = getSampleStyleSheet()
        
        # Custom styles
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#0072B2'),
            spaceAfter=30,
            alignment=TA_CENTER,
            fontName='Helvetica-Bold'
        )
        
        heading_style = ParagraphStyle(
            'CustomHeading',
            parent=styles['Heading2'],
            fontSize=14,
            textColor=colors.HexColor('#0072B2'),
            spaceAfter=6,
            spaceBefore=6,
            fontName='Helvetica-Bold'
        )
        
        description_style = ParagraphStyle(
            'Description',
            parent=styles['Normal'],
            fontSize=10,
            textColor=colors.HexColor('#333333'),
            spaceAfter=12,
            alignment=TA_LEFT,
            fontName='Helvetica'
        )
        
        # Title page
        story.append(Spacer(1, 1.5*inch))
        story.append(Paragraph("Illumina Sequencing Run", title_style))
        story.append(Paragraph("Quality Control Report", title_style))
        story.append(Spacer(1, 0.7*inch))
        
        # Run information overview table
        run_info_data = [
            ['Run Folder:', self.run_folder.name],
            ['Analysis Date:', datetime.now().strftime('%Y-%m-%d %H:%M:%S')],
        ]
        
        try:
            run_info = self.run_metrics.run_info()
            flowcell_barcode = run_info.flowcell().barcode()
            total_cycles = run_info.total_cycles()
            num_reads = run_info.reads().size()
            
            run_info_data.extend([
                ['Flowcell Barcode:', flowcell_barcode],
                ['Total Cycles:', str(total_cycles)],
                ['Number of Reads:', str(num_reads)],
            ])
            
            for i in range(num_reads):
                read = run_info.reads()[i]
                read_info = f"{read.total_cycles()} cycles"
                if read.is_index():
                    read_info += " (Index)"
                run_info_data.append([f'Read {i+1}:', read_info])
        except Exception as e:
            print(f"[WARN] Could not load full run info: {e}")
        
        try:
            tile_df = self._get_tile_metrics()
            if tile_df is not None and len(tile_df) > 0:
                total_clusters = tile_df['Cluster Count'].sum()
                clusters_pf = tile_df['Cluster Count PF'].sum()
                
                run_info_data.append(['Total Clusters:', f'{total_clusters:.0f}'])
                run_info_data.append(['Clusters Passing Filter:', f'{clusters_pf:.0f}'])
                
                if total_clusters > 0:
                    pf_pct = (clusters_pf / total_clusters * 100)
                    run_info_data.append(['% Clusters Passing Filter:', f'{pf_pct:.1f}%'])
        except Exception as e:
            print(f"[WARN] Could not add cluster info: {e}")
        
        info_table = Table(run_info_data, colWidths=[2.5*inch, 5*inch])
        info_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#e8f4f8')),
            ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('ALIGN', (0, 0), (0, -1), 'RIGHT'),
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTNAME', (1, 0), (1, -1), 'Courier'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ]))
        
        story.append(info_table)
        story.append(PageBreak())
        
        # Add all plots - one per page
        for plot_title, plot_description, plot_data in self.plot_images:
            story.append(Paragraph(plot_title, heading_style))
            story.append(Paragraph(plot_description, description_style))
            img = Image(BytesIO(plot_data), width=9*inch, height=5.4*inch, kind='proportional')
            story.append(img)
            story.append(PageBreak())
        
        # Build PDF
        doc.build(story)
        print(f"  ✓ Created PDF report: {output_pdf}")
    
    def generate_all_plots(self, output_pdf):
        """Generate all available QC plots and create PDF"""
        print("\n" + "=" * 80)
        print("GENERATING QC PLOTS")
        print("=" * 80 + "\n")
        
        self.plot_summary_page()
        self.plot_cluster_distribution()
        self.plot_overall_cluster_pie()
    #    self.plot_lane_balance()
        self.plot_percent_pf_by_lane()
        self.plot_cluster_density_heatmap()
        self.plot_qscore_by_cycle()
        self.plot_intensity_distribution()
        self.plot_phasing_prephasing()
        self.plot_fwhm_metrics()
        self.plot_base_composition()
        self.plot_q40_metrics()
        self.generate_pdf_report(output_pdf)
        
        print("\n" + "=" * 80)
        print("QC PLOT GENERATION COMPLETE")
        print("=" * 80)
        print(f"\nPDF report saved to: {output_pdf}")
        print(f"Total plots generated: {len(self.plot_images)}")
        print()


def main():
    if len(sys.argv) != 3:
        print("Usage: python interop_more.py <run_folder> <output_pdf>")
        print("\nExample:")
        print("  python interop_more.py /data/210101_A00000_0001_AHXXXX ./QC_Report.pdf")
        sys.exit(1)
    
    run_folder = sys.argv[1]
    output_pdf = sys.argv[2]
    
    if not os.path.exists(run_folder):
        print(f"[ERROR] Run folder does not exist: {run_folder}")
        sys.exit(1)
    
    if not os.path.exists(os.path.join(run_folder, "InterOp")):
        print(f"[ERROR] InterOp folder not found in: {run_folder}")
        sys.exit(1)
    
    output_dir = os.path.dirname(output_pdf)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    try:
        plotter = InterOpQCPlotter(run_folder)
        plotter.generate_all_plots(output_pdf)
    except Exception as e:
        print(f"\n[ERROR] Failed to generate plots: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()