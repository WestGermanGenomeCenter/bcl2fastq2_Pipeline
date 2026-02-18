#!/usr/bin/env python3
"""
Sequencing Run Report Generator

This script parses Illumina sequencing run output files and generates
a comprehensive PDF report with summary statistics, tables, and plots.

Usage:
    python sequencing_report_generator.py --input <input_folder> --output <output_pdf>
"""

import os
import sys
import argparse
import warnings
from pathlib import Path
import xml.etree.ElementTree as ET

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

warnings.filterwarnings('ignore')

# Set consistent style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("Set2")

# Define consistent colors for metrics
COLORS = {
    'primary': '#2E86AB',
    'secondary': '#A23B72', 
    'success': '#06A77D',
    'warning': '#F18F01',
    'danger': '#C73E1D',
    'q30': '#06A77D',
    'yield': '#2E86AB',
    'adapter': '#F18F01',
}


def format_axis_ticks(ax, max_ticks=10, apply_x=True, apply_y=True):
    """Apply consistent tick formatting to axis - max 10 ticks on each axis.
    
    Args:
        ax: matplotlib axis object
        max_ticks: maximum number of ticks (default 10)
        apply_x: whether to apply limit to x-axis (default True)
        apply_y: whether to apply limit to y-axis (default True)
    """
    if apply_x:
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=max_ticks))
    if apply_y:
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=max_ticks))


class SequencingRunParser:
    """Parser for Illumina sequencing run output files."""
    
    def __init__(self, input_folder):
        self.input_folder = Path(input_folder)
        self.data = {}
        
    def parse_all(self):
        """Parse all available files in the input folder."""
        print("Parsing sequencing run files...")
        
        # Parse each file type
        self.parse_run_info()
        self.parse_sample_sheet()
        self.parse_demultiplex_stats()
        self.parse_quality_metrics()
        self.parse_adapter_metrics()
        self.parse_top_unknown_barcodes()
        self.parse_index_hopping()
        self.parse_tile_metrics()
        
        print("Parsing complete!")
        return self.data
    
    def parse_run_info(self):
        """Parse RunInfo.xml file."""
        file_path = self.input_folder / 'RunInfo.xml'
        if not file_path.exists():
            print(f"Warning: {file_path} not found")
            return
        
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        run = root.find('Run')
        self.data['run_info'] = {
            'run_id': run.get('Id'),
            'run_number': run.get('Number'),
            'flowcell': run.find('Flowcell').text,
            'instrument': run.find('Instrument').text,
            'date': run.find('Date').text,
        }
        
        # Parse reads info
        reads = []
        for read in run.findall('.//Read'):
            reads.append({
                'number': read.get('Number'),
                'cycles': read.get('NumCycles'),
                'is_index': read.get('IsIndexedRead') == 'Y'
            })
        self.data['run_info']['reads'] = reads
        
        # Parse flowcell layout
        layout = run.find('.//FlowcellLayout')
        self.data['run_info']['layout'] = {
            'lanes': int(layout.get('LaneCount')),
            'surfaces': int(layout.get('SurfaceCount')),
            'swaths': int(layout.get('SwathCount')),
            'tiles': int(layout.get('TileCount'))
        }
    
    def parse_sample_sheet(self):
        """Parse SampleSheet.csv file."""
        file_path = self.input_folder / 'SampleSheet.csv'
        if not file_path.exists():
            print(f"Warning: {file_path} not found")
            return
        
        # Read the sample sheet sections
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Parse header section
        header_info = {}
        for line in lines:
            if line.startswith('['):
                break
            parts = line.strip().split(',')
            if len(parts) >= 2 and parts[0]:
                header_info[parts[0]] = parts[1]
        
        self.data['sample_sheet'] = header_info
    
    def parse_demultiplex_stats(self):
        """Parse Demultiplex_Stats.csv file."""
        file_path = self.input_folder / 'Demultiplex_Stats.csv'
        if not file_path.exists():
            print(f"Warning: {file_path} not found")
            return
        
        df = pd.read_csv(file_path)
        self.data['demux_stats'] = df
    
    def parse_quality_metrics(self):
        """Parse Quality_Metrics.csv file."""
        file_path = self.input_folder / 'Quality_Metrics.csv'
        if not file_path.exists():
            print(f"Warning: {file_path} not found")
            return
        
        df = pd.read_csv(file_path)
        self.data['quality_metrics'] = df
    
    def parse_adapter_metrics(self):
        """Parse Adapter_Metrics.csv file."""
        file_path = self.input_folder / 'Adapter_Metrics.csv'
        if not file_path.exists():
            print(f"Warning: {file_path} not found")
            return
        
        df = pd.read_csv(file_path)
        self.data['adapter_metrics'] = df
    
    def parse_top_unknown_barcodes(self):
        """Parse Top_Unknown_Barcodes.csv file."""
        file_path = self.input_folder / 'Top_Unknown_Barcodes.csv'
        if not file_path.exists():
            print(f"Warning: {file_path} not found")
            return
        
        df = pd.read_csv(file_path)
        self.data['unknown_barcodes'] = df
    
    def parse_index_hopping(self):
        """Parse Index_Hopping_Counts.csv file."""
        file_path = self.input_folder / 'Index_Hopping_Counts.csv'
        if not file_path.exists():
            print(f"Warning: {file_path} not found")
            return
        
        df = pd.read_csv(file_path)
        self.data['index_hopping'] = df
    
    def parse_tile_metrics(self):
        """Parse Quality_Tile_Metrics.csv file."""
        file_path = self.input_folder / 'Quality_Tile_Metrics.csv'
        if not file_path.exists():
            print(f"Warning: {file_path} not found")
            return
        
        df = pd.read_csv(file_path)
        self.data['tile_metrics'] = df


class ReportGenerator:
    """Generate PDF report from parsed sequencing data."""
    
    def __init__(self, data):
        self.data = data
        
    def generate(self, output_path):
        """Generate the complete PDF report."""
        print(f"Generating PDF report: {output_path}")
        sys.stdout.flush()
        
        with PdfPages(output_path) as pdf:
            # Title page
            print("  - Creating title page...", end='', flush=True)
            self.add_title_page(pdf)
            print(" Done")
            
            # Run summary
            print("  - Creating run summary...", end='', flush=True)
            self.add_run_summary(pdf)
            print(" Done")
            
            # Demultiplexing statistics
            print("  - Creating demultiplexing statistics...", end='', flush=True)
            self.add_demux_summary(pdf)
            print(" Done")
            
            # Quality metrics
            print("  - Creating quality metrics...", end='', flush=True)
            self.add_quality_metrics(pdf)
            print(" Done")
            
            # Adapter metrics
            print("  - Creating adapter metrics...", end='', flush=True)
            self.add_adapter_metrics(pdf)
            print(" Done")
            
            # Unknown barcodes
            print("  - Creating unknown barcodes analysis...", end='', flush=True)
            self.add_unknown_barcodes(pdf)
            print(" Done")
            
            # Tile metrics
            print("  - Creating tile metrics...", end='', flush=True)
            self.add_tile_metrics(pdf)
            print(" Done")
            
            # Sample details table
            print("  - Creating sample details tables...", end='', flush=True)
            self.add_sample_details_table(pdf)
            print(" Done")
        
        print(f"Report generated successfully: {output_path}")
        sys.stdout.flush()
    
    def add_title_page(self, pdf):
        """Add title page to the report."""
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.7, 'Sequencing Run Report', 
                ha='center', fontsize=28, weight='bold')
        
        if 'run_info' in self.data:
            info = self.data['run_info']
            y_pos = 0.55
            
            details = [
                f"Run ID: {info.get('run_id', 'N/A')}",
                f"Flowcell: {info.get('flowcell', 'N/A')}",
                f"Instrument: {info.get('instrument', 'N/A')}",
                f"Date: {info.get('date', 'N/A')[:10]}",
            ]
            
            for detail in details:
                fig.text(0.5, y_pos, detail, ha='center', fontsize=14)
                y_pos -= 0.05
        
        fig.text(0.5, 0.15, f'Generated: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")}',
                ha='center', fontsize=10, style='italic', color='gray')
        
        plt.axis('off')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        plt.close('all')
    
    def add_run_summary(self, pdf):
        """Add run summary page."""
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Run Summary', fontsize=16, weight='bold', y=0.98)
        
        # Create summary statistics
        summary_text = []
        
        if 'run_info' in self.data:
            info = self.data['run_info']
            summary_text.append("Run Configuration")
            summary_text.append("-" * 50)
            
            if 'reads' in info:
                data_reads = [r for r in info['reads'] if not r['is_index']]
                index_reads = [r for r in info['reads'] if r['is_index']]
                
                # Determine if single or paired-end
                read_type = "Paired-End" if len(data_reads) > 1 else "Single-Read"
                summary_text.append(f"  Read Type: {read_type}")
                
                for read in info['reads']:
                    read_type_label = "Index" if read['is_index'] else "Data"
                    summary_text.append(f"  Read {read['number']}: {read['cycles']} cycles ({read_type_label})")
            
            if 'layout' in info:
                layout = info['layout']
                summary_text.append(f"  Lanes: {layout['lanes']}")
                summary_text.append(f"  Tiles per lane: {layout['surfaces'] * layout['swaths'] * layout['tiles']}")
            summary_text.append("")
        
        if 'demux_stats' in self.data:
            df = self.data['demux_stats']
            # Exclude Undetermined for cleaner stats
            df_determined = df[df['SampleID'] != 'Undetermined']
            df_undetermined = df[df['SampleID'] == 'Undetermined']
            
            total_reads = df['# Reads'].sum()
            determined_reads = df_determined['# Reads'].sum()
            undetermined_reads = df_undetermined['# Reads'].sum() if len(df_undetermined) > 0 else 0
            num_samples = len(df_determined)
            
            summary_text.append("Demultiplexing Summary")
            summary_text.append("-" * 50)
            summary_text.append(f"  Total Samples: {num_samples:,}")
            summary_text.append(f"  Total Reads: {total_reads:,}")
            summary_text.append(f"  Assigned Reads: {determined_reads:,} ({determined_reads/total_reads:.2%})")
            summary_text.append(f"  Undetermined Reads: {undetermined_reads:,} ({undetermined_reads/total_reads:.2%})")
            summary_text.append(f"  Avg Reads/Sample: {determined_reads/num_samples:,.0f}")
            summary_text.append(f"  Perfect Index Rate: {df_determined['% Perfect Index Reads'].mean():.2%}")
            summary_text.append("")
        
        if 'quality_metrics' in self.data:
            df = self.data['quality_metrics']
            # Exclude Undetermined
            df_determined = df[df['SampleID'] != 'Undetermined']
            
            total_yield = df_determined['Yield'].sum()
            avg_q30 = df_determined['% Q30'].mean()
            avg_quality = df_determined['Mean Quality Score (PF)'].mean()
            
            # Check if single or paired-end
            data_reads = [r for r in self.data['run_info']['reads'] if not r['is_index']] if 'run_info' in self.data else []
            read_note = " (single-read)" if len(data_reads) == 1 else " (both reads)" if len(data_reads) == 2 else ""
            
            summary_text.append("Quality Summary")
            summary_text.append("-" * 50)
            summary_text.append(f"  Total Yield{read_note}: {total_yield/1e9:.2f} Gb")
            summary_text.append(f"  Average Q30: {avg_q30:.2%}")
            summary_text.append(f"  Average Quality Score: {avg_quality:.2f}")
            summary_text.append("")
        
        if 'adapter_metrics' in self.data:
            df = self.data['adapter_metrics']
            # Exclude Undetermined
            df_determined = df[df['Sample_ID'] != 'Undetermined']
            avg_adapter = df_determined['% Adapter Bases'].mean()
            
            summary_text.append("Adapter Summary")
            summary_text.append("-" * 50)
            summary_text.append(f"  Average Adapter Content: {avg_adapter:.3%}")
            summary_text.append("")
        
        # Display summary text
        ax = fig.add_subplot(111)
        ax.text(0.1, 0.95, '\n'.join(summary_text), 
               transform=ax.transAxes,
               fontsize=11,
               verticalalignment='top',
               fontfamily='monospace')
        ax.axis('off')
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        plt.close('all')
    
    def add_demux_summary(self, pdf):
        """Add demultiplexing summary with plots."""
        if 'demux_stats' not in self.data:
            return
        
        df = self.data['demux_stats']
        # Exclude Undetermined from plots
        df = df[df['SampleID'] != 'Undetermined'].copy()
        
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Demultiplexing Statistics (Undetermined Excluded)', fontsize=16, weight='bold', y=0.98)
        
        # Sort by reads
        df_sorted = df.sort_values('# Reads', ascending=False)
        
        # Plot 1: Reads per sample (top 20)
        ax1 = plt.subplot(2, 2, 1)
        top_samples = df_sorted.head(20)
        ax1.barh(range(len(top_samples)), top_samples['# Reads'], color=COLORS['primary'])
        ax1.set_yticks(range(len(top_samples)))
        ax1.set_yticklabels(top_samples['SampleID'], fontsize=8)
        ax1.set_xlabel('Number of Reads', fontsize=10)
        ax1.set_title('Reads per Sample (Top 20)', fontsize=11, weight='bold')
        ax1.invert_yaxis()
        format_axis_ticks(ax1, apply_y=False)
        ax1.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
        ax1.text(0.02, 0.02, 'Distribution of sequencing depth across top samples', 
                transform=ax1.transAxes, fontsize=7, style='italic', color='gray')
        
        # Plot 2: Distribution of reads
        ax2 = plt.subplot(2, 2, 2)
        ax2.hist(df['# Reads'], bins=30, color=COLORS['secondary'], alpha=0.7, edgecolor='black')
        ax2.set_xlabel('Number of Reads', fontsize=10)
        ax2.set_ylabel('Number of Samples', fontsize=10)
        ax2.set_title('Distribution of Reads Across Samples', fontsize=11, weight='bold')
        format_axis_ticks(ax2)
        ax2.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
        ax2.tick_params(axis='x', labelrotation=45)
        mean_reads = df['# Reads'].mean()
        ax2.axvline(mean_reads, color='red', linestyle='--', linewidth=2, label=f'Mean')
        ax2.legend(fontsize=8)
        ax2.text(0.02, 0.98, 'Histogram showing read count variability', 
                transform=ax2.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 3: Perfect index rate
        ax3 = plt.subplot(2, 2, 3)
        perfect_rate = df['% Perfect Index Reads'].values * 100
        ax3.hist(perfect_rate, bins=20, color=COLORS['success'], alpha=0.7, edgecolor='black')
        ax3.set_xlabel('Perfect Index Rate (%)', fontsize=10)
        ax3.set_ylabel('Number of Samples', fontsize=10)
        ax3.set_title('Perfect Index Read Rate Distribution', fontsize=11, weight='bold')
        ax3.set_xlim(0, 100)
        format_axis_ticks(ax3)
        ax3.axvline(perfect_rate.mean(), color='red', linestyle='--', label=f'Mean: {perfect_rate.mean():.1f}%')
        ax3.axvline(90, color='orange', linestyle=':', linewidth=2, label='Good (≥90%)')
        ax3.legend(fontsize=8)
        ax3.text(0.02, 0.98, 'Percentage of reads with perfect barcode matches', 
                transform=ax3.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 4: Summary statistics box
        ax4 = plt.subplot(2, 2, 4)
        stats_text = f"""
Demultiplexing Statistics

Total Samples: {len(df):,}
Total Reads: {df['# Reads'].sum():,}

Reads per Sample:
  Min: {df['# Reads'].min():,}
  Max: {df['# Reads'].max():,}
  Mean: {df['# Reads'].mean():,.0f}
  Median: {df['# Reads'].median():,.0f}

Perfect Index Rate:
  Mean: {df['% Perfect Index Reads'].mean():.2%}
  Min: {df['% Perfect Index Reads'].min():.2%}
  Max: {df['% Perfect Index Reads'].max():.2%}
        """
        ax4.text(0.1, 0.9, stats_text.strip(), transform=ax4.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace')
        ax4.axis('off')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        plt.close('all')
    
    def add_quality_metrics(self, pdf):
        """Add quality metrics with plots."""
        if 'quality_metrics' not in self.data:
            return
        
        df = self.data['quality_metrics']
        # Exclude Undetermined
        df = df[df['SampleID'] != 'Undetermined'].copy()
        
        
        # Aggregate by sample (sum yield, average quality scores)
        sample_metrics = df.groupby('SampleID').agg({
            'Yield': 'sum',
            'YieldQ30': 'sum',
            '% Q30': 'mean',
            'Mean Quality Score (PF)': 'mean'
        }).reset_index()
        
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Quality Metrics (Undetermined Excluded)', fontsize=16, weight='bold', y=0.98)
        
        # Plot 1: Yield per sample (top 20)
        ax1 = plt.subplot(2, 2, 1)
        top_yield = sample_metrics.nlargest(20, 'Yield')
        bars = ax1.barh(range(len(top_yield)), top_yield['Yield']/1e9, color=COLORS['yield'])
        ax1.set_yticks(range(len(top_yield)))
        ax1.set_yticklabels(top_yield['SampleID'], fontsize=8)
        ax1.set_xlabel('Yield (Gb)', fontsize=10)
        ax1.set_title('Yield per Sample (Top 20)', fontsize=11, weight='bold')
        ax1.invert_yaxis()
        format_axis_ticks(ax1, apply_y=False)
        ax1.text(0.02, 0.02, 'Total sequencing output per sample', 
                transform=ax1.transAxes, fontsize=7, style='italic', color='gray')
        
        # Plot 2: Q30 distribution
        ax2 = plt.subplot(2, 2, 2)
        q30_values = sample_metrics['% Q30'].values * 100
        ax2.hist(q30_values, bins=20, color=COLORS['q30'], alpha=0.7, edgecolor='black')
        ax2.set_xlabel('% Q30', fontsize=10)
        ax2.set_ylabel('Number of Samples', fontsize=10)
        ax2.set_title('Q30 Score Distribution', fontsize=11, weight='bold')
        format_axis_ticks(ax2)
        ax2.axvline(q30_values.mean(), color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {q30_values.mean():.1f}%')
        ax2.axvline(30, color='orange', linestyle=':', linewidth=2, label='Threshold (≥30%)')
        ax2.legend(fontsize=8)
        ax2.text(0.02, 0.98, 'Percentage of bases with quality score ≥30 (99.9%% accuracy)', 
                transform=ax2.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 3: Mean quality score distribution
        ax3 = plt.subplot(2, 2, 3)
        qual_scores = sample_metrics['Mean Quality Score (PF)'].values
        ax3.hist(qual_scores, bins=20, color=COLORS['secondary'], alpha=0.7, edgecolor='black')
        ax3.set_xlabel('Mean Quality Score', fontsize=10)
        ax3.set_ylabel('Number of Samples', fontsize=10)
        ax3.set_title('Mean Quality Score Distribution', fontsize=11, weight='bold')
        format_axis_ticks(ax3)
        ax3.axvline(qual_scores.mean(), color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {qual_scores.mean():.1f}')
        ax3.axvline(30, color='orange', linestyle=':', linewidth=2, label='Good (≥Q30)')
        ax3.legend(fontsize=8)
        ax3.text(0.02, 0.98, 'Average Phred quality score across all bases', 
                transform=ax3.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 4: Q30 vs Mean Quality Score scatter
        ax4 = plt.subplot(2, 2, 4)
        scatter = ax4.scatter(sample_metrics['% Q30']*100, 
                            sample_metrics['Mean Quality Score (PF)'],
                            c=sample_metrics['Yield']/1e9, 
                            cmap='viridis', alpha=0.6, s=50)
        ax4.set_xlabel('% Q30', fontsize=10)
        ax4.set_ylabel('Mean Quality Score', fontsize=10)
        ax4.set_title('Q30 vs Mean Quality Score', fontsize=11, weight='bold')
        format_axis_ticks(ax4)
        cbar = plt.colorbar(scatter, ax=ax4)
        cbar.set_label('Yield (Gb)', fontsize=8)
        ax4.text(0.02, 0.98, 'Relationship between quality metrics (color = yield)', 
                transform=ax4.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def add_adapter_metrics(self, pdf):
        """Add adapter metrics with plots."""
        if 'adapter_metrics' not in self.data:
            return
        
        df = self.data['adapter_metrics']
        # Exclude Undetermined
        df = df[df['Sample_ID'] != 'Undetermined'].copy()
        
        # Aggregate by sample
        sample_adapter = df.groupby('Sample_ID').agg({
            'AdapterBases': 'sum',
            'SampleBases': 'sum',
            '% Adapter Bases': 'mean'
        }).reset_index()
        
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Adapter Contamination Metrics (Undetermined Excluded)', fontsize=16, weight='bold', y=0.98)
        
        # Plot 1: Adapter percentage distribution
        ax1 = plt.subplot(2, 2, 1)
        adapter_pct = sample_adapter['% Adapter Bases'].values * 100
        ax1.hist(adapter_pct, bins=20, color=COLORS['adapter'], alpha=0.7, edgecolor='black')
        ax1.set_xlabel('% Adapter Bases', fontsize=10)
        ax1.set_ylabel('Number of Samples', fontsize=10)
        ax1.set_title('Adapter Content Distribution', fontsize=11, weight='bold')
        format_axis_ticks(ax1)
        ax1.axvline(adapter_pct.mean(), color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {adapter_pct.mean():.3f}%')
        ax1.axvline(1.0, color='orange', linestyle=':', linewidth=2, label='Warning (≥1%)')
        ax1.legend(fontsize=8)
        ax1.text(0.02, 0.98, 'Percentage of adapter contamination in sequencing reads', 
                transform=ax1.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 2: Top samples by adapter content
        ax2 = plt.subplot(2, 2, 2)
        top_adapter = sample_adapter.nlargest(20, '% Adapter Bases')
        ax2.barh(range(len(top_adapter)), top_adapter['% Adapter Bases']*100, 
                color=COLORS['warning'])
        ax2.set_yticks(range(len(top_adapter)))
        ax2.set_yticklabels(top_adapter['Sample_ID'], fontsize=8)
        ax2.set_xlabel('% Adapter Bases', fontsize=10)
        ax2.set_title('Highest Adapter Content (Top 20)', fontsize=11, weight='bold')
        ax2.invert_yaxis()
        format_axis_ticks(ax2, apply_y=False)  # Don't limit y-axis, keep all sample names
        
        # Plot 3: Adapter bases vs total bases
        ax3 = plt.subplot(2, 2, 3)
        ax3.scatter(sample_adapter['SampleBases']/1e6, 
                   sample_adapter['AdapterBases']/1e3,
                   alpha=0.6, s=50, color=COLORS['adapter'])
        ax3.set_xlabel('Sample Bases (Mb)', fontsize=10)
        ax3.set_ylabel('Adapter Bases (Kb)', fontsize=10)
        ax3.set_title('Adapter Bases vs Sample Bases', fontsize=11, weight='bold')
        format_axis_ticks(ax3)
        ax3.text(0.02, 0.98, 'Correlation between sequencing depth and adapter contamination', 
                transform=ax3.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 4: Summary statistics
        ax4 = plt.subplot(2, 2, 4)
        stats_text = f"""
Adapter Metrics Summary

Total Samples: {len(sample_adapter):,}

Adapter Content (%):
  Mean: {sample_adapter['% Adapter Bases'].mean():.4%}
  Median: {sample_adapter['% Adapter Bases'].median():.4%}
  Min: {sample_adapter['% Adapter Bases'].min():.4%}
  Max: {sample_adapter['% Adapter Bases'].max():.4%}

Total Adapter Bases: {sample_adapter['AdapterBases'].sum():,}
Total Sample Bases: {sample_adapter['SampleBases'].sum():,}

Samples > 0.01%: {(sample_adapter['% Adapter Bases'] > 0.0001).sum()}
Samples > 0.1%: {(sample_adapter['% Adapter Bases'] > 0.001).sum()}
        """
        ax4.text(0.1, 0.9, stats_text.strip(), transform=ax4.transAxes,
                fontsize=9, verticalalignment='top', fontfamily='monospace')
        ax4.axis('off')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def add_unknown_barcodes(self, pdf):
        """Add unknown barcodes analysis."""
        if 'unknown_barcodes' not in self.data:
            return
        
        df = self.data['unknown_barcodes']
        
        if len(df) == 0:
            return
        
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Top Unknown Barcodes', fontsize=16, weight='bold', y=0.98)
        
        # Take top 30 unknown barcodes
        top_n = min(30, len(df))
        top_unknown = df.head(top_n).copy()
        
        # Create comprehensive table with all columns
        ax = plt.subplot(111)
        ax.axis('tight')
        ax.axis('off')
        
        # Add summary at top
        total_unknown = df['# Reads'].sum()
        pct_of_all = df['% of All Reads'].sum()
        summary = f"Total Unknown Reads: {total_unknown:,} ({pct_of_all:.3%} of all reads)\n\n"
        ax.text(0.5, 0.98, summary, transform=ax.transAxes, fontsize=11, 
               ha='center', va='top', weight='bold')
        
        # Prepare table data
        table_data = []
        for idx, row in top_unknown.iterrows():
            table_data.append([
                f"{idx+1}",
                row['index'],
                row['index2'],
                f"{row['# Reads']:,}",
                f"{row['% of Unknown Barcodes']:.2%}",
                f"{row['% of All Reads']:.4%}"
            ])
        
        # Create table
        col_labels = ['Rank', 'Index1', 'Index2', '# Reads', '% of Unknown', '% of All Reads']
        table = ax.table(cellText=table_data,
                        colLabels=col_labels,
                        cellLoc='left',
                        loc='upper center',
                        bbox=[0.05, 0.05, 0.9, 0.87],
                        colWidths=[0.08, 0.18, 0.18, 0.18, 0.18, 0.18])
        
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.5)
        
        # Style header
        for i in range(len(col_labels)):
            cell = table[(0, i)]
            cell.set_facecolor(COLORS['danger'])
            cell.set_text_props(weight='bold', color='white')
        
        # Alternate row colors and add mini bars in the reads column
        max_reads = top_unknown['# Reads'].max()
        for i in range(1, len(table_data) + 1):
            for j in range(len(col_labels)):
                cell = table[(i, j)]
                if i % 2 == 0:
                    cell.set_facecolor('#f0f0f0')
                
                # Highlight top 5
                if i <= 5:
                    cell.set_facecolor('#ffebee')
        
        # Add explanation
        ax.text(0.5, 0.01, 'Reads that could not be assigned to any sample due to barcode mismatches', 
               transform=ax.transAxes, fontsize=8, ha='center', style='italic', color='gray')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        plt.close('all')
    
    def add_tile_metrics(self, pdf):
        """Add tile-level quality metrics."""
        if 'tile_metrics' not in self.data:
            return
        
        df = self.data['tile_metrics']
        # Exclude Undetermined
        df = df[df['SampleID'] != 'Undetermined'].copy()
        
        # Aggregate by tile
        tile_summary = df.groupby('Tile').agg({
            'Yield': 'sum',
            '% Q30': 'mean',
            'Mean Quality Score (PF)': 'mean'
        }).reset_index()
        
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Tile Quality Metrics (Undetermined Excluded)', fontsize=16, weight='bold', y=0.98)
        
        # Plot 1: Yield by tile
        ax1 = plt.subplot(2, 2, 1)
        tiles = tile_summary['Tile'].astype(str)
        ax1.bar(range(len(tiles)), tile_summary['Yield']/1e6, color=COLORS['yield'])
        # Only show subset of tile labels for readability
        num_ticks = min(10, len(tiles))
        tick_positions = np.linspace(0, len(tiles)-1, num_ticks, dtype=int)
        ax1.set_xticks(tick_positions)
        ax1.set_xticklabels([tiles.iloc[i] for i in tick_positions], rotation=90, fontsize=7)
        ax1.set_ylabel('Yield (Mb)', fontsize=10)
        ax1.set_title('Yield per Tile', fontsize=11, weight='bold')
        ax1.yaxis.set_major_locator(plt.MaxNLocator(nbins=10))
        ax1.text(0.02, 0.98, 'Sequencing output across flowcell tiles (aggregated across samples)', 
                transform=ax1.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 2: Q30 by tile
        ax2 = plt.subplot(2, 2, 2)
        ax2.bar(range(len(tiles)), tile_summary['% Q30']*100, color=COLORS['q30'])
        # Only show subset of tile labels for readability
        ax2.set_xticks(tick_positions)
        ax2.set_xticklabels([tiles.iloc[i] for i in tick_positions], rotation=90, fontsize=7)
        ax2.set_ylabel('% Q30', fontsize=10)
        ax2.set_title('Q30 Score per Tile', fontsize=11, weight='bold')
        ax2.yaxis.set_major_locator(plt.MaxNLocator(nbins=10))
        ax2.axhline(tile_summary['% Q30'].mean()*100, color='red', linestyle='--',
                   label=f"Mean: {tile_summary['% Q30'].mean()*100:.1f}%")
        ax2.legend(fontsize=8)
        ax2.text(0.02, 0.98, 'Quality score consistency across tiles (high variance indicates issues)', 
                transform=ax2.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 3: Mean quality score by tile
        ax3 = plt.subplot(2, 2, 3)
        ax3.bar(range(len(tiles)), tile_summary['Mean Quality Score (PF)'], 
               color=COLORS['secondary'])
        # Only show subset of tile labels for readability
        ax3.set_xticks(tick_positions)
        ax3.set_xticklabels([tiles.iloc[i] for i in tick_positions], rotation=90, fontsize=7)
        ax3.set_ylabel('Mean Quality Score', fontsize=10)
        ax3.set_title('Mean Quality Score per Tile', fontsize=11, weight='bold')
        ax3.yaxis.set_major_locator(plt.MaxNLocator(nbins=10))
        ax3.axhline(tile_summary['Mean Quality Score (PF)'].mean(), color='red', 
                   linestyle='--', label=f"Mean: {tile_summary['Mean Quality Score (PF)'].mean():.1f}")
        ax3.legend(fontsize=8)
        ax3.text(0.02, 0.98, 'Mean quality by tile (outliers may indicate instrument problems)', 
                transform=ax3.transAxes, fontsize=7, va='top', style='italic', color='gray')
        
        # Plot 4: Tile statistics summary
        ax4 = plt.subplot(2, 2, 4)
        stats_text = f"""
Tile Metrics Summary

Total Tiles: {len(tile_summary)}

Yield (Mb):
  Min: {tile_summary['Yield'].min()/1e6:.1f}
  Max: {tile_summary['Yield'].max()/1e6:.1f}
  Mean: {tile_summary['Yield'].mean()/1e6:.1f}
  Std: {tile_summary['Yield'].std()/1e6:.1f}

% Q30:
  Min: {tile_summary['% Q30'].min():.2%}
  Max: {tile_summary['% Q30'].max():.2%}
  Mean: {tile_summary['% Q30'].mean():.2%}
  Std: {tile_summary['% Q30'].std():.2%}

Mean Quality Score:
  Min: {tile_summary['Mean Quality Score (PF)'].min():.2f}
  Max: {tile_summary['Mean Quality Score (PF)'].max():.2f}
  Mean: {tile_summary['Mean Quality Score (PF)'].mean():.2f}
  Std: {tile_summary['Mean Quality Score (PF)'].std():.2f}
        """
        ax4.text(0.1, 0.9, stats_text.strip(), transform=ax4.transAxes,
                fontsize=9, verticalalignment='top', fontfamily='monospace')
        ax4.axis('off')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def add_sample_details_table(self, pdf):
        """Add detailed sample table spanning multiple pages if needed."""
        if 'demux_stats' not in self.data or 'quality_metrics' not in self.data:
            return
        
        # Merge demux and quality data
        demux = self.data['demux_stats'].copy()
        quality = self.data['quality_metrics'].groupby('SampleID').agg({
            'Yield': 'sum',
            '% Q30': 'mean',
            'Mean Quality Score (PF)': 'mean'
        }).reset_index()
        
        # Split the Index column into Index1 and Index2
        if 'Index' in demux.columns:
            index_split = demux['Index'].str.split('-', expand=True)
            demux['Index1'] = index_split[0] if len(index_split.columns) > 0 else ''
            demux['Index2'] = index_split[1] if len(index_split.columns) > 1 else ''
        
        merged = pd.merge(demux, quality, on='SampleID', how='left')
        
        # Select and format columns
        display_cols = ['SampleID', 'Index1', 'Index2', '# Reads', 'Yield', '% Q30', 
                       'Mean Quality Score (PF)', '% Perfect Index Reads']
        merged_subset = merged[display_cols].copy()
        
        # Format the data
        merged_subset['# Reads'] = merged_subset['# Reads'].apply(lambda x: f'{x:,}')
        merged_subset['Yield'] = merged_subset['Yield'].apply(lambda x: f'{x/1e9:.2f} Gb' if pd.notna(x) else 'N/A')
        merged_subset['% Q30'] = merged_subset['% Q30'].apply(lambda x: f'{x:.2%}' if pd.notna(x) else 'N/A')
        merged_subset['Mean Quality Score (PF)'] = merged_subset['Mean Quality Score (PF)'].apply(
            lambda x: f'{x:.2f}' if pd.notna(x) else 'N/A')
        merged_subset['% Perfect Index Reads'] = merged_subset['% Perfect Index Reads'].apply(
            lambda x: f'{x:.2%}')
        
        # Rename columns for display
        merged_subset.columns = ['Sample ID', 'Index1', 'Index2', 'Reads', 'Yield', 'Q30%', 
                                'Mean QS', 'Perfect Index%']
        
        # Split into pages (20 rows per page to avoid cutoff)
        rows_per_page = 20
        total_pages = (len(merged_subset) + rows_per_page - 1) // rows_per_page
        
        for page_num in range(total_pages):
            start_idx = page_num * rows_per_page
            end_idx = min((page_num + 1) * rows_per_page, len(merged_subset))
            page_data = merged_subset.iloc[start_idx:end_idx]
            
            fig = plt.figure(figsize=(11, 8.5))
            title = f'Sample Details (Page {page_num + 1}/{total_pages})'
            fig.suptitle(title, fontsize=16, weight='bold', y=0.98)
            
            ax = plt.subplot(111)
            ax.axis('tight')
            ax.axis('off')
            
            # Create table
            table = ax.table(cellText=page_data.values,
                           colLabels=page_data.columns,
                           cellLoc='left',
                           loc='upper left',
                           colWidths=[0.14, 0.11, 0.11, 0.11, 0.11, 0.09, 0.10, 0.13])
            
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.8)
            
            # Style header
            for i in range(len(page_data.columns)):
                cell = table[(0, i)]
                cell.set_facecolor(COLORS['primary'])
                cell.set_text_props(weight='bold', color='white')
            
            # Alternate row colors
            for i in range(1, len(page_data) + 1):
                for j in range(len(page_data.columns)):
                    cell = table[(i, j)]
                    if i % 2 == 0:
                        cell.set_facecolor('#f0f0f0')
            
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()


def main():
    """Main function to parse arguments and generate report."""
    parser = argparse.ArgumentParser(
        description='Generate sequencing run QC report from BCLConvert output files',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input folder containing sequencing run files')
    parser.add_argument('-o', '--output', required=True,
                       help='Output PDF report file path')
    
    args = parser.parse_args()
    
    # Validate input folder
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input folder does not exist: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if not input_path.is_dir():
        print(f"Error: Input path is not a directory: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Ensure output directory exists
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Parse data
        print("Starting report generation...")
        sys.stdout.flush()
        
        parser_obj = SequencingRunParser(input_path)
        data = parser_obj.parse_all()
        
        # Generate report
        report = ReportGenerator(data)
        report.generate(args.output)
        
        print("\n✓ Report generation complete!")
        print(f"✓ Output: {args.output}")
        sys.stdout.flush()
        
    except Exception as e:
        print(f"\n✗ Error during report generation:", file=sys.stderr)
        print(f"  {type(e).__name__}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()