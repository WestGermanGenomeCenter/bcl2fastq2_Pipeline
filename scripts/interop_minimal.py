#!/usr/bin/env python3
"""
Sequencing Run Report Generator

Parses Illumina BCLConvert output files and generates a concise, professional
4-page PDF QC report suitable for sharing with non-specialists.

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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch

warnings.filterwarnings('ignore')

# ── Design system ────────────────────────────────────────────────────────────
NAVY    = '#1B2A4A'
TEAL    = '#0D7377'
TEAL_LT = '#14A085'
AMBER   = '#E8A838'
RED     = '#C0392B'
GREY_BG = '#F4F6F9'
GREY_MID= '#BDC3C7'
WHITE   = '#FFFFFF'
TEXT_DK = '#1B2A4A'
TEXT_MD = '#4A5568'

plt.rcParams.update({
    'font.family':      'DejaVu Sans',
    'axes.spines.top':  False,
    'axes.spines.right':False,
    'axes.edgecolor':   GREY_MID,
    'axes.labelcolor':  TEXT_MD,
    'xtick.color':      TEXT_MD,
    'ytick.color':      TEXT_MD,
    'text.color':       TEXT_DK,
    'figure.facecolor': WHITE,
    'axes.facecolor':   GREY_BG,
    'grid.color':       WHITE,
    'grid.linewidth':   1.0,
})


# ── Helpers ──────────────────────────────────────────────────────────────────

def fmt_num(n, decimals=0):
    """Format a number with thousand separators."""
    if pd.isna(n):
        return 'N/A'
    return f'{n:,.{decimals}f}'

def fmt_pct(v, decimals=1):
    """Format a fraction (0-1) as percentage string."""
    if pd.isna(v):
        return 'N/A'
    return f'{v * 100:.{decimals}f}%'

def pass_fail_color(value, threshold, invert=False):
    """Return TEAL_LT (pass) or RED (fail) based on threshold comparison."""
    if invert:
        return RED if value > threshold else TEAL_LT
    return TEAL_LT if value >= threshold else RED

def draw_metric_card(ax, title, value, subtitle='', value_color=TEAL, bg=GREY_BG):
    """Draw a KPI card on a given axes."""
    ax.set_facecolor(bg)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    # Title
    ax.text(0.5, 0.78, title, ha='center', va='center',
            fontsize=9, color=TEXT_MD, fontweight='normal')
    # Value
    ax.text(0.5, 0.44, value, ha='center', va='center',
            fontsize=20, color=value_color, fontweight='bold')
    # Subtitle
    if subtitle:
        ax.text(0.5, 0.15, subtitle, ha='center', va='center',
                fontsize=8, color=TEXT_MD, style='italic')
    # Border
    for spine in ax.spines.values():
        spine.set_visible(False)
    rect = FancyBboxPatch((0.02, 0.02), 0.96, 0.96,
                           boxstyle='round,pad=0.02',
                           linewidth=1.5, edgecolor=GREY_MID,
                           facecolor=WHITE, transform=ax.transAxes, zorder=0)
    ax.add_patch(rect)


# ── Parser ───────────────────────────────────────────────────────────────────

class SequencingRunParser:
    def __init__(self, input_folder):
        self.input_folder = Path(input_folder)
        self.data = {}

    def parse_all(self):
        print("Parsing sequencing run files...")
        self._parse_run_info()
        self._parse_demultiplex_stats()
        self._parse_quality_metrics()
        self._parse_adapter_metrics()
        self._parse_top_unknown_barcodes()
        print("Parsing complete.")
        return self.data

    def _parse_run_info(self):
        fp = self.input_folder / 'RunInfo.xml'
        if not fp.exists():
            print(f"  Warning: {fp.name} not found")
            return
        root = ET.parse(fp).getroot()
        run = root.find('Run')
        reads = []
        for r in run.findall('.//Read'):
            reads.append({
                'number':   r.get('Number'),
                'cycles':   int(r.get('NumCycles')),
                'is_index': r.get('IsIndexedRead') == 'Y'
            })
        layout = run.find('.//FlowcellLayout')
        self.data['run_info'] = {
            'run_id':    run.get('Id'),
            'flowcell':  run.find('Flowcell').text,
            'instrument':run.find('Instrument').text,
            'date':      run.find('Date').text,
            'reads':     reads,
            'layout': {
                'lanes':    int(layout.get('LaneCount')),
                'surfaces': int(layout.get('SurfaceCount')),
                'swaths':   int(layout.get('SwathCount')),
                'tiles':    int(layout.get('TileCount')),
            }
        }

    def _parse_demultiplex_stats(self):
        fp = self.input_folder / 'Demultiplex_Stats.csv'
        if not fp.exists():
            print(f"  Warning: {fp.name} not found"); return
        self.data['demux_stats'] = pd.read_csv(fp)

    def _parse_quality_metrics(self):
        fp = self.input_folder / 'Quality_Metrics.csv'
        if not fp.exists():
            print(f"  Warning: {fp.name} not found"); return
        self.data['quality_metrics'] = pd.read_csv(fp)

    def _parse_adapter_metrics(self):
        fp = self.input_folder / 'Adapter_Metrics.csv'
        if not fp.exists():
            print(f"  Warning: {fp.name} not found"); return
        self.data['adapter_metrics'] = pd.read_csv(fp)

    def _parse_top_unknown_barcodes(self):
        fp = self.input_folder / 'Top_Unknown_Barcodes.csv'
        if not fp.exists():
            print(f"  Warning: {fp.name} not found"); return
        self.data['unknown_barcodes'] = pd.read_csv(fp)


# ── Report Generator ─────────────────────────────────────────────────────────

class ReportGenerator:
    def __init__(self, data):
        self.data = data
        self._precompute_summaries()

    # ── Pre-computation ───────────────────────────────────────────────────────

    def _precompute_summaries(self):
        """Derive all key numbers once so pages can reuse them."""
        d = self.data
        s = {}

        # Run info
        if 'run_info' in d:
            ri = d['run_info']
            data_reads  = [r for r in ri['reads'] if not r['is_index']]
            index_reads = [r for r in ri['reads'] if r['is_index']]
            s['run_id']     = ri.get('run_id', 'N/A')
            s['flowcell']   = ri.get('flowcell', 'N/A')
            s['instrument'] = ri.get('instrument', 'N/A')
            s['date']       = ri.get('date', 'N/A')[:10] if ri.get('date') else 'N/A'
            s['read_type']  = 'Paired-End' if len(data_reads) > 1 else 'Single-Read'
            s['read_cycles']= ' + '.join(str(r['cycles']) for r in data_reads)
            s['index_cycles']= ' + '.join(str(r['cycles']) for r in index_reads) if index_reads else '—'
            layout = ri.get('layout', {})
            s['lanes']      = layout.get('lanes', 'N/A')
            s['tiles_per_lane'] = (layout.get('surfaces', 0) *
                                   layout.get('swaths', 0) *
                                   layout.get('tiles', 0))

        # Demux stats
        if 'demux_stats' in d:
            df = d['demux_stats']
            df_s = df[df['SampleID'] != 'Undetermined']
            df_u = df[df['SampleID'] == 'Undetermined']
            s['num_samples']      = len(df_s)
            s['total_reads']      = int(df['# Reads'].sum())
            s['assigned_reads']   = int(df_s['# Reads'].sum())
            s['undetermined_reads']= int(df_u['# Reads'].sum()) if len(df_u) > 0 else 0
            s['assigned_pct']     = s['assigned_reads'] / s['total_reads'] if s['total_reads'] else 0
            s['undetermined_pct'] = s['undetermined_reads'] / s['total_reads'] if s['total_reads'] else 0
            s['avg_reads']        = s['assigned_reads'] / s['num_samples'] if s['num_samples'] else 0
            s['perfect_index_mean']= df_s['% Perfect Index Reads'].mean()
            s['demux_df']         = df_s.copy()

        # Quality metrics
        if 'quality_metrics' in d:
            df = d['quality_metrics']
            df_s = df[df['SampleID'] != 'Undetermined']
            sample_q = df_s.groupby('SampleID').agg(
                Yield=('Yield', 'sum'),
                YieldQ30=('YieldQ30', 'sum'),
                Q30=('% Q30', 'mean'),
                MeanQS=('Mean Quality Score (PF)', 'mean')
            ).reset_index()
            s['total_yield_gb']   = sample_q['Yield'].sum() / 1e9
            s['avg_q30']          = sample_q['Q30'].mean()
            s['avg_mean_qs']      = sample_q['MeanQS'].mean()
            s['quality_df']       = sample_q

        # Adapter metrics
        if 'adapter_metrics' in d:
            df = d['adapter_metrics']
            df_s = df[df['Sample_ID'] != 'Undetermined']
            sample_a = df_s.groupby('Sample_ID').agg(
                AdapterBases=('AdapterBases', 'sum'),
                SampleBases=('SampleBases', 'sum'),
                AdapterPct=('% Adapter Bases', 'mean')
            ).reset_index()
            s['avg_adapter_pct']  = sample_a['AdapterPct'].mean()
            s['adapter_df']       = sample_a

        # Unknown barcodes
        if 'unknown_barcodes' in d:
            s['unknown_df'] = d['unknown_barcodes']

        self.summary = s

    # ── PDF entry point ───────────────────────────────────────────────────────

    def generate(self, output_path):
        print(f"Generating report → {output_path}")
        with PdfPages(output_path) as pdf:
            print("  Page 1: Run Overview…", end='', flush=True)
            self._page_overview(pdf);     print(" ✓")
            print("  Page 2: Quality Summary…", end='', flush=True)
            self._page_quality(pdf);      print(" ✓")
            print("  Page 3: Sample Distribution…", end='', flush=True)
            self._page_samples(pdf);      print(" ✓")
            print("  Page 4: QC Flags & Barcodes…", end='', flush=True)
            self._page_flags(pdf);        print(" ✓")
            print("  Page 5: Glossary & Further Reading…", end='', flush=True)
            self._page_glossary(pdf);     print(" ✓")

            # PDF metadata
            info = pdf.infodict()
            info['Title']   = f'Sequencing QC Report – {self.summary.get("flowcell","")}'
            info['Author']  = 'Sequencing Report Generator'
            info['Subject'] = 'Illumina Run Quality Control'

        print(f"\n✓ Report saved: {output_path}")

    # ── Shared header/footer helpers ──────────────────────────────────────────

    def _add_header(self, fig, title, page_num, total_pages=5):
        """Draw a coloured header band and footer line."""
        # Header band
        fig.add_axes([0, 0.945, 1, 0.055]).set_axis_off()
        header_ax = fig.axes[-1]
        header_ax.set_facecolor(NAVY)
        header_ax.set_xlim(0, 1); header_ax.set_ylim(0, 1)
        header_ax.text(0.02, 0.5, '● Sequencing QC Report',
                       color=WHITE, fontsize=9, va='center', alpha=0.7)
        header_ax.text(0.5, 0.5, title,
                       color=WHITE, fontsize=13, fontweight='bold',
                       ha='center', va='center')
        fc = self.summary.get('flowcell', '')
        header_ax.text(0.98, 0.5, f'{fc}  |  {page_num}/{total_pages}',
                       color=WHITE, fontsize=8, va='center', ha='right', alpha=0.7)
        # Footer
        footer = fig.add_axes([0.04, 0.005, 0.92, 0.012])
        footer.set_axis_off()
        footer.set_facecolor(NAVY)
        footer.set_xlim(0,1); footer.set_ylim(0,1)
        ts = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')
        footer.text(0.5, 0.5, f'Generated {ts}   ·   Illumina Sequencing QC',
                    color=WHITE, fontsize=7, ha='center', va='center', alpha=0.8)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PAGE 1 – Run Overview (key numbers at a glance)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    def _page_overview(self, pdf):
        s = self.summary
        fig = plt.figure(figsize=(11, 8.5), facecolor=WHITE)
        self._add_header(fig, 'Run Overview', 1)

        # ── Run identity block ─────────────────────────────────────────────
        # Left column: instrument details
        id_ax = fig.add_axes([0.03, 0.75, 0.40, 0.18])
        id_ax.set_axis_off()
        id_ax.set_facecolor(WHITE)
        id_ax.set_xlim(0,1); id_ax.set_ylim(0,1)

        lines = [
            ('Run ID',      s.get('run_id','N/A')),
            ('Flowcell ID', s.get('flowcell','N/A')),
            ('Instrument',  s.get('instrument','N/A')),
            ('Run Date',    s.get('date','N/A')),
            ('Read Type',   f"{s.get('read_type','N/A')}  ({s.get('read_cycles','?')} bp data  /  {s.get('index_cycles','?')} bp index)"),
            ('Lanes',       str(s.get('lanes','N/A'))),
        ]
        y = 0.92
        id_ax.text(0, 1.02, 'Run Configuration', fontsize=11, fontweight='bold', color=NAVY)
        for label, val in lines:
            id_ax.text(0.00, y, f'{label}:', fontsize=9, color=TEXT_MD)
            id_ax.text(0.40, y, val,         fontsize=9, color=TEXT_DK, fontweight='bold')
            y -= 0.17

        # ── KPI cards row ──────────────────────────────────────────────────
        # 5 cards across the full width
        card_w, card_h = 0.17, 0.19
        card_y = 0.52
        gap = 0.016
        starts = [0.03 + i*(card_w+gap) for i in range(5)]

        # Card 1: Total clusters (= total reads)
        total_reads = s.get('total_reads', 0)
        val_str = (f'{total_reads/1e6:.1f} M' if total_reads < 1e9
                   else f'{total_reads/1e9:.2f} B')
        cax = fig.add_axes([starts[0], card_y, card_w, card_h])
        draw_metric_card(cax, 'Total Clusters', val_str,
                         subtitle='Reads generated this run',
                         value_color=NAVY)

        # Card 2: Assigned reads %
        ap = s.get('assigned_pct', 0)
        ap_col = pass_fail_color(ap, 0.80)
        cax = fig.add_axes([starts[1], card_y, card_w, card_h])
        draw_metric_card(cax, 'Assigned Reads', fmt_pct(ap),
                         subtitle='Successfully demultiplexed',
                         value_color=ap_col)

        # Card 3: Q30 %
        q30 = s.get('avg_q30', 0)
        q30_col = pass_fail_color(q30, 0.80)
        cax = fig.add_axes([starts[2], card_y, card_w, card_h])
        draw_metric_card(cax, '% Bases ≥ Q30', fmt_pct(q30),
                         subtitle='Avg. across all samples (≥80% good)',
                         value_color=q30_col)

        # Card 4: Total yield
        yld = s.get('total_yield_gb', 0)
        cax = fig.add_axes([starts[3], card_y, card_w, card_h])
        draw_metric_card(cax, 'Total Yield', f'{yld:.1f} Gb',
                         subtitle='Gigabases of sequencing data',
                         value_color=TEAL)

        # Card 5: Perfect index rate
        pir = s.get('perfect_index_mean', 0)
        pir_col = pass_fail_color(pir, 0.85)
        cax = fig.add_axes([starts[4], card_y, card_w, card_h])
        draw_metric_card(cax, 'Perfect Index Rate', fmt_pct(pir),
                         subtitle='Exact barcode match (≥85% good)',
                         value_color=pir_col)

        # ── What do these metrics mean? ────────────────────────────────────
        expl_ax = fig.add_axes([0.03, 0.36, 0.93, 0.14])
        expl_ax.set_axis_off()
        expl_ax.set_facecolor(GREY_BG)
        expl_ax.set_xlim(0,1); expl_ax.set_ylim(0,1)

        rect = FancyBboxPatch((0,0), 1, 1, boxstyle='round,pad=0.01',
                               linewidth=1, edgecolor=GREY_MID,
                               facecolor=GREY_BG, transform=expl_ax.transAxes)
        expl_ax.add_patch(rect)
        expl_ax.text(0.01, 0.90, 'How to read this report',
                     fontsize=9, fontweight='bold', color=NAVY, va='top')
        explanations = (
            "Total Clusters — the number of DNA fragments detected and sequenced on the flowcell.  "
            "Assigned Reads — fraction of clusters matched to a sample barcode; low values suggest barcode or sample-loading issues.  "
            "% Bases ≥ Q30 — a Q30 base has a 99.9% accuracy; this is the primary quality benchmark for Illumina runs.  "
            "Total Yield — total gigabases (Gb) of sequence produced, excluding unassigned reads.  "
            "Perfect Index Rate — fraction of reads where the barcode matched exactly (no mismatches); "
            "a low rate can increase undetermined reads."
        )
        expl_ax.text(0.01, 0.65, explanations,
                     fontsize=8, color=TEXT_MD, va='top', wrap=True,
                     linespacing=1.55)

        # ── Overall QC verdict ─────────────────────────────────────────────
        verdict_ax = fig.add_axes([0.03, 0.18, 0.93, 0.16])
        verdict_ax.set_axis_off()
        verdict_ax.set_xlim(0,1); verdict_ax.set_ylim(0,1)

        checks = [
            ('Demultiplexing ≥ 80%', s.get('assigned_pct', 0),  0.80, False),
            ('Q30 ≥ 80%',            s.get('avg_q30', 0),        0.80, False),
            ('Perfect Index ≥ 85%',  s.get('perfect_index_mean',0), 0.85, False),
        ]
        if 'avg_adapter_pct' in s:
            checks.append(('Adapter < 1%', s['avg_adapter_pct'], 0.01, True))

        all_pass = all(
            (v < thr if inv else v >= thr)
            for _, v, thr, inv in checks
        )
        verdict_col   = TEAL_LT if all_pass else AMBER
        verdict_label = '✔  All key metrics PASS' if all_pass else '⚠  One or more metrics need attention'

        verdict_ax.text(0.0, 0.92, 'QC Verdict', fontsize=10, fontweight='bold', color=NAVY)
        verdict_ax.text(0.0, 0.65, verdict_label,
                        fontsize=13, fontweight='bold', color=verdict_col)

        x = 0.0
        for label, val, thr, inv in checks:
            passed = (val < thr) if inv else (val >= thr)
            icon   = '✔' if passed else '✗'
            col    = TEAL_LT if passed else RED
            verdict_ax.text(x, 0.28, f'{icon} {label}',
                            fontsize=8.5, color=col, va='center')
            x += 0.26

        # Secondary stats
        verdict_ax.text(0.0, 0.02,
            f"Samples: {s.get('num_samples','N/A')}   ·   "
            f"Avg reads/sample: {fmt_num(s.get('avg_reads',0))}   ·   "
            f"Undetermined: {fmt_pct(s.get('undetermined_pct',0))}   ·   "
            f"Mean quality score: {s.get('avg_mean_qs',0):.1f}",
            fontsize=8, color=TEXT_MD)

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PAGE 2 – Quality Metrics
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    def _page_quality(self, pdf):
        s = self.summary
        fig = plt.figure(figsize=(11, 8.5), facecolor=WHITE)
        self._add_header(fig, 'Quality Metrics', 2)

        has_quality  = 'quality_df' in s
        has_adapter  = 'adapter_df' in s

        # ── Annotation boxes ───────────────────────────────────────────────
        def note(fig, x, y, w, h, text):
            ax = fig.add_axes([x, y, w, h])
            ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
            rect = FancyBboxPatch((0,0), 1, 1, boxstyle='round,pad=0.02',
                                   linewidth=0.8, edgecolor=GREY_MID,
                                   facecolor=GREY_BG, transform=ax.transAxes)
            ax.add_patch(rect)
            ax.text(0.05, 0.5, text, fontsize=7.5, color=TEXT_MD,
                    va='center', linespacing=1.5)

        # ── Plot 1: Q30 distribution ───────────────────────────────────────
        if has_quality:
            ax1 = fig.add_axes([0.07, 0.57, 0.38, 0.32])
            q30 = s['quality_df']['Q30'].values * 100
            n, bins, patches = ax1.hist(q30, bins=25, color=TEAL, edgecolor=WHITE, linewidth=0.4)
            for p, b in zip(patches, bins):
                p.set_facecolor(TEAL_LT if b >= 80 else RED)
            ax1.axvline(q30.mean(), color=NAVY, linestyle='--', linewidth=1.5,
                        label=f'Mean: {q30.mean():.1f}%')
            ax1.axvline(80, color=AMBER, linestyle=':', linewidth=1.5,
                        label='Target ≥ 80%')
            ax1.set_xlabel('% Bases with Quality ≥ Q30', fontsize=9)
            ax1.set_ylabel('Number of Samples', fontsize=9)
            ax1.set_title('Q30 Score Distribution', fontsize=11, fontweight='bold',
                          color=NAVY, pad=8)
            ax1.legend(fontsize=8)
            ax1.xaxis.set_major_locator(plt.MaxNLocator(8))
            ax1.yaxis.set_major_locator(plt.MaxNLocator(6))
            note(fig, 0.07, 0.49, 0.38, 0.07,
                 "Q30 means ≥99.9% base accuracy. Samples left of the dashed line may "
                 "need closer inspection. Most runs target ≥80%.")

        # ── Plot 2: Mean quality score distribution ────────────────────────
        if has_quality:
            ax2 = fig.add_axes([0.57, 0.57, 0.38, 0.32])
            qs = s['quality_df']['MeanQS'].values
            ax2.hist(qs, bins=25, color=TEAL, edgecolor=WHITE, linewidth=0.4)
            ax2.axvline(qs.mean(), color=NAVY, linestyle='--', linewidth=1.5,
                        label=f'Mean: {qs.mean():.1f}')
            ax2.axvline(30, color=AMBER, linestyle=':', linewidth=1.5,
                        label='Target ≥ Q30')
            ax2.set_xlabel('Mean Phred Quality Score', fontsize=9)
            ax2.set_ylabel('Number of Samples', fontsize=9)
            ax2.set_title('Mean Quality Score Distribution', fontsize=11,
                          fontweight='bold', color=NAVY, pad=8)
            ax2.legend(fontsize=8)
            ax2.xaxis.set_major_locator(plt.MaxNLocator(8))
            ax2.yaxis.set_major_locator(plt.MaxNLocator(6))
            note(fig, 0.57, 0.49, 0.38, 0.07,
                 "The Phred score (Q-score) measures base-call accuracy. "
                 "Q30 = 99.9% accurate; Q20 = 99%. Higher is better.")

        # ── Plot 3: Adapter contamination ────────────────────────────────
        if has_adapter:
            ax3 = fig.add_axes([0.07, 0.13, 0.38, 0.30])
            adp = s['adapter_df']['AdapterPct'].values * 100
            n, bins, patches = ax3.hist(adp, bins=25, color=AMBER, edgecolor=WHITE, linewidth=0.4)
            for p, b in zip(patches, bins):
                p.set_facecolor(RED if b >= 1.0 else AMBER)
            ax3.axvline(adp.mean(), color=NAVY, linestyle='--', linewidth=1.5,
                        label=f'Mean: {adp.mean():.3f}%')
            ax3.axvline(1.0, color=RED, linestyle=':', linewidth=1.5,
                        label='Warning ≥ 1%')
            ax3.set_xlabel('% Adapter Bases per Sample', fontsize=9)
            ax3.set_ylabel('Number of Samples', fontsize=9)
            ax3.set_title('Adapter Contamination Distribution', fontsize=11,
                          fontweight='bold', color=NAVY, pad=8)
            ax3.legend(fontsize=8)
            ax3.xaxis.set_major_locator(plt.MaxNLocator(8))
            ax3.yaxis.set_major_locator(plt.MaxNLocator(6))
            note(fig, 0.07, 0.05, 0.38, 0.07,
                 "Adapter sequences are sequencing kit artefacts, not biological sequence. "
                 "Values >1% suggest short insert sizes or library prep issues.")

        # ── KPI summary panel ──────────────────────────────────────────────
        kpi_ax = fig.add_axes([0.57, 0.08, 0.38, 0.38])
        kpi_ax.set_axis_off(); kpi_ax.set_xlim(0,1); kpi_ax.set_ylim(0,1)
        rect = FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.02',
                               linewidth=1,edgecolor=GREY_MID,
                               facecolor=GREY_BG,transform=kpi_ax.transAxes)
        kpi_ax.add_patch(rect)
        kpi_ax.text(0.5, 0.94, 'Quality Summary', fontsize=11, fontweight='bold',
                    color=NAVY, ha='center', va='top')

        rows = []
        if has_quality:
            qd = s['quality_df']
            rows += [
                ('Avg Q30 (%)',         fmt_pct(s['avg_q30']),
                 pass_fail_color(s['avg_q30'], 0.80)),
                ('Min Q30 – any sample',fmt_pct(qd['Q30'].min()),
                 pass_fail_color(qd['Q30'].min(), 0.75)),
                ('Avg Mean Quality Score', f"{s['avg_mean_qs']:.2f}",
                 pass_fail_color(s['avg_mean_qs'], 30)),
                ('Total Yield',         f"{s['total_yield_gb']:.1f} Gb", TEAL),
                ('Avg Yield / Sample',  f"{s['total_yield_gb']/s['num_samples']*1000:.1f} Mb"
                 if s.get('num_samples') else 'N/A', TEAL),
            ]
        if has_adapter:
            avg_a = s['avg_adapter_pct']
            rows += [
                ('Avg Adapter Content', f"{avg_a*100:.3f}%",
                 pass_fail_color(avg_a, 0.01, invert=True)),
                ('Max Adapter – any sample',
                 f"{s['adapter_df']['AdapterPct'].max()*100:.3f}%",
                 pass_fail_color(s['adapter_df']['AdapterPct'].max(), 0.01, invert=True)),
            ]

        y = 0.82
        for label, val, col in rows:
            kpi_ax.text(0.06, y, label, fontsize=9, color=TEXT_MD)
            kpi_ax.text(0.94, y, val,   fontsize=9, color=col,
                        fontweight='bold', ha='right')
            kpi_ax.axhline(y - 0.035, xmin=0.04, xmax=0.96,
                           color=GREY_MID, linewidth=0.5,
                           transform=kpi_ax.transAxes)
            y -= 0.10

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PAGE 3 – Sample Distribution
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    def _page_samples(self, pdf):
        s = self.summary
        if 'demux_df' not in s:
            return
        fig = plt.figure(figsize=(11, 8.5), facecolor=WHITE)
        self._add_header(fig, 'Sample Distribution', 3)

        demux = s['demux_df'].sort_values('# Reads', ascending=False).reset_index(drop=True)
        top_n  = min(25, len(demux))
        top    = demux.head(top_n)

        # ── Plot 1: Reads per sample bar chart ─────────────────────────────
        ax1 = fig.add_axes([0.06, 0.50, 0.55, 0.38])
        colors = [TEAL_LT if v >= s.get('avg_reads',0)*0.5 else RED
                  for v in top['# Reads']]
        bars = ax1.barh(range(len(top)), top['# Reads']/1e6,
                        color=colors, edgecolor=WHITE, linewidth=0.3)
        ax1.set_yticks(range(len(top)))
        ax1.set_yticklabels(top['SampleID'], fontsize=7.5)
        ax1.invert_yaxis()
        ax1.set_xlabel('Reads (millions)', fontsize=9)
        ax1.set_title(f'Reads per Sample  (top {top_n} shown)', fontsize=11,
                      fontweight='bold', color=NAVY, pad=8)
        ax1.xaxis.set_major_locator(plt.MaxNLocator(7))
        avg_m = s.get('avg_reads', 0)/1e6
        ax1.axvline(avg_m, color=AMBER, linestyle='--', linewidth=1.5,
                    label=f'Average: {avg_m:.1f} M')
        ax1.legend(fontsize=8)
        note_text = (f"Red bars = samples with <50% of average read depth "
                     f"({fmt_num(s.get('avg_reads',0))} reads avg). "
                     f"Unevenness may reflect library prep variation.")
        self._note_box(fig, 0.06, 0.43, 0.55, 0.06, note_text)

        # ── Plot 2: Read-depth distribution histogram ─────────────────────
        ax2 = fig.add_axes([0.70, 0.50, 0.26, 0.38])
        ax2.hist(demux['# Reads']/1e6, bins=20, color=TEAL, edgecolor=WHITE, linewidth=0.4,
                 orientation='horizontal')
        ax2.axhline(demux['# Reads'].mean()/1e6, color=NAVY, linestyle='--',
                    linewidth=1.5, label='Mean')
        ax2.set_ylabel('Reads per Sample (M)', fontsize=9)
        ax2.set_xlabel('Sample Count', fontsize=9)
        ax2.set_title('Read Depth\nDistribution', fontsize=11,
                      fontweight='bold', color=NAVY, pad=8)
        ax2.legend(fontsize=8)
        ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax2.yaxis.set_major_locator(plt.MaxNLocator(7))

        # ── Stats table ────────────────────────────────────────────────────
        stats_ax = fig.add_axes([0.06, 0.08, 0.88, 0.33])
        stats_ax.set_axis_off(); stats_ax.set_xlim(0,1); stats_ax.set_ylim(0,1)
        rect = FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.02',
                               linewidth=1,edgecolor=GREY_MID,
                               facecolor=GREY_BG,transform=stats_ax.transAxes)
        stats_ax.add_patch(rect)
        stats_ax.text(0.5, 0.94, 'Sample Statistics', fontsize=11, fontweight='bold',
                      color=NAVY, ha='center', va='top')

        # 2-column layout
        col_data = [
            ('Total samples', fmt_num(s.get('num_samples',0))),
            ('Total reads (all)', fmt_num(s.get('total_reads',0))),
            ('Assigned reads', f"{fmt_num(s.get('assigned_reads',0))}  ({fmt_pct(s.get('assigned_pct',0))})"),
            ('Undetermined reads', f"{fmt_num(s.get('undetermined_reads',0))}  ({fmt_pct(s.get('undetermined_pct',0))})"),
        ]
        col_data2 = [
            ('Reads / sample – mean',   fmt_num(demux['# Reads'].mean())),
            ('Reads / sample – median', fmt_num(demux['# Reads'].median())),
            ('Reads / sample – min',    fmt_num(demux['# Reads'].min())),
            ('Reads / sample – max',    fmt_num(demux['# Reads'].max())),
        ]

        y = 0.78
        for (l1,v1),(l2,v2) in zip(col_data, col_data2):
            stats_ax.text(0.03, y, l1, fontsize=9, color=TEXT_MD)
            stats_ax.text(0.25, y, v1, fontsize=9, color=TEXT_DK, fontweight='bold')
            stats_ax.text(0.53, y, l2, fontsize=9, color=TEXT_MD)
            stats_ax.text(0.76, y, v2, fontsize=9, color=TEXT_DK, fontweight='bold')
            stats_ax.axhline(y-0.05, xmin=0.02, xmax=0.98, color=GREY_MID,
                             linewidth=0.5, transform=stats_ax.transAxes)
            y -= 0.16

        # Low-depth warning
        low_n = (demux['# Reads'] < s.get('avg_reads',0)*0.5).sum()
        if low_n > 0:
            stats_ax.text(0.03, 0.12,
                f"⚠  {low_n} sample(s) have fewer than 50% of average reads — "
                "consider re-sequencing or checking library QC.",
                fontsize=8.5, color=RED, va='center')

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PAGE 4 – QC Flags & Top Unknown Barcodes
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    def _page_flags(self, pdf):
        s = self.summary
        fig = plt.figure(figsize=(11, 8.5), facecolor=WHITE)
        self._add_header(fig, 'QC Flags & Barcode Analysis', 4)

        # ── QC checklist ───────────────────────────────────────────────────
        chk_ax = fig.add_axes([0.04, 0.62, 0.42, 0.30])
        chk_ax.set_axis_off(); chk_ax.set_xlim(0,1); chk_ax.set_ylim(0,1)
        rect = FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.02',
                               linewidth=1,edgecolor=GREY_MID,
                               facecolor=GREY_BG,transform=chk_ax.transAxes)
        chk_ax.add_patch(rect)
        chk_ax.text(0.05, 0.93, 'QC Checklist', fontsize=11, fontweight='bold',
                    color=NAVY, va='top')

        checks = [
            ('Assigned reads ≥ 80%',     s.get('assigned_pct',0),      0.80, False,
             'Low: check barcode design or sample loading'),
            ('Undetermined reads < 20%',  s.get('undetermined_pct',0),  0.20, True,
             'High: barcode mismatch or index contamination'),
            ('Avg Q30 ≥ 80%',             s.get('avg_q30',0),           0.80, False,
             'Low: sequencing chemistry or cluster density issue'),
            ('Perfect index rate ≥ 85%',  s.get('perfect_index_mean',0),0.85, False,
             'Low: index read errors; may increase undetermined reads'),
        ]
        if 'avg_adapter_pct' in s:
            checks.append(
                ('Avg adapter content < 1%', s['avg_adapter_pct'],      0.01, True,
                 'High: short inserts or library prep problem')
            )

        y = 0.76
        for label, val, thr, inv, reason in checks:
            passed = (val < thr) if inv else (val >= thr)
            icon   = '✔' if passed else '✗'
            col    = TEAL_LT if passed else RED
            chk_ax.text(0.04, y, icon,   fontsize=12, color=col, va='center')
            chk_ax.text(0.12, y, label,  fontsize=9,  color=TEXT_DK, va='center')
            if not passed:
                chk_ax.text(0.12, y-0.065, f'↳ {reason}', fontsize=7.5,
                            color=RED, va='center', style='italic')
            y -= 0.155 if not passed else 0.13

        # ── Perfect index rate per-sample bars ────────────────────────────
        if 'demux_df' in s:
            pir_ax = fig.add_axes([0.55, 0.62, 0.40, 0.30])
            df = s['demux_df'].copy()
            df_sorted = df.sort_values('% Perfect Index Reads').reset_index(drop=True)
            top_show = min(20, len(df_sorted))
            df_show  = df_sorted.head(top_show)
            cols = [RED if v < 0.85 else TEAL_LT
                    for v in df_show['% Perfect Index Reads']]
            pir_ax.barh(range(len(df_show)),
                        df_show['% Perfect Index Reads']*100,
                        color=cols, edgecolor=WHITE, linewidth=0.3)
            pir_ax.axvline(85, color=AMBER, linestyle=':', linewidth=1.5,
                           label='Target 85%')
            pir_ax.set_yticks(range(len(df_show)))
            pir_ax.set_yticklabels(df_show['SampleID'], fontsize=7)
            pir_ax.set_xlabel('Perfect Index Rate (%)', fontsize=9)
            pir_ax.set_title(f'Perfect Index Rate — Lowest {top_show} Samples',
                             fontsize=10, fontweight='bold', color=NAVY, pad=6)
            pir_ax.legend(fontsize=8)
            pir_ax.xaxis.set_major_locator(plt.MaxNLocator(7))

        # ── Unknown barcodes table ─────────────────────────────────────────
        if 'unknown_df' in s:
            unk = s['unknown_df']
            unk_ax = fig.add_axes([0.04, 0.08, 0.91, 0.50])
            unk_ax.set_axis_off(); unk_ax.set_xlim(0,1); unk_ax.set_ylim(0,1)
            rect = FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.02',
                                   linewidth=1,edgecolor=GREY_MID,
                                   facecolor=GREY_BG,transform=unk_ax.transAxes)
            unk_ax.add_patch(rect)

            total_u  = unk['# Reads'].sum()
            total_a  = s.get('total_reads', total_u)
            unk_pct  = total_u / total_a if total_a else 0

            unk_ax.text(0.03, 0.94,
                        f'Top Unknown Barcodes — {fmt_num(total_u)} undetermined reads '
                        f'({fmt_pct(unk_pct)} of total)',
                        fontsize=10, fontweight='bold', color=NAVY, va='top')
            unk_ax.text(0.03, 0.86,
                        'These barcodes were detected but could not be matched to any sample. '
                        'Dominant unknown barcodes may indicate a missing sample, index swap, '
                        'or a sample-sheet error.',
                        fontsize=8, color=TEXT_MD, va='top')

            # Column headers
            cols_def = [
                ('Rank',          0.03, 0.06),
                ('Index 1',       0.08, 0.18),
                ('Index 2',       0.27, 0.18),
                ('Reads',         0.46, 0.14),
                ('% of Unknown',  0.61, 0.15),
                ('% of All Reads',0.77, 0.20),
            ]
            hdr_y = 0.75
            for hdr, x, _ in cols_def:
                unk_ax.add_patch(FancyBboxPatch(
                    (x, hdr_y-0.04), _-0.005, 0.08,
                    boxstyle='round,pad=0.005', linewidth=0,
                    facecolor=NAVY, transform=unk_ax.transAxes))
                unk_ax.text(x+0.005, hdr_y, hdr, fontsize=8,
                            color=WHITE, fontweight='bold', va='center')

            top_show = min(12, len(unk))
            row_h = 0.055
            for i, (_, row) in enumerate(unk.head(top_show).iterrows()):
                ry = hdr_y - 0.04 - (i+1)*row_h
                bg = '#EDF2F7' if i % 2 == 0 else WHITE
                unk_ax.add_patch(FancyBboxPatch(
                    (0.02, ry-0.005), 0.96, row_h-0.004,
                    boxstyle='round,pad=0.002', linewidth=0,
                    facecolor=bg, transform=unk_ax.transAxes))
                vals = [
                    str(i+1),
                    row.get('index',''),
                    row.get('index2',''),
                    fmt_num(row['# Reads']),
                    fmt_pct(row['% of Unknown Barcodes']),
                    fmt_pct(row['% of All Reads'], decimals=3),
                ]
                for val, (_, x, w) in zip(vals, cols_def):
                    unk_ax.text(x+0.005, ry + row_h*0.38, val,
                                fontsize=8, color=TEXT_DK, va='center')
        else:
            # No unknown barcodes data
            no_ax = fig.add_axes([0.04, 0.08, 0.91, 0.50])
            no_ax.set_axis_off()
            no_ax.text(0.5, 0.5, 'Top_Unknown_Barcodes.csv not found.',
                       ha='center', va='center', fontsize=11, color=TEXT_MD)

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PAGE 5 – Glossary & Further Reading
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    def _page_glossary(self, pdf):
        fig = plt.figure(figsize=(11, 8.5), facecolor=WHITE)
        self._add_header(fig, 'Understanding Your Sequencing Report', 5)

        # ── Intro sentence ────────────────────────────────────────────────
        intro_ax = fig.add_axes([0.04, 0.87, 0.92, 0.055])
        intro_ax.set_axis_off(); intro_ax.set_xlim(0,1); intro_ax.set_ylim(0,1)
        intro_ax.text(0.0, 0.5,
            'This page explains the key terms used in this report. '
            'No prior knowledge of sequencing is needed — if anything is still unclear, '
            'please contact us and we will be happy to help.',
            fontsize=9, color=TEXT_MD, va='center', linespacing=1.6)

        # ── Helper to draw a glossary card ────────────────────────────────
        def gcard(fig, x, y, w, h, term, emoji, definition, threshold=None):
            ax = fig.add_axes([x, y, w, h])
            ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
            rect = FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.03',
                                   linewidth=1.2, edgecolor=TEAL,
                                   facecolor=WHITE, transform=ax.transAxes)
            ax.add_patch(rect)
            # Accent bar on left
            ax.add_patch(FancyBboxPatch((0,0),0.045,1,boxstyle='square,pad=0',
                                         linewidth=0, facecolor=TEAL,
                                         transform=ax.transAxes))
            ax.text(0.07, 0.80, f'{emoji}  {term}',
                    fontsize=10, fontweight='bold', color=NAVY, va='top')
            ax.text(0.07, 0.54, definition,
                    fontsize=8, color=TEXT_MD, va='top', linespacing=1.55,
                    wrap=True)
            if threshold:
                ax.text(0.07, 0.10, f'✔  {threshold}',
                        fontsize=7.5, color=TEAL_LT, va='bottom', fontweight='bold')

        # ── 6 glossary cards in a 3 × 2 grid ─────────────────────────────
        cw, ch = 0.285, 0.195
        gap_x, gap_y = 0.025, 0.018
        x0, y0 = 0.04, 0.645

        cards = [
            ('Total Clusters', '',
             'A "cluster" is a single DNA fragment that has been amplified into a bright '
             'spot on the flowcell. Each cluster is sequenced independently. '
             'More clusters = more data. Typical runs produce tens to hundreds of millions.',
             'More is generally better — depends on your experiment'),
            ('% Bases ≥ Q30', '',
             'Every base (A, T, C, G) the sequencer calls is given a quality score (Q-score). '
             'Q30 means there is a 1-in-1000 chance of error — i.e. 99.9% accuracy. '
             'This is the single most important quality metric for most analyses.',
             'Target: ≥ 80% of all bases should reach Q30'),
            ('Assigned Reads', '',
             'Each sample in a sequencing run is labelled with one or two short unique DNA tags called '
             'barcode or index. After sequencing, reads are matched to their sample by '
             'this barcode or barcode combination. "Assigned" reads are those successfully matched.',
             'Target: ≥ 80% of reads assigned to a sample'),
            ('Undetermined Reads', '',
             'Reads whose barcode could not be matched to any sample. A small fraction '
             '(< 5%) is normal. Higher values may mean a sample was accidentally left off '
             'the sample sheet, or barcodes were confused during library preparation.',
             'Concern: > 20% undetermined warrants investigation'),
            ('Adapter Content', '',
             'Adapters are short synthetic sequences added during library preparation to '
             'allow DNA to attach to the flowcell. If a DNA fragment is very short, the '
             'sequencer may read into the adapter. These bases are artefacts, not real biology.',
             'Target: < 1% adapter bases per sample'),
            ('Perfect Index Rate', '',
             'The percentage of reads where the barcode matched its expected sequence '
             'with zero errors. Even one mismatched base can cause a read to be assigned '
             'to the wrong sample or discarded. Higher rates mean cleaner sample separation.',
             'Target: ≥ 85% perfect index matches'),
        ]

        positions = [
            (x0,               y0),
            (x0 + cw + gap_x,  y0),
            (x0 + 2*(cw+gap_x),y0),
            (x0,               y0 - ch - gap_y),
            (x0 + cw + gap_x,  y0 - ch - gap_y),
            (x0 + 2*(cw+gap_x),y0 - ch - gap_y),
        ]
        for (term, emoji, defn, thr), (px, py) in zip(cards, positions):
            gcard(fig, px, py, cw, ch, term, emoji, defn, thr)

        # ── Further reading section ────────────────────────────────────────
        link_ax = fig.add_axes([0.04, 0.07, 0.92, 0.22])
        link_ax.set_axis_off(); link_ax.set_xlim(0,1); link_ax.set_ylim(0,1)
        rect = FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.02',
                               linewidth=1, edgecolor=GREY_MID,
                               facecolor=GREY_BG, transform=link_ax.transAxes)
        link_ax.add_patch(rect)

        link_ax.text(0.02, 0.90, '  Further Reading & Resources',
                     fontsize=10, fontweight='bold', color=NAVY, va='top')
        link_ax.text(0.02, 0.75,
            'The following Illumina resources explain the technology and quality metrics '
            'in more depth — all are freely available online:',
            fontsize=8.5, color=TEXT_MD, va='top')

        links = [
            ('How Illumina Sequencing Works  (short video + overview)',
             'illumina.com/science/technology/next-generation-sequencing/sequencing-technology.html'),
            ('Understanding Sequencing Data Quality',
             'illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html'),
            ('Sequencing Coverage & Depth — how much data do you need?',
             'illumina.com/science/technology/next-generation-sequencing/plan-experiments/coverage-depth-replicates.html'),
            ('BCLConvert Output Files — what each file contains',
             'support.illumina.com/sequencing/sequencing_software/bcl-convert/documentation.html'),
            ('Interpreting Run Quality Metrics (BaseSpace guide)',
             'support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreFilter_swBS.htm'),
        ]

        y = 0.55
        for title_txt, url in links:
            link_ax.text(0.025, y, '→', fontsize=9, color=TEAL, va='center', fontweight='bold')
            link_ax.text(0.05,  y, title_txt, fontsize=8.5, color=TEXT_DK, va='center',
                         fontweight='bold')
            link_ax.text(0.05,  y - 0.095, url, fontsize=7.5, color=TEAL, va='center',
                         style='italic')
            y -= 0.175

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # ── Utility ───────────────────────────────────────────────────────────────
    def _note_box(self, fig, x, y, w, h, text):
        ax = fig.add_axes([x, y, w, h])
        ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
        rect = FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.02',
                               linewidth=0.8,edgecolor=GREY_MID,
                               facecolor=GREY_BG,transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(0.02, 0.5, text, fontsize=7.5, color=TEXT_MD, va='center',
                linespacing=1.4)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Generate a concise Illumina sequencing QC report (≤4 pages).',
    )
    parser.add_argument('-i', '--input',  required=True,
                        help='Folder containing BCLConvert output files')
    parser.add_argument('-o', '--output', required=True,
                        help='Output PDF path')
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_dir():
        print(f"Error: {args.input} is not a directory", file=sys.stderr)
        sys.exit(1)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    try:
        run_parser = SequencingRunParser(input_path)
        data = run_parser.parse_all()
        ReportGenerator(data).generate(args.output)
        print(f"\n✓ Done — {args.output}")
    except Exception as exc:
        import traceback
        print(f"\n✗ Failed: {exc}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()