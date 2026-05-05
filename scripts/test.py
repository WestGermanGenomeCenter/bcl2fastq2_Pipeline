#!/usr/bin/env python3
"""
Illumina InterOp QC Report Generator
Compatible with NovaSeqX+, NextSeq, NextSeq2000, and MiSeq
Produces a concise 4-page PDF for sequencing groups.

Usage:
    python interop_more.py <run_folder> <output_pdf>

Requirements:
    conda install -c bioconda illumina-interop
    conda install pandas matplotlib seaborn numpy reportlab
"""

import sys
import os
from pathlib import Path
from io import BytesIO
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from interop import py_interop_run_metrics
    import interop.core as ic
    from reportlab.lib.pagesizes import letter, landscape
    from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                    Image, PageBreak, Table, TableStyle,
                                    HRFlowable)
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
except ImportError as e:
    print(f"Error: Missing required package: {e}")
    print("  conda install -c bioconda illumina-interop")
    print("  conda install pandas matplotlib seaborn numpy reportlab")
    sys.exit(1)

# ── Style ────────────────────────────────────────────────────────────────────
sns.set_style("whitegrid")
plt.rcParams.update({
    'figure.dpi': 150, 'savefig.dpi': 150,
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 13,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
})

LANE_COLORS = {1:'#0072B2', 2:'#E69F00', 3:'#009E73', 4:'#D55E00',
               5:'#CC79A7', 6:'#F0E442', 7:'#56B4E9', 8:'#999999'}

CHANNEL_COLORS = {'A':'#009E73', 'C':'#0072B2', 'G':'#444444', 'T':'#D55E00'}

# ── Useful reference links shown on the title page ───────────────────────────
REFERENCE_LINKS = [
    ("Illumina SAV (Sequencing Analysis Viewer)",
     "https://support.illumina.com/sequencing/sequencing_software/sequencing_analysis_viewer_sav.html"),
    ("InterOp File Format Documentation",
     "https://illumina.github.io/interop/index.html"),
    ("Understanding Sequencing Quality Scores (Phred / Q-scores)",
     "https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html"),

]


# ── Helpers ───────────────────────────────────────────────────────────────────
def _has_data(series):
    """True if series has at least some non-zero, non-NaN, non-constant values."""
    if series is None:
        return False
    s = series.dropna() if isinstance(series, pd.Series) else pd.Series(series).dropna()
    return len(s) > 0 and not (s == 0).all() and s.std() > 1e-10


def _fig_to_bytes(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    buf.seek(0)
    plt.close(fig)
    return buf.getvalue()


def _placeholder_plot(title, message):
    """Grey placeholder image used when real data is unavailable."""
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.set_facecolor('#f5f5f5')
    fig.patch.set_facecolor('#f5f5f5')
    ax.text(0.5, 0.58, title, transform=ax.transAxes,
            ha='center', va='center', fontsize=15, fontweight='bold',
            color='#888888')
    ax.text(0.5, 0.42, message, transform=ax.transAxes,
            ha='center', va='center', fontsize=10, color='#aaaaaa',
            style='italic', wrap=True)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor('#cccccc')
    fig.tight_layout()
    return _fig_to_bytes(fig)


# ── Data loaders ──────────────────────────────────────────────────────────────
def load_run(run_folder):
    run_metrics = py_interop_run_metrics.run_metrics()
    run_metrics.read(str(run_folder))
    imaging_df = pd.DataFrame(ic.imaging(str(run_folder)))
    for col in ['Lane', 'Tile', 'Read', 'Cycle', 'Cycle Within Read']:
        if col in imaging_df.columns:
            imaging_df[col] = imaging_df[col].astype(int)
    return run_metrics, imaging_df


def get_tile_df(run_metrics):
    tms = run_metrics.tile_metric_set()
    rows = []
    for i in range(tms.size()):
        m = tms.at(i)
        rows.append({
            'Lane': m.lane(),
            'Tile': m.tile(),
            'Cluster Count':    m.cluster_count(),
            'Cluster Count PF': m.cluster_count_pf(),
            'Percent PF':       m.percent_pf(),
        })
    return pd.DataFrame(rows) if rows else None


# ── Plot 1: % Clusters Passing Filter by Lane ────────────────────────────────
def plot_pf_by_lane(tile_df):
    if tile_df is None or len(tile_df) == 0:
        return _placeholder_plot(
            'Plot 1 — % Clusters Passing Filter by Lane',
            'Tile metrics not available for this run.')

    lg = tile_df.groupby('Lane').agg(
        {'Cluster Count': 'sum', 'Cluster Count PF': 'sum'}).reset_index()
    lg = lg[lg['Cluster Count'] > 0]

    if len(lg) == 0:
        return _placeholder_plot(
            'Plot 1 — % Clusters Passing Filter by Lane',
            'All lanes reported zero clusters.')

    pf_pct = (lg['Cluster Count PF'] / lg['Cluster Count'] * 100).tolist()
    lanes  = lg['Lane'].values.astype(int)
    x      = np.arange(len(lanes))

    bar_colors = ['#2ecc71' if p >= 80 else '#f39c12' if p >= 70 else '#e74c3c'
                  for p in pf_pct]

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.bar(x, pf_pct, color=bar_colors, alpha=0.85,
           edgecolor='black', linewidth=1.2, width=0.6)
    ax.axhline(80, color='#27ae60', linestyle='--', linewidth=1.8,
               alpha=0.8, label='Excellent  >= 80%')
    ax.axhline(70, color='#e67e22', linestyle='--', linewidth=1.8,
               alpha=0.8, label='Acceptable >= 70%')
    for i, v in enumerate(pf_pct):
        ax.text(i, min(v + 1.5, 104), f'{v:.1f}%',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f'Lane {l}' for l in lanes], fontsize=10)
    ax.set_ylabel('% Clusters Passing Filter', fontsize=11)
    ax.set_ylim(0, 110)
    ax.set_title('Plot 1 - % Clusters Passing Filter by Lane',
                 fontweight='bold', pad=10)
    ax.legend(loc='lower right', framealpha=0.95)
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    ax.set_axisbelow(True)
    fig.tight_layout()
    return _fig_to_bytes(fig)


# ── Plot 2: % Bases >= Q30 by Cycle ──────────────────────────────────────────
def plot_qscore_by_cycle(imaging_df):
    # Try several column-name variants produced by different interop versions
    metric = None
    for candidate in ['%>= Q30', '% >= Q30', '%>=Q30', 'Pct Q30', '% Q30',
                      '%>=Q30 (R1)', 'Q30']:
        if candidate in imaging_df.columns:
            metric = candidate
            break

    if metric is None:
        return _placeholder_plot(
            'Plot 2 - % Bases >= Q30 by Cycle',
            'Q30 metric column not found in imaging table.\n'
            'This can happen on older run configurations.')

    if not _has_data(imaging_df[metric]):
        return _placeholder_plot(
            'Plot 2 - % Bases >= Q30 by Cycle',
            'Q30 column found but all values are zero or missing.')

    fig, ax = plt.subplots(figsize=(10, 5.5))
    plotted = False
    for lane in sorted(imaging_df['Lane'].unique()):
        cm = (imaging_df[imaging_df['Lane'] == lane]
              .groupby('Cycle')[metric].mean())
        if not _has_data(cm):
            continue
        ax.plot(cm.index, cm.values,
                label=f'Lane {lane}', linewidth=2, alpha=0.85,
                color=LANE_COLORS.get(lane, f'C{lane-1}'),
                marker='o', markersize=3)
        plotted = True

    if not plotted:
        plt.close(fig)
        return _placeholder_plot(
            'Plot 2 - % Bases >= Q30 by Cycle',
            'No valid per-lane Q30 data could be extracted.')

    ax.axhline(80, color='#27ae60', linestyle='--', linewidth=1.8,
               alpha=0.7, label='Target >= 80%')
    ax.set_xlabel('Cycle', fontsize=11)
    ax.set_ylabel('% Bases >= Q30', fontsize=11)
    ax.set_ylim(0, 105)
    ax.set_title('Plot 2 - % Bases >= Q30 by Cycle', fontweight='bold', pad=10)
    ax.legend(loc='lower left', framealpha=0.95, ncol=4)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    fig.tight_layout()
    return _fig_to_bytes(fig)


# ── Plot 3: Base Composition by Cycle ────────────────────────────────────────
def plot_base_composition(imaging_df):
    base_cols = [c for c in imaging_df.columns
                 if any(pat in c for pat in ['% Base', '%Base', '% base'])]

    if not base_cols:
        return _placeholder_plot(
            'Plot 3 - Base Composition by Cycle',
            'Base-composition columns not found in imaging table.\n'
            'This metric may not be available for this instrument / run type.')

    cycle_data = imaging_df.groupby('Cycle')[base_cols].mean()

    fig, ax = plt.subplots(figsize=(10, 5.5))
    plotted = False
    for col in base_cols:
        if not _has_data(cycle_data[col]):
            continue
        # Robustly extract the single-letter base name
        raw = col.split('/')[-1] if '/' in col else col
        raw = (raw.replace('% Base', '').replace('%Base', '')
                  .replace('% base', '').strip())
        base = raw[-1] if raw and raw[-1] in 'ACGT' else raw
        ax.plot(cycle_data.index, cycle_data[col],
                label=f'Base {base}', linewidth=2, alpha=0.85,
                color=CHANNEL_COLORS.get(base),
                marker='o', markersize=3)
        plotted = True

    if not plotted:
        plt.close(fig)
        return _placeholder_plot(
            'Plot 3 - Base Composition by Cycle',
            'All base-composition values are zero or missing.')

    ax.set_xlabel('Cycle', fontsize=11)
    ax.set_ylabel('% of Bases Called', fontsize=11)
    ax.set_ylim(0, 100)
    ax.set_title('Plot 3 - Base Composition by Cycle', fontweight='bold', pad=10)
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    fig.tight_layout()
    return _fig_to_bytes(fig)


# ── PDF builder ───────────────────────────────────────────────────────────────
def build_pdf(run_folder, run_metrics, tile_df, plots, output_pdf):
    """Assemble a 4-page landscape PDF."""
    doc = SimpleDocTemplate(
        output_pdf, pagesize=landscape(letter),
        rightMargin=0.55*inch, leftMargin=0.55*inch,
        topMargin=0.45*inch, bottomMargin=0.45*inch)

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle('T', parent=styles['Heading1'],
        fontSize=22, textColor=colors.HexColor('#0072B2'),
        spaceAfter=6, alignment=TA_CENTER, fontName='Helvetica-Bold')

    subtitle_style = ParagraphStyle('ST', parent=styles['Normal'],
        fontSize=11, textColor=colors.HexColor('#555555'),
        spaceAfter=10, alignment=TA_CENTER, fontName='Helvetica')

    head_style = ParagraphStyle('H', parent=styles['Heading2'],
        fontSize=13, textColor=colors.HexColor('#0072B2'),
        spaceAfter=3, spaceBefore=3, fontName='Helvetica-Bold')

    section_style = ParagraphStyle('SEC', parent=styles['Heading3'],
        fontSize=10, textColor=colors.HexColor('#0072B2'),
        spaceAfter=3, spaceBefore=8, fontName='Helvetica-Bold')

    note_style = ParagraphStyle('N', parent=styles['Normal'],
        fontSize=9.5, textColor=colors.HexColor('#333333'),
        spaceAfter=6, alignment=TA_LEFT, fontName='Helvetica', leading=13)

    guide_style = ParagraphStyle('G', parent=styles['Normal'],
        fontSize=8.8, textColor=colors.HexColor('#1a1a1a'),
        spaceAfter=2, alignment=TA_LEFT, fontName='Helvetica', leading=12,
        leftIndent=4)

    link_style = ParagraphStyle('L', parent=styles['Normal'],
        fontSize=8, textColor=colors.HexColor('#0057a8'),
        spaceAfter=2, fontName='Helvetica', leading=11)

    footer_style = ParagraphStyle('F', parent=styles['Normal'],
        fontSize=7.5, textColor=colors.HexColor('#999999'),
        alignment=TA_CENTER, fontName='Helvetica')

    story = []

    # ══════════════════════════════════════════════════════════════════════════
    # PAGE 1 — Title + run summary table + reader guide + links
    # ══════════════════════════════════════════════════════════════════════════
    story.append(Spacer(1, 0.2*inch))
    story.append(Paragraph('Illumina Sequencing Run', title_style))
    story.append(Paragraph('Quality Control Report', title_style))
    story.append(HRFlowable(width='100%', thickness=1.5,
                             color=colors.HexColor('#0072B2'), spaceAfter=5))
    story.append(Paragraph(
        f'Run: <b>{Path(run_folder).name}</b> &nbsp;·&nbsp; '
        f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}',
        subtitle_style))

    # Run summary table ──────────────────────────────────────────────────────
    info = [['Parameter', 'Value']]
    try:
        ri = run_metrics.run_info()
        info.append(['Flowcell Barcode', ri.flowcell().barcode()])
        info.append(['Total Cycles',     str(ri.total_cycles())])
        num_reads = ri.reads().size()
        for i in range(num_reads):
            r = ri.reads()[i]
            label = f'Read {i+1}' + (' (Index)' if r.is_index() else '')
            info.append([label, f'{r.total_cycles()} cycles'])
    except Exception:
        pass

    if tile_df is not None and len(tile_df):
        total  = tile_df['Cluster Count'].sum()
        pf     = tile_df['Cluster Count PF'].sum()
        pf_pct = pf / total * 100 if total else 0
        info.append(['Total Clusters',               f'{total:,.0f}'])
        info.append(['Clusters Passing Filter (PF)', f'{pf:,.0f}'])
        info.append(['% Clusters Passing Filter',    f'{pf_pct:.1f}%'])
        qc_text  = ('EXCELLENT  (> 80% PF)' if pf_pct >= 80
                    else 'GOOD  (70-80% PF)' if pf_pct >= 70
                    else 'WARNING  (< 70% PF - please investigate)')
        qc_color = ('#27ae60' if pf_pct >= 80
                    else '#e67e22' if pf_pct >= 70 else '#c0392b')
        info.append(['Overall Quality',
                     Paragraph(f'<font color="{qc_color}"><b>{qc_text}</b></font>',
                               guide_style)])

    summary_tbl = Table(info, colWidths=[2.5*inch, 3.8*inch])
    summary_tbl.setStyle(TableStyle([
        ('BACKGROUND',     (0,0), (-1,0),  colors.HexColor('#0072B2')),
        ('TEXTCOLOR',      (0,0), (-1,0),  colors.white),
        ('FONTNAME',       (0,0), (-1,0),  'Helvetica-Bold'),
        ('FONTSIZE',       (0,0), (-1,-1), 9),
        ('BACKGROUND',     (0,1), (0,-1),  colors.HexColor('#deeaf5')),
        ('FONTNAME',       (0,1), (0,-1),  'Helvetica-Bold'),
        ('FONTNAME',       (1,1), (1,-1),  'Courier'),
        ('GRID',           (0,0), (-1,-1), 0.7, colors.HexColor('#aaaaaa')),
        ('ROWBACKGROUNDS', (0,1), (-1,-1),
            [colors.white, colors.HexColor('#f4f9fd')]),
        ('BOTTOMPADDING',  (0,0), (-1,-1), 5),
        ('TOPPADDING',     (0,0), (-1,-1), 5),
        ('VALIGN',         (0,0), (-1,-1), 'MIDDLE'),
        ('ALIGN',          (0,0), (-1,-1), 'LEFT'),
    ]))

    # Reader guide (right column) ────────────────────────────────────────────
    guide_paragraphs = [
        Paragraph('<b>How to read this report</b>', section_style),
        Paragraph(
            '<b>Plot 1 - % Clusters Passing Filter:</b> Fraction of raw '
            'clusters passing Illumina\'s purity filter. '
            'Green = excellent (>= 80%), amber = acceptable (70-80%), '
            'red = investigate (< 70%). Low values may indicate '
            'over/under-clustering or surface contamination.',
            guide_style),
        Spacer(1, 3),
        Paragraph(
            '<b>Plot 2 - % Bases >= Q30:</b> Q30 means a 1-in-1,000 '
            'chance of an incorrect base call. Most cycles should stay '
            'above 80%. A gradual late-cycle decline is normal; a sudden '
            'mid-run drop is not.',
            guide_style),
        Spacer(1, 3),
        Paragraph(
            '<b>Plot 3 - Base Composition:</b> The four bases should run '
            'roughly parallel after the first few cycles. Crossing lines '
            'suggest library bias or reagent issues. Index-read cycles '
            'are intentionally imbalanced.',
            guide_style),
        Spacer(1, 3),
        Paragraph(
            'Questions? Contact your sequencing facility or consult the '
            'references below.',
            guide_style),
        Spacer(1, 6),
        Paragraph('<b>Useful references and tools</b>', section_style),
    ] + [
        Paragraph(
            f'<link href="{url}"><font color="#0057a8">{label}</font></link>',
            link_style)
        for label, url in REFERENCE_LINKS
    ]

    # Side-by-side layout: summary table | guide + links ────────────────────
    page1_table = Table(
        [[summary_tbl, guide_paragraphs]],
        colWidths=[6.55*inch, 4.1*inch],
        style=TableStyle([
            ('VALIGN',       (0,0), (-1,-1), 'TOP'),
            ('LEFTPADDING',  (1,0), (1,0),   12),
            ('RIGHTPADDING', (0,0), (0,0),   8),
            ('TOPPADDING',   (0,0), (-1,-1), 0),
            ('BOTTOMPADDING',(0,0), (-1,-1), 0),
        ])
    )
    story.append(page1_table)
    story.append(Spacer(1, 0.12*inch))
    story.append(HRFlowable(width='100%', thickness=0.5,
                             color=colors.HexColor('#cccccc'), spaceAfter=4))
    story.append(Paragraph(
        'This report is generated automatically from Illumina InterOp binary metrics and is '
        'intended as a first-pass QC summary. Always verify important findings with Illumina SAV '
        'or your facility\'s primary analysis pipeline.',
        footer_style))
    story.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # PAGES 2-4 — One numbered plot per page with explanation
    # ══════════════════════════════════════════════════════════════════════════
    plot_meta = [
        (
            'Plot 1 - % Clusters Passing Filter by Lane',
            'Each bar shows the fraction of raw clusters in that lane that passed '
            "Illumina's chastity (purity) filter. Only passing-filter (PF) clusters are used "
            'in downstream base calling and FASTQ generation. '
            '<b>Green bars (>= 80%)</b> indicate excellent performance; '
            '<b>amber (70-80%)</b> is generally acceptable; '
            '<b>red (< 70%)</b> should be investigated — common causes include '
            'over- or under-clustering, surface contamination, or fluidics issues. '
            'Strong lane-to-lane differences may indicate a hardware or sample-loading problem.',
        ),
        (
            'Plot 2 - % Bases >= Q30 by Cycle',
            'Phred quality scores measure base-call confidence: Q30 = 99.9% accuracy '
            '(1 error per 1,000 bases), Q20 = 99% (1 per 100). '
            'This plot tracks the fraction of bases meeting Q30 across every cycle, per lane. '
            '<b>A gradual decline in later cycles is normal</b> as signal-to-noise decreases; '
            '<b>a sharp mid-run drop</b> may indicate reagent exhaustion, focus drift, or a bubble event. '
            'Most variant callers and assemblers expect >= 80% of bases above Q30 for reliable results. '
            'Index-read cycles are typically excluded from this calculation.',
        ),
        (
            'Plot 3 - Base Composition by Cycle',
            'In a diverse, well-prepared library all four bases (A, C, G, T) should be represented '
            'roughly equally — four nearly parallel lines near 25% each. '
            '<b>The first 1-3 cycles</b> commonly show imbalance due to priming artefacts, which is expected. '
            '<b>Diverging lines mid-read</b> can indicate low-diversity libraries (e.g. amplicon panels), '
            'adapter contamination, or index-hopping. '
            '<b>Index-read cycles</b> are intentionally non-uniform because each cluster carries a specific '
            'barcode sequence. Contact the sequencing facility if you see persistent imbalance in '
            'data-generating reads.',
        ),
    ]

    for img_bytes, (heading, explanation) in zip(plots, plot_meta):
        # img_bytes is always set (real data or grey placeholder)
        story.append(Paragraph(heading, head_style))
        story.append(HRFlowable(width='100%', thickness=0.8,
                                 color=colors.HexColor('#0072B2'), spaceAfter=4))
        story.append(Paragraph(explanation, note_style))
        story.append(Spacer(1, 0.04*inch))
        story.append(Image(BytesIO(img_bytes),
                           width=9.5*inch, height=4.8*inch,
                           kind='proportional'))
        story.append(PageBreak())

    doc.build(story)
    print(f"[OK] PDF saved: {output_pdf}")


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 3:
        print("Usage: python interop_more.py <run_folder> <output_pdf>")
        sys.exit(1)

    run_folder, output_pdf = sys.argv[1], sys.argv[2]

    if not os.path.exists(run_folder):
        print(f"[ERROR] Run folder not found: {run_folder}")
        sys.exit(1)
    if not os.path.exists(os.path.join(run_folder, 'InterOp')):
        print(f"[ERROR] InterOp folder missing in: {run_folder}")
        sys.exit(1)

    out_dir = os.path.dirname(output_pdf)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"[INFO] Loading run: {run_folder}")
    run_metrics, imaging_df = load_run(run_folder)
    tile_df = get_tile_df(run_metrics)

    print("[INFO] Generating plots...")
    plots = [
        plot_pf_by_lane(tile_df),
        plot_qscore_by_cycle(imaging_df),
        plot_base_composition(imaging_df),
    ]
    n_real = sum(1 for _ in plots)   # all 3 always produce an image
    print(f"[INFO] {n_real}/3 plot images generated")

    print("[INFO] Building PDF...")
    build_pdf(run_folder, run_metrics, tile_df, plots, output_pdf)


if __name__ == '__main__':
    main()
    # interop_more - cluster_file is output
    # script = p + "/scripts/interop_more.py"