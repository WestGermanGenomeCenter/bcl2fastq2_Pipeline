"""
RNA-seq QC Interactive HTML Report
Usage: python rnaseq_html_report.py <featureCounts_file> <output.html>

Sections (logical flow):
  1. Summary          — per-sample stats table + quick visual overview
  2. QC Metrics       — library sizes, gene detection, CPM distributions, outlier Z-scores
  3. Sample Structure — correlation heatmap, Ward dendrogram, PCA, UMAP
  4. Gene Analysis    — mean-CV plot, sortable top-genes table, variable-gene heatmap
  5. Group Estimation — k-means (UMAP), hierarchical cut, consensus clustering
"""

import os, sys, json, time, datetime, warnings, traceback
from contextlib import contextmanager

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.stats import gaussian_kde, median_abs_deviation
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

warnings.filterwarnings("ignore")

# ── 30 visually distinct sample colours (consistent across all plots) ─────────
PALETTE = [
    "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
    "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
    "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5",
    "#c49c94","#f7b6d2","#c7c7c7","#dbdb8d","#9edae5",
    "#393b79","#637939","#8c6d31","#843c39","#7b4173",
    "#5254a3","#8ca252","#bd9e39","#ad494a","#a55194",
]

# matplotlib colour codes returned by scipy dendrogram → hex
_MPL = {"b":"#1f77b4","g":"#2ca02c","r":"#d62728","c":"#17becf",
        "m":"#9467bd","y":"#bcbd22","k":"#333333","w":"#ffffff"}
_MPL.update({f"C{i}":c for i,c in enumerate(
    ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
     "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"])})

def safe_color(c):
    return _MPL.get(c, c) if not (c.startswith("#") or c.startswith("rgb")) else c

def make_colors(names):
    """Consistent sample → colour mapping (cycles after 30)."""
    return {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(names)}

# ── Shared Plotly layout ───────────────────────────────────────────────────────
BASE = dict(
    font=dict(family="Arial, sans-serif", size=12),
    paper_bgcolor="white",
    plot_bgcolor="#f9f9f9",
    margin=dict(l=70, r=40, t=65, b=75),
    hoverlabel=dict(bgcolor="white", font_size=11),
)

def base(fig, title="", height=480):
    fig.update_layout(**BASE, title=dict(text=title, font_size=14), height=height)
    return fig

def to_json(fig):
    return pio.to_json(fig, validate=False)

def hline_mean_sd(fig, vals, row=None, col=None):
    m, s = float(np.mean(vals)), float(np.std(vals))
    kw = {"row": row, "col": col} if row else {}
    for v, dash, lbl in [(m,"solid",f"Mean {m:.2f}"),
                          (m+s,"dash",f"+1 SD {m+s:.2f}"),
                          (m-s,"dash",f"−1 SD {m-s:.2f}")]:
        fig.add_hline(y=v, line_dash=dash, line_color="#666", line_width=1.2,
                      annotation_text=lbl, annotation_font_size=9,
                      annotation_position="top right", **kw)

# ── Logging ────────────────────────────────────────────────────────────────────
@contextmanager
def step(name):
    t = time.time()
    print(f"  {name} ...", end=" ", flush=True)
    yield
    print(f"({time.time()-t:.1f}s)")

def warn(msg):
    print(f"\n  [WARNING] {msg}")

# ── Data loading ───────────────────────────────────────────────────────────────
def load_counts(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    raw = df.iloc[:, 6:].to_numpy(dtype=float)    # genes × samples
    genes = df["Geneid"].astype(str).values
    names = []
    for c in df.columns[6:]:
        n = os.path.basename(str(c))
        for suf in ("_Aligned.sortedByCoord.out.bam", ".bam",
                    "_Aligned.out.bam", ".sortedByCoord.out.bam"):
            n = n.replace(suf, "")
        names.append(n)
    return raw, genes, names

def normalise(raw):
    lib = raw.sum(axis=0)
    lib = np.where(lib == 0, 1, lib)
    cpm = raw / lib[np.newaxis, :] * 1e6
    keep = (cpm > 1).any(axis=1)
    if keep.sum() == 0:
        keep = np.ones(raw.shape[0], dtype=bool)
    cpm_f = cpm[keep]
    lcpm = np.log2(cpm_f + 1).T        # samples × genes
    return lib, cpm, cpm_f, lcpm, keep

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE BUILDERS
# ══════════════════════════════════════════════════════════════════════════════

# ── Section 1: Summary ────────────────────────────────────────────────────────

def fig_overview(names, lib, n_before, n_after, colors):
    """Two-panel overview: library sizes (left) + raw vs filtered gene counts (right)."""
    n   = len(names)
    tfs = max(7, 10 - n // 12)
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Library sizes (total mapped reads)",
                        "Genes detected: raw (faded) vs CPM-filtered (solid)"],
        horizontal_spacing=0.08,
    )
    # Left — library sizes + mean line
    fig.add_trace(go.Bar(
        x=names, y=(lib / 1e6).tolist(),
        marker_color=[colors[s] for s in names], opacity=0.85,
        name="Library size", showlegend=False,
        hovertemplate="%{x}<br><b>%{y:.2f} M reads</b><extra></extra>",
    ), row=1, col=1)
    m = float(np.mean(lib / 1e6))
    fig.add_hline(y=m, line_dash="dash", line_color="#666", line_width=1.2,
                  annotation_text=f"Mean {m:.1f} M", annotation_font_size=9,
                  annotation_position="top right", row=1, col=1)

    # Right — raw genes (faded) + filtered genes (solid), same colour per sample
    fig.add_trace(go.Bar(
        x=names, y=list(n_before),
        marker_color=[colors[s] for s in names], opacity=0.30,
        name="Genes detected (raw)", offsetgroup="a",
        hovertemplate="%{x}<br>Raw: <b>%{y:,} genes</b><extra></extra>",
    ), row=1, col=2)
    fig.add_trace(go.Bar(
        x=names, y=list(n_after),
        marker_color=[colors[s] for s in names], opacity=0.90,
        name="Genes detected (CPM-filtered)", offsetgroup="b",
        hovertemplate="%{x}<br>Filtered: <b>%{y:,} genes</b><extra></extra>",
    ), row=1, col=2)

    base(fig, "Dataset Overview — library sizes and gene detection (faded = raw, solid = CPM-filtered)",
         max(420, n * 5 + 170))
    fig.update_layout(
        barmode="group",
        legend=dict(x=0.52, y=1.10, orientation="h", font_size=10),
    )
    for col in [1, 2]:
        fig.update_xaxes(tickangle=-40, tickfont_size=tfs, row=1, col=col)
    fig.update_yaxes(title_text="Millions of reads", row=1, col=1)
    fig.update_yaxes(title_text="Number of genes", row=1, col=2)
    return fig

# ── Section 2: QC Metrics ─────────────────────────────────────────────────────

def fig_library_sizes(names, lib, colors):
    n = len(names)
    fig = go.Figure(go.Bar(
        x=names, y=(lib / 1e6).tolist(),
        marker_color=[colors[s] for s in names],
        hovertemplate="%{x}<br><b>%{y:.2f} M reads</b><extra></extra>",
    ))
    hline_mean_sd(fig, lib / 1e6)
    base(fig, "Library Sizes — total mapped reads per sample (all raw counts)",
         max(420, n * 5 + 150))
    fig.update_layout(showlegend=False,
                      xaxis=dict(title="Sample", tickangle=-40,
                                 tickfont_size=max(7, 10 - n // 12)),
                      yaxis_title="Millions of reads")
    return fig

def fig_detected_raw(names, n_before, colors):
    n = len(names)
    fig = go.Figure(go.Bar(
        x=names, y=list(n_before),
        marker_color=[colors[s] for s in names],
        hovertemplate="%{x}<br><b>%{y:,} genes</b><extra></extra>",
    ))
    hline_mean_sd(fig, np.array(n_before, dtype=float))
    base(fig, "Genes Detected per Sample — raw counts (count > 0)",
         max(420, n * 5 + 150))
    fig.update_layout(showlegend=False,
                      xaxis=dict(title="Sample", tickangle=-40,
                                 tickfont_size=max(7, 10 - n // 12)),
                      yaxis_title="Number of genes (raw)")
    return fig

def fig_detected_filt(names, n_after, colors):
    n = len(names)
    fig = go.Figure(go.Bar(
        x=names, y=list(n_after),
        marker_color=[colors[s] for s in names],
        hovertemplate="%{x}<br><b>%{y:,} genes</b><extra></extra>",
    ))
    hline_mean_sd(fig, np.array(n_after, dtype=float))
    base(fig, "Genes Detected per Sample — after CPM > 1 filter",
         max(420, n * 5 + 150))
    fig.update_layout(showlegend=False,
                      xaxis=dict(title="Sample", tickangle=-40,
                                 tickfont_size=max(7, 10 - n // 12)),
                      yaxis_title="Number of genes (CPM-filtered)")
    return fig

def fig_cpm_boxplot(names, lcpm, colors):
    n = len(names)
    fig = go.Figure()
    for i, s in enumerate(names):
        fig.add_trace(go.Box(
            y=lcpm[i].tolist(), name=s,
            marker_color=colors[s], line_color=colors[s],
            boxpoints=False,
            hovertemplate=f"<b>{s}</b><br>%{{y:.2f}}<extra></extra>",
        ))
    base(fig, "CPM Distributions — log₂(CPM+1) per sample (CPM-filtered genes)",
         max(460, n * 6 + 150))
    fig.update_layout(showlegend=False,
                      xaxis=dict(title="Sample", tickangle=-40,
                                 tickfont_size=max(7, 10 - n // 12)),
                      yaxis_title="log₂(CPM+1)")
    return fig

def fig_cpm_density(names, lcpm, colors):
    fig = go.Figure()
    for i, s in enumerate(names):
        v = lcpm[i]; v = v[np.isfinite(v)]
        try:
            x = np.linspace(float(v.min()), float(v.max()), 300)
            y = gaussian_kde(v, bw_method=0.3)(x)
            fig.add_trace(go.Scatter(
                x=x.tolist(), y=y.tolist(), mode="lines", name=s,
                line=dict(color=colors[s], width=1.5),
                hovertemplate=f"<b>{s}</b>: %{{y:.4f}}<extra></extra>",
            ))
        except Exception:
            pass
    base(fig, "CPM Density Curves — log₂(CPM+1); all samples should broadly overlap")
    fig.update_layout(xaxis_title="log₂(CPM+1)", yaxis_title="Density",
                      legend=dict(font_size=9, itemclick="toggleothers"))
    return fig

def fig_outlier_heatmap(names, lcpm, lib, n_after, n_filt):
    n = len(names)
    # ── Eight QC metrics ──────────────────────────────────────────────────────
    cor_mat  = np.corrcoef(lcpm)
    np.fill_diagonal(cor_mat, np.nan)          # exclude self-correlation
    metrics = {
        "Library size\n(raw)":               lib.astype(float),
        "Genes detected\n(CPM-filtered)":    np.array(n_after, dtype=float),
        "Median log₂ CPM\n(CPM-filtered)":   np.median(lcpm, axis=1),
        "IQR log₂ CPM\n(CPM-filtered)":      np.percentile(lcpm, 75, axis=1)
                                              - np.percentile(lcpm, 25, axis=1),
        "CV (std/mean)\n(CPM-filtered)":     np.array([lcpm[i].std() / (lcpm[i].mean() + 1e-9)
                                                        for i in range(n)]),
        "% zeros\n(CPM-filtered)":           np.array([(lcpm[i] == 0).mean() * 100
                                                        for i in range(n)]),
        "Mean sample r\n(CPM-filtered)":     np.nanmean(cor_mat, axis=1),
        "Min sample r\n(CPM-filtered)":      np.nanmin(cor_mat, axis=1),
    }
    Z = np.zeros((n, len(metrics)))
    for j, v in enumerate(metrics.values()):
        mad        = median_abs_deviation(v) + 1e-9
        Z[:, j]   = (v - np.median(v)) / mad

    ann    = [[f"{Z[i,j]:+.1f}" for j in range(len(metrics))] for i in range(n)]
    shapes = []
    for i in range(n):
        if np.any(np.abs(Z[i]) >= 2):
            shapes.append(dict(type="rect", x0=-0.5, x1=len(metrics) - 0.5,
                               y0=i - 0.5, y1=i + 0.5,
                               line=dict(color="red", width=2),
                               fillcolor="rgba(0,0,0,0)", layer="above"))
    tfs = max(7, 11 - n // 8)
    fig = go.Figure(go.Heatmap(
        z=Z.tolist(), x=list(metrics.keys()), y=list(names),
        text=ann, texttemplate="%{text}", textfont_size=max(7, tfs - 1),
        colorscale="RdBu", reversescale=True, zmid=0,
        colorbar=dict(title="Robust Z", thickness=14),
        hovertemplate="<b>%{y}</b> — %{x}<br>Z = %{z:.3f}<extra></extra>",
    ))
    base(fig, "Multi-Metric Outlier Overview — 8 QC metrics, robust Z-score per sample "
              "(red border = |Z| ≥ 2 on any metric; hover for exact values)",
         max(420, n * 28 + 140))
    fig.update_layout(
        shapes=shapes,
        xaxis=dict(tickangle=-30, tickfont_size=9, automargin=True),
        yaxis=dict(tickfont_size=tfs, autorange="reversed", automargin=True),
        margin=dict(l=90, r=40, t=70, b=120),
    )
    return fig

# ── Section 3: Sample Structure ───────────────────────────────────────────────

def fig_correlation(names, lcpm, n_filt, group_labels=None):
    """
    Two-panel correlation heatmap:
      Left  — samples ordered by hierarchical clustering
      Right — samples ordered by suspected group (if available), else by mean r descending
    """
    n   = len(names)
    cor = np.corrcoef(lcpm)
    tfs = max(6, 10 - n // 8)

    # Order 1: hierarchical clustering
    try:
        d      = np.clip(squareform(1 - cor, checks=False), 0, None)
        order1 = dendrogram(linkage(d, method="average"), no_plot=True)["leaves"]
    except Exception:
        order1 = list(range(n))

    # Order 2: by group label (if provided), else by mean correlation descending
    if group_labels is not None and len(group_labels) == n:
        order2 = list(np.argsort(group_labels))
        order2_title = "Ordered by estimated group"
    else:
        mean_r = cor.mean(axis=1)
        order2 = list(np.argsort(mean_r)[::-1])
        order2_title = "Ordered by mean inter-sample r (high→low)"

    ann_thresh = 40   # only annotate cells when n is small
    vmin = max(float(cor.min()) - 0.01, 0.5)

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[f"Ordered by hierarchical clustering",
                        order2_title],
        horizontal_spacing=0.08,
    )

    for col_idx, order in enumerate([order1, order2], start=1):
        co  = cor[np.ix_(order, order)]
        no  = [names[i] for i in order]
        ann = [[f"{co[i,j]:.2f}" for j in range(n)] for i in range(n)] if n <= ann_thresh else None
        fig.add_trace(go.Heatmap(
            z=co.tolist(), x=no, y=no,
            text=ann, texttemplate="%{text}" if ann else None,
            textfont_size=max(5, tfs - 2),
            colorscale="RdBu", reversescale=True,
            zmin=vmin, zmax=1.0,
            colorbar=dict(title="Pearson r", thickness=12,
                          x=0.46 if col_idx == 1 else 1.01),
            hovertemplate="<b>%{x}</b> vs <b>%{y}</b><br>r = %{z:.4f}<extra></extra>",
            showscale=True,
        ), row=1, col=col_idx)
        fig.update_xaxes(tickangle=-45, tickfont_size=tfs, row=1, col=col_idx,
                         automargin=True)
        fig.update_yaxes(tickfont_size=tfs, autorange="reversed", row=1, col=col_idx,
                         automargin=True)

    base(fig, f"Sample-to-Sample Pearson Correlation — log₂(CPM+1), "
              f"CPM-filtered genes (n={n_filt})",
         max(520, n * 22 + 130))
    fig.update_layout(margin=dict(l=80, r=80, t=100, b=80))
    return fig

def fig_ward_dendrogram(names, lcpm, n_filt):
    n = len(names)
    try:
        link = linkage(pdist(lcpm, metric="euclidean"), method="ward")
        dend = dendrogram(link, labels=list(names), no_plot=True,
                          color_threshold=0.7 * float(link[:, 2].max()))
    except Exception as e:
        warn(f"Dendrogram: {e}"); return None
    fig = go.Figure()
    for xs, ys, col in zip(dend["icoord"], dend["dcoord"], dend["color_list"]):
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines",
                                 line=dict(color=safe_color(col), width=1.5),
                                 hoverinfo="skip", showlegend=False))
    lo = dend["leaves"]
    tfs = max(7, 11 - n // 10)
    base(fig, f"Hierarchical Sample Clustering — Euclidean distance, Ward linkage "
              f"(CPM-filtered genes n={n_filt}; longer branches = more dissimilar)", 450)
    fig.update_layout(
        xaxis=dict(tickvals=[5 + 10 * i for i in range(n)],
                   ticktext=[names[i] for i in lo],
                   tickangle=-45, tickfont_size=tfs, showgrid=False),
        yaxis_title="Euclidean distance (log₂ CPM+1)",
    )
    return fig

def _pca_coords(lcpm, n_top, n):
    X = lcpm[:, np.argsort(lcpm.var(axis=0))[-n_top:]]
    nc = min(n, X.shape[1], 10)
    if nc < 2:
        return None, None, 0
    p = PCA(n_components=nc)
    C = p.fit_transform(X)
    return C, p.explained_variance_ratio_ * 100, nc

def fig_pca_3d(names, lcpm, n_top, colors):
    """Interactive 3D PCA scatter — drag to rotate, hover for sample name."""
    n = len(names)
    C, pct, nc = _pca_coords(lcpm, n_top, n)
    if C is None or nc < 2:
        return None, C, pct
    # Use 3 components if available, else pad with zeros
    x = C[:, 0].tolist()
    y = C[:, 1].tolist()
    z = (C[:, 2].tolist() if nc >= 3
         else [0.0] * n)
    pct_z = pct[2] if nc >= 3 else 0.0
    fig = go.Figure()
    for i, s in enumerate(names):
        fig.add_trace(go.Scatter3d(
            x=[x[i]], y=[y[i]], z=[z[i]],
            mode="markers+text" if n <= 20 else "markers",
            name=s,
            text=[s] if n <= 20 else None,
            textposition="top center",
            textfont=dict(size=9),
            marker=dict(color=colors[s], size=7,
                        line=dict(color="white", width=0.5)),
            hovertemplate=(f"<b>{s}</b><br>"
                           f"PC1=%{{x:.2f}}<br>PC2=%{{y:.2f}}<br>PC3=%{{z:.2f}}"
                           f"<extra></extra>"),
            showlegend=(n <= 25),
        ))
    # Use tighter margins to avoid label/axis overlap
    fig.update_layout(
        **{k: v for k, v in BASE.items() if k != "margin"},
        margin=dict(l=0, r=0, t=80, b=0),
        title=dict(
            text=(f"PCA 3D — PC1 ({pct[0]:.1f}%) × PC2 ({pct[1]:.1f}%) × "
                  f"PC3 ({pct_z:.1f}%) — top {n_top} variable genes; drag to rotate"),
            font_size=13,
        ),
        height=max(560, n * 7 + 180),
        scene=dict(
            xaxis=dict(title=dict(text=f"PC1 ({pct[0]:.1f}%)", font=dict(size=11)),
                       tickfont=dict(size=9), showbackground=True,
                       backgroundcolor="#f0f0f0"),
            yaxis=dict(title=dict(text=f"PC2 ({pct[1]:.1f}%)", font=dict(size=11)),
                       tickfont=dict(size=9), showbackground=True,
                       backgroundcolor="#f0f0f0"),
            zaxis=dict(title=dict(text=f"PC3 ({pct_z:.1f}%)" if nc >= 3 else "PC3 (n/a)",
                                  font=dict(size=11)),
                       tickfont=dict(size=9), showbackground=True,
                       backgroundcolor="#f0f0f0"),
            bgcolor="white",
        ),
        legend=dict(font_size=9, itemsizing="constant"),
    )
    return fig, C, pct

def fig_scree(lcpm, n_top):
    n = lcpm.shape[0]
    X = lcpm[:, np.argsort(lcpm.var(axis=0))[-n_top:]]
    nc = min(n, X.shape[1], 10)
    if nc < 2:
        return None
    p = PCA(n_components=nc); p.fit(X)
    pct = p.explained_variance_ratio_ * 100
    cum = np.cumsum(pct)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=list(range(1, nc + 1)), y=pct.tolist(), name="Per-PC variance",
                         hovertemplate="PC%{x}<br>%{y:.1f}%<extra></extra>"))
    fig.add_trace(go.Scatter(x=list(range(1, nc + 1)), y=cum.tolist(),
                             mode="lines+markers", name="Cumulative %",
                             hovertemplate="PC%{x}<br>cumulative %{y:.1f}%<extra></extra>"))
    fig.add_hline(y=80, line_dash="dot", line_color="grey", line_width=1,
                  annotation_text="80%", annotation_font_size=9)
    base(fig, f"PCA Scree Plot — variance explained per PC (top {n_top} variable genes)", 390)
    fig.update_layout(xaxis=dict(title="Principal Component",
                                 tickvals=list(range(1, nc + 1))),
                      yaxis_title="% Variance Explained")
    return fig

def fig_pca_loadings(lcpm, genes_f, n_top):
    """Top 15 genes contributing to PC1 and PC2 (loading bar chart)."""
    n = lcpm.shape[0]
    X = lcpm[:, np.argsort(lcpm.var(axis=0))[-n_top:]]
    sel_genes = genes_f[np.argsort(lcpm.var(axis=0))[-n_top:]]
    nc = min(n, X.shape[1], 10)
    if nc < 2:
        return None
    p = PCA(n_components=nc); p.fit(X)
    comp = p.components_          # nc × n_top
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=["Top 15 genes driving PC1",
                                        "Top 15 genes driving PC2"],
                        horizontal_spacing=0.12)
    for pc_idx, col in enumerate([1, 2], start=1):
        loadings = comp[pc_idx - 1]              # n_top values
        top15    = np.argsort(np.abs(loadings))[-15:][::-1]
        gnames   = [sel_genes[i] for i in top15]
        lvals    = [float(loadings[i]) for i in top15]
        colors_bar = ["#d62728" if v > 0 else "#1f77b4" for v in lvals]
        fig.add_trace(go.Bar(
            x=lvals[::-1], y=gnames[::-1], orientation="h",
            marker_color=colors_bar[::-1],
            hovertemplate="<b>%{y}</b><br>Loading = %{x:.4f}<extra></extra>",
            showlegend=False,
        ), row=1, col=col)
        fig.update_xaxes(title_text="Loading", row=1, col=col)
        fig.update_yaxes(tickfont_size=9, row=1, col=col)
    base(fig, f"PCA Loadings — top 15 genes contributing to PC1 and PC2 "
              f"(top {n_top} variable genes; red = positive, blue = negative loading)", 480)
    return fig


def fig_sample_scatter_pairs(names, lcpm, colors):
    """Log₂ CPM scatter for all pairs — only shown for n ≤ 6 samples."""
    n = len(names)
    if n > 6:
        return None
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if not pairs:
        return None
    nc = min(3, len(pairs))
    nr = int(np.ceil(len(pairs) / nc))
    # Compute spacing so it never exceeds Plotly's 1/(rows-1) limit
    v_space = min(0.08, 0.9 / max(nr - 1, 1))
    h_space = min(0.08, 0.9 / max(nc - 1, 1))
    titles  = [f"{names[i]} vs {names[j]}" for i, j in pairs]
    titles += [""] * (nr * nc - len(pairs))
    fig = make_subplots(rows=nr, cols=nc, subplot_titles=titles,
                        horizontal_spacing=h_space, vertical_spacing=v_space)
    for idx, (i, j) in enumerate(pairs):
        r, c = divmod(idx, nc)
        xi = lcpm[i].tolist(); xj = lcpm[j].tolist()
        r_val = float(np.corrcoef(lcpm[i], lcpm[j])[0, 1])
        fig.add_trace(go.Scatter(
            x=xi, y=xj, mode="markers",
            marker=dict(size=3, opacity=0.4, color="#555"),
            hovertemplate=f"{names[i]}=%{{x:.2f}}<br>{names[j]}=%{{y:.2f}}<extra></extra>",
            showlegend=False,
        ), row=r + 1, col=c + 1)
        mn = min(min(xi), min(xj)); mx = max(max(xi), max(xj))
        fig.add_trace(go.Scatter(x=[mn, mx], y=[mn, mx], mode="lines",
                                 line=dict(color="red", width=1, dash="dot"),
                                 showlegend=False), row=r + 1, col=c + 1)
        # r annotation inside the panel (no paper coords to avoid crowding)
        fig.add_annotation(
            text=f"r={r_val:.3f}", x=mx, y=mn,
            xref=f"x{idx+1}" if idx > 0 else "x",
            yref=f"y{idx+1}" if idx > 0 else "y",
            showarrow=False, font_size=9, xanchor="right", yanchor="bottom",
            bgcolor="rgba(255,255,255,0.7)", borderpad=2,
        )
        fig.update_xaxes(title_text="log₂ CPM+1", row=r + 1, col=c + 1)
        fig.update_yaxes(title_text="log₂ CPM+1", row=r + 1, col=c + 1)
    base(fig, f"Pairwise Sample Scatter — log₂(CPM+1), all {len(pairs)} pairs "
              f"(red dashed = perfect agreement; r = Pearson correlation; "
              f"CPM-filtered genes; shown only for n ≤ 6 samples)",
         max(350, nr * 300))
    return fig


def html_silhouette_table(names, lcpm, n_top):
    """HTML table of silhouette scores for k=2..min(6,n-1) using Ward groups."""
    from sklearn.metrics import silhouette_score
    n     = len(names)
    max_k = min(6, n - 1)
    if max_k < 2:
        return ""
    X     = lcpm[:, np.argsort(lcpm.var(axis=0))[-n_top:]]
    rows  = ""
    best  = -1; best_k = 2
    try:
        link = linkage(pdist(X, metric="euclidean"), method="ward")
    except Exception:
        return ""
    for k in range(2, max_k + 1):
        try:
            labs = fcluster(link, k, criterion="maxclust") - 1
            if len(np.unique(labs)) < 2:
                continue
            sil  = float(silhouette_score(X, labs))
            if sil > best:
                best = sil; best_k = k
        except Exception:
            sil = float("nan")
        rows += f"<tr{'style=\"background:#e8f4e8\"' if k==best_k else ''}>" \
                f"<td style='text-align:center'>{k}</td>" \
                f"<td style='text-align:center'>{sil:.3f}</td>" \
                f"<td style='text-align:center'>{'★ best' if k==best_k else ''}</td></tr>"
    if not rows:
        return ""
    return f"""<div style="margin-top:10px">
<p style="font-size:11px;color:#555;margin-bottom:6px">
  Silhouette score (−1..1) measures cluster cohesion vs separation.
  Higher = better-separated groups. Computed on Ward hierarchical groups
  (top {n_top} variable genes).
</p>
<table style="border-collapse:collapse;font-size:12px;width:280px">
  <thead><tr style="background:#336699;color:white">
    <th style="padding:6px 16px">k</th>
    <th style="padding:6px 16px">Silhouette score</th>
    <th style="padding:6px 16px">Note</th>
  </tr></thead>
  <tbody>{rows}</tbody>
</table></div>"""


def fig_umap_sample(names, lcpm, n_top, colors, run_warnings):
    """UMAP (or PCA fallback) coloured by sample identity."""
    n = len(names)
    X = lcpm[:, np.argsort(lcpm.var(axis=0))[-n_top:]]
    try:
        from umap import UMAP
        emb = UMAP(n_components=2, n_neighbors=max(2, min(15, n - 1)),
                   min_dist=0.3, random_state=42, verbose=False).fit_transform(X)
        elbl = "UMAP"
    except Exception as e:
        run_warnings.append(f"UMAP unavailable ({e}), using PCA.")
        emb = PCA(n_components=2).fit_transform(X)
        elbl = "PC"
    mode = "markers+text" if n <= 30 else "markers"
    fig = go.Figure()
    for i, s in enumerate(names):
        fig.add_trace(go.Scatter(
            x=[float(emb[i, 0])], y=[float(emb[i, 1])], mode=mode, name=s,
            text=[s] if n <= 30 else None, textposition="top center", textfont_size=8,
            marker=dict(color=colors[s], size=9, line=dict(color="white", width=1)),
            hovertemplate=f"<b>{s}</b><br>{elbl}1=%{{x:.2f}}<br>{elbl}2=%{{y:.2f}}<extra></extra>",
            showlegend=(n <= 25),
        ))
    base(fig, f"{elbl} Embedding — coloured by sample "
              f"(top {n_top} variable genes, log₂ CPM+1, CPM-filtered)",
         max(500, n * 7 + 160))
    fig.update_layout(xaxis_title=f"{elbl} 1", yaxis_title=f"{elbl} 2")
    return fig, emb, elbl

def fig_ma(names, lcpm):
    """MA plots vs pseudo-bulk mean. Only for n ≤ 12."""
    n, ref = len(names), lcpm.mean(axis=0)
    nc, nr = min(3, n), int(np.ceil(n / min(3, n)))
    titles = list(names) + [""] * (nr * nc - n)
    fig = make_subplots(rows=nr, cols=nc, subplot_titles=titles,
                        horizontal_spacing=0.06, vertical_spacing=0.12)
    for idx, s in enumerate(names):
        r, c = divmod(idx, nc)
        A = ((lcpm[idx] + ref) / 2).tolist()
        M = (lcpm[idx] - ref)
        mm, ms = float(M.mean()), float(M.std())
        try:
            stk = np.vstack([A, M.tolist()])
            z = gaussian_kde(stk)(stk); o = z.argsort()
            fig.add_trace(go.Scatter(
                x=np.array(A)[o].tolist(), y=M[o].tolist(), mode="markers",
                marker=dict(color=z[o].tolist(), colorscale="Viridis", size=3,
                            opacity=0.6, showscale=False),
                showlegend=False,
                hovertemplate="A=%{x:.2f}<br>M=%{y:.2f}<extra></extra>",
            ), row=r + 1, col=c + 1)
        except Exception:
            fig.add_trace(go.Scatter(x=A, y=M.tolist(), mode="markers",
                                     marker=dict(size=3, opacity=0.4),
                                     showlegend=False), row=r + 1, col=c + 1)
        for val, cl, dash in [(0, "red", "dash"), (mm, "green", "solid"),
                               (mm + ms, "green", "dot"), (mm - ms, "green", "dot")]:
            fig.add_hline(y=val, line_color=cl, line_dash=dash,
                          line_width=1.2, row=r + 1, col=c + 1)
    base(fig, "MA Plots — each sample vs pseudo-bulk mean "
              "(red dashed = M=0, green = mean M ±1 SD; CPM-filtered genes)", nr * 290)
    return fig

# ── Section 4: Gene Analysis ──────────────────────────────────────────────────

def fig_mean_cv(lcpm, genes_f, n_filt):
    G = lcpm.T; mu = G.mean(axis=1); cv = G.std(axis=1) / (mu + 1e-9)
    bins = np.array_split(np.argsort(mu), 30)
    bx = np.array([mu[b].mean() for b in bins])
    by = np.array([np.median(cv[b]) for b in bins])
    exp = mu > 1
    top = (np.where(exp)[0][np.argsort(cv[exp])[-20:]]
           if exp.sum() >= 20 else np.argsort(cv)[-20:])
    try:
        z = gaussian_kde(np.vstack([mu, cv]))(np.vstack([mu, cv]))
    except Exception:
        z = np.ones(len(mu))
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=mu.tolist(), y=cv.tolist(), mode="markers",
        marker=dict(color=z.tolist(), colorscale="Viridis", size=3, opacity=0.5,
                    showscale=True, colorbar=dict(title="Density", thickness=12)),
        text=list(genes_f), showlegend=False,
        hovertemplate="<b>%{text}</b><br>Mean log₂ CPM=%{x:.2f}<br>CV=%{y:.3f}<extra></extra>",
    ))
    fig.add_trace(go.Scatter(x=bx.tolist(), y=by.tolist(), mode="lines",
                             name="Running median CV", line_width=2))
    fig.add_trace(go.Scatter(
        x=mu[top].tolist(), y=cv[top].tolist(), mode="markers+text",
        text=list(genes_f[top]), textposition="top right", textfont_size=8,
        marker=dict(size=7, symbol="diamond", color="red"),
        name="Top 20 high-CV genes",
        hovertemplate="<b>%{text}</b><br>Mean=%{x:.2f}<br>CV=%{y:.3f}<extra></extra>",
    ))
    base(fig, f"Mean Expression vs Coefficient of Variation — "
              f"CPM-filtered genes (n={n_filt}); top 20 high-CV genes labelled (mean > 1)", 530)
    fig.update_layout(xaxis_title="Mean log₂(CPM+1)", yaxis_title="CV (std/mean)")
    return fig


def html_top_genes_table(names, cpm_f, genes_f, n_top=20):
    """
    Sortable, searchable HTML table of top expressed genes — all samples shown.
    Returns a raw HTML string (embedded directly, not a Plotly figure).
    """
    mu    = cpm_f.mean(axis=1)
    idx   = np.argsort(mu)[::-1][:n_top]
    tg    = genes_f[idx]; tm = mu[idx]
    tp    = tm / (mu.sum() + 1e-9) * 100
    tc    = cpm_f[idx, :]

    hdrs = ["#", "Gene ID", "Mean CPM", "% of total"] + list(names)
    th   = "".join(f'<th onclick="sortTbl({i})">{h}</th>'
                   for i, h in enumerate(hdrs))

    rows = ""
    for rank, (g, m, p, row) in enumerate(zip(tg, tm, tp, tc), 1):
        bg = "background:#ffe8e8" if p > 5 else ("background:#f8f8f8" if rank % 2 == 0 else "")
        st = f' style="{bg}"' if bg else ""
        cells = (f'<td data-val="{rank}">{rank}</td>'
                 f'<td style="font-weight:bold">{g}</td>'
                 f'<td data-val="{m:.6f}">{m:.1f}</td>'
                 f'<td data-val="{p:.8f}">{p:.2f}%</td>')
        for v in row:
            cells += f'<td data-val="{v:.4f}">{v:.0f}</td>'
        rows += f'<tr{st}>{cells}</tr>\n'

    return f"""
<div style="margin-bottom:8px;display:flex;align-items:center;gap:12px">
  <input id="gene-search" type="text" placeholder="Search gene or sample..."
    style="padding:5px 10px;border:1px solid #ccc;border-radius:3px;
           width:260px;font-size:12px"
    oninput="searchGeneTable(this.value)">
  <span style="color:#888;font-size:11px">
    Click any column header to sort &nbsp;|&nbsp;
    <span style="background:#ffe8e8;padding:2px 6px;border-radius:3px">red row</span>
    = gene >5% of total mean CPM
  </span>
</div>
<div style="overflow-x:auto;overflow-y:auto;max-height:540px;
            border:1px solid #ddd;border-radius:3px">
  <table id="top-genes-tbl"
    style="border-collapse:collapse;font-size:11px;min-width:100%;white-space:nowrap">
    <thead>
      <tr style="background:#336699;color:white;cursor:pointer;
                 position:sticky;top:0;z-index:2">{th}</tr>
    </thead>
    <tbody>{rows}</tbody>
  </table>
</div>"""


def fig_heatmap_with_dendro(names, lcpm, genes_f, n_filt):
    """
    Heatmap of top 500 variable genes with:
      - sample dendrogram on the left (horizontal)
      - genes ordered by clustering on the x-axis
      - samples ordered by clustering on the y-axis
    """
    n     = len(names)
    n_top = min(500, n_filt)

    var_idx   = np.argsort(lcpm.var(axis=0))[-n_top:]
    mat       = lcpm[:, var_idx]
    sub_genes = genes_f[var_idx]

    # Cluster genes
    try:
        gene_link  = linkage(pdist(mat.T, metric="euclidean"), method="ward")
        gene_order = leaves_list(gene_link)
    except Exception:
        gene_order = np.arange(n_top)

    # Cluster samples — keep dendrogram for drawing
    try:
        samp_link  = linkage(pdist(mat, metric="euclidean"), method="ward")
        samp_order = leaves_list(samp_link)
        samp_dend  = dendrogram(samp_link, labels=list(names), no_plot=True)
    except Exception:
        samp_order = np.arange(n)
        samp_dend  = None

    mat_ord   = mat[np.ix_(samp_order, gene_order)]
    names_ord = [names[i] for i in samp_order]
    genes_ord = [sub_genes[i] for i in gene_order]
    mat_c     = mat_ord - mat_ord.mean(axis=0)   # centre per gene

    # Two-column layout: narrow dendrogram | wide heatmap
    fig = make_subplots(rows=1, cols=2,
                        column_widths=[0.09, 0.91],
                        horizontal_spacing=0.004)

    # ── Sample dendrogram (horizontal: x=height, y=normalised position) ──────
    # scipy positions leaves at 5, 15, 25 … → normalise to 0, 1, 2 …
    if samp_dend:
        for xs_raw, ys_h, col in zip(samp_dend["icoord"],
                                      samp_dend["dcoord"],
                                      samp_dend["color_list"]):
            y_norm = [(p - 5) / 10 for p in xs_raw]
            fig.add_trace(go.Scatter(
                x=ys_h, y=y_norm, mode="lines",
                line=dict(color=safe_color(col), width=1),
                hoverinfo="skip", showlegend=False,
            ), row=1, col=1)

    # ── Heatmap ───────────────────────────────────────────────────────────────
    tfs_y = max(6, 10 - n // 8)
    tfs_x = max(4, 7 - n_top // 100)
    fig.add_trace(go.Heatmap(
        z=mat_c.tolist(),
        x=genes_ord if n_top <= 60 else None,
        y=names_ord,
        colorscale="RdBu", reversescale=True, zmid=0,
        colorbar=dict(title="Centred<br>log₂(CPM+1)", thickness=14),
        hovertemplate="Sample: %{y}<br>Value: %{z:.2f}<extra></extra>",
    ), row=1, col=2)

    base(fig,
         f"Heatmap — top {n_top} most variable genes "
         f"(centred log₂ CPM+1; rows & columns ordered by Ward clustering; "
         f"left panel = sample dendrogram)",
         max(520, n * 20 + 130))

    # Style dendrogram: reversed x (root at left, leaves touching heatmap),
    # y range matches heatmap categorical positions 0 … n-1
    fig.update_xaxes(autorange="reversed", showticklabels=False,
                     showgrid=False, zeroline=False, row=1, col=1)
    fig.update_yaxes(range=[-0.5, n - 0.5], showticklabels=False,
                     showgrid=False, zeroline=False, row=1, col=1)

    # Style heatmap
    fig.update_xaxes(
        title=dict(text="Genes (ordered by clustering)" if n_top > 60 else "Gene"),
        showticklabels=(n_top <= 60),
        tickangle=-45, tickfont_size=tfs_x,
        row=1, col=2,
    )
    fig.update_yaxes(tickfont_size=tfs_y, row=1, col=2)

    return fig

# ── Section 5: Group Estimation ───────────────────────────────────────────────

def _kmeans_on_embedding(names, lcpm, n_top, emb, elbl):
    """k-means (k=2) on top-variable gene space, shown on existing embedding."""
    n  = len(names)
    X  = lcpm[:, np.argsort(lcpm.var(axis=0))[-n_top:]]
    kl = KMeans(n_clusters=2, random_state=42, n_init=10).fit_predict(X)
    GP = ["#1f77b4", "#d62728"]
    mode = "markers+text" if n <= 30 else "markers"
    fig = go.Figure()
    for g in range(2):
        idx = [i for i in range(n) if kl[i] == g]
        fig.add_trace(go.Scatter(
            x=[float(emb[i, 0]) for i in idx],
            y=[float(emb[i, 1]) for i in idx],
            mode=mode,
            text=[names[i] for i in idx] if n <= 30 else None,
            textposition="top center", textfont_size=8,
            marker=dict(color=GP[g], size=9, line=dict(color="white", width=1)),
            name=f"k-means cluster {g + 1}",
            hovertemplate=("<b>%{text}</b><br>" if n <= 30 else "") +
                          f"Cluster {g+1}<extra></extra>",
        ))
    base(fig, f"{elbl} Embedding — k-means clusters (k=2, "
              f"top {n_top} variable genes; hover to identify samples)",
         max(500, n * 7 + 160))
    fig.update_layout(xaxis_title=f"{elbl} 1", yaxis_title=f"{elbl} 2")
    return fig, kl


def _hierarchical_k(link, max_k):
    """Find optimal k from Ward dendrogram by the largest gap in merge heights."""
    dists = link[:, 2]
    # look at the last max_k merges (from the top of the tree)
    top_dists = dists[-(max_k - 1):][::-1]   # decreasing order
    if len(top_dists) < 2:
        return 2
    gaps = np.diff(top_dists) * -1            # positive where height drops fast
    k    = int(np.argmax(gaps)) + 2           # gap[0] → k=2, gap[1] → k=3, …
    return max(2, min(k, max_k))


def figs_hierarchical_groups(names, lcpm, n_top):
    """PCA + group table from hierarchical dendrogram auto-cut."""
    n  = len(names)
    if n < 3:
        return None, None
    X  = lcpm[:, np.argsort(lcpm.var(axis=0))[-n_top:]]
    mk = min(6, n - 1)
    try:
        link = linkage(pdist(X, metric="euclidean"), method="ward")
        k    = _hierarchical_k(link, mk)
        kl   = fcluster(link, k, criterion="maxclust") - 1   # 0-indexed
    except Exception as e:
        warn(f"Hierarchical groups: {e}"); return None, None

    GP = ["#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b"]
    try:
        emb  = PCA(n_components=2).fit_transform(X)
        pct2 = PCA(n_components=2).fit(X).explained_variance_ratio_ * 100
    except Exception:
        return None, None

    mode = "markers+text" if n <= 30 else "markers"
    fig_p = go.Figure()
    for g in range(k):
        idx = [i for i in range(n) if kl[i] == g]
        fig_p.add_trace(go.Scatter(
            x=[float(emb[i, 0]) for i in idx],
            y=[float(emb[i, 1]) for i in idx],
            mode=mode,
            text=[names[i] for i in idx] if n <= 30 else None,
            textposition="top center", textfont_size=8,
            marker=dict(color=GP[g % len(GP)], size=10,
                        line=dict(color="white", width=1)),
            name=f"Group {g + 1}",
            hovertemplate=("<b>%{text}</b><br>" if n <= 30 else "") +
                          f"Group {g+1}<extra></extra>",
        ))
    base(fig_p, f"PCA — hierarchical groups (k={k}, automatic dendrogram cut at largest gap; "
                f"top {n_top} variable genes, log₂ CPM+1)",
         max(480, n * 7 + 160))
    fig_p.update_layout(xaxis_title=f"PC1 ({pct2[0]:.1f}%)",
                        yaxis_title=f"PC2 ({pct2[1]:.1f}%)")

    grp = [[f"Group {g+1}", f"{int(sum(kl==g))}",
             ", ".join(names[i] for i in range(n) if kl[i] == g)]
            for g in range(k)]
    fig_t = go.Figure(go.Table(
        header=dict(values=["<b>Group</b>", "<b>n samples</b>", "<b>Members</b>"],
                    fill_color="#336699", font=dict(color="white", size=11),
                    align="center", height=30),
        cells=dict(values=[[r[0] for r in grp], [r[1] for r in grp],
                            [r[2] for r in grp]],
                   fill_color=[["#f5f5f5" if i % 2 == 0 else "white"
                                 for i in range(k)]] * 3,
                   font_size=10, height=26, align=["center", "center", "left"]),
    ))
    base(fig_t, f"Hierarchical Group Assignments — k={k}, Ward linkage, "
                f"automatic cut at largest inter-cluster distance gap",
         max(200, k * 44 + 130))
    fig_t.update_layout(margin=dict(l=20, r=20, t=60, b=20))
    return fig_p, fig_t


def figs_consensus(names, lcpm, n_top, run_warnings):
    n     = len(names)
    max_k = min(6, n - 1)
    if max_k < 2:
        run_warnings.append("Consensus clustering needs ≥ 3 samples.")
        return {}

    X       = lcpm[:, np.argsort(lcpm.var(axis=0))[-n_top:]]
    N_ITER  = 100
    seeds   = np.random.default_rng(42).integers(0, 100_000, N_ITER)
    K_RANGE = range(2, max_k + 1)

    cooc, labs = {}, {}
    for k in K_RANGE:
        M = np.zeros((n, n))
        for seed in seeds:
            try:
                l = KMeans(n_clusters=k, random_state=int(seed),
                           n_init=1, max_iter=100).fit_predict(X)
                for i in range(n):
                    for j in range(n):
                        if l[i] == l[j]: M[i, j] += 1
            except Exception:
                pass
        M /= N_ITER; cooc[k] = M
        try:
            labs[k] = KMeans(n_clusters=k, random_state=42, n_init=20).fit_predict(M)
        except Exception:
            labs[k] = np.zeros(n, dtype=int)

    areas  = {k: float(np.mean(cooc[k][np.triu_indices(n, k=1)])) for k in K_RANGE}
    klist  = sorted(K_RANGE)
    deltas = {klist[i]: areas[klist[i]] - areas[klist[i - 1]] for i in range(1, len(klist))}
    best_k = max(deltas, key=deltas.get) if deltas else 2
    best_l = labs[best_k]

    order    = np.argsort(best_l)
    nms      = [names[i] for i in order]
    co_ord   = cooc[best_k][np.ix_(order, order)]
    tfs      = max(6, 10 - n // 8)
    ann      = [[f"{co_ord[i,j]:.2f}" for j in range(n)] for i in range(n)] if n <= 40 else None
    shapes   = []
    for b in np.where(np.diff(best_l[order]))[0] + 1:
        for kw in [{"x0": b-.5,"x1": b-.5,"y0": -.5,"y1": n-.5},
                   {"y0": b-.5,"y1": b-.5,"x0": -.5,"x1": n-.5}]:
            shapes.append(dict(type="line", line=dict(color="red", width=2),
                               layer="above", **kw))

    # Co-occurrence heatmap
    fig_c = go.Figure(go.Heatmap(
        z=co_ord.tolist(), x=nms, y=nms,
        text=ann, texttemplate="%{text}" if ann else None,
        textfont_size=max(5, tfs - 2),
        colorscale="Blues", zmin=0, zmax=1,
        colorbar=dict(title="Co-occurrence", thickness=14),
        hovertemplate="<b>%{x}</b> vs <b>%{y}</b><br>co-occurrence = %{z:.3f}<extra></extra>",
    ))
    base(fig_c, f"Consensus Co-occurrence Matrix — best k={best_k} "
                f"({N_ITER} iterations, top {n_top} variable genes; "
                f"red lines = cluster boundaries)", max(480, n * 22 + 120))
    fig_c.update_layout(shapes=shapes,
                        xaxis=dict(tickangle=-45, tickfont_size=tfs),
                        yaxis=dict(tickfont_size=tfs, autorange="reversed"))

    # Stability chart
    ks     = sorted(areas.keys())
    fig_s  = make_subplots(specs=[[{"secondary_y": True}]])
    fig_s.add_trace(go.Bar(
        x=ks, y=[areas[k] for k in ks],
        marker_color=["red" if k == best_k else "#336699" for k in ks],
        name="CDF area",
        hovertemplate="k=%{x}<br>CDF area=%{y:.3f}<extra></extra>",
    ))
    if deltas:
        dk = sorted(deltas)
        fig_s.add_trace(go.Scatter(
            x=dk, y=[deltas[k] for k in dk], mode="lines+markers",
            name="Delta", line_color="green", marker_size=7,
            hovertemplate="k=%{x}<br>Δ=%{y:.4f}<extra></extra>",
        ), secondary_y=True)
    base(fig_s, f"Consensus Clustering Stability — CDF area per k "
                f"(red = chosen k={best_k}; green = delta between adjacent k)", 390)
    fig_s.update_layout(
        xaxis=dict(title="k (number of clusters)", tickvals=ks),
        yaxis=dict(title="Mean co-occurrence (CDF area)"),
        yaxis2=dict(
            title=dict(text="Delta CDF area", font=dict(color="green")),
            overlaying="y", side="right",
            tickfont=dict(color="green"),
        ),
    )

    # PCA coloured by consensus group
    GP = ["#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b"]
    fig_p = None
    try:
        emb  = PCA(n_components=2).fit_transform(X)
        pct2 = PCA(n_components=2).fit(X).explained_variance_ratio_ * 100
        mode = "markers+text" if n <= 30 else "markers"
        fig_p = go.Figure()
        for g in range(best_k):
            idx = [i for i in range(n) if best_l[i] == g]
            fig_p.add_trace(go.Scatter(
                x=[float(emb[i, 0]) for i in idx],
                y=[float(emb[i, 1]) for i in idx],
                mode=mode,
                text=[names[i] for i in idx] if n <= 30 else None,
                textposition="top center", textfont_size=8,
                marker=dict(color=GP[g % len(GP)], size=10,
                            line=dict(color="white", width=1)),
                name=f"Consensus group {g + 1}",
                hovertemplate=("<b>%{text}</b><br>" if n <= 30 else "") +
                              f"Group {g+1}<extra></extra>",
            ))
        base(fig_p, f"PCA — samples coloured by consensus group (k={best_k})",
             max(480, n * 7 + 160))
        fig_p.update_layout(xaxis_title=f"PC1 ({pct2[0]:.1f}%)",
                            yaxis_title=f"PC2 ({pct2[1]:.1f}%)")
    except Exception as e:
        run_warnings.append(f"Consensus PCA: {e}")

    # Group assignment table
    grp = [[f"Group {g+1}", f"{int(sum(best_l==g))}",
             ", ".join(names[i] for i in range(n) if best_l[i] == g)]
            for g in range(best_k)]
    fig_g = go.Figure(go.Table(
        header=dict(values=["<b>Group</b>", "<b>n samples</b>", "<b>Members</b>"],
                    fill_color="#336699", font=dict(color="white", size=11),
                    align="center", height=30),
        cells=dict(values=[[r[0] for r in grp], [r[1] for r in grp],
                            [r[2] for r in grp]],
                   fill_color=[["#f5f5f5" if i % 2 == 0 else "white"
                                 for i in range(best_k)]] * 3,
                   font_size=10, height=26, align=["center", "center", "left"]),
    ))
    base(fig_g, f"Estimated Group Assignments — consensus k={best_k}, "
                f"no prior label information used",
         max(200, best_k * 44 + 130))
    fig_g.update_layout(margin=dict(l=20, r=20, t=60, b=20))

    return {"cooc": fig_c, "stab": fig_s, "pca": fig_p, "grp": fig_g,
            "best_k": best_k, "grp_labels": best_l}

# ══════════════════════════════════════════════════════════════════════════════
# HTML TEMPLATE
# ══════════════════════════════════════════════════════════════════════════════

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RNA-seq QC Report</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:Arial,sans-serif;font-size:13px;color:#222;background:#f4f4f4}
.sidebar{position:fixed;top:0;left:0;bottom:0;width:200px;background:#fff;
  border-right:1px solid #ddd;overflow-y:auto;padding:14px 0}
.sidebar h1{font-size:13px;font-weight:bold;padding:0 14px 10px;
  border-bottom:1px solid #eee;line-height:1.4}
.nav-group{margin-top:10px}
.nav-label{font-size:10px;font-weight:bold;text-transform:uppercase;
  letter-spacing:.05em;color:#888;padding:0 14px 4px}
.nav-item{display:block;padding:5px 14px;font-size:12px;color:#444;
  text-decoration:none;border-left:3px solid transparent}
.nav-item:hover{background:#f0f0f0;color:#000}
.nav-item.active{border-left-color:#336699;color:#336699;font-weight:bold}
.sidebar-footer{padding:12px 14px;font-size:10px;color:#aaa;
  border-top:1px solid #eee;margin-top:14px;line-height:1.7}
.main{margin-left:200px;padding:20px 28px 60px;max-width:1400px}
.section{padding-top:36px}
.sec-title{font-size:16px;font-weight:bold;border-bottom:2px solid #336699;
  padding-bottom:6px;margin-bottom:4px}
.sec-desc{color:#555;font-size:12px;margin-bottom:14px;max-width:840px;line-height:1.5}
.card{background:#fff;border:1px solid #ddd;border-radius:4px;
  margin-bottom:16px;padding:14px}
.card-title{font-size:12px;font-weight:bold;color:#336699;margin-bottom:3px}
.card-note{font-size:11px;color:#777;margin-bottom:10px}
.stats-row{display:flex;flex-wrap:wrap;gap:10px;margin-bottom:16px}
.stat{background:#fff;border:1px solid #ddd;border-radius:4px;
  padding:8px 14px;min-width:120px}
.stat-label{font-size:10px;text-transform:uppercase;color:#888;margin-bottom:2px}
.stat-val{font-size:18px;font-weight:bold;color:#336699}
.stat-sub{font-size:10px;color:#aaa;margin-top:1px}
table.qc{width:100%;border-collapse:collapse;font-size:11px}
table.qc th{background:#336699;color:white;padding:7px 10px;text-align:left}
table.qc td{padding:6px 10px;border-bottom:1px solid #eee}
table.qc tr:nth-child(even) td{background:#f7f7f7}
table.qc tr:hover td{background:#eef4ff}
.warn{background:#fff8f0;border:1px solid #f5c06a;border-radius:4px;
  padding:10px 14px;margin-top:10px;font-size:11px;color:#8a5000}
.warn strong{display:block;margin-bottom:4px}
#top-genes-tbl th{cursor:pointer;user-select:none}
#top-genes-tbl th:hover{background:#2a5580}
hr.div{border:none;border-top:1px solid #ddd;margin:40px 0 0}
</style>
</head>
<body>
<aside class="sidebar">
  <h1>RNA-seq QC<br>Report</h1>
  __NAV__
  <div class="sidebar-footer">
    __TS__<br>__INPUT__<br>
    __NS__ samples<br>__NGR__ genes (raw)<br>__NGF__ genes (filtered)
  </div>
</aside>
<main class="main">__CONTENT__</main>
<script>
// ── Plotly rendering ──
const plots = __PLOTS__;
Object.entries(plots).forEach(([id,js]) => {
  const el = document.getElementById(id);
  if (!el) return;
  const s = JSON.parse(js);
  Plotly.newPlot(el, s.data, s.layout, {
    responsive:true, displayModeBar:true,
    modeBarButtonsToRemove:['lasso2d','select2d','autoScale2d'],
    toImageButtonOptions:{format:'svg',filename:id},
  });
});

// ── Sidebar scroll-spy ──
const secs = document.querySelectorAll('.section[id]');
const navs = document.querySelectorAll('.nav-item');
new IntersectionObserver(es => {
  es.forEach(e => { if (e.isIntersecting)
    navs.forEach(a => a.classList.toggle('active',
      a.getAttribute('href')==='#'+e.target.id));
  });
},{ rootMargin:'-20% 0px -70% 0px' }).observe
&& secs.forEach(s =>
  new IntersectionObserver(es => {
    if (es[0].isIntersecting)
      navs.forEach(a => a.classList.toggle('active',
        a.getAttribute('href')==='#'+es[0].target.id));
  },{ rootMargin:'-20% 0px -70% 0px' }).observe(s));

// ── Sortable top-genes table ──
let _gs = {};
function sortTbl(col) {
  const tbl = document.getElementById('top-genes-tbl');
  if (!tbl) return;
  const asc = _gs[col] !== 'asc'; _gs[col] = asc ? 'asc' : 'desc';
  const tbody = tbl.tBodies[0];
  [...tbody.rows].sort((a,b) => {
    const av = a.cells[col].dataset.val ?? a.cells[col].textContent.trim();
    const bv = b.cells[col].dataset.val ?? b.cells[col].textContent.trim();
    const an = parseFloat(av), bn = parseFloat(bv);
    if (!isNaN(an)&&!isNaN(bn)) return asc ? an-bn : bn-an;
    return asc ? av.localeCompare(bv) : bv.localeCompare(av);
  }).forEach(r => tbody.appendChild(r));
  [...tbl.querySelectorAll('thead th')].forEach((th,i) => {
    th.textContent = th.textContent.replace(/ [▲▼]$/,'');
    if (i===col) th.textContent += asc ? ' ▲' : ' ▼';
  });
}
function searchGeneTable(q) {
  const tbl = document.getElementById('top-genes-tbl');
  if (!tbl) return;
  const ql = q.toLowerCase();
  [...tbl.tBodies[0].rows].forEach(r => {
    r.style.display = r.textContent.toLowerCase().includes(ql) ? '' : 'none';
  });
}
</script>
</body>
</html>"""

# ── HTML assembly helpers ──────────────────────────────────────────────────────

def card(pid, title, note):
    return (f'<div class="card"><div class="card-title">{title}</div>'
            f'<div class="card-note">{note}</div>'
            f'<div id="{pid}" style="width:100%;min-height:300px"></div></div>')

def html_card(html_content, title, note):
    return (f'<div class="card"><div class="card-title">{title}</div>'
            f'<div class="card-note">{note}</div>'
            f'{html_content}</div>')

def section(sid, title, desc, body):
    return (f'<section class="section" id="{sid}">'
            f'<div class="sec-title">{title}</div>'
            f'<p class="sec-desc">{desc}</p>{body}'
            f'<hr class="div"></section>')

def nav_group(label, items):
    links = "".join(f'<a class="nav-item" href="#{s}">{n}</a>' for s, n in items)
    return f'<div class="nav-group"><div class="nav-label">{label}</div>{links}</div>'

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    if len(sys.argv) < 3:
        print("Usage: python rnaseq_html_report.py <counts_file> <output.html>")
        sys.exit(1)

    counts_f, out_html = sys.argv[1], sys.argv[2]
    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)

    warns_list = []
    plots      = {}    # plot_id → JSON string
    html_blks  = {}    # block_id → raw HTML string
    content    = []
    navs       = []

    def store(k, fn, *a, **kw):
        try:
            f = fn(*a, **kw)
            if f is not None:
                plots[k] = to_json(f)
        except Exception as e:
            warns_list.append(f"'{k}': {e}")
            warn(f"'{k}' failed:\n    {traceback.format_exc(limit=2)}")

    print("\nRNA-seq QC HTML Report")
    print("=" * 50)

    # ── Load ─────────────────────────────────────────────────────────────────
    with step("Loading"):
        try:
            raw, genes, names = load_counts(counts_f)
        except Exception as e:
            print(f"\n[ERROR] {e}"); sys.exit(1)
        n_samp, n_raw = raw.shape[1], raw.shape[0]
        print(f"\n    Samples: {n_samp}  |  Genes: {n_raw}")
        if n_samp < 2:
            warns_list.append("Only 1 sample — many plots will be skipped.")

    # ── Normalise ─────────────────────────────────────────────────────────────
    with step("Normalising"):
        lib, cpm_all, cpm_f, lcpm, keep = normalise(raw)
        genes_f  = genes[keep]
        n_filt   = lcpm.shape[1]
        n_before = (raw > 0).sum(axis=0).astype(int)
        n_after  = (lcpm > 0).sum(axis=1).astype(int)
        n_top    = min(2000, n_filt)
        print(f"\n    Filtered: {n_filt} genes  |  Top-var: {n_top}")

    colors = make_colors(names)

    # ── 1. Summary ────────────────────────────────────────────────────────────
    with step("Summary"):
        store("overview", fig_overview, names, lib, n_before, n_after, colors)

        # Stats mini-table
        lib_m = lib / 1e6
        stats_html = "".join(
            f'<div class="stat"><div class="stat-label">{l}</div>'
            f'<div class="stat-val">{v}</div>'
            f'<div class="stat-sub">{s}</div></div>'
            for l, v, s in [
                ("Samples",        str(n_samp),                ""),
                ("Genes (raw)",    f"{n_raw:,}",               ""),
                ("Genes (filtered)",f"{n_filt:,}",             "CPM > 1"),
                ("Lib size — mean",f"{lib_m.mean():.1f} M",    f"SD {lib_m.std():.1f} M"),
                ("Lib size — min", f"{lib_m.min():.1f} M",     ""),
                ("Lib size — max", f"{lib_m.max():.1f} M",     ""),
            ])

        rows = "".join(
            f"<tr><td>{s}</td><td>{lib[i]:,.0f}</td>"
            f"<td>{n_before[i]:,}</td><td>{n_after[i]:,}</td>"
            f"<td>{n_after[i]/max(n_before[i],1)*100:.1f}%</td></tr>"
            for i, s in enumerate(names))
        warn_html = ""
        if warns_list:
            items = "".join(f"• {w}<br>" for w in warns_list)
            warn_html = f'<div class="warn"><strong>Warnings</strong>{items}</div>'

        tbl_html = (
            f'<div class="card"><div class="card-title">Per-sample counts</div>'
            f'<div class="card-note">Library size = total mapped reads (raw counts). '
            f'Filtered = genes with CPM > 1 in ≥ 1 sample.</div>'
            f'<table class="qc"><thead><tr><th>Sample</th><th>Library size</th>'
            f'<th>Genes (raw)</th><th>Genes (CPM-filtered)</th><th>Kept %</th>'
            f'</tr></thead><tbody>{rows}</tbody></table>{warn_html}</div>')

        content.append(section(
            "s-summary", "Summary",
            "Quick dataset overview. The chart shows library sizes (bars) and "
            "CPM-filtered gene counts (diamonds) per sample. "
            "The table provides exact numbers. "
            "Outliers in either metric warrant investigation.",
            f'<div class="stats-row">{stats_html}</div>'
            + card("overview",
                   "Dataset Overview — library sizes and gene detection",
                   "Left: total mapped reads per sample; dashed line = mean. "
                   "Right: genes detected before (faded) and after (solid) CPM > 1 filtering — "
                   "the drop shows how many lowly-expressed genes were removed. "
                   "Same colour per sample across both panels.")
            + tbl_html))
        navs.append(nav_group("Overview", [("s-summary", "Summary")]))

    # ── 2. QC Metrics ─────────────────────────────────────────────────────────
    with step("QC metrics"):
        store("lib_sizes",  fig_library_sizes,  names, lib, colors)
        store("det_raw",    fig_detected_raw,   names, n_before, colors)
        store("det_filt",   fig_detected_filt,  names, n_after, colors)
        store("cpm_box",    fig_cpm_boxplot,    names, lcpm, colors)
        store("cpm_dens",   fig_cpm_density,    names, lcpm, colors)
        store("outlier_z",  fig_outlier_heatmap,names, lcpm, lib, n_after, n_filt)

        body = (
            card("lib_sizes", "Library Sizes",
                 "Total mapped reads per sample (raw counts). "
                 "Lines: mean (solid) and ±1 SD (dashed). Outliers deviate strongly from the mean.") +
            card("det_raw", "Genes Detected — Raw Counts",
                 "Genes with at least one read per sample. "
                 "Mean ±1 SD shown. Data: all raw counts.") +
            card("det_filt", "Genes Detected — After CPM > 1 Filter",
                 "Genes passing CPM > 1 in ≥ 1 sample. "
                 "A sample with far fewer genes than peers may be degraded or swapped.") +
            card("cpm_box", "CPM Distribution — Boxplots",
                 "log₂(CPM+1) per sample (CPM-filtered genes). "
                 "Medians should align across all samples; systematic shifts indicate bias.") +
            card("cpm_dens", "CPM Distribution — Density Curves",
                 "log₂(CPM+1) density per sample. "
                 "All samples should broadly overlap. "
                 "Toggle individual samples via the legend.") +
            card("outlier_z", "Multi-Metric Outlier Overview",
                 "Robust Z-score = (value − median) / MAD per QC metric. "
                 "Red border = |Z| ≥ 2 on ≥ 1 metric. "
                 "Hover for exact values. "
                 "Library size uses raw counts; all other metrics use CPM-filtered genes."))
        content.append(section("s-qc", "QC Metrics",
            "Per-sample quality control. Check library sizes for under-sequenced samples, "
            "CPM distributions for systematic shifts, and the outlier heatmap for "
            "samples that deviate on multiple metrics simultaneously.",
            body))
        navs.append(nav_group("QC", [("s-qc", "QC Metrics")]))

    # ── 3. Sample Structure ───────────────────────────────────────────────────
    with step("Sample structure"):
        # Run consensus early so we can pass group labels to the correlation plot
        try:
            _cres_early = figs_consensus(names, lcpm, n_top, [])
            _early_labels = (_cres_early["grp_labels"]
                             if _cres_early and "grp_labels" in _cres_early else None)
        except Exception:
            _early_labels = None

        store("corr", fig_correlation, names, lcpm, n_filt, _early_labels)
        store("dend", fig_ward_dendrogram, names, lcpm, n_filt)

        f3d, C, pct = fig_pca_3d(names, lcpm, n_top, colors)
        if f3d: plots["pca3d"] = to_json(f3d)
        store("pca_scree",    fig_scree,        lcpm, n_top)
        store("pca_loadings", fig_pca_loadings, lcpm, genes_f, n_top)
        store("scatter_pair", fig_sample_scatter_pairs, names, lcpm, colors)

        if n_samp <= 12:
            store("embed", fig_ma, names, lcpm)
            emb_card = card("embed", "MA Plots — each sample vs pseudo-bulk mean",
                            "M = log₂(sample) − log₂(mean). Red dashed = M=0. "
                            "Green = mean M ±1 SD. Systematic shift = sample-level bias. "
                            "Data: CPM-filtered genes.")
            emb_lbl = "MA Plots"
            emb, elbl = None, "PC"
        else:
            try:
                f_es, emb, elbl = fig_umap_sample(names, lcpm, n_top, colors, warns_list)
                plots["embed"] = to_json(f_es)
            except Exception as e:
                warns_list.append(f"UMAP: {e}"); emb, elbl = None, "PC"
            emb_card = card("embed",
                            f"{elbl} Embedding — coloured by sample identity",
                            f"Top {n_top} variable genes (log₂ CPM+1, CPM-filtered). "
                            "Hover to identify samples. Isolated points = candidate outliers.")
            emb_lbl = elbl

        body = (
            card("corr", "Sample-to-Sample Pearson Correlation — two orderings",
                 f"Left: ordered by hierarchical clustering. "
                 f"Right: ordered by estimated group (if available) or by mean inter-sample r. "
                 f"Hover for exact r values. "
                 f"Data: log₂(CPM+1), CPM-filtered genes (n={n_filt}).") +
            card("dend", "Hierarchical Sample Clustering — Ward Linkage",
                 f"Euclidean distance on log₂(CPM+1), Ward linkage. "
                 f"Longer branches = more dissimilar. "
                 f"Isolated samples on long branches are outlier candidates.") +
            (card("pca3d", "PCA — Interactive 3D (PC1 × PC2 × PC3)",
                  f"Top {n_top} most variable genes. Drag to rotate. "
                  "Hover for exact coordinates. "
                  "Isolated points or tight clusters suggest outliers or biological groups.") if f3d else "") +
            card("pca_scree", "PCA Scree Plot",
                 "% variance explained per PC and cumulative total. "
                 "Many PCs needed to reach 80% (dashed) suggests complex structure.") +
            card("pca_loadings", "PCA Loadings — Top 15 Genes per PC",
                 "Genes with the largest absolute loading on PC1 (left) and PC2 (right). "
                 "Red = positive contribution, blue = negative. "
                 "These genes drive the separation seen in the 3D PCA.") +
            (card("scatter_pair", "Pairwise Sample Scatter — all sample pairs",
                  "Each panel = one sample pair. Red dashed = y=x (perfect agreement). "
                  "r = Pearson correlation. Points below the line = lower expression in that sample. "
                  "Only shown for n ≤ 6 samples.")
             if "scatter_pair" in plots else "") +
            emb_card)

        content.append(section("s-struct", "Sample Structure",
            "How similar are samples to each other? "
            "The correlation heatmap and dendrogram show pairwise similarity. "
            "PCA and the embedding reveal global low-dimensional structure. "
            "Outliers appear as isolated points or samples with low average correlation.",
            body))
        navs.append(nav_group("Structure", [
            ("s-struct", "Sample Structure"),
        ]))

    # ── 4. Gene Analysis ──────────────────────────────────────────────────────
    with step("Gene analysis"):
        store("mean_cv", fig_mean_cv, lcpm, genes_f, n_filt)

        try:
            tg_html = html_top_genes_table(names, cpm_f, genes_f, n_top=20)
            html_blks["top_genes"] = tg_html
        except Exception as e:
            warns_list.append(f"Top genes table: {e}")
            html_blks["top_genes"] = f"<p style='color:red'>Table failed: {e}</p>"

        store("heatmap", fig_heatmap_with_dendro, names, lcpm, genes_f, n_filt)

        body = (
            card("mean_cv", "Mean Expression vs Coefficient of Variation",
                 f"Each dot = one CPM-filtered gene (n={n_filt}). "
                 "Colour = local density. Red curve = running median CV (30 bins). "
                 "Red diamonds = top 20 high-CV genes (mean log₂ CPM > 1, labelled). "
                 "High-CV genes drive PCA/UMAP structure.") +
            html_card(html_blks["top_genes"],
                      "Top 20 Most Highly Expressed Genes — all samples (linear CPM)",
                      "Ranked by mean CPM across all samples. "
                      "Click any column header to sort. "
                      "Use the search box to filter by gene ID or sample name. "
                      "Light red rows = gene accounts for >5% of total mean CPM "
                      "(may indicate rRNA contamination or a dominant transcript). "
                      "Data: linear CPM, CPM-filtered genes.") +
            card("heatmap",
                 f"Heatmap — Top {min(500,n_filt)} Most Variable Genes",
                 "Centred log₂(CPM+1). "
                 "Rows (samples) and columns (genes) ordered by Ward hierarchical clustering. "
                 "Left panel = sample dendrogram. "
                 "Red = above mean, blue = below mean. "
                 "Use zoom/pan to inspect gene clusters."))
        content.append(section("s-genes", "Gene Analysis",
            "Gene-level diagnostics: the CV plot reveals the mean–variance relationship "
            "and flags the most variable genes. The top-genes table shows the "
            "most abundant transcripts across all samples (useful for spotting "
            "contamination). The heatmap shows how the most variable genes "
            "separate samples.",
            body))
        navs.append(nav_group("Genes", [("s-genes", "Gene Analysis")]))

    # ── 5. Group Estimation ───────────────────────────────────────────────────
    with step("Group estimation"):
        ge_body = ""

        # Silhouette score table
        try:
            sil_html = html_silhouette_table(names, lcpm, n_top)
            if sil_html:
                ge_body += html_card(
                    sil_html,
                    "Silhouette Scores — Ward Hierarchical Groups (k = 2 … min(6, n−1))",
                    "Silhouette score ranges from −1 to 1. "
                    "Values > 0.5 indicate well-separated groups; < 0.2 suggests overlap. "
                    "Use alongside the consensus matrix to judge group confidence. "
                    "Green row = best k by silhouette."
                )
        except Exception as e:
            warns_list.append(f"Silhouette table: {e}")

        # k-means on embedding
        if emb is not None:
            try:
                f_ec, kl = _kmeans_on_embedding(names, lcpm, n_top, emb, elbl)
                plots["km_embed"] = to_json(f_ec)
                ge_body += card("km_embed",
                                f"{elbl} Embedding — k-means clusters (k=2)",
                                f"Same embedding as Sample Structure section, "
                                f"coloured by k-means assignment (k=2, "
                                f"top {n_top} variable genes). "
                                f"Purely exploratory — no ground-truth used.")
            except Exception as e:
                warns_list.append(f"k-means embed: {e}")

        # Hierarchical cut
        try:
            f_hp, f_ht = figs_hierarchical_groups(names, lcpm, n_top)
            if f_hp:
                plots["hc_pca"] = to_json(f_hp)
                plots["hc_tbl"] = to_json(f_ht)
                ge_body += (
                    card("hc_pca",
                         "PCA — Hierarchical Groups (automatic dendrogram cut)",
                         "Ward dendrogram cut at the k where the largest gap in "
                         "merge heights occurs (k = 2..6 considered). "
                         "Complements consensus clustering with a faster, "
                         "deterministic method.") +
                    card("hc_tbl", "Hierarchical Group Assignments",
                         "Sample-to-group assignments from the hierarchical cut method."))
        except Exception as e:
            warns_list.append(f"Hierarchical groups: {e}")

        # Consensus clustering
        try:
            cres = figs_consensus(names, lcpm, n_top, warns_list)
            if cres and "cooc" in cres:
                bk = cres["best_k"]
                plots["cons_cooc"] = to_json(cres["cooc"])
                plots["cons_stab"] = to_json(cres["stab"])
                if cres.get("pca"): plots["cons_pca"] = to_json(cres["pca"])
                plots["cons_grp"]  = to_json(cres["grp"])
                ge_body += (
                    card("cons_cooc",
                         f"Consensus Co-occurrence Matrix (best k={bk})",
                         f"Fraction of 100 k-means runs (k=2..{min(6,n_samp-1)}) "
                         f"where each pair of samples shared a cluster. "
                         f"Values near 1 (dark) = consistently grouped. "
                         f"Values near 0.5 = ambiguous. "
                         f"Red lines = cluster boundaries at best k={bk}. "
                         f"Data: top {n_top} variable genes.") +
                    card("cons_stab",
                         "Consensus Stability — CDF Area per k",
                         f"Higher bar = more decisive clustering at that k. "
                         f"Red bar = chosen k={bk} (largest jump in delta). "
                         f"Green line = delta between adjacent k values.") +
                    (card("cons_pca",
                          f"PCA — Consensus Group Colours (k={bk})",
                          "PCA coloured by the consensus group assignment.")
                     if "cons_pca" in plots else "") +
                    card("cons_grp",
                         f"Consensus Group Assignments (k={bk})",
                         "Final group assignments. No prior sample labels were used."))
        except Exception as e:
            warns_list.append(f"Consensus: {e}")
            warn(f"Consensus: {traceback.format_exc(limit=2)}")

        if ge_body:
            content.append(section("s-groups", "Group Estimation",
                "Three complementary unsupervised methods estimate potential sample groups "
                "without using any prior labels. "
                "(1) k-means (k=2) on the UMAP/PCA embedding — fast, simple. "
                "(2) Hierarchical dendrogram cut at the largest distance gap — deterministic. "
                "(3) Consensus clustering over 100 random restarts — more robust, "
                "quantifies uncertainty via co-occurrence. "
                "Agreement across methods increases confidence in a grouping.",
                ge_body))
            navs.append(nav_group("Groups", [("s-groups", "Group Estimation")]))

    # ── Write HTML ─────────────────────────────────────────────────────────────
    with step("Writing"):
        ts   = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        html = (HTML
                .replace("__NAV__",     "".join(navs))
                .replace("__TS__",      ts)
                .replace("__INPUT__",   os.path.basename(counts_f))
                .replace("__NS__",      str(n_samp))
                .replace("__NGR__",     f"{n_raw:,}")
                .replace("__NGF__",     f"{n_filt:,}")
                .replace("__CONTENT__", "".join(content))
                .replace("__PLOTS__",   json.dumps(plots)))
        with open(out_html, "w", encoding="utf-8") as f:
            f.write(html)
        kb = os.path.getsize(out_html) // 1024
        print(f"\n    {out_html}  ({kb} KB)")

    print("\n" + "=" * 50)
    if warns_list:
        print("Warnings:")
        for w in warns_list:
            print(f"  * {w}")
    print(f"\nDone → {out_html}\n")

if __name__ == "__main__":
    main()