"""
RNA-seq QC Report Generator
Produces a single PDF with QC and exploratory plots from a featureCounts count matrix.
Usage: python comprehensive_report_python.py <counts_file> <output_pdf>
"""

import os
import sys
import time
import warnings
import traceback
from contextlib import contextmanager

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
from scipy.stats import gaussian_kde, median_abs_deviation

warnings.filterwarnings("ignore")

# ── Palette ───────────────────────────────────────────────────────────────────
BLUE    = "#2166AC"
RED     = "#D6604D"
GREEN   = "#4DAC26"
GREY    = "#878787"
BG      = "#F7F7F7"
DIVIDER = "#CCCCCC"

sns.set_theme(style="whitegrid", font_scale=1.05)
plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor":   BG,
    "axes.edgecolor":   DIVIDER,
    "grid.color":       "white",
    "grid.linewidth":   1.0,
    "font.family":      "DejaVu Sans",
})

# ── Scaling thresholds ────────────────────────────────────────────────────────
# Above this many samples, inline scatter labels are replaced by a legend
LABEL_INLINE_MAX   = 30
# Above this, tick labels on heatmaps are hidden entirely (too dense)
HEATMAP_TICK_MAX   = 60
# Font size floor for x-tick labels on bar charts
BAR_FONT_MIN       = 5


# ── Helpers ───────────────────────────────────────────────────────────────────

@contextmanager
def log_step(name):
    t = time.time()
    print(f"  [{time.strftime('%H:%M:%S')}] {name} ...", end=" ", flush=True)
    yield
    print(f"done ({time.time()-t:.1f}s)")


def warn(msg):
    print(f"\n  [WARNING] {msg}", flush=True)


def assign_colors(labels):
    palette = sns.color_palette("tab20", n_colors=max(len(labels), 1))
    return {s: palette[i % len(palette)] for i, s in enumerate(labels)}


def add_page_title(fig, title, subtitle=""):
    fig.text(0.5, 0.97, title, ha="center", va="top",
             fontsize=14, fontweight="bold", color="#222222")
    if subtitle:
        fig.text(0.5, 0.945, subtitle, ha="center", va="top",
                 fontsize=9, color=GREY)


def add_data_note(ax, note):
    """Small italic label in bottom-right corner of an axes."""
    ax.text(1.0, -0.18, note, transform=ax.transAxes,
            fontsize=7, color=GREY, style="italic", ha="right", va="top")


def bar_tick_fontsize(n):
    """Return a reasonable x-tick font size that scales down with sample count."""
    return max(BAR_FONT_MIN, min(9, int(9 - (n - 20) * 0.12)))


def heatmap_tick_fontsize(n):
    """Return tick font size for heatmaps; 0 = hide ticks entirely."""
    if n > HEATMAP_TICK_MAX:
        return 0
    return max(4, min(9, int(9 - (n - 20) * 0.10)))


def mean_sd_lines(ax, values, orientation="h"):
    """Draw mean ± 1 SD reference lines."""
    m  = np.mean(values)
    sd = np.std(values)
    if orientation == "h":
        ax.axhline(m,      color=BLUE, linestyle="-",  linewidth=1.2,
                   label=f"Mean: {m:.2f}", zorder=4)
        ax.axhline(m + sd, color=BLUE, linestyle="--", linewidth=0.8,
                   label=f"+1 SD: {m+sd:.2f}", zorder=4)
        ax.axhline(m - sd, color=BLUE, linestyle="--", linewidth=0.8,
                   label=f"-1 SD: {m-sd:.2f}", zorder=4)
    else:
        ax.axvline(m,      color=BLUE, linestyle="-",  linewidth=1.2, zorder=4)
        ax.axvline(m + sd, color=BLUE, linestyle="--", linewidth=0.8, zorder=4)
        ax.axvline(m - sd, color=BLUE, linestyle="--", linewidth=0.8, zorder=4)


def adjust_text_labels(ax, texts):
    """Iterative repulsion to prevent label overlap in scatter plots."""
    if len(texts) < 2:
        return
    fig = ax.get_figure()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    MAX_ITER = 60
    PAD = 2.0
    for _ in range(MAX_ITER):
        moved  = False
        bboxes = [t.get_window_extent(renderer=renderer) for t in texts]
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                bi, bj = bboxes[i], bboxes[j]
                if (bi.x0 < bj.x1 + PAD and bi.x1 + PAD > bj.x0 and
                        bi.y0 < bj.y1 + PAD and bi.y1 + PAD > bj.y0):
                    ox = (bi.x0 + bi.x1) / 2 - (bj.x0 + bj.x1) / 2
                    oy = (bi.y0 + bi.y1) / 2 - (bj.y0 + bj.y1) / 2
                    dx = (bi.x1 - bi.x0 + bj.x1 - bj.x0) / 2 + PAD - abs(ox)
                    dy = (bi.y1 - bi.y0 + bj.y1 - bj.y0) / 2 + PAD - abs(oy)
                    inv = ax.transData.inverted()
                    if dy <= dx:
                        shift = inv.transform((0, dy / 2)) - inv.transform((0, 0))
                        sign  = 1 if oy >= 0 else -1
                        xi, yi = texts[i].get_position()
                        xj, yj = texts[j].get_position()
                        texts[i].set_position((xi, yi + sign * abs(shift[1])))
                        texts[j].set_position((xj, yj - sign * abs(shift[1])))
                    else:
                        shift = inv.transform((dx / 2, 0)) - inv.transform((0, 0))
                        sign  = 1 if ox >= 0 else -1
                        xi, yi = texts[i].get_position()
                        xj, yj = texts[j].get_position()
                        texts[i].set_position((xi + sign * abs(shift[0]), yi))
                        texts[j].set_position((xj - sign * abs(shift[0]), yj))
                    moved  = True
                    bboxes = [t.get_window_extent(renderer=renderer) for t in texts]
        if not moved:
            break


def scatter_with_labels(ax, xd, yd, colnames, sample_colors, s=70):
    """
    Scatter plot with sample labels.
    - n <= LABEL_INLINE_MAX : labels placed next to each point, overlap-adjusted
    - n >  LABEL_INLINE_MAX : legend box instead (stays readable at any n)
    """
    n = len(colnames)
    # Scale marker size down for large n
    s_scaled = max(20, s - max(0, n - 20) * 1.5)

    for i, name in enumerate(colnames):
        ax.scatter(xd[i], yd[i], color=sample_colors[name],
                   s=s_scaled, zorder=5, edgecolors="white", linewidths=0.5,
                   label=name if n > LABEL_INLINE_MAX else None)

    if n <= LABEL_INLINE_MAX:
        y_pad  = (yd.max() - yd.min()) * 0.03 if (yd.max() - yd.min()) > 0 else 0.01
        # Font size scales with n: 7 at n=5, 5 at n=30
        fs     = max(4, 7 - int((n - 5) * 0.10))
        texts  = []
        for i, name in enumerate(colnames):
            t = ax.text(xd[i], yd[i] + y_pad, name,
                        fontsize=fs, ha="center", va="bottom", color="#333333")
            texts.append(t)
        adjust_text_labels(ax, texts)
    else:
        # Legend — two columns if more than 20 entries
        ncol = 2 if n > 20 else 1
        ax.legend(fontsize=6, loc="best", framealpha=0.7,
                  markerscale=0.8, ncol=ncol,
                  handlelength=1, borderpad=0.4, labelspacing=0.3)


def set_bar_xlabels(ax, colnames, rotation=45):
    """Set x-tick labels on bar charts with auto-scaled font size."""
    n  = len(colnames)
    fs = bar_tick_fontsize(n)
    ax.set_xticks(range(n))
    ax.set_xticklabels(colnames, rotation=rotation, ha="right", fontsize=fs)


def apply_heatmap_ticks(ax, names_x, names_y=None):
    """
    Apply tick labels to a heatmap axes.
    Hides labels entirely when n > HEATMAP_TICK_MAX.
    names_y defaults to names_x if not provided.
    """
    if names_y is None:
        names_y = names_x
    nx = len(names_x)
    ny = len(names_y)
    fs_x = heatmap_tick_fontsize(nx)
    fs_y = heatmap_tick_fontsize(ny)

    ax.set_xticks(range(nx))
    ax.set_yticks(range(ny))
    if fs_x > 0:
        ax.set_xticklabels(names_x, rotation=45, ha="right", fontsize=fs_x)
    else:
        ax.set_xticklabels([])
    if fs_y > 0:
        ax.set_yticklabels(names_y, fontsize=fs_y)
    else:
        ax.set_yticklabels([])


def save_fig(pdf, fig):
    pdf.savefig(fig, bbox_inches="tight", dpi=150)
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# PLOT FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def plot_summary_table(pdf, colnames, lib_sizes, n_genes_before, n_genes_after,
                       warnings_list):
    """Page 1 — dataset overview table. Splits into multiple pages if n is large."""
    ROWS_PER_PAGE = 40
    col_labels = ["Sample", "Library Size", "Genes (all counts)",
                  "Genes (CPM > 1 filtered)", "Kept %"]

    all_rows = []
    for i, s in enumerate(colnames):
        all_rows.append([
            s,
            f"{lib_sizes[i]:,.0f}",
            f"{n_genes_before[i]:,}",
            f"{n_genes_after[i]:,}",
            f"{n_genes_after[i] / max(n_genes_before[i], 1) * 100:.1f}%",
        ])

    chunks = [all_rows[i:i + ROWS_PER_PAGE]
              for i in range(0, len(all_rows), ROWS_PER_PAGE)]

    for page_idx, rows in enumerate(chunks):
        fig_h = max(4, len(rows) * 0.30 + 2.5)
        fig, ax = plt.subplots(figsize=(11, fig_h))
        ax.axis("off")
        suffix = f" (page {page_idx+1}/{len(chunks)})" if len(chunks) > 1 else ""
        add_page_title(fig, f"RNA-seq QC Report — Dataset Summary{suffix}")

        table = ax.table(cellText=rows, colLabels=col_labels,
                         loc="center", cellLoc="center")
        table.auto_set_font_size(False)
        fs = max(6, 9 - max(0, len(rows) - 20) // 5)
        table.set_fontsize(fs)
        table.scale(1, max(1.2, 1.6 - len(rows) * 0.01))

        for (r, c), cell in table.get_celld().items():
            if r == 0:
                cell.set_facecolor(BLUE)
                cell.set_text_props(color="white", fontweight="bold")
            elif r % 2 == 0:
                cell.set_facecolor("#E8EEF5")
            else:
                cell.set_facecolor("white")
            cell.set_edgecolor(DIVIDER)

        if warnings_list and page_idx == 0:
            warn_text = "Warnings:\n" + "\n".join(f"  * {w}" for w in warnings_list)
            fig.text(0.05, 0.01, warn_text, fontsize=7, color="#CC4400",
                     va="bottom", family="monospace")

        save_fig(pdf, fig)


def plot_library_sizes(pdf, colnames, lib_sizes, sample_colors):
    """Page 2 — library sizes bar chart."""
    n   = len(colnames)
    fig_w = max(8, n * 0.35 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 6))
    add_page_title(fig, "Library Sizes",
                   "Total mapped reads per sample (all counts, pre-filter). "
                   "Blue line = mean, dashed = mean +/- 1 SD.")

    ax.bar(colnames, lib_sizes / 1e6,
           color=[sample_colors[s] for s in colnames],
           edgecolor="white", linewidth=0.5, zorder=3)
    mean_sd_lines(ax, lib_sizes / 1e6)
    ax.set_ylabel("Library size (millions of reads)")
    set_bar_xlabels(ax, colnames)
    ax.legend(fontsize=8)
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    add_data_note(ax, "Data: all raw counts")

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(pdf, fig)


def plot_cpm_distributions(pdf, colnames, lcpm, sample_colors, n_genes_filt):
    """Page 3 — density + boxplot of log2 CPM+1."""
    n         = len(colnames)
    data_note = f"Data: CPM-filtered genes (n={n_genes_filt}, CPM > 1 in >= 1 sample)"
    fig_w     = max(14, n * 0.25 + 8)
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, 6))
    add_page_title(fig, "CPM Distributions (log2 CPM+1)",
                   "Left: density curves — all samples should overlap. "
                   "Right: boxplots — medians should align. "
                   "Systematic shifts indicate sample-level bias.")

    # Density
    ax = axes[0]
    for i, s in enumerate(colnames):
        vals = lcpm[i, :]
        vals = vals[np.isfinite(vals)]
        try:
            kde_x = np.linspace(vals.min(), vals.max(), 300)
            kde   = gaussian_kde(vals, bw_method=0.3)
            ax.plot(kde_x, kde(kde_x), color=sample_colors[s],
                    linewidth=1.2, alpha=0.8, label=s)
        except Exception:
            warn(f"KDE failed for {s}, skipping density curve.")
    ax.set_xlabel("log2 CPM+1")
    ax.set_ylabel("Density")
    ax.set_title("Density per sample")
    if n <= 20:
        ax.legend(fontsize=7, loc="upper right", framealpha=0.7)
    add_data_note(ax, data_note)

    # Boxplot
    ax2 = axes[1]
    bp  = ax2.boxplot(
        [lcpm[i, :] for i in range(n)],
        patch_artist=True, notch=False,
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(linewidth=0.8),
        capprops=dict(linewidth=0.8),
        flierprops=dict(marker=".", markersize=1, alpha=0.2, color=GREY),
        zorder=3,
    )
    for patch, s in zip(bp["boxes"], colnames):
        patch.set_facecolor(sample_colors[s])
        patch.set_alpha(0.7)
    ax2.set_xticks(range(1, n + 1))
    ax2.set_xticklabels(colnames, rotation=45, ha="right",
                        fontsize=bar_tick_fontsize(n))
    ax2.set_ylabel("log2 CPM+1")
    ax2.set_title("Boxplot per sample")
    ax2.yaxis.grid(True, zorder=0)
    ax2.set_axisbelow(True)
    add_data_note(ax2, data_note)

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(pdf, fig)


def plot_detected_genes(pdf, colnames, n_genes_before, n_genes_after, sample_colors):
    """Page 4 — detected genes before/after filter."""
    n     = len(colnames)
    fig_w = max(14, n * 0.35 + 6)
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, 6))
    add_page_title(fig, "Detected Genes per Sample",
                   "Left: all genes with any count. Right: after CPM > 1 filter. "
                   "Blue line = mean, dashed = mean +/- 1 SD.")

    for ax, vals, title, note in [
        (axes[0], n_genes_before, "Genes detected (all counts, raw)",
         "Data: all raw counts, genes with count > 0"),
        (axes[1], n_genes_after,  "Genes detected (CPM > 1 filtered)",
         "Data: CPM-filtered counts (CPM > 1 in >= 1 sample)"),
    ]:
        ax.bar(colnames, vals,
               color=[sample_colors[s] for s in colnames],
               edgecolor="white", linewidth=0.5, zorder=3)
        mean_sd_lines(ax, np.array(vals, dtype=float))
        ax.set_ylabel("Number of detected genes")
        set_bar_xlabels(ax, colnames)
        ax.set_title(title, fontweight="bold")
        ax.legend(fontsize=8)
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        add_data_note(ax, note)

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(pdf, fig)


def plot_ma_plots(pdf, colnames, lcpm, n_genes_filt):
    """Page 5a — MA plots (only for n <= 10 samples)."""
    n         = len(colnames)
    ref       = lcpm.mean(axis=0)
    data_note = f"Data: CPM-filtered genes (n={n_genes_filt})"
    ncols     = min(4, n)
    nrows     = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 3.5, nrows * 3.2),
                             squeeze=False)
    add_page_title(fig, "MA Plots — Each Sample vs Pseudo-bulk Mean",
                   "M = log2(sample CPM+1) - log2(mean CPM+1).  "
                   "Green = sample mean M +/- 1 SD.  Points should scatter around M=0.")

    for idx, s in enumerate(colnames):
        r, c  = divmod(idx, ncols)
        ax    = axes[r][c]
        A     = (lcpm[idx, :] + ref) / 2
        M     = lcpm[idx, :] - ref
        try:
            xy    = np.vstack([A, M])
            z     = gaussian_kde(xy)(xy)
            order = z.argsort()
            ax.scatter(A[order], M[order], c=z[order], s=2, cmap="plasma",
                       alpha=0.6, linewidths=0, rasterized=True)
        except Exception:
            ax.scatter(A, M, s=2, alpha=0.3, color=BLUE, rasterized=True)

        mean_m = np.mean(M)
        sd_m   = np.std(M)
        ax.axhline(0,             color=RED,   linewidth=1.0, linestyle="--", zorder=4)
        ax.axhline(mean_m,        color=GREEN, linewidth=1.0, linestyle="-",
                   label=f"mean={mean_m:.2f}", zorder=4)
        ax.axhline(mean_m + sd_m, color=GREEN, linewidth=0.7, linestyle="--", zorder=4)
        ax.axhline(mean_m - sd_m, color=GREEN, linewidth=0.7, linestyle="--", zorder=4)
        ax.set_title(s, fontsize=8, fontweight="bold")
        ax.set_xlabel("A (mean log2 CPM+1)", fontsize=7)
        ax.set_ylabel("M (log2 ratio)", fontsize=7)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc="upper right")
        add_data_note(ax, data_note)

    for idx in range(n, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(pdf, fig)


def plot_umap_kmeans(pdf, colnames, lcpm, sample_colors, n_top, run_warnings):
    """Page 5b — UMAP + k=2 k-means (n > 10 samples)."""
    gene_vars = lcpm.var(axis=0)
    top_idx   = np.argsort(gene_vars)[-n_top:]
    X         = lcpm[:, top_idx]
    n         = len(colnames)
    data_note = f"Data: top {n_top} most variable genes (log2 CPM+1, CPM-filtered)"

    try:
        from umap import UMAP
        reducer   = UMAP(n_components=2,
                         n_neighbors=max(2, min(15, n - 1)),
                         min_dist=0.3, random_state=42, verbose=False)
        embedding = reducer.fit_transform(X)
        embed_lbl = "UMAP"
    except Exception as e:
        run_warnings.append(f"UMAP failed ({e}); using PCA 2D fallback.")
        embedding = PCA(n_components=2).fit_transform(X)
        embed_lbl = "PC"

    km             = KMeans(n_clusters=2, random_state=42, n_init=10)
    cluster_labels = km.fit_predict(X)
    cluster_pal    = [BLUE, RED]

    fig = plt.figure(figsize=(14, 6))
    add_page_title(fig, f"{embed_lbl} Embedding + k-means (k=2)",
                   f"Left: coloured by sample name. "
                   f"Right: coloured by k-means cluster. "
                   f"Labels shown inline (n<={LABEL_INLINE_MAX}) or as legend (n>{LABEL_INLINE_MAX}).")

    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.38)

    ax1 = fig.add_subplot(gs[0])
    scatter_with_labels(ax1, embedding[:, 0], embedding[:, 1],
                        colnames, sample_colors)
    ax1.set_xlabel(f"{embed_lbl} 1", fontsize=9)
    ax1.set_ylabel(f"{embed_lbl} 2", fontsize=9)
    ax1.set_title(f"{embed_lbl} — by sample", fontweight="bold")
    add_data_note(ax1, data_note)

    ax2 = fig.add_subplot(gs[1])
    cl_colors = {s: cluster_pal[cluster_labels[i]] for i, s in enumerate(colnames)}
    scatter_with_labels(ax2, embedding[:, 0], embedding[:, 1],
                        colnames, cl_colors)
    ax2.legend(handles=[Patch(facecolor=cluster_pal[k], label=f"Cluster {k+1}")
                         for k in range(2)],
               fontsize=9, loc="best", framealpha=0.8)
    ax2.set_xlabel(f"{embed_lbl} 1", fontsize=9)
    ax2.set_ylabel(f"{embed_lbl} 2", fontsize=9)
    ax2.set_title(f"{embed_lbl} — k-means (k=2)", fontweight="bold")
    add_data_note(ax2, data_note)

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    save_fig(pdf, fig)


def plot_sample_correlation(pdf, colnames, lcpm, sample_colors, n_genes_filt):
    """Page 6 — sample-to-sample Pearson correlation heatmap."""
    n         = len(colnames)
    cor       = np.corrcoef(lcpm)
    data_note = f"Data: CPM-filtered genes (n={n_genes_filt}, log2 CPM+1)"

    try:
        dist  = squareform(1 - cor, checks=False)
        dist  = np.clip(dist, 0, None)
        link  = linkage(dist, method="average")
        order = dendrogram(link, no_plot=True)["leaves"]
    except Exception:
        warn("Hierarchical clustering of correlation matrix failed; using original order.")
        order = list(range(n))

    cor_ord   = cor[np.ix_(order, order)]
    names_ord = [colnames[i] for i in order]

    # Scale figure so every cell is at least 0.18 inches
    cell_sz = max(0.18, min(0.45, 7.0 / n))
    fig_sz  = max(6, n * cell_sz + 2.5)
    fig     = plt.figure(figsize=(fig_sz, fig_sz))
    add_page_title(fig, "Sample-to-Sample Pearson Correlation",
                   "Ordered by hierarchical clustering. "
                   "Values close to 1 = high similarity. "
                   "Unexpected low-correlation samples are outlier candidates.")

    ax   = fig.add_axes([0.18, 0.14, 0.68, 0.70])
    vmin = max(float(cor_ord.min()) - 0.01, 0.5)
    im   = ax.imshow(cor_ord, cmap="RdYlBu_r", vmin=vmin, vmax=1.0, aspect="auto")
    apply_heatmap_ticks(ax, names_ord)

    # Annotate cells only if they are large enough to read
    if n <= 20:
        for i in range(n):
            for j in range(n):
                v = cor_ord[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=5, color="black" if v > 0.7 else "white")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Pearson r", fontsize=8)
    fig.text(0.99, 0.01, data_note, ha="right", va="bottom",
             fontsize=7, color=GREY, style="italic")
    save_fig(pdf, fig)


def plot_pca(pdf, colnames, lcpm, sample_colors, n_top):
    """Page 7 — PCA: PC1/PC2, PC1/PC3, scree plot."""
    gene_vars = lcpm.var(axis=0)
    X         = lcpm[:, np.argsort(gene_vars)[-n_top:]]
    data_note = f"Data: top {n_top} most variable genes (log2 CPM+1, CPM-filtered)"
    n_comp    = min(len(colnames), X.shape[1], 10)

    if n_comp < 2:
        warn("Not enough samples/genes for PCA.")
        return

    coords = PCA(n_components=n_comp).fit(X)
    pct    = coords.explained_variance_ratio_ * 100
    emb    = coords.transform(X)

    fig = plt.figure(figsize=(14, 5.5))
    add_page_title(fig, "Principal Component Analysis",
                   f"Top {n_top} most variable genes (log2 CPM+1, CPM-filtered). "
                   "Isolated samples = outliers. Natural groupings suggest batch or biology.")

    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.42)

    def pc_panel(ax, xd, yd, xl, yl, title):
        scatter_with_labels(ax, xd, yd, colnames, sample_colors)
        ax.set_xlabel(xl, fontsize=9)
        ax.set_ylabel(yl, fontsize=9)
        ax.set_title(title, fontweight="bold")
        ax.axhline(0, color=DIVIDER, linewidth=0.5)
        ax.axvline(0, color=DIVIDER, linewidth=0.5)
        add_data_note(ax, data_note)

    pc_panel(fig.add_subplot(gs[0]),
             emb[:, 0], emb[:, 1],
             f"PC1 ({pct[0]:.1f}%)", f"PC2 ({pct[1]:.1f}%)", "PC1 vs PC2")

    ax2 = fig.add_subplot(gs[1])
    if n_comp >= 3:
        pc_panel(ax2, emb[:, 0], emb[:, 2],
                 f"PC1 ({pct[0]:.1f}%)", f"PC3 ({pct[2]:.1f}%)", "PC1 vs PC3")
    else:
        ax2.set_visible(False)

    ax3 = fig.add_subplot(gs[2])
    cumvar = np.cumsum(pct)
    ax3.bar(range(1, n_comp + 1), pct, color=BLUE, alpha=0.75, label="Per-PC")
    ax3.plot(range(1, n_comp + 1), cumvar, color=RED, marker="o",
             markersize=4, linewidth=1.4, label="Cumulative")
    ax3.axhline(80, color=GREY, linestyle="--", linewidth=0.8, label="80% line")
    ax3.set_xlabel("PC")
    ax3.set_ylabel("% variance")
    ax3.set_title("Scree Plot", fontweight="bold")
    ax3.legend(fontsize=8)
    ax3.set_xticks(range(1, n_comp + 1))
    add_data_note(ax3, data_note)

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    save_fig(pdf, fig)


def plot_outlier_summary(pdf, colnames, lcpm, lib_sizes, n_genes_after, n_genes_filt):
    """Page 8 — robust Z-score heatmap across QC metrics."""
    n       = len(colnames)
    metrics = {
        "Library size\n(all counts)":          lib_sizes.astype(float),
        "Detected genes\n(CPM filtered)":       np.array(n_genes_after, dtype=float),
        "Median log2 CPM\n(CPM filtered)":      np.median(lcpm, axis=1),
        "CV (std/mean)\n(CPM filtered)":        np.array([
            np.std(lcpm[i]) / (np.mean(lcpm[i]) + 1e-9) for i in range(n)]),
        "Mean inter-sample r\n(CPM filtered)":  np.array([
            np.corrcoef(lcpm)[i].mean() for i in range(n)]),
    }

    zmat = np.zeros((n, len(metrics)))
    for j, vals in enumerate(metrics.values()):
        med        = np.median(vals)
        mad        = median_abs_deviation(vals) + 1e-9
        zmat[:, j] = (vals - med) / mad

    row_h = max(0.25, min(0.55, 8.0 / n))
    fig_h = max(5, n * row_h + 2.5)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    add_page_title(fig, "Multi-Metric Outlier Overview (Robust Z-scores)",
                   "Each cell = (value - median) / MAD. "
                   "Red border = |Z| >= 2 on at least one metric. "
                   "Data sources labelled in column headers.")

    vmax = max(3.0, float(np.abs(zmat).max()))
    im   = ax.imshow(zmat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(list(metrics.keys()), rotation=30, ha="right", fontsize=8)
    fs_y = heatmap_tick_fontsize(n)
    ax.set_yticks(range(n))
    ax.set_yticklabels(colnames if fs_y > 0 else [], fontsize=fs_y)

    # Cell annotations — skip if rows are too thin to read
    if row_h >= 0.22:
        ann_fs = max(5, min(8, int(row_h * 14)))
        for i in range(n):
            for j in range(len(metrics)):
                v = zmat[i, j]
                ax.text(j, i, f"{v:+.1f}", ha="center", va="center",
                        fontsize=ann_fs,
                        color="black" if abs(v) < 1.5 else "white",
                        fontweight="bold" if abs(v) >= 2 else "normal")

    for i in range(n):
        if np.any(np.abs(zmat[i]) >= 2):
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(metrics), 1,
                                        fill=False, edgecolor=RED, linewidth=1.8))

    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label(
        "Robust Z-score", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    save_fig(pdf, fig)


def plot_hierarchical_clustering(pdf, colnames, lcpm, n_genes_filt):
    """Page 9 — Ward linkage dendrogram."""
    data_note = f"Data: CPM-filtered genes (n={n_genes_filt}, log2 CPM+1)"
    try:
        link = linkage(pdist(lcpm, metric="euclidean"), method="ward")
    except Exception as e:
        warn(f"Hierarchical clustering failed: {e}")
        return

    n     = len(colnames)
    fig_w = max(8, n * 0.45 + 2)
    # Scale leaf font size
    lfs   = max(5, min(9, int(9 - (n - 20) * 0.12)))

    fig, ax = plt.subplots(figsize=(fig_w, 6))
    add_page_title(fig, "Hierarchical Sample Clustering (Euclidean / Ward)",
                   "Samples that cluster tightly are similar. "
                   "Long branch lengths indicate dissimilarity — potential outliers or distinct groups.")

    dendrogram(link, labels=colnames, ax=ax,
               leaf_rotation=45, leaf_font_size=lfs,
               color_threshold=0.7 * max(link[:, 2]),
               above_threshold_color=GREY)

    ax.set_ylabel("Euclidean distance (log2 CPM+1)")
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    add_data_note(ax, data_note)

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    save_fig(pdf, fig)


def plot_mean_cv(pdf, colnames, lcpm, cpm, geneids, n_genes_filt):
    """Page 10 — mean expression vs CV scatter with running-median floor."""
    data_note = f"Data: CPM-filtered genes (n={n_genes_filt}, log2 CPM+1)"
    G         = lcpm.T                             # genes x samples
    mean_expr = G.mean(axis=1)
    cv        = G.std(axis=1) / (G.mean(axis=1) + 1e-9)

    n_bins  = 30
    bin_idx = np.array_split(np.argsort(mean_expr), n_bins)
    bin_x   = np.array([mean_expr[b].mean() for b in bin_idx])
    bin_y   = np.array([np.median(cv[b])     for b in bin_idx])

    try:
        z     = gaussian_kde(np.vstack([mean_expr, cv]))(np.vstack([mean_expr, cv]))
        order = z.argsort()
    except Exception:
        order = np.arange(len(mean_expr))
        z     = np.ones(len(mean_expr))

    fig, ax = plt.subplots(figsize=(9, 6))
    add_page_title(fig, "Mean Expression vs Coefficient of Variation",
                   "Each point = one gene (CPM-filtered). "
                   "Colour = local density. "
                   "Curve = running median CV. "
                   "Labelled = top 20 high-CV genes (mean log2 CPM > 1).")

    sc = ax.scatter(mean_expr[order], cv[order], c=z[order],
                    s=3, cmap="viridis", alpha=0.5, linewidths=0, rasterized=True)
    ax.plot(bin_x, bin_y, color=RED, linewidth=1.8, label="Running median CV", zorder=5)

    expressed  = mean_expr > 1.0
    top_cv_idx = (np.where(expressed)[0][np.argsort(cv[expressed])[-20:]]
                  if expressed.sum() >= 20 else np.argsort(cv)[-20:])

    texts = []
    for idx in top_cv_idx:
        t = ax.text(mean_expr[idx], cv[idx], geneids[idx],
                    fontsize=6, color="#222222", ha="left", va="bottom")
        texts.append(t)
    adjust_text_labels(ax, texts)

    ax.set_xlabel("Mean log2 CPM+1")
    ax.set_ylabel("Coefficient of Variation (CV)")
    ax.legend(fontsize=8)
    fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02).set_label("Local density", fontsize=8)
    add_data_note(ax, data_note)

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(pdf, fig)


def plot_top_genes_table(pdf, colnames, cpm, geneids, n_top_genes=20):
    """Page 11 — top N most highly expressed genes."""
    n         = len(colnames)
    data_note = f"Data: CPM-filtered genes, linear CPM (top {n_top_genes} by mean CPM)"
    mean_cpm  = cpm.mean(axis=1)
    top_idx   = np.argsort(mean_cpm)[::-1][:n_top_genes]
    top_genes = geneids[top_idx]
    top_mean  = mean_cpm[top_idx]
    top_pct   = top_mean / (mean_cpm.sum() + 1e-9) * 100
    top_cpm_m = cpm[top_idx, :]

    # Per-sample columns only if they fit
    MAX_SAMPLE_COLS = min(n, 8)
    show_cols       = list(colnames[:MAX_SAMPLE_COLS])
    truncated       = n > MAX_SAMPLE_COLS

    col_labels = ["Rank", "Gene", "Mean CPM", "% of total"] + show_cols
    if truncated:
        col_labels[-1] = col_labels[-1] + f" ... (+{n - MAX_SAMPLE_COLS} more)"

    cell_data = []
    for rank, (gene, mn, pct, row) in enumerate(
            zip(top_genes, top_mean, top_pct, top_cpm_m), start=1):
        cell_data.append(
            [str(rank), gene, f"{mn:.1f}", f"{pct:.2f}%"]
            + [f"{v:.0f}" for v in row[:MAX_SAMPLE_COLS]]
        )

    n_rows  = len(cell_data)
    fig_h   = max(5, n_rows * 0.32 + 2.5)
    n_cols_t = len(col_labels)
    fig_w   = min(20, max(8, n_cols_t * 1.4))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    add_page_title(fig, f"Top {n_top_genes} Most Highly Expressed Genes",
                   "Ranked by mean CPM across all samples (CPM-filtered). "
                   "High % of total from a single gene may indicate contamination. "
                   + (f"Per-sample CPM shown for first {MAX_SAMPLE_COLS} samples." if truncated else ""))

    table = ax.table(cellText=cell_data, colLabels=col_labels,
                     loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    fs = max(5, min(8, int(72 / max(n_cols_t, 6))))
    table.set_fontsize(fs)
    table.scale(1, max(1.2, 1.5 - n_rows * 0.01))
    table.auto_set_column_width(list(range(n_cols_t)))

    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor(BLUE)
            cell.set_text_props(color="white", fontweight="bold")
            cell.set_edgecolor(DIVIDER)
        else:
            gene_pct = top_pct[r - 1] if r <= n_rows else 0
            if gene_pct > 5.0:
                cell.set_facecolor("#FDECEA")
            elif r % 2 == 0:
                cell.set_facecolor("#E8EEF5")
            else:
                cell.set_facecolor("white")
            cell.set_edgecolor(DIVIDER)

    fig.text(0.05, 0.01,
             "Light red rows: gene accounts for > 5% of total mean CPM.",
             fontsize=7, color="#AA2200", style="italic")
    fig.text(0.99, 0.01, data_note, ha="right", va="bottom",
             fontsize=7, color=GREY, style="italic")

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(pdf, fig)


def plot_consensus_clustering(pdf, colnames, lcpm, sample_colors, n_top, run_warnings):
    """Page 12 — unsupervised group estimation via consensus clustering."""
    n     = len(colnames)
    max_k = min(6, n - 1)

    if max_k < 2:
        run_warnings.append(
            "Consensus clustering requires at least 3 samples — skipping.")
        return

    gene_vars = lcpm.var(axis=0)
    X         = lcpm[:, np.argsort(gene_vars)[-n_top:]]
    data_note = f"Data: top {n_top} most variable genes (log2 CPM+1, CPM-filtered)"
    N_ITER    = 100
    K_RANGE   = range(2, max_k + 1)
    seeds     = np.random.default_rng(42).integers(0, 100_000, size=N_ITER)

    # ── Build co-occurrence matrices ─────────────────────────────────────────
    cooc   = {}
    labels = {}
    for k in K_RANGE:
        M = np.zeros((n, n), dtype=float)
        for seed in seeds:
            try:
                lab = KMeans(n_clusters=k, random_state=int(seed),
                             n_init=1, max_iter=100).fit_predict(X)
                for i in range(n):
                    for j in range(n):
                        if lab[i] == lab[j]:
                            M[i, j] += 1
            except Exception:
                pass
        M       /= N_ITER
        cooc[k]  = M
        try:
            labels[k] = KMeans(n_clusters=k, random_state=42,
                                n_init=20).fit_predict(M)
        except Exception:
            labels[k] = np.zeros(n, dtype=int)

    # ── Choose best k via CDF-area delta ────────────────────────────────────
    cdf_areas = {k: float(np.mean(cooc[k][np.triu_indices(n, k=1)]))
                 for k in K_RANGE}
    k_list    = sorted(K_RANGE)
    deltas    = {k_list[i]: cdf_areas[k_list[i]] - cdf_areas[k_list[i - 1]]
                 for i in range(1, len(k_list))}
    best_k    = max(deltas, key=deltas.get) if deltas else 2
    best_lbl  = labels[best_k]

    order      = np.argsort(best_lbl)
    names_ord  = [colnames[i] for i in order]
    cooc_ord   = cooc[best_k][np.ix_(order, order)]
    cp         = sns.color_palette("Set2", n_colors=best_k)

    # ── Tick font sizes for heatmaps ─────────────────────────────────────────
    fs_heat      = heatmap_tick_fontsize(n)
    fs_heat_sm   = max(0, fs_heat - 1)   # smaller for the bottom-row thumbnails

    # ── Figure ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 10))
    add_page_title(
        fig,
        f"Consensus Clustering — Estimated Groups (best k={best_k})",
        f"k-means run {N_ITER}x per k=2..{max_k}.  "
        f"Co-occurrence = fraction of runs two samples share a cluster.  "
        f"Best k = largest CDF-area increase.  No prior group info used.")

    gs = gridspec.GridSpec(2, 3, figure=fig,
                           hspace=0.50, wspace=0.40,
                           top=0.88, bottom=0.10, left=0.07, right=0.97)

    # ── Panel A: co-occurrence heatmap (best k) ──────────────────────────────
    ax_h = fig.add_subplot(gs[0, 0])
    im   = ax_h.imshow(cooc_ord, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    apply_heatmap_ticks(ax_h, names_ord)
    ax_h.set_title(f"Co-occurrence matrix (k={best_k})", fontweight="bold", fontsize=9)
    # Cell annotations only if large enough
    if n <= 25:
        for i in range(n):
            for j in range(n):
                v = cooc_ord[i, j]
                ax_h.text(j, i, f"{v:.2f}", ha="center", va="center",
                          fontsize=max(4, 5 - n // 10),
                          color="white" if v > 0.6 else "black")
    # Cluster boundaries
    for b in np.where(np.diff(best_lbl[order]))[0] + 1:
        ax_h.axhline(b - 0.5, color=RED, linewidth=1.2)
        ax_h.axvline(b - 0.5, color=RED, linewidth=1.2)
    fig.colorbar(im, ax=ax_h, fraction=0.046, pad=0.04).set_label(
        "Co-occurrence", fontsize=7)
    add_data_note(ax_h, data_note)

    # ── Panel B: stability bar chart ─────────────────────────────────────────
    ax_s = fig.add_subplot(gs[0, 1])
    ks   = sorted(cdf_areas.keys())
    ax_s.bar(ks, [cdf_areas[k] for k in ks],
             color=[RED if k == best_k else BLUE for k in ks],
             edgecolor="white", zorder=3, alpha=0.85)
    ax_s.set_xticks(ks)
    ax_s.set_xlabel("k", fontsize=9)
    ax_s.set_ylabel("Mean co-occurrence (CDF area)", fontsize=8)
    ax_s.set_title("Stability per k  (red = chosen)", fontweight="bold", fontsize=9)
    ax_s.yaxis.grid(True, zorder=0)
    ax_s.set_axisbelow(True)
    if deltas:
        ax_d = ax_s.twinx()
        dk   = sorted(deltas.keys())
        ax_d.plot(dk, [deltas[k] for k in dk], color=GREEN, marker="o",
                  markersize=5, linewidth=1.4, label="Delta")
        ax_d.set_ylabel("Delta CDF area", fontsize=8, color=GREEN)
        ax_d.tick_params(axis="y", colors=GREEN, labelsize=7)

    # ── Panel C: PCA coloured by consensus group ─────────────────────────────
    ax_p = fig.add_subplot(gs[0, 2])
    try:
        emb    = PCA(n_components=2).fit_transform(X)
        pct_p  = PCA(n_components=2).fit(X).explained_variance_ratio_ * 100
        c_cols = {s: cp[best_lbl[i]] for i, s in enumerate(colnames)}
        scatter_with_labels(ax_p, emb[:, 0], emb[:, 1], colnames, c_cols, s=60)
        ax_p.legend(handles=[Patch(facecolor=cp[k], label=f"Group {k+1}")
                              for k in range(best_k)],
                    fontsize=8, framealpha=0.8)
        ax_p.set_xlabel(f"PC1 ({pct_p[0]:.1f}%)", fontsize=9)
        ax_p.set_ylabel(f"PC2 ({pct_p[1]:.1f}%)", fontsize=9)
        ax_p.set_title("PCA — consensus groups", fontweight="bold", fontsize=9)
        ax_p.axhline(0, color=DIVIDER, linewidth=0.5)
        ax_p.axvline(0, color=DIVIDER, linewidth=0.5)
    except Exception as e:
        run_warnings.append(f"PCA panel in consensus plot failed: {e}")
        ax_p.set_visible(False)
    add_data_note(ax_p, data_note)

    # ── Bottom row: co-occurrence thumbnails for other k values ──────────────
    other_ks = [k for k in K_RANGE if k != best_k]
    for idx in range(3):
        ax_t = fig.add_subplot(gs[1, idx])
        if idx < len(other_ks):
            k       = other_ks[idx]
            ord_k   = np.argsort(labels[k])
            cooc_k  = cooc[k][np.ix_(ord_k, ord_k)]
            nms_k   = [colnames[i] for i in ord_k]
            ax_t.imshow(cooc_k, cmap="Blues", vmin=0, vmax=1, aspect="auto")
            apply_heatmap_ticks(ax_t, nms_k)
            ax_t.set_title(f"Co-occurrence k={k}", fontweight="bold", fontsize=9)
            for b in np.where(np.diff(labels[k][ord_k]))[0] + 1:
                ax_t.axhline(b - 0.5, color=RED, linewidth=1.0)
                ax_t.axvline(b - 0.5, color=RED, linewidth=1.0)
        else:
            ax_t.set_visible(False)

    # ── Group assignment summary text ─────────────────────────────────────────
    # Wrap long member lists so they don't run off the page
    max_chars   = 110
    group_lines = [f"Estimated groups (best k={best_k}):"]
    for g in range(best_k):
        members = [colnames[i] for i in range(n) if best_lbl[i] == g]
        line    = f"  Group {g+1}: {', '.join(members)}"
        # Hard-wrap if too long
        if len(line) > max_chars:
            line = f"  Group {g+1}: {', '.join(members[:len(members)//2])},\n" \
                   f"          {', '.join(members[len(members)//2:])}"
        group_lines.append(line)

    fig.text(0.01, 0.005, "\n".join(group_lines),
             fontsize=7, family="monospace", color="#222222", va="bottom")

    save_fig(pdf, fig)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    if len(sys.argv) < 3:
        print("Usage: python comprehensive_report_python.py <counts_file> <output_pdf>")
        sys.exit(1)

    counts_f   = sys.argv[1]
    output_pdf = sys.argv[2]
    os.makedirs(os.path.dirname(os.path.abspath(output_pdf)), exist_ok=True)

    run_warnings = []

    print("\nRNA-seq QC Report Generator")
    print("=" * 50)

    # 1. Load
    with log_step("Loading count matrix"):
        try:
            raw = pd.read_csv(counts_f, sep="\t", comment="#")
        except Exception as e:
            print(f"\n[ERROR] Could not read input file: {e}")
            sys.exit(1)

        raw_counts  = raw.iloc[:, 6:].to_numpy(dtype=float)
        geneids     = raw["Geneid"].values
        colnames    = [c.split("/")[-1].split("_Aligned")[0].rstrip(".bam")
                       for c in raw.columns[6:]]
        n_samples   = raw_counts.shape[1]
        n_genes_raw = raw_counts.shape[0]
        print(f"\n    Samples: {n_samples}  |  Genes: {n_genes_raw}")
        if n_samples < 2:
            run_warnings.append("Only 1 sample — several multi-sample plots will be skipped.")

    # 2. Pre-filter stats
    with log_step("Pre-filter stats"):
        n_genes_before = (raw_counts > 0).sum(axis=0).astype(int)

    # 3. CPM + filter
    with log_step("CPM normalisation and filtering"):
        lib_sizes_raw = raw_counts.sum(axis=0)
        if np.any(lib_sizes_raw == 0):
            run_warnings.append("One or more samples have zero total counts.")
        lib_sizes_raw = np.where(lib_sizes_raw == 0, 1, lib_sizes_raw)

        cpm_raw = raw_counts / lib_sizes_raw[np.newaxis, :] * 1e6
        keep    = (cpm_raw > 1).any(axis=1)
        if keep.sum() == 0:
            run_warnings.append("No genes passed CPM > 1 filter — skipping gene filter.")
            keep = np.ones(n_genes_raw, dtype=bool)

        counts_filt  = raw_counts[keep, :]
        geneids      = geneids[keep]
        lib_sizes    = counts_filt.sum(axis=0)
        lib_sizes    = np.where(lib_sizes == 0, 1, lib_sizes)
        cpm          = counts_filt / lib_sizes[np.newaxis, :] * 1e6
        lcpm         = np.log2(cpm + 1).T                     # samples x genes
        n_genes_filt = lcpm.shape[1]
        n_genes_after = (lcpm > 0).sum(axis=1).astype(int)
        n_top         = min(2000, n_genes_filt)

        print(f"\n    Genes after filter: {n_genes_filt}  |  Top-var used: {n_top}")
        if n_genes_filt < 10:
            run_warnings.append(
                f"Very few genes passed filter ({n_genes_filt}). Results may be unreliable.")

    # 4. Sample colours
    sample_colors = assign_colors(colnames)

    # 5. Build PDF
    print(f"\n  Writing PDF: {output_pdf}")
    with PdfPages(output_pdf) as pdf:
        meta            = pdf.infodict()
        meta["Title"]   = "RNA-seq QC Report"
        meta["Subject"] = counts_f

        embedding_step = (
            ("MA plots",
             lambda: plot_ma_plots(pdf, colnames, lcpm, n_genes_filt))
            if n_samples <= 10 else
            ("UMAP + k-means (k=2)",
             lambda: plot_umap_kmeans(
                 pdf, colnames, lcpm, sample_colors, n_top, run_warnings))
        )

        steps = [
            ("Summary table",
             lambda: plot_summary_table(
                 pdf, colnames, lib_sizes, n_genes_before, n_genes_after, run_warnings)),
            ("Library sizes",
             lambda: plot_library_sizes(pdf, colnames, lib_sizes, sample_colors)),
            ("CPM distributions",
             lambda: plot_cpm_distributions(
                 pdf, colnames, lcpm, sample_colors, n_genes_filt)),
            ("Detected genes",
             lambda: plot_detected_genes(
                 pdf, colnames, n_genes_before, n_genes_after, sample_colors)),
            embedding_step,
            ("Sample correlation",
             lambda: plot_sample_correlation(
                 pdf, colnames, lcpm, sample_colors, n_genes_filt)),
            ("PCA",
             lambda: plot_pca(pdf, colnames, lcpm, sample_colors, n_top)),
            ("Outlier overview",
             lambda: plot_outlier_summary(
                 pdf, colnames, lcpm, lib_sizes, n_genes_after, n_genes_filt)),
            ("Hierarchical clustering",
             lambda: plot_hierarchical_clustering(pdf, colnames, lcpm, n_genes_filt)),
            ("Mean vs CV",
             lambda: plot_mean_cv(pdf, colnames, lcpm, cpm, geneids, n_genes_filt)),
           # ("Top expressed genes",
           #  lambda: plot_top_genes_table(pdf, colnames, cpm, geneids)),
            ("Consensus clustering",
             lambda: plot_consensus_clustering(
                 pdf, colnames, lcpm, sample_colors, n_top, run_warnings)),
        ]

        for step_name, fn in steps:
            with log_step(step_name):
                try:
                    fn()
                except Exception as e:
                    run_warnings.append(f"Plot '{step_name}' failed: {e}")
                    warn(f"'{step_name}' failed — skipping.\n"
                         f"    {traceback.format_exc(limit=3)}")

    print("\n" + "=" * 50)
    if run_warnings:
        print("Warnings encountered:")
        for w in run_warnings:
            print(f"  * {w}")
    print(f"\nDone. Report saved to: {output_pdf}\n")


if __name__ == "__main__":
    main()