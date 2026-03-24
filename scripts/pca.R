# scripts/pca_umap.R
# Dependencies: r-matrixstats, r-ggplot2, r-ggrepel, r-uwot  (all CRAN, no Bioconductor)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(uwot)
})
options(warn=-1)
args     <- commandArgs(trailingOnly = TRUE)
counts_f <- args[1]
outdir   <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load and clean featureCounts output ────────────────────────────────
raw    <- read.delim(counts_f, comment.char = "#", check.names = FALSE)
counts <- as.matrix(raw[, 7:ncol(raw)])
rownames(counts) <- raw$Geneid

# Strip full BAM paths → clean sample names
colnames(counts) <- gsub(".*/", "", colnames(counts))
colnames(counts) <- gsub("_Aligned.*\\.bam$", "", colnames(counts))
colnames(counts) <- gsub("\\.bam$", "", colnames(counts))

n_samples <- ncol(counts)
message("Samples: ", n_samples)
message("Sample names: ", paste(colnames(counts), collapse = ", "))

# ── 2. Filter lowly expressed genes, then log2(CPM+1) normalise ───────────
# Keep genes with CPM > 1 in at least 2 samples — removes noise without
# requiring any group info
lib_sizes <- colSums(counts)
cpm       <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
keep      <- rowSums(cpm > 1) >= min(2, n_samples)
counts    <- counts[keep, ]
message("Genes after CPM filter: ", nrow(counts))

# Recalculate CPM on filtered matrix and log-transform
lib_sizes <- colSums(counts)
cpm       <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
lcpm      <- log2(cpm + 1)               # genes × samples
mat       <- t(lcpm)                      # samples × genes for prcomp / umap

# ── 3. Helper: consistent ggplot theme ────────────────────────────────────
base_theme <- function() {
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(size = 13, face = "plain"))
}

save_plot <- function(p, name, w = 10, h = 8) {
  ggsave(file.path(outdir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(outdir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
}

# ── 4. PCA ────────────────────────────────────────────────────────────────
# Select top 500 most variable genes (by row variance) before PCA
gene_vars <- apply(lcpm, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(500, length(gene_vars))])
pca_res   <- prcomp(t(lcpm[top_genes, ]), center = TRUE, scale. = FALSE)

pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2))
pca_df  <- data.frame(PC1    = pca_res$x[, 1],
                      PC2    = pca_res$x[, 2],
                      sample = rownames(pca_res$x))

p_pca <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3, colour = "#3B82C4") +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 20, seed = 42) +
  labs(title = "PCA — top 500 variable genes (log2 CPM+1)",
       x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance")) +
  base_theme()

save_plot(p_pca, "pca")
message("PCA written.")

# ── 5. UMAP ───────────────────────────────────────────────────────────────
if (n_samples < 4) {
  message("WARNING: Only ", n_samples, " samples — skipping UMAP (need >= 4).")
} else {
  nn     <- min(15L, n_samples - 1L)
  n_pcs  <- min(50L, n_samples - 1L)
  pc_mat <- pca_res$x[, 1:n_pcs, drop = FALSE]

  set.seed(42)
  umap_coords <- uwot::umap(pc_mat,
                             n_neighbors = nn,
                             min_dist    = 0.3,
                             metric      = "euclidean",
                             n_epochs    = 500,
                             verbose     = FALSE)

  umap_df <- data.frame(UMAP1  = umap_coords[, 1],
                        UMAP2  = umap_coords[, 2],
                        sample = colnames(counts))

  p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2)) +
    geom_point(size = 3, colour = "#E07B39") +
    geom_text_repel(aes(label = sample), size = 3, max.overlaps = 20, seed = 42) +
    labs(title = paste0("UMAP — log2 CPM+1 (n_neighbors=", nn, ", min_dist=0.3)"),
         x = "UMAP 1", y = "UMAP 2") +
    base_theme()

  save_plot(p_umap, "umap")
  message("UMAP written.")
}

message("Done. Output in: ", outdir)