# Local 10x PBMC 3k (filtered_gene_bc_matrices) through the app's own
# cell-level path, with the QC violin and the expression heatmaps.
suppressMessages(devtools::load_all("c:/Users/Lenovo/Desktop/PERCEPTION", quiet = TRUE))
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

root <- "C:/Users/Lenovo/Desktop/python_class/practical_course/class8/filtered_gene_bc_matrices/hg19"
counts <- as(readMM(file.path(root, "matrix.mtx")), "CsparseMatrix")
genes <- read.delim(file.path(root, "genes.tsv"), header = FALSE)
rownames(counts) <- genes$V2
colnames(counts) <- readLines(file.path(root, "barcodes.tsv"))
counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
cat("input:", nrow(counts), "genes x", ncol(counts), "cells\n")

patient <- stats::setNames(rep("PBMC", ncol(counts)), colnames(counts))
t0 <- proc.time()[3]
invisible(capture.output(
  lab <- suppressMessages(suppressWarnings(PERCEPTIONx::preprocess_sc(
    counts = counts, patient_labels = patient,
    nmads = 3, min_genes = 200, max_mito = 25,
    nfeatures = 2000, dims = 20, resolution = 0.8, method = "umap"
  )))
))
cat(sprintf("preprocess_sc: %.1fs\n", proc.time()[3] - t0))

cl <- lab$cell_labels
qm <- lab$qc_metrics
qc <- lab$qc
kept <- cl$cell_id
cat(sprintf("kept %d / %d cells (%.1f%%)\n", length(kept), ncol(counts),
            100 * length(kept) / ncol(counts)))
cat(sprintf("mito fence %.1f%% | gene ceiling %.0f | clusters %d\n",
            qc$mito_threshold[1], qc$genes_max[1], length(unique(cl$clone_id))))
print(sort(table(cl$clone_id), decreasing = TRUE))

out <- "C:/Users/Lenovo/Desktop/data-PERCEPTION/qc_bench/pbmc_local"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- 1. QC violin (pooled, the app's own plot) -----------------------------
p_qc <- PERCEPTIONx::plot_qc_report(
  qm, kept_cells = kept, qc_report = qc, metric = "all"
)
ggsave(file.path(out, "01_qc_violin.png"), p_qc, width = 6.4, height = 8.4, dpi = 200)

p_gene <- PERCEPTIONx::plot_qc_report(
  qm, kept_cells = kept, qc_report = qc, metric = "Detected genes"
)
ggsave(file.path(out, "02_qc_violin_genes.png"), p_gene, width = 6, height = 4, dpi = 200)

# ---- 2. marker x cluster heatmap ------------------------------------------
markers <- c(
  "CD3D", "CD3E", "CD8A", "MS4A1", "CD79A", "LST1", "CD14",
  "FCGR3A", "NKG7", "GNLY", "PPBP", "HBB"
)
markers <- intersect(markers, rownames(counts))
logc <- as.matrix(log1p(
  t(t(counts[markers, , drop = FALSE]) / Matrix::colSums(counts) * 1e4)
))
logc <- logc[, kept, drop = FALSE]
cell_cl <- cl$clone_id[match(kept, cl$cell_id)]

cl_mean <- t(vapply(
  split(seq_along(kept), cell_cl),
  function(i) rowMeans(logc[, i, drop = FALSE]),
  numeric(length(markers))
))
cl_n <- table(cell_cl)
df_mean <- data.frame(
  cluster = rep(rownames(cl_mean), times = length(markers)),
  marker = rep(colnames(cl_mean), each = nrow(cl_mean)),
  value = as.vector(cl_mean)
)
df_mean$cluster <- factor(df_mean$cluster, levels = names(sort(cl_n, decreasing = TRUE)))
df_mean$marker <- factor(df_mean$marker, levels = markers)

p_mean <- ggplot(df_mean, aes(cluster, marker, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", value)), size = 2.6, colour = "grey15") +
  scale_fill_gradientn(
    colours = c("#F7F9FC", "#93B5D0", "#165C91"),
    name = "log expr"
  ) +
  labs(x = "cluster (size in cells)", y = NULL,
       title = "PBMC 3k - cluster mean marker expression") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(face = "bold"))
ggsave(file.path(out, "03_cluster_marker_heatmap.png"), p_mean,
       width = 8, height = 5, dpi = 220)

# ---- 3. cell x marker heatmap (ordered by cluster) ------------------------
set.seed(3)
per_cl <- 70
idx <- unlist(lapply(names(sort(cl_n, decreasing = TRUE)), function(k) {
  ii <- which(cell_cl == k)
  head(sample(ii, min(per_cl, length(ii))), per_cl)
}))
sub <- logc[, idx, drop = FALSE]
df_cell <- data.frame(
  x = rep(seq_along(idx), each = length(markers)),
  marker = factor(rep(markers, times = length(idx)), levels = rev(markers)),
  value = as.vector(sub),
  cluster = rep(cell_cl[idx], each = length(markers))
)
df_cell$cluster <- factor(df_cell$cluster,
                          levels = names(sort(cl_n, decreasing = TRUE)))

p_cell <- ggplot(df_cell, aes(x, marker, fill = value)) +
  geom_raster() +
  facet_grid(~cluster, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#F7F9FC", "#93B5D0", "#165C91"),
    limits = c(0, 4), oob = scales::squish, name = "log expr"
  ) +
  labs(x = NULL, y = NULL,
       title = "PBMC 3k - per-cell marker expression (<=70 cells/cluster)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid = element_blank(), panel.spacing.x = unit(0.12, "lines"),
    strip.text = element_text(size = 8), plot.title = element_text(face = "bold")
  )
ggsave(file.path(out, "04_cell_marker_heatmap.png"), p_cell,
       width = 13, height = 5, dpi = 220)

# ---- 4. UMAP: clusters and mitochondrial fraction -------------------------
dim_df <- data.frame(
  cell_id = cl$cell_id, dim1 = cl$dim_1, dim2 = cl$dim_2,
  cluster = factor(cl$clone_id)
)
mt <- qm$pct_mito[match(dim_df$cell_id, qm$cell_id)]
dim_df$mito <- mt

p_umap <- ggplot(dim_df, aes(dim1, dim2, colour = cluster)) +
  geom_point(size = 0.55, alpha = 0.85) +
  labs(title = "PBMC 3k - clusters after QC", x = "UMAP 1", y = "UMAP 2",
       colour = "cluster") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out, "05_umap_clusters.png"), p_umap, width = 7, height = 5.2, dpi = 220)

p_mito <- ggplot(dim_df[order(dim_df$mito), ], aes(dim1, dim2, colour = mito)) +
  geom_point(size = 0.55, alpha = 0.85) +
  scale_colour_gradientn(colours = c("#E9EEF5", "#93B5D0", "#c13232"),
                         name = "mito %") +
  labs(title = "PBMC 3k - mitochondrial fraction", x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out, "06_umap_mito.png"), p_mito, width = 7, height = 5.2, dpi = 220)

# ---- 5. tables -------------------------------------------------------------
rep <- qc
rep$n_clusters <- length(unique(cl$clone_id))
write.csv(rep, file.path(out, "qc_report.csv"), row.names = FALSE)
write.csv(
  data.frame(cluster = rownames(cl_mean), n_cells = as.integer(cl_n),
             round(cl_mean, 3), check.names = FALSE),
  file.path(out, "cluster_marker_means.csv"), row.names = FALSE
)
cat("\nsaved to", out, "\n")
print(list.files(out))
