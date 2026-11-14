# Raw 10x PBMC 3k (no preprocessing) through PERCEPTIONx::preprocess_sc().
# Reports what QC removes, how the mito fence behaves, and whether the
# surviving cells form recognisable PBMC populations.
suppressMessages(devtools::load_all("c:/Users/Lenovo/Desktop/PERCEPTION", quiet = TRUE))
suppressPackageStartupMessages(library(Matrix))

root <- "C:/Users/Lenovo/AppData/Local/Temp/pbmc3k_raw/raw_gene_bc_matrices/hg19"
if (!file.exists(file.path(root, "matrix.mtx"))) {
  dest <- "C:/Users/Lenovo/AppData/Local/Temp/pbmc3k_raw.tar.gz"
  if (!file.exists(dest)) {
    old <- getOption("timeout")
    options(timeout = 1200)
    on.exit(options(timeout = old), add = TRUE)
    download.file(
      "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_raw_gene_bc_matrices.tar.gz",
      dest, quiet = TRUE, mode = "wb"
    )
  }
  utils::untar(dest, exdir = "C:/Users/Lenovo/AppData/Local/Temp/pbmc3k_raw")
}
counts <- readMM(file.path(root, "matrix.mtx"))
genes <- read.delim(file.path(root, "genes.tsv"), header = FALSE)
barcodes <- readLines(file.path(root, "barcodes.tsv"))
counts <- as(counts, "CsparseMatrix")
rownames(counts) <- genes$V2
colnames(counts) <- barcodes
# duplicated gene symbols are common in 10x references: keep the first
counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
# Stays sparse: the raw matrix is 32k x 737k barcodes, and densifying it is
# 179 GB. Only the ~2.7k cells that survive QC are ever embedded.
cat("raw matrix:", nrow(counts), "genes x", ncol(counts), "barcodes\n")
cat("total UMI:", format(sum(counts), big.mark = ","), "\n")
n_umi_bc <- Matrix::colSums(counts)
cat(sprintf(
  "UMI per barcode: <10: %.1f%%  <200: %.1f%%  >=200: %.1f%%\n",
  100 * mean(n_umi_bc < 10), 100 * mean(n_umi_bc < 200), 100 * mean(n_umi_bc >= 200)
))

patient <- rep("PBMC", ncol(counts))
names(patient) <- colnames(counts)

qc <- PERCEPTIONx:::qc_filter_cells(counts, patient, nmads = 3, min_genes = 200, max_mito = 25)
kept <- colnames(counts) %in% qc$keep
m <- qc$metrics

q <- function(x, p) round(stats::quantile(x, p, names = FALSE), 1)
cat("\n== QC (3 MAD, min 200 genes, mito cap 25%) ==\n")
cat(sprintf("input %d -> kept %d (%.1f%% removed)\n", ncol(counts), sum(kept), 100 * mean(!kept)))
cat("removed by reason:\n")
print(qc$report[, c("n_input", "n_kept", "pct_removed", "main_reason", "mito_threshold", "genes_max")])
cat(sprintf(
  "barcode-level genes/cell: median %s  p99 %s  max %s   (empty droplets are %.1f%% of input)\n",
  q(m$n_genes, .5), q(m$n_genes, .99), max(m$n_genes), 100 * mean(m$n_genes < 200)
))
cat(sprintf(
  "kept cells: genes median %s IQR %s-%s | counts median %s | mito%% median %s p95 %s\n",
  q(m$n_genes[kept], .5), q(m$n_genes[kept], .25), q(m$n_genes[kept], .75),
  format(q(m$n_counts[kept], .5), big.mark = ","),
  q(m$pct_mito[kept], .5), q(m$pct_mito[kept], .95)
))
cat(sprintf(
  "removed:    genes median %s | mito%% median %s p95 %s\n",
  q(m$n_genes[!kept], .5), q(m$pct_mito[!kept], .5), q(m$pct_mito[!kept], .95)
))

# ---- three mito policies on the same matrix --------------------------------
policies <- list(
  hybrid = list(nmads = 3, max_mito = 25),
  mad_only = list(nmads = 3, max_mito = 100),
  fixed_only = list(nmads = 50, max_mito = 25)
)
cat("\n== mito policy comparison ==\n")
for (nm in names(policies)) {
  p <- policies[[nm]]
  r <- PERCEPTIONx:::qc_filter_cells(counts, patient, nmads = p$nmads,
                                     min_genes = 200, max_mito = p$max_mito)
  k <- colnames(counts) %in% r$keep
  cat(sprintf(
    "%-10s kept %4d (%.1f%%)  of kept mito>25%%: %4d  of removed mito>25%%: %4d  mito fence %s%%\n",
    nm, sum(k), 100 * mean(k), sum(k & m$pct_mito > 25), sum(!k & m$pct_mito > 25),
    round(r$report$mito_threshold, 1)
  ))
}

# ---- full pipeline, filtered -------------------------------------------------
run_embed <- function(mat, pat, label, nmads, min_genes, max_mito) {
  t0 <- proc.time()[3]
  invisible(capture.output(
    lab <- suppressMessages(suppressWarnings(PERCEPTIONx::preprocess_sc(
      counts = mat, patient_labels = pat,
      nmads = nmads, min_genes = min_genes, max_mito = max_mito,
      nfeatures = 2000, dims = 20, resolution = 0.8, method = "umap"
    )))
  ))
  cl <- lab$cell_labels
  sizes <- sort(table(cl$clone_id), decreasing = TRUE)
  cat(sprintf(
    "\n[%s] %.1fs | %d cells | %d clusters | biggest %.0f%% of cells | singletons %d\n",
    label, proc.time()[3] - t0, nrow(cl), length(sizes),
    100 * max(sizes) / nrow(cl), sum(sizes == 1)
  ))
  print(sizes)
  lab
}

filtered <- run_embed(counts, patient, "filtered", 3, 200, 25)

# "Before QC": the 737k raw barcodes are mostly empty droplets, so they cannot
# all be embedded. Take every barcode with >= 100 UMI (cells AND dying cells)
# and embed it with the QC rules effectively switched off.
keep_raw <- Matrix::colSums(counts) >= 100
counts_raw <- counts[, keep_raw, drop = FALSE]
patient_raw <- stats::setNames(rep("PBMC", ncol(counts_raw)), colnames(counts_raw))
cat(sprintf("\nno-QC subset (>=100 UMI): %d barcodes\n", ncol(counts_raw)))
raw_embed <- run_embed(counts_raw, patient_raw, "no QC", 50, 0, 101)

# ---- biological sanity of the clusters --------------------------------------
markers <- c(
  "CD3D" = "T", "CD8A" = "CD8 T", "MS4A1" = "B", "CD79A" = "B",
  "LST1" = "mono/DC", "CD14" = "mono", "FCGR3A" = "NK/mono", "NKG7" = "NK",
  "GNLY" = "NK", "PPBP" = "platelet", "HBB" = "erythroid"
)
present <- intersect(names(markers), rownames(counts))
logcounts <- log1p(as.matrix(
  t(t(counts[present, , drop = FALSE]) / Matrix::colSums(counts) * 10000)
))
cl <- filtered$cell_labels
cat("\n== marker mean expression per cluster (log-normalised, filtered data) ==\n")
tab <- sapply(split(seq_len(nrow(cl)), cl$clone_id), function(idx) {
  ids <- cl$cell_id[idx]
  rowMeans(logcounts[, ids, drop = FALSE])
})
tab <- t(tab)
colnames(tab) <- paste0(present, "(", markers[present], ")")
print(round(tab, 2))
cat("\ncluster size / mito%:\n")
mt <- filtered$qc_metrics
print(data.frame(
  cluster = as.integer(names(sizes <- table(cl$clone_id))),
  n = as.integer(sizes),
  mito_median = round(tapply(mt$pct_mito[match(cl$cell_id, mt$cell_id)], cl$clone_id, median), 1)
))

# ---- figures ---------------------------------------------------------------
contam <- function(lab, label) {
  ids <- lab$cell_labels$cell_id
  mq <- m[match(ids, m$cell_id), ]
  cat(sprintf(
    "[%s] embedded %d cells | mito>25%%: %.1f%% | mito>10%%: %.1f%% | median mito %.1f%%\n",
    label, length(ids), 100 * mean(mq$pct_mito > 25),
    100 * mean(mq$pct_mito > 10), stats::median(mq$pct_mito)
  ))
}
cat("\n== high-mito contamination in each embedding ==\n")
contam(filtered, "filtered")
contam(raw_embed, "no QC")

out <- "C:/Users/Lenovo/Desktop/data-PERCEPTION/qc_bench"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(ggplot2)
plot_umap <- function(lab, title, file) {
  d <- lab$cell_labels
  d$cluster <- factor(d$clone_id)
  p <- ggplot(d, aes(dim_1, dim_2, colour = cluster)) +
    geom_point(size = 0.5, alpha = 0.8) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = "cluster") +
    theme_bw(base_size = 11) + theme(legend.position = "right")
  ggsave(file.path(out, file), p, width = 7.5, height = 5.5, dpi = 200)
}
plot_umap(raw_embed, "Raw PBMC 3k — no QC", "10_raw_umap.png")
plot_umap(filtered, "PBMC 3k — after preprocess_sc QC", "11_filtered_umap.png")
write.csv(qc$report, file.path(out, "10_qc_report.csv"), row.names = FALSE)
saveRDS(list(filtered = filtered, raw = raw_embed, qc = qc), file.path(out, "10_pbmc_result.rds"))
cat("\nsaved to", out, "\n")
