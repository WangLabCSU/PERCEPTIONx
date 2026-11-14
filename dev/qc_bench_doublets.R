# Head-to-head: does the built-in "minimum viable" QC agree with a dedicated
# doublet detector? PBMC 3k raw -> preprocess_sc QC -> scDblFinder, then compare
# which cells each removes and whether the clustering changes.
suppressMessages(devtools::load_all("c:/Users/Lenovo/Desktop/PERCEPTION", quiet = TRUE))
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

root <- "C:/Users/Lenovo/AppData/Local/Temp/pbmc3k_raw/raw_gene_bc_matrices/hg19"
counts <- as(readMM(file.path(root, "matrix.mtx")), "CsparseMatrix")
genes <- read.delim(file.path(root, "genes.tsv"), header = FALSE)
rownames(counts) <- genes$V2
colnames(counts) <- readLines(file.path(root, "barcodes.tsv"))
counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
counts <- counts[, Matrix::colSums(counts) >= 100, drop = FALSE]
patient <- stats::setNames(rep("PBMC", ncol(counts)), colnames(counts))
cat("candidate cells:", ncol(counts), "\n")

# --- A. built-in QC ---------------------------------------------------------
qc <- PERCEPTIONx:::qc_filter_cells(counts, patient, nmads = 3, min_genes = 200, max_mito = 25)
qc_keep <- colnames(counts) %in% qc$keep
cat(sprintf("PERCEPTIONx QC kept %d (%.1f%%)\n", sum(qc_keep), 100 * mean(qc_keep)))

# --- B. scDblFinder on the same candidates ---------------------------------
sce <- SingleCellExperiment(assays = list(counts = counts))
set.seed(1)
sce <- suppressMessages(scDblFinder::scDblFinder(sce, verbose = FALSE))
is_dbl <- sce$scDblFinder.class == "doublet"
cat(sprintf("scDblFinder doublets: %d (%.1f%%)\n", sum(is_dbl), 100 * mean(is_dbl)))

tab <- table(
  `PERCEPTIONx QC` = ifelse(qc_keep, "kept", "removed"),
  `scDblFinder` = ifelse(is_dbl, "doublet", "singlet")
)
cat("\n== overlap ==\n")
print(tab)
cat(sprintf(
  "doublet rate among removed cells: %.1f%% | among kept cells: %.1f%%\n",
  100 * mean(is_dbl[!qc_keep]), 100 * mean(is_dbl[qc_keep])
))
cat(sprintf(
  "of all doublets, %.1f%% were already removed by the built-in QC\n",
  100 * mean(!qc_keep[is_dbl])
))

# --- C. does removing doublets change the clustering? ----------------------
run_cl <- function(mat, pat, label) {
  invisible(capture.output(
    lab <- suppressMessages(suppressWarnings(PERCEPTIONx::preprocess_sc(
      counts = mat, patient_labels = pat,
      nmads = 3, min_genes = 200, max_mito = 25,
      nfeatures = 2000, dims = 20, resolution = 0.8, method = "umap"
    )))
  ))
  cl <- lab$cell_labels
  cat(sprintf("[%s] %d cells, %d clusters\n", label, nrow(cl), length(unique(cl$clone_id))))
  cl
}
cl_base <- run_cl(counts, patient, "built-in QC only")

mat2 <- counts[, !is_dbl, drop = FALSE]
pat2 <- stats::setNames(rep("PBMC", ncol(mat2)), colnames(mat2))
cl_dbl <- run_cl(mat2, pat2, "built-in QC + scDblFinder")

# ARI on the cells both runs share
ari <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  ij <- sum(choose(tab, 2))
  ai <- sum(choose(rowSums(tab), 2))
  bj <- sum(choose(colSums(tab), 2))
  expct <- ai * bj / choose(n, 2)
  maxi <- (ai + bj) / 2
  (ij - expct) / (maxi - expct)
}
common <- intersect(cl_base$cell_id, cl_dbl$cell_id)
a <- cl_base$clone_id[match(common, cl_base$cell_id)]
b <- cl_dbl$clone_id[match(common, cl_dbl$cell_id)]
cat(sprintf(
  "\ncluster agreement on %d shared cells: ARI = %.3f | same label %.1f%%\n",
  length(common), ari(a, b), 100 * mean(a == b)
))

# --- D. what do the doublets look like in the built-in metrics? ------------
mm <- qc$metrics[match(colnames(counts), qc$metrics$cell_id), ]
cat("\n== QC metrics by scDblFinder class ==\n")
for (k in c("singlet", "doublet")) {
  idx <- is_dbl == (k == "doublet")
  cat(sprintf(
    "%-8s n=%4d | genes median %4.0f p95 %4.0f | mito%% median %.2f\n",
    k, sum(idx), stats::median(mm$n_genes[idx]),
    stats::quantile(mm$n_genes[idx], .95), stats::median(mm$pct_mito[idx])
  ))
}
