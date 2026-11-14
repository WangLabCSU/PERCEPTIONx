# Stress test: a sample whose mitochondrial distribution is shifted far above
# the PBMC baseline (what necrotic tumour tissue looks like). The counts are
# real PBMC 3k; only the mitochondrial fraction is perturbed, so the comparison
# isolates the threshold policy rather than the biology.
suppressMessages(devtools::load_all("c:/Users/Lenovo/Desktop/PERCEPTION", quiet = TRUE))
suppressPackageStartupMessages(library(Matrix))

root <- "C:/Users/Lenovo/AppData/Local/Temp/pbmc3k_raw/raw_gene_bc_matrices/hg19"
counts <- as(readMM(file.path(root, "matrix.mtx")), "CsparseMatrix")
genes <- read.delim(file.path(root, "genes.tsv"), header = FALSE)
rownames(counts) <- genes$V2
colnames(counts) <- readLines(file.path(root, "barcodes.tsv"))
counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
counts <- counts[, Matrix::colSums(counts) >= 100, drop = FALSE]
cat("starting cells (>=100 UMI):", ncol(counts), "\n")

is_mito <- grepl("^MT-", toupper(rownames(counts)))
mt_idx <- which(is_mito)

# ---- shift the mitochondrial fraction to a necrotic-tumour profile ----------
set.seed(11)
target <- pmin(stats::rgamma(ncol(counts), shape = 2.2, scale = 9), 70)
lib <- Matrix::colSums(counts)
new_mito <- round(lib * target / 100)

counts_stress <- counts
mt_prob <- rowSums(as.matrix(counts[mt_idx, , drop = FALSE])) + 1
for (j in seq_len(ncol(counts_stress))) {
  counts_stress[mt_idx, j] <- stats::rmultinom(1, new_mito[j], prob = mt_prob)
}
counts_stress <- as(counts_stress, "CsparseMatrix")

pct <- function(m) {
  mc <- Matrix::colSums(m[is_mito, , drop = FALSE])
  as.numeric(mc) / pmax(as.numeric(Matrix::colSums(m)), 1) * 100
}
mp <- pct(counts_stress)
cat(sprintf(
  "mito%% after perturbation: median %.1f  IQR %.1f-%.1f  p90 %.1f  max %.1f\n",
  stats::median(mp), stats::quantile(mp, .25), stats::quantile(mp, .75),
  stats::quantile(mp, .9), max(mp)
))

patient <- stats::setNames(rep("TUMOR", ncol(counts_stress)), colnames(counts_stress))
policies <- list(
  hybrid = list(nmads = 3, max_mito = 25),
  mad_only = list(nmads = 3, max_mito = 100),
  fixed_only = list(nmads = 50, max_mito = 25)
)
cat("\n== policy behaviour on the stressed sample ==\n")
for (nm in names(policies)) {
  p <- policies[[nm]]
  r <- PERCEPTIONx:::qc_filter_cells(
    counts_stress, patient,
    nmads = p$nmads, min_genes = 200, max_mito = p$max_mito
  )
  k <- colnames(counts_stress) %in% r$keep
  cat(sprintf(
    "%-10s kept %4d (%.1f%%) | mito fence %5s%% | kept mito>25%%: %3d | dropped with mito<10%%: %3d | kept mito median %.1f%%\n",
    nm, sum(k), 100 * mean(k), round(r$report$mito_threshold, 1),
    sum(k & mp > 25), sum(!k & mp < 10), stats::median(mp[k])
  ))
}
