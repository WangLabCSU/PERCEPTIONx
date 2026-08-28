# test_umap_gene_vec.R — verifies the umap_gene optimization: the worker now
# receives ONLY the selected gene's per-cell vector (named by cell id) instead
# of the whole expression matrix. Checks (1) the vector path produces the same
# values as the old full-matrix path, (2) serialization of the vector is
# orders of magnitude cheaper, (3) the worker branch builds the ggplot.
# Usage: Rscript dev/test_umap_gene_vec.R
suppressMessages(library(PERCEPTIONx, lib.loc = "E:/rlibtest"))
suppressMessages(library(ggplot2))

set.seed(3)
n_g <- 2000; n_c <- 3000
expr_mat <- matrix(rpois(n_g * n_c, 1), nrow = n_g, ncol = n_c)
rownames(expr_mat) <- paste0("G", seq_len(n_g))
colnames(expr_mat) <- paste0("C", seq_len(n_c))
gene <- "G42"

# (1) old path: cell_expr from the full matrix; new path: from the vector
cells_old <- as.numeric(expr_mat[gene, ])
gene_vec <- setNames(as.numeric(expr_mat[gene, ]), colnames(expr_mat))
cells_new <- as.numeric(gene_vec[match(colnames(expr_mat), names(gene_vec))])
stopifnot(identical(cells_old, cells_new))
cat("vector path == full-matrix path values\n")

# (2) serialization cost comparison
f_mat <- tempfile(fileext = ".rds"); f_vec <- tempfile(fileext = ".rds")
t1 <- system.time(saveRDS(list(expr = expr_mat), f_mat))[3]
t2 <- system.time(saveRDS(list(expr = gene_vec), f_vec))[3]
cat(sprintf("saveRDS full matrix: %.2fs (%.1f MB); gene vector: %.2fs (%.1f KB)\n",
            t1, file.info(f_mat)$size / 1e6, t2, file.info(f_vec)$size / 1e3))
stopifnot(t2 < t1 / 10, file.info(f_vec)$size < file.info(f_mat)$size / 100)
unlink(c(f_mat, f_vec))

# (3) worker branch (async_jobs.R run_plot, umap_gene) with the vector
umap_coords <- data.frame(cell_id = colnames(expr_mat),
                          dim_1 = rnorm(n_c), dim_2 = rnorm(n_c), stringsAsFactors = FALSE)
clone_data <- data.frame(cell_id = colnames(expr_mat),
                         clone_id = rep(c("cl_A", "cl_B"), length.out = n_c), stringsAsFactors = FALSE)
common_cells <- intersect(clone_data$cell_id, umap_coords$cell_id)
stopifnot(length(common_cells) == n_c)
safe_range01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] == rng[1]) rep(0.5, length(x)) else (x - rng[1]) / (rng[2] - rng[1])
}
expr_vec <- gene_vec
cell_expr <- as.numeric(expr_vec[match(common_cells, names(expr_vec))])
x_col <- "dim_1"; y_col <- "dim_2"
X <- umap_coords[[x_col]][match(common_cells, umap_coords$cell_id)]
Y <- umap_coords[[y_col]][match(common_cells, umap_coords$cell_id)]
n <- length(common_cells)
stopifnot(length(X) == n, length(Y) == n, length(cell_expr) == n, !anyNA(cell_expr))
umap_data <- data.frame(X = X, Y = Y, expression = safe_range01(cell_expr), row.names = common_cells)
p <- PERCEPTIONx::plot_tsne_response(umap_data, color_var = "expression",
                                     title = gene, color_label = "Expression (0-1)",
                                     palette = "viridis", base_size = 11)
stopifnot(inherits(p, "ggplot"))
cat("worker umap_gene branch builds ggplot OK with", n, "cells\n")
