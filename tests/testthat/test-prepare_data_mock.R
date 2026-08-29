# =============================================================================
# PERCEPTIONx prepare_data smoke test on a tiny mock expression matrix
# -----------------------------------------------------------------------------
# Runs the REAL Seurat clustering + clone annotation + rank normalization path
# of prepare_data() on a 60-gene x 12-cell mock matrix, so the annotation
# module executes on CI (the real patient scRNA data cannot ship with tests).
# =============================================================================

test_that("prepare_data runs the full clustering pipeline on mock data", {
  skip_if_not_installed("Seurat")

  set.seed(7)
  genes <- sprintf("GENE%03d", seq_len(60))
  cells <- paste0("C", seq_len(30))
  expr <- matrix(rpois(60 * 30, lambda = 5), nrow = 60,
                 dimnames = list(genes, cells))

  patient_mapping <- data.frame(
    cell_id    = cells,
    patient_id = c(rep("P1", 15), rep("P2", 15)),
    stringsAsFactors = FALSE
  )

  prepared <- prepare_data(
    method = "umap",
    expression_matrix = expr,
    patient_mapping = patient_mapping,
    seurat_resolution = 0.5,
    seurat_dims = 5
  )

  expect_type(prepared, "list")
  expect_true("clone_expression_rnorm" %in% names(prepared))
  expect_true("clone_counts" %in% names(prepared))
  expect_true("reduction_method" %in% names(prepared))
  expect_equal(prepared$reduction_method, "umap")

  # Rank-normalized clone expression: values in [0, 1], genes as rows.
  norm <- prepared$clone_expression_rnorm
  expect_true(is.matrix(norm))
  expect_true(all(norm >= 0 & norm <= 1))
  expect_true(ncol(norm) > 0)
})

test_that("prepare_data skip_clustering treats each column as a clone", {
  set.seed(7)
  genes <- sprintf("GENE%03d", seq_len(60))
  clones <- c("P1_c1", "P1_c2", "P2_c1")
  expr <- matrix(rpois(60 * 3, lambda = 5), nrow = 60,
                 dimnames = list(genes, clones))

  prepared <- prepare_data(
    expression_matrix = expr,
    skip_clustering = TRUE
  )

  expect_type(prepared, "list")
  expect_true("clone_expression_rnorm" %in% names(prepared))
  expect_equal(ncol(prepared$clone_expression_rnorm), 3)
})
