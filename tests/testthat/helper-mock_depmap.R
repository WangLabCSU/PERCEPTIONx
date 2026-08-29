# =============================================================================
# Minimal mock DepMap for CI smoke tests
# -----------------------------------------------------------------------------
# The real DepMap bundle is a ~567 MB download, so it cannot ship with the
# package or be fetched in CI. These tests build a tiny structurally faithful
# fake DepMap so the CORE pipeline (feature ranking -> train_models ->
# predict_drugs -> predict_patients) actually executes on CI instead of being
# skipped. Values are random noise - we only assert that the code runs and the
# output shapes are right, never the numbers.
# =============================================================================

mock_depmap <- function(seed = 42,
                        n_genes = 60,
                        n_train = 8,
                        n_test = 4,
                        drugs = c("erlotinib", "gefitinib")) {
  set.seed(seed)
  genes   <- sprintf("GENE%03d", seq_len(n_genes))
  cl_all  <- paste0("CL", sprintf("%02d", seq_len(n_train + n_test)))
  cl_train <- cl_all[seq_len(n_train)]
  cl_test  <- cl_all[seq_len(n_test) + n_train]

  # Single-cell sample names: 2 cells per test cell line.
  cell_names <- paste0("SC_", rep(cl_test, each = 2), "_", rep(c("a", "b"), n_test))

  # Bulk expression: genes x ALL cell lines (rows must match scRNA_subset_rnorm).
  expr <- matrix(runif(n_genes * length(cl_all), 0, 10), nrow = n_genes,
                 dimnames = list(genes, cl_all))
  # Single-cell expression per cell line (columns = cell-line IDs).
  sc_by_line <- matrix(runif(n_genes * n_test, 0, 10), nrow = n_genes,
                       dimnames = list(genes, cl_test))
  # Single-cell expression per single cell.
  sc_by_cell <- matrix(runif(n_genes * length(cell_names), 0, 10), nrow = n_genes,
                       dimnames = list(genes, cell_names))
  # Drug response: one row per drug, columns = cell lines.
  resp <- matrix(runif(length(drugs) * length(cl_all), -2, 2),
                 nrow = length(drugs), byrow = TRUE,
                 dimnames = list(drugs, cl_all))

  list(
    expression_rnorm        = expr,
    expression_20Q4         = expr,
    scRNA_complete          = sc_by_line,
    scRNA_subset_rnorm      = sc_by_line,
    CPM_scRNA_CCLE_rnorm    = sc_by_cell,
    annotation_20Q4         = data.frame(DepMap_ID = cl_all,
                                         lineage   = rep("Lung", length(cl_all)),
                                         stringsAsFactors = FALSE),
    metadata_CPM_scRNA      = data.frame(NAME      = cell_names,
                                         DepMap_ID = rep(cl_test, each = 2),
                                         stringsAsFactors = FALSE),
    secondary_prism         = resp,
    secondary_screen_drugAnnotation = data.frame(
      CommonName = drugs,
      Screen_id  = paste0("MTS", seq_along(drugs)),
      stringsAsFactors = FALSE
    )
  )
}
