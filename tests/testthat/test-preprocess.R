# =============================================================================
# Automatic preprocessing (path A): raw counts -> clone labels
# =============================================================================

make_mock_counts <- function(n_genes = 300, n_cells = 80, n_patients = 2) {
  set.seed(11)
  counts <- matrix(
    stats::rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes,
    dimnames = list(
      sprintf("GENE%03d", seq_len(n_genes)),
      sprintf("CELL%03d", seq_len(n_cells))
    )
  )
  # A few mitochondrial genes so the QC metric is not constant.
  rownames(counts)[seq_len(3)] <- paste0("MT-G", seq_len(3))

  patient <- rep(
    sprintf("P%d", seq_len(n_patients)),
    each = n_cells %/% n_patients
  )
  names(patient) <- colnames(counts)

  list(counts = counts, patient = patient)
}


test_that("qc_upper_threshold ignores extreme outliers", {
  x <- c(1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 45, 60)
  threshold <- qc_upper_threshold(x, nmads = 3)

  # median = 4.5, MAD = 3.71 -> fence 15.6; both outliers fall outside.
  expect_gt(threshold, 10)
  expect_lt(threshold, 20)
})


test_that("qc_upper_threshold returns Inf when every value is identical", {
  expect_equal(qc_upper_threshold(rep(5, 10)), Inf)
})


test_that("qc_filter_cells reports one row per patient", {
  sim <- make_mock_counts(n_cells = 40)
  res <- qc_filter_cells(sim$counts, sim$patient, nmads = 5)

  expect_named(res, c("keep", "report", "metrics"))
  expect_equal(nrow(res$report), 2L)
  expect_true(
    all(
      c(
        "patient",
        "n_input",
        "n_kept",
        "pct_removed",
        "main_reason",
        "mito_threshold"
      ) %in%
        names(res$report)
    )
  )
  expect_true(all(res$keep %in% colnames(sim$counts)))
})


test_that("detect_input_type aborts on a transposed matrix", {
  counts <- matrix(
    1,
    nrow = 5,
    ncol = 3,
    dimnames = list(
      c(
        "AAACCTGAGCGTCTAT",
        "AAACCTGAGTACGTTC",
        "AAACCTGGTCAGAGCC",
        "AAACCTGGTCTTGTCC",
        "AAACCTGTCTCGATGA"
      ),
      c("G1", "G2", "G3")
    )
  )

  expect_error(detect_input_type(counts), "transposed")
})


test_that("detect_input_type warns on non-count data", {
  counts <- matrix(
    c(1.5, 2.2, 3.1, 0.4),
    nrow = 2,
    dimnames = list(c("G1", "G2"), c("C1", "C2"))
  )

  expect_warning(detect_input_type(counts), "raw counts")
})


test_that("detect_input_type accepts an integer matrix", {
  counts <- matrix(
    c(1L, 2L, 3L, 0L),
    nrow = 2,
    dimnames = list(c("G1", "G2"), c("C1", "C2"))
  )

  info <- detect_input_type(counts)

  expect_equal(info$kind, "matrix")
  expect_true(info$is_counts)
})


test_that("preprocess_sc returns patient-local clones", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  sim <- make_mock_counts(n_genes = 300, n_cells = 80)
  res <- preprocess_sc(
    sim$counts,
    sim$patient,
    nmads = 5,
    dims = 10,
    resolution = 0.5
  )

  expect_named(
    res,
    c(
      "cell_labels",
      "qc",
      "qc_metrics",
      "excluded_cells",
      "used_harmony",
      "label_source"
    )
  )
  expect_equal(res$label_source, "auto_cluster")
  expect_true(
    all(
      c(
        "cell_id",
        "clone_id",
        "clone_key",
        "patient",
        "cluster_size_class",
        "dim_1",
        "dim_2"
      ) %in%
        names(res$cell_labels)
    )
  )

  # Cluster ids restart at 0 inside every patient, so the pair
  # (patient, clone_id) must be unique even though clone_id alone is not.
  pairs <- unique(res$cell_labels[, c("patient", "clone_id")])
  expect_equal(nrow(pairs), length(unique(res$cell_labels$clone_key)))

  expect_true(all(res$cell_labels$cluster_size_class %in% c("large", "small")))
  expect_true(all(res$cell_labels$patient %in% c("P1", "P2")))
  expect_equal(nrow(res$excluded_cells), 0L)
})


test_that("exclude_markers removes cells that score high for a marker set", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  sim <- make_mock_counts(n_genes = 300, n_cells = 60)
  rownames(sim$counts)[1] <- "PTPRC"
  sim$counts[1, ] <- 0
  sim$counts[1, seq_len(20)] <- 200

  res <- preprocess_sc(
    sim$counts,
    sim$patient,
    nmads = 5,
    dims = 10,
    resolution = 0.5,
    exclude_markers = list(immune = "PTPRC"),
    exclude_min_expr = 1
  )

  expect_gt(nrow(res$excluded_cells), 0L)
  expect_equal(unique(res$excluded_cells$set), "immune")
  expect_false(
    any(res$excluded_cells$cell_id %in% res$cell_labels$cell_id)
  )
})


test_that("filter_marker_sets keeps every cell when markers are absent", {
  mat <- matrix(
    1,
    nrow = 3,
    ncol = 4,
    dimnames = list(c("A", "B", "C"), c("c1", "c2", "c3", "c4"))
  )
  res <- filter_marker_sets(mat, list(immune = c("NOT_A_GENE")))

  expect_equal(res$keep, colnames(mat))
  expect_equal(nrow(res$removed), 0L)
})


test_that("filter_marker_sets requires a named list", {
  mat <- matrix(
    1,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("A", "B"), c("c1", "c2"))
  )
  expect_error(
    filter_marker_sets(mat, c("A")),
    "named list"
  )
})


test_that("prepare_data accepts clone_labels and records their source", {
  sim <- make_mock_counts(n_genes = 100, n_cells = 20)

  labels <- data.frame(
    cell_id = colnames(sim$counts),
    clone_id = rep(c("c1", "c2"), length.out = ncol(sim$counts)),
    stringsAsFactors = FALSE
  )
  mapping <- data.frame(
    cell_id = colnames(sim$counts),
    patient_id = as.character(sim$patient),
    stringsAsFactors = FALSE
  )

  prepared <- prepare_data(
    expression_matrix = sim$counts,
    patient_mapping = mapping,
    clone_labels = labels
  )

  expect_equal(prepared$label_source, "user_labels")
  expect_true(
    all(c("clone_expression_rnorm", "clone_counts") %in% names(prepared))
  )
})


test_that("prepare_data extracts labels from a Seurat object", {
  skip_if_not_installed("Seurat")

  sim <- make_mock_counts(n_genes = 100, n_cells = 20)

  so <- Seurat::CreateSeuratObject(
    counts = sim$counts,
    meta.data = data.frame(
      patient = as.character(sim$patient),
      row.names = colnames(sim$counts)
    )
  )
  so <- Seurat::NormalizeData(so, verbose = FALSE)
  so$clone <- rep(c("c1", "c2"), length.out = ncol(so))
  Seurat::Idents(so) <- so$clone

  prepared <- prepare_data(data = so, patient_col = "patient")

  expect_equal(prepared$label_source, "seurat_idents")
  expect_true("clone_counts" %in% names(prepared))
})
