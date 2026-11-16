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


test_that("compute_qc_metrics detects mitochondrial genes given Ensembl IDs", {
  # 10x GRCh38 references and many GEO uploads use Ensembl IDs as row names, so
  # the "MT-" prefix never matches. If the fallback list were missing the mito
  # fraction would be silently 0 and the mito fence would stop working.
  counts <- matrix(
    c(10, 10, 5, 5, 0, 0, 10, 5),
    nrow = 4,
    dimnames = list(
      c(
        "ENSG00000198888",
        "ENSG00000198763",
        "ENSG00000000003",
        "ENSG00000000005"
      ),
      c("C1", "C2")
    )
  )
  m <- compute_qc_metrics(counts)

  expect_equal(m$pct_mito[1], 20 / 30 * 100, tolerance = 1e-8)
  expect_equal(m$pct_mito[2], 0)
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


test_that("preprocess_sc labels every cell with a (patient, clone) pair", {
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
      "doublet_cells",
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

  # clone_key = "patient@@clone", so it is unique even when one clone id is
  # shared by several patients (the default global clustering).
  pairs <- unique(res$cell_labels[, c("patient", "clone_id")])
  expect_equal(nrow(pairs), length(unique(res$cell_labels$clone_key)))

  expect_true(all(res$cell_labels$cluster_size_class %in% c("large", "small")))
  expect_true(all(res$cell_labels$patient %in% c("P1", "P2")))
  expect_equal(nrow(res$excluded_cells), 0L)
})


test_that("cluster_scope controls whether a clone can span patients", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  sim <- make_mock_counts(n_genes = 300, n_cells = 80)
  global <- preprocess_sc(
    sim$counts,
    sim$patient,
    nmads = 5,
    dims = 10,
    resolution = 0.5
  )
  local <- preprocess_sc(
    sim$counts,
    sim$patient,
    nmads = 5,
    dims = 10,
    resolution = 0.5,
    cluster_scope = "patient"
  )

  # global (the default): at least one clone id contains cells of both patients
  by_patient <- split(global$cell_labels$clone_id, global$cell_labels$patient)
  expect_gt(length(intersect(by_patient$P1, by_patient$P2)), 0L)

  # patient: numbers restart at 0 inside every patient, which is what makes the
  # labels patient-local (the same id means different clusters in P1 and P2).
  for (p in c("P1", "P2")) {
    ids <- sort(unique(local$cell_labels$clone_id[
      local$cell_labels$patient == p
    ]))
    expect_equal(ids, as.character(seq_along(ids) - 1L))
  }
  expect_equal(nrow(local$cell_labels), nrow(global$cell_labels))
})


test_that("skip_qc keeps every cell and marks the fences as off", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  sim <- make_mock_counts(n_genes = 300, n_cells = 60)
  # A gene floor above the panel size would normally remove every cell.
  res <- preprocess_sc(
    sim$counts,
    sim$patient,
    nmads = 5,
    min_genes = 1e6,
    dims = 10,
    resolution = 0.5,
    skip_qc = TRUE
  )

  expect_equal(nrow(res$cell_labels), ncol(sim$counts))
  expect_equal(sum(res$qc$n_input), ncol(sim$counts))
  expect_equal(sum(res$qc$n_kept), ncol(sim$counts))
  expect_true(all(is.na(res$qc$genes_max)))
  expect_true(all(res$qc$main_reason == "-"))
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


test_that("remove_doublets reports what it removed and keeps the counts consistent", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")
  skip_if_not_installed("scDblFinder")
  skip_if_not_installed("SingleCellExperiment")

  sim <- make_mock_counts(n_genes = 300, n_cells = 120)
  # Artificial doublets: sums of two random cells, the signal scDblFinder
  # is designed to detect.
  set.seed(7)
  for (i in seq_len(10)) {
    sim$counts[, i] <- sim$counts[, sample(ncol(sim$counts), 1)] +
      sim$counts[, sample(ncol(sim$counts), 1)]
  }

  res <- suppressWarnings(preprocess_sc(
    sim$counts,
    sim$patient,
    nmads = 5,
    dims = 10,
    resolution = 0.5,
    remove_doublets = TRUE
  ))

  expect_true("n_doublet_removed" %in% names(res$qc))
  expect_named(res$doublet_cells, c("cell_id", "patient", "score"))
  expect_equal(
    nrow(res$cell_labels),
    sum(res$qc$n_kept) - sum(res$qc$n_doublet_removed)
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


# Two clearly separated transcriptional groups inside every patient, used to
# check that clustering actually separates populations rather than always
# returning a single cluster.
make_two_group_counts <- function(n_genes = 400, n_per = 25) {
  set.seed(23)
  block <- function(rows) {
    m <- matrix(stats::rpois(n_genes * n_per, 2), n_genes, n_per)
    m[rows, ] <- m[rows, ] + 60
    m
  }
  # Two blocks per patient: the second block over-expresses a different gene
  # set, so each patient contains two transcriptionally distinct populations.
  counts <- cbind(
    block(seq_len(120)),
    block(121:240),
    block(seq_len(120)),
    block(121:240)
  )
  colnames(counts) <- sprintf("C%03d", seq_len(ncol(counts)))
  rownames(counts) <- sprintf("G%03d", seq_len(n_genes))

  patient <- rep(c("P1", "P2"), each = n_per * 2)
  names(patient) <- colnames(counts)
  stopifnot(ncol(counts) == length(patient))
  list(counts = counts, patient = patient)
}


test_that("every retained cell keeps its patient assignment", {
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

  expect_false(anyNA(res$cell_labels$patient))
  expect_equal(
    unname(res$cell_labels$patient),
    unname(sim$patient[res$cell_labels$cell_id])
  )
})


test_that("cluster labels are one per cell and are not deparsed vectors", {
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

  expect_equal(nrow(res$cell_labels), length(unique(res$cell_labels$cell_id)))
  expect_false(anyDuplicated(res$cell_labels$cell_id) > 0)

  # Regression guard for the Seurat 5 FindClusters() data.frame return: a
  # whole column used to be deparsed into "c(1, 1, 1, ...)" and recycled.
  expect_false(any(grepl("c\\(", res$cell_labels$clone_id)))
  expect_true(all(nchar(res$cell_labels$clone_id) < 10L))
})


test_that("clustering separates two distinct populations per patient", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  sim <- make_two_group_counts()
  res <- preprocess_sc(
    sim$counts,
    sim$patient,
    nmads = 5,
    dims = 15,
    resolution = 0.5
  )

  clusters_per_patient <- tapply(
    res$cell_labels$clone_id,
    res$cell_labels$patient,
    function(x) length(unique(x))
  )

  expect_true(all(clusters_per_patient >= 2L))
})


test_that("a patient too small to cluster becomes a single clone", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  sim <- make_mock_counts(n_genes = 300, n_cells = 60)
  patient <- sim$patient
  # Two cells cannot support a neighbour graph; they must not abort the run.
  patient[colnames(sim$counts)[1:2]] <- "P_tiny"

  # The single-clone fallback is part of the per-patient clustering path.
  res <- preprocess_sc(
    sim$counts,
    patient,
    nmads = 5,
    dims = 10,
    resolution = 0.5,
    cluster_scope = "patient"
  )

  tiny <- res$cell_labels[res$cell_labels$patient == "P_tiny", ]
  expect_equal(nrow(tiny), 2L)
  expect_equal(unique(tiny$clone_id), "1")
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


test_that("clone_labels supply the patient grouping when no mapping is given", {
  sim <- make_mock_counts(n_genes = 100, n_cells = 20)

  labels <- data.frame(
    cell_id = colnames(sim$counts),
    clone_id = rep(c("c1", "c2"), length.out = ncol(sim$counts)),
    patient = as.character(sim$patient),
    stringsAsFactors = FALSE
  )

  prepared <- prepare_data(
    expression_matrix = sim$counts,
    clone_labels = labels
  )

  expect_equal(
    sort(unique(prepared$clone_counts$patients)),
    sort(unique(as.character(sim$patient)))
  )
  expect_equal(nrow(prepared$cell_clone_map), ncol(sim$counts))
})


test_that("cells without a clone label are dropped instead of failing", {
  sim <- make_mock_counts(n_genes = 100, n_cells = 20)

  # Simulates the QC case: labels exist for a subset of the matrix columns, and
  # the unfiltered matrix is passed as is.
  labels <- data.frame(
    cell_id = colnames(sim$counts)[-1],
    clone_id = rep(c("c1", "c2"), length.out = ncol(sim$counts) - 1L),
    patient = as.character(sim$patient)[-1],
    stringsAsFactors = FALSE
  )

  prepared <- prepare_data(
    expression_matrix = sim$counts,
    clone_labels = labels
  )

  expect_equal(nrow(prepared$cell_clone_map), ncol(sim$counts) - 1L)
  expect_false(colnames(sim$counts)[1] %in% prepared$cell_clone_map$cell_id)
})


test_that("Seurat input prefers the counts layer and can be overridden", {
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

  auto <- PERCEPTIONx:::seurat_expression_matrix(so)
  counts <- SeuratObject::LayerData(so, assay = "RNA", layer = "counts")
  normalised <- SeuratObject::LayerData(so, assay = "RNA", layer = "data")

  # Default = counts, so the Seurat path and the raw-counts path average the
  # same scale. Passing layer = "data" still gets the normalised values.
  expect_equal(auto, counts)
  expect_equal(
    PERCEPTIONx:::seurat_expression_matrix(so, layer = "data"),
    normalised
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
