# ===========================================================================
# Automatic preprocessing: raw counts -> clone labels
#
# These functions implement the "path A" workflow: a user who only has a raw
# count matrix (no QC, no clustering, no clone labels) can still run the full
# PERCEPTION pipeline. Every step records what it did, so the result can be
# audited afterwards -- see the `label_source` field returned by
# prepare_data().
#
# Design notes:
#   * QC thresholds are MAD-based, i.e. relative to each patient's own
#     distribution. No tissue-specific threshold table is needed.
#   * HVG selection and PCA are computed on all cells, so the embedding is
#     shared; clustering then runs inside each patient, so clones never merge
#     across patients.
#   * Harmony is only used when a single patient has more than one sample
#     (e.g. longitudinal biopsies), because that is the only case where
#     technical batch effects split a real clone in two.
# ===========================================================================


# Upper fence of a robust (MAD-based) filter: median + nmads * MAD.
# Returns Inf when MAD is zero (all values identical) so the filter keeps
# everything instead of discarding half the cells on a degenerate sample.
qc_upper_threshold <- function(x, nmads = 3) {
  med <- stats::median(x, na.rm = TRUE)
  spread <- stats::mad(x, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) {
    return(Inf)
  }
  med + nmads * spread
}


# Per-cell QC metrics: detected genes, total counts, mitochondrial fraction.
# Mitochondrial genes are matched case-insensitively ("MT-" in human, "mt-"
# in mouse), so no species argument is required.
compute_qc_metrics <- function(counts) {
  n_genes <- Matrix::colSums(counts > 0)
  n_counts <- Matrix::colSums(counts)
  is_mito <- grepl("^MT-", toupper(rownames(counts)))

  pct_mito <- if (any(is_mito)) {
    mito_counts <- Matrix::colSums(counts[is_mito, , drop = FALSE])
    as.numeric(mito_counts) / pmax(as.numeric(n_counts), 1) * 100
  } else {
    rep(0, ncol(counts))
  }

  data.frame(
    cell_id = colnames(counts),
    n_genes = as.numeric(n_genes),
    n_counts = as.numeric(n_counts),
    pct_mito = as.numeric(pct_mito),
    stringsAsFactors = FALSE
  )
}


# Filter low-quality cells, one patient at a time.
#
# Thresholds are derived from each patient's own distribution
# (median +/- nmads * MAD) and then capped by absolute limits, so a sample
# with >50% bad cells cannot silently widen the fence.
qc_filter_cells <- function(counts,
                            patient,
                            nmads = 3,
                            min_genes = 200,
                            max_mito = 25,
                            max_genes = NULL) {
  metrics <- compute_qc_metrics(counts)
  metrics$patient <- as.character(patient)

  patients <- unique(metrics$patient)
  keep <- rep(TRUE, nrow(metrics))
  rows <- vector("list", length(patients))

  for (i in seq_along(patients)) {
    p <- patients[[i]]
    idx <- which(metrics$patient == p)
    n_genes <- metrics$n_genes[idx]
    pct_mito <- metrics$pct_mito[idx]

    genes_upper <- if (is.null(max_genes)) {
      qc_upper_threshold(n_genes, nmads)
    } else {
      max_genes
    }
    mito_upper <- min(qc_upper_threshold(pct_mito, nmads), max_mito)

    fail_low_genes <- n_genes < min_genes
    fail_high_genes <- n_genes > genes_upper
    fail_mito <- pct_mito > mito_upper
    ok <- !(fail_low_genes | fail_high_genes | fail_mito)
    keep[idx] <- ok

    reasons <- c(
      "low gene count" = sum(fail_low_genes),
      "high gene count" = sum(fail_high_genes),
      "high mitochondrial fraction" = sum(fail_mito)
    )

    rows[[i]] <- data.frame(
      patient = p,
      n_input = length(idx),
      n_kept = sum(ok),
      pct_removed = round(100 * mean(!ok), 1),
      main_reason = if (all(reasons == 0L)) {
        "-"
      } else {
        names(which.max(reasons))
      },
      mito_threshold = round(mito_upper, 1),
      genes_min = min_genes,
      genes_max = round(genes_upper, 1),
      stringsAsFactors = FALSE
    )
  }

  list(
    keep = metrics$cell_id[keep],
    report = do.call(rbind, rows),
    metrics = metrics
  )
}


# Check that the matrix is genes x cells and looks like raw counts.
# Returns a list describing what was detected; aborts on a transposed matrix
# because that error is silent and expensive to debug downstream.
detect_input_type <- function(counts) {
  if (inherits(counts, "Seurat")) {
    return(list(kind = "seurat", is_counts = NULL))
  }

  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    cli::cli_abort(c(
      "{.arg counts} must be a matrix or a Seurat object.",
      "x" = "Got {.cls {class(counts)}}."
    ))
  }

  rn <- rownames(counts)
  if (!is.null(rn) && mean(grepl("^[ACGTN]{12,}$", rn)) > 0.5) {
    cli::cli_abort(c(
      "{.arg counts} looks transposed: its rows are cell barcodes.",
      "i" = "Genes must be rows and cells must be columns.",
      "i" = "Try {.code t(counts)}."
    ))
  }

  nr <- min(200L, nrow(counts))
  nc <- min(200L, ncol(counts))
  probe <- as.numeric(counts[seq_len(nr), seq_len(nc), drop = FALSE])
  is_counts <- all(probe >= 0) && all(abs(probe - round(probe)) < 1e-8)

  if (!is_counts) {
    cli::cli_warn(c(
      "{.arg counts} does not look like raw counts.",
      "i" = "Non-integer values found; they will be log-normalised again.",
      "i" = "Pass raw counts when possible."
    ))
  }

  list(kind = "matrix", is_counts = is_counts)
}


# t-SNE counterpart of umap_from_pca(), used when method = "tsne".
tsne_from_pca <- function(so, dims, seed = 42) {
  so <- Seurat::RunTSNE(
    so,
    dims = seq_len(dims),
    check_duplicates = FALSE,
    seed.use = seed
  )
  emb <- Seurat::Embeddings(so, "tsne")[, c(1, 2), drop = FALSE]
  colnames(emb) <- c("TSNE_1", "TSNE_2")
  emb
}


# Copy the uwot UMAP settings used by run_seurat_pipeline(): the default
# RunUMAP settings cost seconds of fixed overhead on small datasets, and this
# embedding is only used for visualisation.
umap_from_pca <- function(so, dims, seed = 42, n_threads = 1L) {
  pca <- Seurat::Embeddings(so, "pca")[, seq_len(dims), drop = FALSE]

  n_neighbors <- min(15L, max(2L, nrow(pca) - 1L))
  coords <- uwot::umap(
    X = pca,
    n_neighbors = n_neighbors,
    n_components = 2,
    n_epochs = 200,
    n_trees = 10,
    metric = "cosine",
    learning_rate = 1,
    min_dist = 0.3,
    spread = 1,
    set_op_mix_ratio = 1,
    local_connectivity = 1,
    repulsion_strength = 1,
    negative_sample_rate = 5,
    n_threads = n_threads,
    seed = seed,
    verbose = FALSE
  )
  rownames(coords) <- colnames(so)
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  coords
}


# Cluster cells inside each patient on the shared PCA embedding.
#
# The embedding comes from all cells (so distances mean the same thing
# everywhere), but neighbours are only searched within a patient, which makes
# it impossible for a clone to span two patients.
cluster_within_patient <- function(so,
                                   patient_col = "patient",
                                   sample_col = NULL,
                                   dims = 20,
                                   resolution = 0.8,
                                   integrate = c("auto", "none", "harmony"),
                                   seed = 42,
                                   progress_cb = NULL) {
  integrate <- match.arg(integrate)
  set.seed(seed)

  emb <- Seurat::Embeddings(so, "pca")
  dims <- min(dims, ncol(emb))
  meta <- so[[]]
  patient <- as.character(meta[[patient_col]])
  patients <- unique(patient)

  out <- vector("list", length(patients))
  used_harmony <- FALSE

  for (i in seq_along(patients)) {
    p <- patients[[i]]
    idx <- which(patient == p)
    emb_p <- emb[idx, seq_len(dims), drop = FALSE]

    n_samples <- if (is.null(sample_col)) {
      1L
    } else {
      length(unique(meta[[sample_col]][idx]))
    }

    if (n_samples > 1L && integrate != "none") {
      if (!rlang::is_installed("harmony")) {
        cli::cli_warn(c(
          "Patient {.val {p}} has {n_samples} samples, but {.pkg harmony} is not installed.",
          "i" = "Clones from different samples may be split apart.",
          "i" = "Install with {.run install.packages(\"harmony\")}."
        ))
      } else {
        emb_p <- harmony::RunHarmony(
          data_mat = emb_p,
          meta_data = meta[idx, , drop = FALSE],
          vars_use = sample_col,
          verbose = FALSE
        )
        used_harmony <- TRUE
      }
    }

    nn <- Seurat::FindNeighbors(
      emb_p,
      dims = seq_len(ncol(emb_p)),
      verbose = FALSE
    )
    cl <- Seurat::FindClusters(
      nn$snn,
      resolution = resolution,
      verbose = FALSE
    )

    out[[i]] <- data.frame(
      cell_id = rownames(emb_p),
      clone_id = as.character(cl),
      patient = p,
      stringsAsFactors = FALSE
    )

    if (is.function(progress_cb)) {
      progress_cb(p, i, length(patients))
    }
  }

  list(labels = do.call(rbind, out), used_harmony = used_harmony)
}


#' Preprocess raw single-cell data into clone labels
#'
#' End-to-end preprocessing for users who only have a raw count matrix: quality
#' control, normalisation, variable-feature selection, PCA, and per-patient
#' clustering, followed by a QC report. No tissue-specific thresholds are
#' required; the QC fences adapt to each patient's own distribution.
#'
#' Clustering is performed **within each patient** on a shared PCA embedding,
#' so transcriptional clones can never span two patients. When a patient has
#' more than one sample (e.g. longitudinal biopsies), harmony is applied inside
#' that patient to prevent one biological clone from being split in two by a
#' technical batch effect. \pkg{harmony} is an optional dependency: without it,
#' the step is skipped with a warning.
#'
#' @param counts Matrix (genes as rows, cells as columns) or a \code{Seurat}
#'        object. Raw counts are expected; normalised values trigger a warning
#'        because they would be normalised twice.
#' @param patient_labels Character vector, one entry per cell, giving the patient
#'        a cell belongs to. Names are matched to \code{colnames(counts)}.
#' @param sample_labels Character vector, one entry per cell, giving the sample
#'        or batch a cell belongs to. Required only when a patient has multiple
#'        samples. Default = \code{NULL}.
#' @param nmads Numeric. Width of the MAD-based QC fence, in median absolute
#'        deviations. Default = 3.
#' @param min_genes Integer. Minimum detected genes per cell. Default = 200.
#' @param max_mito Numeric. Hard cap on the mitochondrial fraction (%). The
#'        effective threshold is the smaller of this and the MAD fence.
#'        Default = 25.
#' @param nfeatures Integer. Number of variable features used for PCA.
#'        Default = 2000.
#' @param dims Integer. Number of principal components. Capped automatically
#'        for small datasets. Default = 20.
#' @param method Character. 2D embedding for visualisation: \code{"umap"}
#'        (default) or \code{"tsne"}. Clustering itself always runs on the PCA
#'        embedding; this only changes the coordinates returned in
#'        \code{cell_labels}.
#' @param resolution Numeric. Clustering resolution. Default = 0.8.
#' @param min_clone_cells Integer. Clones with fewer cells are kept but marked
#'        with \code{confidence = "low"}. Default = 10.
#' @param integrate Character. Sample integration strategy within a patient:
#'        \code{"auto"} (default) applies harmony when a patient has multiple
#'        samples, \code{"none"} never integrates, \code{"harmony"} behaves
#'        like \code{"auto"}.
#' @param seed Integer. Random seed. Default = 42.
#' @param progress_cb Optional function called as \code{progress_cb(phase, i, n)}
#'        between stages. Used to drive progress bars in the Shiny app.
#'
#' @return A list with:
#' \describe{
#'   \item{cell_labels}{Data frame with \code{cell_id}, \code{clone_id},
#'         \code{clone_key} (unique \code{"patient@@clone"} key),
#'         \code{patient}, \code{confidence}, and \code{dim_1}/\code{dim_2}
#'         (UMAP coordinates). This is the \code{clone_labels} argument of
#'         \code{prepare_data()}.}
#'   \item{qc}{Data frame. Per-patient QC summary.}
#'   \item{qc_metrics}{Data frame. Per-cell QC metrics before filtering.}
#'   \item{used_harmony}{Logical. Whether harmony was applied.}
#'   \item{label_source}{Character. Always \code{"auto_cluster"}.}
#' }
#'
#' @examples
#' \dontrun{
#'   labels <- preprocess_sc(counts, patient_labels = cell_patients)
#'   prepared <- prepare_data(
#'     expression_matrix = counts,
#'     patient_mapping = labels$cell_labels[, c("cell_id", "patient")],
#'     clone_labels = labels$cell_labels
#'   )
#' }
#'
#' @export
preprocess_sc <- function(counts,
                          patient_labels,
                          sample_labels = NULL,
                          nmads = 3,
                          min_genes = 200,
                          max_mito = 25,
                          nfeatures = 2000,
                          dims = 20,
                          method = c("umap", "tsne"),
                          resolution = 0.8,
                          min_clone_cells = 10,
                          integrate = c("auto", "none", "harmony"),
                          seed = 42,
                          progress_cb = NULL) {
  integrate <- match.arg(integrate)
  method <- match.arg(method)
  set.seed(seed)

  notify <- function(phase, i = NULL, n = NULL) {
    if (is.function(progress_cb)) {
      progress_cb(phase, i, n)
    }
  }

  # --- 1. Input validation -------------------------------------------------
  notify("validating")
  detect_input_type(counts)

  if (inherits(counts, "Seurat")) {
    counts <- Seurat::GetAssayData(counts, layer = "counts")
  }
  if (length(patient_labels) != ncol(counts)) {
    cli::cli_abort(c(
      "{.arg patient_labels} must have one entry per cell.",
      "x" = "Got {length(patient_labels)} labels for {ncol(counts)} cells."
    ))
  }
  names(patient_labels) <- colnames(counts)
  patient_labels <- as.character(patient_labels)

  if (!is.null(sample_labels)) {
    if (length(sample_labels) != ncol(counts)) {
      cli::cli_abort(c(
        "{.arg sample_labels} must have one entry per cell.",
        "x" = "Got {length(sample_labels)} labels for {ncol(counts)} cells."
      ))
    }
    names(sample_labels) <- colnames(counts)
    sample_labels <- as.character(sample_labels)
  }

  # --- 2. Quality control --------------------------------------------------
  notify("qc")
  qc <- qc_filter_cells(
    counts,
    patient_labels,
    nmads = nmads,
    min_genes = min_genes,
    max_mito = max_mito
  )

  keep <- qc$keep
  counts <- counts[, keep, drop = FALSE]
  patient_labels <- patient_labels[keep]
  if (!is.null(sample_labels)) {
    sample_labels <- sample_labels[keep]
  }

  if (ncol(counts) < 10L) {
    cli::cli_abort(c(
      "Only {ncol(counts)} cell{?s} left after quality control.",
      "i" = "Relax {.arg nmads} or {.arg min_genes} to keep more cells."
    ))
  }

  # --- 3. Normalise and embed on ALL cells ---------------------------------
  notify("normalizing")
  meta <- data.frame(
    patient = patient_labels,
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  if (!is.null(sample_labels)) {
    meta$sample <- sample_labels
  }

  # min.cells / min.features = 0: QC already ran above and downstream users
  # may pass an already-filtered matrix, so no second silent filter here.
  so <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    min.cells = 0,
    min.features = 0,
    project = "PERCEPTIONx"
  )

  so <- Seurat::NormalizeData(
    so,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

  notify("variable-features")
  nfeatures <- min(nfeatures, nrow(so))
  so <- Seurat::FindVariableFeatures(
    so,
    selection.method = "vst",
    nfeatures = nfeatures
  )

  notify("scaling")
  so <- Seurat::ScaleData(so, features = Seurat::VariableFeatures(so))

  npcs <- min(
    dims,
    ncol(so) - 1L,
    length(Seurat::VariableFeatures(so)) - 1L
  )
  npcs <- max(1L, as.integer(npcs))

  notify("pca")
  so <- Seurat::RunPCA(
    so,
    features = Seurat::VariableFeatures(so),
    npcs = npcs
  )

  # --- 4. Cluster within each patient --------------------------------------
  notify("clustering")
  clustered <- cluster_within_patient(
    so,
    patient_col = "patient",
    sample_col = if (is.null(sample_labels)) NULL else "sample",
    dims = npcs,
    resolution = resolution,
    integrate = integrate,
    seed = seed,
    progress_cb = function(p, i, n) notify("clustering", i, n)
  )

  # --- 5. Shared 2D embedding for visualisation ----------------------------
  notify("embedding")
  n_threads <- max(1L, min(4L, as.integer(parallel::detectCores()), na.rm = TRUE))
  coords <- if (method == "umap") {
    umap_from_pca(so, npcs, seed = seed, n_threads = n_threads)
  } else {
    tsne_from_pca(so, npcs, seed = seed)
  }

  # --- 6. Assemble labels --------------------------------------------------
  notify("assembling")
  labels <- clustered$labels
  labels$dim_1 <- coords[labels$cell_id, 1]
  labels$dim_2 <- coords[labels$cell_id, 2]

  # Cluster ids restart at 0 for every patient, so they are only unique within
  # a patient. `clone_key` is the globally unique form used downstream.
  labels$clone_key <- paste(labels$patient, labels$clone_id, sep = "@@")

  clone_size <- table(labels$clone_key)[labels$clone_key]
  labels$confidence <- ifelse(
    as.integer(clone_size) >= min_clone_cells,
    "high",
    "low"
  )
  rownames(labels) <- NULL

  qc$report$n_clones <- vapply(
    qc$report$patient,
    function(p) length(unique(labels$clone_id[labels$patient == p])),
    integer(1)
  )

  list(
    cell_labels = labels,
    qc = qc$report,
    qc_metrics = qc$metrics,
    used_harmony = clustered$used_harmony,
    label_source = "auto_cluster"
  )
}


#' Plot QC metrics before and after filtering
#'
#' Draws one violin per metric, split by whether the cell was kept, with the
#' per-patient thresholds drawn as horizontal lines. This is what makes the
#' automatic QC auditable: a user does not need to understand the MAD rule to
#' see whether their data were cut too aggressively.
#'
#' @param qc_metrics Data frame returned as \code{qc_metrics} by
#'        \code{preprocess_sc()} (columns \code{cell_id}, \code{n_genes},
#'        \code{n_counts}, \code{pct_mito}, \code{patient}).
#' @param kept_cells Character vector of cell IDs that passed QC.
#' @param qc_report Data frame returned as \code{qc} by \code{preprocess_sc()}.
#' @param base_size Numeric. Base font size. Default = 12.
#'
#' @return A \code{ggplot} object, faceted by metric.
#'
#' @export
plot_qc_report <- function(qc_metrics, kept_cells, qc_report,
                           base_size = 12) {
  metrics <- qc_metrics
  metrics$status <- ifelse(
    metrics$cell_id %in% kept_cells,
    "kept",
    "removed"
  )

  long <- data.frame(
    patient = rep(metrics$patient, 3L),
    metric = rep(
      c("Detected genes", "Total counts", "Mitochondrial %"),
      each = nrow(metrics)
    ),
    value = c(metrics$n_genes, metrics$n_counts, metrics$pct_mito),
    status = rep(metrics$status, 3L),
    stringsAsFactors = FALSE
  )

  thresholds <- data.frame(
    patient = rep(qc_report$patient, 2L),
    metric = rep(
      c("Detected genes", "Mitochondrial %"),
      each = nrow(qc_report)
    ),
    value = c(qc_report$genes_min, qc_report$mito_threshold),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$status, y = .data$value, fill = .data$status)
  ) +
    ggplot2::geom_violin(scale = "width", linewidth = 0.25) +
    ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white") +
    ggplot2::geom_hline(
      data = thresholds,
      ggplot2::aes(yintercept = .data$value),
      linetype = "dashed",
      colour = "grey35",
      linewidth = 0.3
    ) +
    ggplot2::facet_grid(metric ~ patient, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(kept = "#009E73", removed = "#D55E00")) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.22, colour = "grey91"),
      panel.border = ggplot2::element_rect(linewidth = 0.4, colour = "grey35"),
      strip.background = ggplot2::element_rect(fill = "grey95", linewidth = 0.3),
      legend.position = "top"
    )
}
