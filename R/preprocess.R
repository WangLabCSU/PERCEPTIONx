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

# ---------------------------------------------------------------------------
# Reading 10x Genomics count matrices
#
# The entry point of path A. The Shiny upload uses the same function, so there
# is one implementation of "10x files in, sparse matrix out" for both users.
# ---------------------------------------------------------------------------

# Read a text file, transparently handling gzip.
read_10x_lines <- function(path) {
  if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    readLines(con, warn = FALSE)
  } else {
    readLines(path, warn = FALSE)
  }
}

# Turn the three files of a 10x triplet into a sparse genes x cells matrix.
# Stays sparse end to end: a raw matrix has hundreds of thousands of barcodes
# and densifying it is not an option.
read_10x_triplet_files <- function(mtx, barcodes, features) {
  feat <- utils::read.delim(
    text = read_10x_lines(features),
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = ""
  )
  bc <- read_10x_lines(barcodes)
  # Column 2 is the gene symbol; column 1 is the Ensembl ID. Fall back to
  # column 1 when a file carries only one column.
  genes <- if (ncol(feat) >= 2L) feat[[2]] else feat[[1]]

  mat <- if (grepl("[.]gz$", mtx, ignore.case = TRUE)) {
    con <- gzfile(mtx, "rt")
    on.exit(close(con), add = TRUE)
    Matrix::readMM(con)
  } else {
    Matrix::readMM(mtx)
  }
  mat <- methods::as(mat, "CsparseMatrix")

  if (nrow(mat) != length(genes) || ncol(mat) != length(bc)) {
    cli::cli_abort(c(
      "The 10x files do not belong together.",
      "x" = "The matrix is {nrow(mat)} genes x {ncol(mat)} cells, but features/barcodes give {length(genes)} x {length(bc)}."
    ))
  }

  # Feature files repeat symbols occasionally; Seurat needs unique row names.
  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- bc
  mat
}

# Locate the triplet inside a directory, e.g. the `filtered_feature_bc_matrix`
# folder CellRanger writes.
find_10x_triplet <- function(dir) {
  files <- list.files(dir, recursive = TRUE, full.names = TRUE)
  nm <- tolower(basename(files))

  pick <- function(pattern) {
    hit <- which(grepl(pattern, nm))
    if (length(hit) == 1L) files[[hit]] else NULL
  }
  found <- list(
    mtx = pick("^matrix[.]mtx([.]gz)?$"),
    barcodes = pick("^barcodes[.]tsv([.]gz)?$"),
    features = pick("^(genes|features)[.]tsv([.]gz)?$")
  )
  missing <- names(found)[vapply(found, is.null, logical(1))]
  if (length(missing) > 0L) {
    cli::cli_abort(c(
      "No complete 10x triplet in {.path {dir}}.",
      "x" = "Missing: {.file {missing}}.",
      "i" = "Expected matrix.mtx[.gz], barcodes.tsv[.gz] and genes/features.tsv[.gz]."
    ))
  }
  found
}

#' Read a 10x Genomics count matrix
#'
#' One reader for the inputs a user actually has: a CellRanger output folder, a
#' `.h5` file, or the three files passed explicitly. The result is a sparse
#' genes x cells matrix with gene symbols as row names, which is what
#' `preprocess_sc()` and `prepare_data()` expect.
#'
#' @param path Directory holding `matrix.mtx[.gz]`, `barcodes.tsv[.gz]` and
#'   `genes.tsv[.gz]` / `features.tsv[.gz]`, or a `.h5` file written by
#'   CellRanger. Ignored when the three files are given explicitly.
#' @param mtx Path to `matrix.mtx` or `matrix.mtx.gz`.
#' @param barcodes Path to `barcodes.tsv` or `barcodes.tsv.gz`.
#' @param features Path to `genes.tsv` or `features.tsv` (either `.gz` or not).
#'
#' @return A `dgCMatrix` with genes as rows and cells as columns.
#'
#' @examples
#' \dontrun{
#' counts <- read_10x("path/to/filtered_feature_bc_matrix")
#' counts <- read_10x("path/to/filtered_feature_bc_matrix.h5")
#' counts <- read_10x(
#'   mtx = "matrix.mtx.gz",
#'   barcodes = "barcodes.tsv.gz",
#'   features = "features.tsv.gz"
#' )
#' }
#'
#' @export
read_10x <- function(
  path = NULL,
  mtx = NULL,
  barcodes = NULL,
  features = NULL
) {
  if (!is.null(path)) {
    if (dir.exists(path)) {
      found <- find_10x_triplet(path)
      return(read_10x_triplet_files(
        found$mtx,
        found$barcodes,
        found$features
      ))
    }
    if (grepl("[.]h5$", path, ignore.case = TRUE)) {
      rlang::check_installed("hdf5r", reason = "to read 10x .h5 matrices")
      mat <- Seurat::Read10X_h5(path)
      # Multi-modal files return a named list; the gene expression block is the
      # one preprocess_sc() wants.
      if (is.list(mat)) {
        hit <- grep("gene expression", names(mat), ignore.case = TRUE)
        mat <- mat[[if (length(hit) > 0L) hit[[1]] else 1L]]
      }
      return(methods::as(mat, "CsparseMatrix"))
    }
    cli::cli_abort(c(
      "{.path {path}} is neither a directory nor a .h5 file.",
      "i" = "Point {.arg path} at the CellRanger output folder, or pass {.arg mtx}, {.arg barcodes} and {.arg features}."
    ))
  }

  if (is.null(mtx) || is.null(barcodes) || is.null(features)) {
    cli::cli_abort(
      "Supply either {.arg path}, or all three of {.arg mtx}, {.arg barcodes} and {.arg features}."
    )
  }
  read_10x_triplet_files(mtx, barcodes, features)
}

# Upper fence of a robust (MAD-based) filter: median + nmads * MAD.
# Returns Inf when MAD is zero (all values identical) so the filter keeps
# everything instead of discarding half the cells on a degenerate sample.
qc_upper_threshold <- function(x, nmads = 5) {
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

  # Some references name genes by Ensembl ID (10x GRCh38 with `features.tsv`
  # column 1, many GEO uploads), so "MT-" never matches and the mitochondrial
  # fraction would silently be 0 for every cell -- the mito fence would do
  # nothing. Fall back to the fixed Ensembl IDs of the 13 protein-coding
  # mitochondrial genes plus the two rRNA genes where a reference carries them.
  if (!any(is_mito)) {
    mt_ids <- c(
      "ENSG00000198888",
      "ENSG00000198763",
      "ENSG00000198804",
      "ENSG00000198712",
      "ENSG00000228253",
      "ENSG00000198899",
      "ENSG00000198938",
      "ENSG00000198840",
      "ENSG00000212907",
      "ENSG00000198886",
      "ENSG00000198786",
      "ENSG00000198695",
      "ENSG00000198727",
      "ENSG00000210049",
      "ENSG00000210082"
    )
    is_mito <- sub("[.][0-9]+$", "", rownames(counts)) %in% mt_ids
  }

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
qc_filter_cells <- function(
  counts,
  patient,
  nmads = 5,
  min_genes = 200,
  max_mito = 25,
  max_genes = NULL,
  skip = FALSE
) {
  metrics <- compute_qc_metrics(counts)
  metrics$patient <- as.character(patient)

  # "No QC" still needs a report: the app and the vignette both read the same
  # table, so every cell is simply marked as kept and the fences are NA.
  if (isTRUE(skip)) {
    patients <- unique(metrics$patient)
    report <- data.frame(
      patient = patients,
      n_input = as.integer(table(metrics$patient)[patients]),
      n_kept = as.integer(table(metrics$patient)[patients]),
      pct_removed = 0,
      main_reason = "-",
      mito_threshold = NA_real_,
      genes_min = NA_integer_,
      genes_max = NA_real_,
      stringsAsFactors = FALSE
    )
    rownames(report) <- NULL
    return(list(keep = metrics$cell_id, report = report, metrics = metrics))
  }

  patients <- unique(metrics$patient)
  keep <- rep(TRUE, nrow(metrics))
  rows <- vector("list", length(patients))

  # Pooled fences across every sample. A sample whose plausible-cell pool is
  # dominated by low-quality barcodes gets a collapsed median/MAD -- the fence
  # lands on the junk tail -- and would then remove its real cells. When that
  # happens we fall back to the pooled estimate instead of silently losing the
  # sample.
  plausible <- metrics$n_genes >= min_genes
  pooled <- if (sum(plausible) >= 10L) {
    list(
      genes = if (is.null(max_genes)) {
        qc_upper_threshold(metrics$n_genes[plausible], nmads)
      } else {
        max_genes
      },
      mito = min(
        qc_upper_threshold(metrics$pct_mito[plausible], nmads),
        max_mito
      )
    )
  } else {
    NULL
  }

  for (i in seq_along(patients)) {
    p <- patients[[i]]
    idx <- which(metrics$patient == p)
    n_genes <- metrics$n_genes[idx]
    pct_mito <- metrics$pct_mito[idx]

    # Estimate the fence on plausible cells only. A raw CellRanger matrix is
    # >99% empty droplets, and those zero-gene barcodes drag the median and
    # MAD down to 0 — which turns both upper fences into +/-Inf, so the "MAD"
    # half of the rule silently does nothing on exactly the input that needs
    # it most. The absolute gene floor is a safe pre-filter: it is applied to
    # every cell anyway, so it cannot change which cells survive.
    est <- idx[n_genes >= min_genes]
    if (length(est) < 10L) {
      est <- idx
    }
    n_genes_est <- metrics$n_genes[est]
    pct_mito_est <- metrics$pct_mito[est]

    genes_upper <- if (is.null(max_genes)) {
      qc_upper_threshold(n_genes_est, nmads)
    } else {
      max_genes
    }
    mito_upper <- min(qc_upper_threshold(pct_mito_est, nmads), max_mito)

    # Two ways the local fence can be nonsense: it collapses onto the
    # low-quality tail (below ~2x the gene floor), or it keeps nothing at all.
    # Both mean "the plausible-cell pool was not representative", so use the
    # pooled fence rather than deleting the sample.
    if (!is.null(pooled)) {
      collapsed <- is.null(max_genes) && genes_upper < max(500, 2 * min_genes)
      keeps_none <- !any(
        !(n_genes < min_genes | n_genes > genes_upper | pct_mito > mito_upper)
      )
      if (collapsed || keeps_none) {
        if (is.null(max_genes)) {
          genes_upper <- pooled$genes
        }
        mito_upper <- min(pooled$mito, max_mito)
      }
    }

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


# Score each cell against user-supplied marker sets and drop the cells whose
# strongest set exceeds `min_expr`.
#
# This is an EXCLUSION rule, not a tumour-cell classifier: "epithelial" is not
# the same as "malignant", and no marker rule can establish malignancy. The
# only promise is that a cell strongly expressing an excluded programme is left
# out of the clustering, so immune or stromal populations do not form clusters
# that later get read as resistance clones.
filter_marker_sets <- function(log_counts, exclude_markers, min_expr = 1) {
  empty <- list(
    keep = colnames(log_counts),
    removed = data.frame(
      cell_id = character(0),
      set = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    )
  )

  if (is.null(exclude_markers) || length(exclude_markers) == 0L) {
    return(empty)
  }
  if (!is.list(exclude_markers) || is.null(names(exclude_markers))) {
    cli::cli_abort(
      "{.arg exclude_markers} must be a named list of gene vectors."
    )
  }

  sets <- lapply(exclude_markers, function(g) {
    intersect(g, rownames(log_counts))
  })
  sets <- sets[lengths(sets) > 0L]
  if (length(sets) == 0L) {
    cli::cli_warn(
      "None of the genes in {.arg exclude_markers} were found in the matrix."
    )
    return(empty)
  }

  scores <- vapply(
    sets,
    function(g) Matrix::colMeans(log_counts[g, , drop = FALSE]),
    numeric(ncol(log_counts))
  )
  if (is.null(dim(scores))) {
    scores <- matrix(
      scores,
      ncol = 1L,
      dimnames = list(colnames(log_counts), names(sets))
    )
  }

  best <- max.col(scores, ties.method = "first")
  best_score <- scores[cbind(seq_len(nrow(scores)), best)]
  hit <- best_score >= min_expr

  removed <- data.frame(
    cell_id = colnames(log_counts)[hit],
    set = names(sets)[best[hit]],
    score = round(best_score[hit], 3),
    stringsAsFactors = FALSE
  )

  list(keep = colnames(log_counts)[!hit], removed = removed)
}


# Detect doublets with scDblFinder, one patient at a time.
#
# Doublet detection is not a QC threshold: two cells captured in one droplet can
# have perfectly normal gene counts and mitochondrial fractions, so the MAD
# fences cannot see them. They inflate the upper gene-count tail, and a doublet
# of two different clones gets read as a "resistant clone" of its own.
#
# Returns the cells classified as doublets (empty when scDblFinder is absent).
detect_doublets <- function(counts, patient, seed = 42L, min_cells = 30L) {
  none <- data.frame(
    cell_id = character(0),
    patient = character(0),
    score = numeric(0),
    stringsAsFactors = FALSE
  )

  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    cli::cli_warn(c(
      "{.pkg scDblFinder} is not installed; skipping doublet detection.",
      "i" = "Install it with {.code{BiocManager::install(\"scDblFinder\")}}.",
      "i" = "The upper gene-count fence still removes most high-complexity doublets."
    ))
    return(none)
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    cli::cli_warn(c(
      "{.pkg SingleCellExperiment} is not installed; skipping doublet detection."
    ))
    return(none)
  }

  patient <- as.character(patient)
  out <- vector("list", 0L)
  for (p in unique(patient)) {
    idx <- which(patient == p)
    # scDblFinder needs enough neighbours to build its artificial doublets.
    if (length(idx) < min_cells) {
      cli::cli_warn(
        "Patient {.val {p}} has {length(idx)} cell{?s}; skipping doublet detection."
      )
      next
    }
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = counts[, idx, drop = FALSE])
    )
    set.seed(as.integer(seed))
    sce <- suppressMessages(scDblFinder::scDblFinder(sce, verbose = FALSE))
    is_dbl <- sce$scDblFinder.class == "doublet"
    if (any(is_dbl)) {
      out[[length(out) + 1L]] <- data.frame(
        cell_id = colnames(sce)[is_dbl],
        patient = p,
        score = round(sce$scDblFinder.score[is_dbl], 3),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0L) {
    return(none)
  }
  do.call(rbind, out)
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


# Copy the uwot UMAP settings used by the original Seurat pipeline: the default
# RunUMAP settings cost seconds of fixed overhead on small datasets, and this
# embedding is only used for visualisation.
umap_from_pca <- function(so, dims, seed = 42, n_threads = 1L) {
  pca <- Seurat::Embeddings(so, "pca")[, seq_len(dims), drop = FALSE]

  n_neighbors <- min(15L, max(2L, nrow(pca) - 1L))
  # UMAP is the single most expensive step of the pipeline, and its cost is
  # almost all SGD epochs. Below ~5k cells the default 200 epochs are cheap;
  # at 20k cells they cost ~60 s on a build without OpenMP, for a 2D layout
  # that is only used for visualisation. The layout is stable well before 200
  # epochs, so the budget is scaled down with cell count.
  n_epochs <- if (nrow(pca) > 15000L) {
    50L
  } else if (nrow(pca) > 5000L) {
    100L
  } else {
    200L
  }
  coords <- uwot::umap(
    X = pca,
    n_neighbors = n_neighbors,
    n_components = 2,
    n_epochs = n_epochs,
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
cluster_within_patient <- function(
  so,
  patient_col = "patient",
  sample_col = NULL,
  dims = 20,
  resolution = 0.8,
  cluster_algorithm = 1L,
  integrate = c("auto", "none", "harmony"),
  seed = 42,
  progress_cb = NULL
) {
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

    # FindNeighbors() needs a few cells: with one or two it produces NA edges
    # and errors out. Tiny biopsies also carry no evidence of heterogeneity, so
    # they are reported as a single clone rather than failing the whole run.
    if (nrow(emb_p) < 3L) {
      out[[i]] <- data.frame(
        cell_id = rownames(emb_p),
        clone_id = rep("1", nrow(emb_p)),
        patient = p,
        stringsAsFactors = FALSE
      )
      if (is.function(progress_cb)) {
        progress_cb(p, i, length(patients))
      }
      next
    }

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
    cl_raw <- Seurat::FindClusters(
      nn$snn,
      resolution = resolution,
      algorithm = cluster_algorithm,
      verbose = FALSE
    )
    # Seurat 5 returns a one-column data.frame when given a Graph; older
    # versions return a factor. Normalise both to a character vector.
    if (is.data.frame(cl_raw)) {
      cl_names <- rownames(cl_raw)
      cl <- as.character(cl_raw[[1]])
    } else {
      cl_names <- names(cl_raw)
      cl <- as.character(cl_raw)
    }
    # Labels are attached positionally, so a differently ordered (but
    # same-length) return value would silently assign clones to the wrong
    # cells -- check the cell names, not just the count.
    default_names <- as.character(seq_len(nrow(emb_p)))
    if (
      !is.null(cl_names) &&
        !identical(cl_names, rownames(emb_p)) &&
        !identical(cl_names, default_names)
    ) {
      cli::cli_abort(c(
        "Cluster labels are not aligned with the input cells.",
        "i" = "They must follow the row order of the PCA embedding."
      ))
    }
    if (length(cl) != nrow(emb_p)) {
      cli::cli_abort(
        "Clustering returned {length(cl)} label{?s} for {nrow(emb_p)} cell{?s}."
      )
    }

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


# Shared core: variable features -> scaling -> PCA -> clustering -> 2D embedding.
#
# One implementation for every entry point. preprocess_sc() calls it after QC,
# while annotate_clones() / plot_seurat_clustering() call it directly on an
# already-normalised object (they keep their own min_cells/min_features
# filters). Keeping a single core is what makes "the new pipeline" the only
# pipeline: cluster ids, embedding settings and the small-dataset PCA cap
# cannot drift between the exported functions.
#
# The input must already be normalised (LogNormalize) and carry a `patient`
# column in its metadata; a single-patient matrix simply gets
# `patient = "patient1"`.
cluster_embedding <- function(
  so,
  nfeatures = 2000L,
  dims = 20L,
  resolution = 0.8,
  method = c("umap", "tsne"),
  seed = 42L,
  cluster_algorithm = 1L,
  variable_selection = c("vst", "dispersion", "mvp"),
  sample_col = NULL,
  integrate = c("auto", "none", "harmony"),
  cluster_scope = c("global", "patient"),
  progress_cb = NULL
) {
  method <- rlang::arg_match(method)
  variable_selection <- rlang::arg_match(variable_selection)
  integrate <- rlang::arg_match(integrate)
  cluster_scope <- rlang::arg_match(cluster_scope)
  set.seed(seed)
  notify <- function(phase, i = NULL, n = NULL) {
    if (is.function(progress_cb)) {
      progress_cb(phase, i, n)
    }
  }

  notify("variable-features")
  nfeatures <- min(nfeatures, nrow(so))
  so <- Seurat::FindVariableFeatures(
    so,
    selection.method = variable_selection,
    nfeatures = nfeatures
  )

  notify("scaling")
  so <- Seurat::ScaleData(so, features = Seurat::VariableFeatures(so))

  # PCA needs npcs < min(cells, variable features); small datasets get capped
  # instead of erroring out.
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

  notify("clustering")
  if (cluster_scope == "global") {
    # Original PERCEPTION behaviour: one clustering over every cell, so a clone
    # is a transcriptional state shared across patients. Harmony is not applied
    # here -- it only makes sense inside a patient, which is the "patient" scope.
    so <- Seurat::FindNeighbors(
      so,
      dims = seq_len(npcs),
      reduction = "pca",
      verbose = FALSE
    )
    so <- Seurat::FindClusters(
      so,
      resolution = resolution,
      algorithm = cluster_algorithm,
      verbose = FALSE
    )
    clustered <- list(
      labels = data.frame(
        cell_id = colnames(so),
        clone_id = as.character(so$seurat_clusters),
        patient = as.character(so$patient),
        stringsAsFactors = FALSE
      ),
      used_harmony = FALSE
    )
  } else {
    clustered <- cluster_within_patient(
      so,
      patient_col = "patient",
      sample_col = sample_col,
      dims = npcs,
      resolution = resolution,
      cluster_algorithm = cluster_algorithm,
      integrate = integrate,
      seed = seed,
      progress_cb = function(p, i, n) notify("clustering", i, n)
    )
  }

  notify(if (method == "umap") "umap" else "tsne")
  coords <- if (method == "umap") {
    n_threads <- max(
      1L,
      min(4L, as.integer(parallel::detectCores()), na.rm = TRUE)
    )
    umap_from_pca(so, npcs, seed = seed, n_threads = n_threads)
  } else {
    tsne_from_pca(so, npcs, seed = seed)
  }

  # Attach the embedding as a real reduction so Seurat::DimPlot()/Idents() work
  # for callers that use the object directly (plot_seurat_clustering()).
  key <- if (method == "umap") "UMAP_" else "TSNE_"
  so[[method]] <- Seurat::CreateDimReducObject(
    embeddings = coords,
    key = key,
    assay = Seurat::DefaultAssay(so)
  )

  labels <- clustered$labels
  # Idents follow the clone labels, not the cell names Seurat would default to.
  idents <- labels$clone_id[match(colnames(so), labels$cell_id)]
  so$seurat_clusters <- factor(idents)
  Seurat::Idents(so) <- "seurat_clusters"

  list(
    so = so,
    labels = labels,
    npcs = npcs,
    used_harmony = clustered$used_harmony
  )
}


#' Preprocess raw single-cell data into clone labels
#'
#' End-to-end preprocessing for users who only have a raw count matrix: quality
#' control, normalisation, variable-feature selection, PCA, and clustering,
#' followed by a QC report. No tissue-specific thresholds are
#' required; the QC fences adapt to each patient's own distribution.
#'
#' Clustering is performed on a shared PCA embedding. By default
#' (\code{cluster_scope = "global"}) every cell is clustered together, which is
#' what the original PERCEPTION method does: a clone is a transcriptional state
#' that can appear in several patients. With \code{cluster_scope = "patient"}
#' neighbours are searched inside each patient only, so a clone can never span
#' two patients -- useful when clone labels are meant to be patient-specific.
#' For a single-patient dataset the two settings are identical.
#'
#' When a patient has more than one sample (e.g. longitudinal biopsies) and
#' \code{cluster_scope = "patient"}, harmony is applied inside that patient to
#' prevent one biological clone from being split in two by a technical batch
#' effect. \pkg{harmony} is an optional dependency: without it, the step is
#' skipped with a warning.
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
#'        deviations. Default = 5. Tumour tissue has genuinely higher
#'        mitochondrial content and larger cells than blood, so a 3-MAD fence
#'        removes real cells: against CellRanger's own calls on a real lung
#'        squamous carcinoma sample, 5 MAD recovers ~88% of its cells versus
#'        ~78% at 3 MAD, with precision unchanged.
#' @param min_genes Integer. Minimum detected genes per cell. Default = 200.
#' @param max_mito Numeric. Hard cap on the mitochondrial fraction (%). The
#'        effective threshold is the smaller of this and the MAD fence.
#'        Default = 25.
#' @param max_genes Integer. Fixed upper bound on detected genes per cell. When
#'        \code{NULL} (default) the bound is estimated per patient as
#'        \code{median + nmads * MAD}, which is what removes empty droplets and
#'        high-complexity doublets. Set it to reproduce a fixed threshold such
#'        as the 2,500 genes used in most tutorials.
#' @param remove_doublets Logical. When \code{TRUE}, cells classified as
#'        doublets by \pkg{scDblFinder} are removed after QC, one patient at a
#'        time. Doublets can pass every QC threshold -- they have normal gene
#'        counts and mitochondrial fractions -- yet they inflate the
#'        high-complexity tail and can form a spurious "resistant clone" made of
#'        two fused clones. Requires the optional \pkg{scDblFinder} package;
#'        without it the step is skipped with a warning. Default = \code{FALSE},
#'        because detection is probabilistic and should be opted into.
#' @param skip_qc Logical. When \code{TRUE}, the MAD fences are not applied and
#'        every cell is kept; the QC table then only describes the input. This
#'        reproduces the original PERCEPTION pipeline, which expects an
#'        already-curated matrix (e.g. malignant cells downloaded from a paper),
#'        where an extra per-patient fence only perturbs the result. Default
#'        \code{FALSE} -- raw 10x data always needs the fences.
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
#'        \code{cluster_size_class = "small"}. Default = 10.
#' @param exclude_markers Optional named list of gene vectors, e.g.
#'        \code{list(immune = c("PTPRC", "CD3D"), stromal = c("COL1A1"))}.
#'        Cells whose strongest marker-set score reaches \code{exclude_min_expr}
#'        are removed before clustering. This is an **exclusion** rule for
#'        reducing non-target cells, not a tumour-cell classifier: epithelial
#'        markers do not imply malignancy. Default = \code{NULL} (keep all
#'        cells).
#' @param exclude_min_expr Numeric. Mean log-normalised expression a cell must
#'        reach in a marker set to be excluded. Default = 1.
#' @param integrate Character. Sample integration strategy within a patient:
#'        \code{"auto"} (default) applies harmony when a patient has multiple
#'        samples, \code{"none"} never integrates, \code{"harmony"} behaves
#'        like \code{"auto"}. Only used when \code{cluster_scope = "patient"}.
#' @param cluster_scope Character. \code{"global"} (default) clusters all cells
#'        together, matching the original PERCEPTION method, so one clone can
#'        span several patients. \code{"patient"} clusters inside each patient
#'        only, so clone labels are patient-specific. Identical for a
#'        single-patient dataset.
#' @param seed Integer. Random seed. Default = 42.
#' @param progress_cb Optional function called as \code{progress_cb(phase, i, n)}
#'        between stages. Used to drive progress bars in the Shiny app.
#'
#' @section What "clone" means here:
#' In the cell-level workflow, clone labels are **operationally defined** by
#' transcriptomic clustering. Cells are grouped by transcriptional similarity
#' within each patient, and the resulting cluster assignments are used as
#' surrogate clone labels for downstream analysis. These labels do not
#' necessarily represent genetically validated clonal populations and should not
#' be interpreted as evidence of shared genomic ancestry. CNV- or genotype-based
#' evidence is required to call true clones.
#'
#' @return A list with:
#' \describe{
#'   \item{cell_labels}{Data frame with \code{cell_id}, \code{clone_id},
#'         \code{clone_key} (unique \code{"patient@@clone"} key),
#'         \code{patient}, \code{cluster_size_class} (\code{"large"} /
#'         \code{"small"}, a cell-count category and nothing more), and
#'         \code{dim_1}/\code{dim_2}. This is the \code{clone_labels} argument
#'         of \code{prepare_data()}.}
#'   \item{qc}{Data frame. Per-patient QC summary, including
#'         \code{n_marker_removed} when \code{exclude_markers} was used and
#'         \code{n_doublet_removed} when \code{remove_doublets = TRUE}.}
#'   \item{qc_metrics}{Data frame. Per-cell QC metrics before filtering.}
#'   \item{excluded_cells}{Data frame. Cells removed by \code{exclude_markers}
#'         (columns \code{cell_id}, \code{set}, \code{score}); empty when no
#'         marker filtering was requested.}
#'   \item{doublet_cells}{Data frame. Cells removed by \code{remove_doublets}
#'         (columns \code{cell_id}, \code{patient}, \code{score}); empty when
#'         the step was not requested or no doublets were found.}
#'   \item{used_harmony}{Logical. Whether harmony was applied.}
#'   \item{label_source}{Character. Always \code{"auto_cluster"}.}
#' }
#'
#' @examples
#' \dontrun{
#'   labels <- preprocess_sc(counts, patient_labels = cell_patients)
#'   prepared <- prepare_data(
#'     expression_matrix = counts,
#'     clone_labels = labels$cell_labels
#'   )
#' }
#'
#' @export
preprocess_sc <- function(
  counts,
  patient_labels,
  sample_labels = NULL,
  nmads = 5,
  min_genes = 200,
  max_mito = 25,
  max_genes = NULL,
  skip_qc = FALSE,
  remove_doublets = FALSE,
  exclude_markers = NULL,
  exclude_min_expr = 1,
  nfeatures = 2000,
  dims = 20,
  method = c("umap", "tsne"),
  resolution = 0.8,
  min_clone_cells = 10,
  integrate = c("auto", "none", "harmony"),
  cluster_scope = c("global", "patient"),
  seed = 42,
  progress_cb = NULL
) {
  integrate <- match.arg(integrate)
  cluster_scope <- match.arg(cluster_scope)
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
  # as.character() drops names, so the names must be (re)assigned after it --
  # the QC step below indexes this vector by cell name.
  patient_labels <- as.character(patient_labels)
  names(patient_labels) <- colnames(counts)

  if (!is.null(sample_labels)) {
    if (length(sample_labels) != ncol(counts)) {
      cli::cli_abort(c(
        "{.arg sample_labels} must have one entry per cell.",
        "x" = "Got {length(sample_labels)} labels for {ncol(counts)} cells."
      ))
    }
    sample_labels <- as.character(sample_labels)
    names(sample_labels) <- colnames(counts)
  }

  # --- 2. Quality control --------------------------------------------------
  notify("qc")
  qc <- qc_filter_cells(
    counts,
    patient_labels,
    nmads = nmads,
    min_genes = min_genes,
    max_mito = max_mito,
    max_genes = max_genes,
    skip = skip_qc
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

  # --- 2b. Optional doublet removal ---------------------------------------
  # Doublets survive every QC fence (normal genes, normal mito), so they need
  # their own detector. Runs on the QC-passed counts, before normalisation.
  doublet_cells <- data.frame(
    cell_id = character(0),
    patient = character(0),
    score = numeric(0),
    stringsAsFactors = FALSE
  )
  if (isTRUE(remove_doublets)) {
    notify("doublets")
    doublet_cells <- detect_doublets(counts, patient_labels, seed = seed)
    if (nrow(doublet_cells) > 0L) {
      keep_dbl <- setdiff(colnames(counts), doublet_cells$cell_id)
      counts <- counts[, keep_dbl, drop = FALSE]
      patient_labels <- patient_labels[keep_dbl]
      if (!is.null(sample_labels)) {
        sample_labels <- sample_labels[keep_dbl]
      }
      cli::cli_alert_info(
        "Removed {nrow(doublet_cells)} doublet{?s} detected by scDblFinder."
      )
      if (ncol(counts) < 10L) {
        cli::cli_abort(c(
          "Only {ncol(counts)} cell{?s} left after doublet removal.",
          "i" = "Set {.arg remove_doublets = FALSE} to keep them."
        ))
      }
    }
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

  # --- 3b. Optional exclusion of non-target cells --------------------------
  # Runs on normalised values so the marker scores are comparable across cells.
  # Excluding (rather than keeping only) is deliberate: "epithelial" is not the
  # same as "malignant" and no marker rule can establish malignancy.
  excluded_cells <- data.frame(
    cell_id = character(0),
    set = character(0),
    score = numeric(0),
    stringsAsFactors = FALSE
  )
  marker_removed <- integer(0)
  if (!is.null(exclude_markers)) {
    notify("filtering")
    logdata <- Seurat::GetAssayData(so, assay = "RNA", layer = "data")
    filtered <- filter_marker_sets(logdata, exclude_markers, exclude_min_expr)
    if (nrow(filtered$removed) > 0L) {
      so <- subset(so, cells = filtered$keep)
      excluded_cells <- filtered$removed
      marker_removed <- table(patient_labels[excluded_cells$cell_id])

      # A marker rule that removes most of the data is almost certainly
      # mis-specified, so surface it instead of silently continuing.
      pct_excluded <- nrow(excluded_cells) / length(patient_labels)
      if (pct_excluded > 0.5) {
        cli::cli_warn(c(
          "{round(100 * pct_excluded)}% of cells were removed by {.arg exclude_markers}.",
          "i" = "Check that the marker sets are specific to the populations you meant to drop.",
          "i" = "Inspect {.code{table(excluded_cells$set)}} and consider a higher {.arg exclude_min_expr}."
        ))
      }

      cli::cli_alert_info(
        "Excluded {nrow(excluded_cells)} cell{?s} expressing non-target marker sets."
      )
    }
    if (ncol(so) < 10L) {
      cli::cli_abort(c(
        "Only {ncol(so)} cell{?s} left after marker filtering.",
        "i" = "Relax {.arg exclude_min_expr} or drop {.arg exclude_markers}."
      ))
    }
  }

  # --- 4. Variable features, PCA, clustering and the 2D embedding ----------
  # Same core as annotate_clones()/plot_seurat_clustering(), so the exported
  # entry points cannot drift apart.
  clustered <- cluster_embedding(
    so,
    nfeatures = nfeatures,
    dims = dims,
    resolution = resolution,
    method = method,
    seed = seed,
    variable_selection = "vst",
    sample_col = if (is.null(sample_labels)) NULL else "sample",
    integrate = integrate,
    cluster_scope = cluster_scope,
    progress_cb = notify
  )
  so <- clustered$so
  labels <- clustered$labels
  coords <- Seurat::Embeddings(so, method)

  # --- 5. Assemble labels --------------------------------------------------
  notify("assembling")
  labels$dim_1 <- coords[labels$cell_id, 1]
  labels$dim_2 <- coords[labels$cell_id, 2]

  # Cluster ids restart at 0 for every patient, so they are only unique within
  # a patient. `clone_key` is the globally unique form used downstream.
  labels$clone_key <- paste(labels$patient, labels$clone_id, sep = "@@")

  clone_size <- table(labels$clone_key)[labels$clone_key]
  # "large"/"small" rather than "high"/"low" on purpose: this only reflects how
  # many cells the cluster contains, NOT how credible the clone is.
  labels$cluster_size_class <- ifelse(
    as.integer(clone_size) >= min_clone_cells,
    "large",
    "small"
  )
  rownames(labels) <- NULL

  # Carry provenance on the object itself, so prepare_data(clone_labels=)
  # can report "auto_cluster" instead of assuming the labels were hand-made.
  attr(labels, "label_source") <- "auto_cluster"

  qc$report$n_clones <- vapply(
    qc$report$patient,
    function(p) length(unique(labels$clone_id[labels$patient == p])),
    integer(1)
  )

  # Per-patient count of cells removed by exclude_markers, so the QC table
  # reports the full n_input -> n_kept -> n_marker_removed -> clones chain.
  qc$report$n_marker_removed <- 0L
  if (length(marker_removed) > 0L) {
    removed <- as.integer(marker_removed[qc$report$patient])
    removed[is.na(removed)] <- 0L
    qc$report$n_marker_removed <- removed
  }

  qc$report$n_doublet_removed <- 0L
  if (nrow(doublet_cells) > 0L) {
    removed <- as.integer(table(doublet_cells$patient)[qc$report$patient])
    removed[is.na(removed)] <- 0L
    qc$report$n_doublet_removed <- removed
  }

  list(
    cell_labels = labels,
    qc = qc$report,
    qc_metrics = qc$metrics,
    excluded_cells = excluded_cells,
    doublet_cells = doublet_cells,
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
#' @param metric Character. \code{"all"} (default) pools every patient and
#'        stacks the three metrics in one column -- the only readable layout
#'        when a cohort has dozens of patients. Pass one of
#'        \code{"Detected genes"}, \code{"Total counts"} or
#'        \code{"Mitochondrial \%"} to get one panel per patient (wrapped into
#'        rows), each with that patient's own threshold line; this is the
#'        layout to use when per-patient QC needs to be inspected.
#' @param ncol Integer. Number of facet columns in the per-patient layout.
#'        Default = 6.
#'
#' @return A \code{ggplot} object, faceted by metric.
#'
#' @export
plot_qc_report <- function(
  qc_metrics,
  kept_cells,
  qc_report,
  base_size = 12,
  metric = c(
    "all",
    "Detected genes",
    "Total counts",
    "Mitochondrial %"
  ),
  ncol = 6L
) {
  metric <- match.arg(metric)
  metrics <- qc_metrics
  metrics$status <- factor(
    ifelse(metrics$cell_id %in% kept_cells, "kept", "removed"),
    levels = c("kept", "removed")
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
  long$metric <- factor(long$metric, levels = unique(long$metric))

  thresholds <- data.frame(
    patient = rep(qc_report$patient, 2L),
    metric = rep(
      c("Detected genes", "Mitochondrial %"),
      each = nrow(qc_report)
    ),
    value = c(qc_report$genes_min, qc_report$mito_threshold),
    stringsAsFactors = FALSE
  )

  if (metric != "all") {
    long <- long[long$metric == metric, , drop = FALSE]
    thresholds <- thresholds[thresholds$metric == metric, , drop = FALSE]
  }

  # A violin needs at least two cells. ggplot silently drops smaller groups
  # (and warns), which is why a sample with a single removed cell showed an
  # empty slot instead of a mark. Those cells are drawn as points instead.
  group_key <- paste(long$patient, long$metric, long$status, sep = "\r")
  group_n <- ave(seq_len(nrow(long)), group_key, FUN = length)
  violin_df <- long[group_n >= 2L, , drop = FALSE]
  point_df <- long[group_n < 2L, , drop = FALSE]

  n_cells <- nrow(metrics)
  n_kept <- sum(metrics$status == "kept")
  n_samples <- length(unique(metrics$patient))
  subtitle <- sprintf(
    "%s cells · %s retained · %s removed (%.1f%%) · %s sample%s",
    format(n_cells, big.mark = ","),
    format(n_kept, big.mark = ","),
    format(n_cells - n_kept, big.mark = ","),
    if (n_cells > 0) 100 * (n_cells - n_kept) / n_cells else 0,
    n_samples,
    if (n_samples == 1L) "" else "s"
  )

  caption <- if (metric == "all") {
    paste0(
      "Violins split cells by QC outcome (blue = retained, red = removed). ",
      "Per-sample cut-offs are listed in the table.\n",
      "Total counts has no cut-off and is shown for context."
    )
  } else {
    paste0(
      "Dashed lines mark each sample's effective cut-off (median ± n×MAD).\n",
      "Total counts has no cut-off and is shown for context."
    )
  }

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$status, y = .data$value, fill = .data$status)
  ) +
    ggplot2::geom_violin(
      data = violin_df,
      scale = "width",
      width = 0.62,
      linewidth = 0.25
    ) +
    # Groups too small for a density estimate (1 cell) still deserve a mark.
    ggplot2::geom_point(
      data = point_df,
      size = 1.1,
      position = ggplot2::position_jitter(width = 0.06, height = 0, seed = 42)
    ) +
    # Boxplots stay slim: with only two x categories a wide violin reads as a
    # flat blob, so the fill is deliberately narrower than ggplot's 0.9 default.
    ggplot2::geom_boxplot(
      data = violin_df,
      width = 0.1,
      outlier.size = 0.4,
      fill = "white"
    ) +
    # Kept = theme blue, removed = theme red; the same pairing used by the
    # response plots (red = "worse", blue = "better").
    ggplot2::scale_fill_manual(
      values = c(kept = "#165C91", removed = "#c13232"),
      labels = c(kept = "Retained", removed = "Removed"),
      drop = FALSE
    ) +
    # A metric that is constant within a panel (e.g. 0% mitochondrial reads
    # when no MT genes were detected) collapses the range to zero, which made
    # every y tick print the same value. The additive margin gives such a
    # panel a real interval without moving the data.
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = 0.05, add = 0.5)
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = NULL,
      title = "Quality control metrics by cell status",
      subtitle = subtitle,
      caption = caption
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.22,
        colour = "grey91"
      ),
      panel.border = ggplot2::element_rect(linewidth = 0.4, colour = "grey35"),
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(0.85)
      ),
      legend.position = "top",
      legend.key.size = grid::unit(0.9, "lines"),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.05),
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "grey30",
        size = ggplot2::rel(0.82),
        hjust = 0.5
      ),
      plot.caption = ggplot2::element_text(
        colour = "grey40",
        size = ggplot2::rel(0.7),
        hjust = 0
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )

  if (metric == "all") {
    # Pooled view: per-patient thresholds would all overlap on one panel, so
    # they are omitted (the QC table reports them per patient instead).
    p + ggplot2::facet_wrap(~metric, ncol = 1L, scales = "free")
  } else {
    # Per-patient view, wrapped into rows like the clone lollipop: one panel
    # per patient with its own threshold line, so nobody ends up as a sliver
    # in a 40-column strip. scales = "free" also removes the empty status
    # slot in a sample where nothing was removed.
    p +
      ggplot2::geom_hline(
        data = thresholds,
        ggplot2::aes(yintercept = .data$value),
        linetype = "dashed",
        colour = "grey35",
        linewidth = 0.3
      ) +
      ggplot2::facet_wrap(~patient, ncol = ncol, scales = "free") +
      ggplot2::labs(y = metric)
  }
}
