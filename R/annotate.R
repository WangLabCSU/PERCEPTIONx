# encoding: UTF-8

#' PERCEPTIONx Patient Data Annotation Functions
#'
#' Functions for annotating single-cell data with clone and patient information,
#' and preparing patient data for the prediction pipeline.
#'
#' @name annotate_perception
#' @keywords internal
NULL


# Internal: normalise a single matrix and build clones with the shared core.
#
# This used to duplicate NormalizeData -> FindVariableFeatures -> ScaleData ->
# PCA -> FindNeighbors -> FindClusters -> UMAP/t-SNE. It now delegates to
# cluster_embedding() in preprocess.R, so annotate_clones(),
# plot_seurat_clustering() and preprocess_sc() share one implementation --
# including the small-dataset PCA cap and the adaptive UMAP epoch budget.
cluster_cells_seurat <- function(
  expression_matrix,
  method,
  min_cells,
  min_features,
  nfeatures,
  dims,
  resolution,
  seed,
  progress_cb = NULL,
  cluster_algorithm = 1L,
  variable_selection = c("vst", "dispersion", "mvp")
) {
  variable_selection <- rlang::arg_match(variable_selection)
  notify <- function(phase) {
    if (is.function(progress_cb)) progress_cb(phase)
  }

  notify("normalizing")
  so <- Seurat::CreateSeuratObject(
    counts = expression_matrix,
    project = "PERCEPTIONx",
    min.cells = min_cells,
    min.features = min_features
  )
  so <- Seurat::NormalizeData(
    so,
    normalization.method = "LogNormalize",
    scale.factor = 10000L
  )
  # Single patient: clone ids only need to be unique inside this matrix.
  so$patient <- "patient1"

  cluster_embedding(
    so,
    nfeatures = nfeatures,
    dims = dims,
    resolution = resolution,
    method = method,
    seed = seed,
    cluster_algorithm = cluster_algorithm,
    variable_selection = variable_selection,
    progress_cb = function(phase, i = NULL, n = NULL) notify(phase)
  )$so
}


#' Annotate cells with clone IDs via Seurat clustering
#'
#' Performs Seurat clustering on a single-cell expression matrix and returns
#' a mapping of each cell to its cluster (clone) ID. This matches the original
#' PERCEPTIONx pipeline where Seurat clusters define transcriptional subclones.
#'
#' @param method Character. Dimensionality reduction method. One of \code{"umap"}
#'        (default) or \code{"tsne"}. UMAP is faster and preserves global structure
#'        better; t-SNE emphasizes local neighborhoods.
#' @param expression_matrix Matrix. Gene expression matrix with genes as rows
#'        and cells as columns. Raw counts or normalized values are both accepted.
#' @param min_cells Integer. Minimum cells per feature. Default = 3.
#' @param min_features Integer. Minimum features per cell. Default = 200.
#' @param nfeatures Integer. Number of variable features. Default = 2000.
#' @param dims Integer. Number of PCA dimensions for clustering. Default = 10.
#' @param resolution Numeric. Clustering resolution. Default = 0.8.
#' @param seed Integer. Random seed for reproducibility. Default = 42.
#' @param progress_cb Optional function called with the current Seurat stage
#'        name (e.g. \code{"pca"}, \code{"clustering"}) between pipeline steps.
#'        Used to show live progress in the Shiny app.
#' @param cluster_algorithm Integer. Community-detection algorithm passed to
#'        \code{Seurat::FindClusters()}. Default = 1 (Louvain, no extra
#'        dependency). Use 4 (Leiden) only if the \pkg{leidenbase} package is
#'        installed.
#' @param variable_selection Character. Variable-feature selection method for
#'        \code{Seurat::FindVariableFeatures()}. One of \code{"vst"} (default,
#'        most accurate), \code{"dispersion"} or \code{"mvp"} (faster on large
#'        matrices at a small cost in feature quality).
#'
#' @return A data frame with columns: \code{cell_id}, \code{clone_id}, and
#'         \code{dim_1}, \code{dim_2} (2D embedding coordinates for visualization).
#'
#' @examples
#' \dontrun{
#'   cell_clone_map <- annotate_clones(patient_expression)
#'   cell_clone_map <- annotate_clones("tsne", patient_expression)
#' }
#'
#' @export
annotate_clones <- function(
  method = c("umap", "tsne"),
  expression_matrix,
  min_cells = 3L,
  min_features = 200L,
  nfeatures = 2000L,
  dims = 10L,
  resolution = 0.8,
  seed = 42L,
  progress_cb = NULL,
  cluster_algorithm = 1L,
  variable_selection = c("vst", "dispersion", "mvp")
) {
  method <- rlang::arg_match(method)
  variable_selection <- rlang::arg_match(variable_selection)

  rlang::check_installed("Seurat", "for clone annotation.")

  # TODO: deprecated this function
  so <- cluster_cells_seurat(
    expression_matrix,
    method,
    min_cells,
    min_features,
    nfeatures,
    dims,
    resolution,
    seed,
    progress_cb,
    cluster_algorithm,
    variable_selection
  )

  cluster_ids <- Seurat::Idents(so)
  emb <- Seurat::Embeddings(so, reduction = method)

  data.table::data.table(
    cell_id = names(cluster_ids),
    clone_id = as.character(cluster_ids),
    dim_1 = emb[, 1L],
    dim_2 = emb[, 2L]
  )
}


#' Build clone abundance table from cell-clone mapping
#'
#' Computes the number of cells per clone per patient, producing the
#' \code{clone_counts} data frame required by \code{predict_patients()}.
#'
#' @param cell_clone_map Data frame with columns \code{cell_id} and \code{clone_id}.
#' @param patient_ids Character vector. Patient ID for each cell, in the same
#'        order as rows in cell_clone_map.
#'
#' @return A data frame with first column \code{patients} and remaining columns
#'         as clone IDs with cell counts as values.
#' @export
build_clone_counts <- function(cell_clone_map, patient_ids) {
  if (!all(c("cell_id", "clone_id") %in% colnames(cell_clone_map))) {
    rlang::abort("cell_clone_map must have columns 'cell_id' and 'clone_id'.")
  }

  if (length(patient_ids) != nrow(cell_clone_map)) {
    rlang::abort(
      "Length of patient_ids must match number of rows in cell_clone_map."
    )
  }

  map_dt <- data.table::data.table(
    patient = patient_ids,
    # Clone names are only unique *within* a patient (e.g. "Kydar01_c1" and
    # "Kydar02_c1" both shorten to "c1"). Keying the abundance columns on
    # "patient@@clone" keeps them aligned with the clone_expression_rnorm
    # columns and prevents abundance from being summed across patients.
    clone_key = paste(patient_ids, cell_clone_map$clone_id, sep = "@@")
  )

  all_patients <- unique(map_dt$patient)
  all_keys <- unique(map_dt$clone_key)

  # Count cells per (patient, clone) pair, then pivot wide so that every clone
  # key becomes a column of per-patient abundances.
  patient <- clone_key <- NULL # * ease checking NOTE
  counts <- map_dt[, list(n_cells = .N), by = list(patient, clone_key)]
  counts <- data.table::dcast(
    counts,
    patient ~ clone_key,
    value.var = "n_cells",
    fill = 0L
  )

  # dcast sorts the clone keys; restore first-appearance column order so the
  # columns line up with the clone_expression_rnorm columns.
  data.table::setnames(counts, "patient", "patients")
  data.table::setcolorder(counts, c("patients", all_keys))

  # Plain data.frame is the documented contract: predict_patients() subsets it
  # as clone_counts[i, clone_cols], which is data.frame (not data.table) syntax.
  as.data.frame(counts)
}


# Pull the expression matrix out of a Seurat object.
#
# `layer = NULL` means "counts when the object has them, otherwise data". That
# is also what the app does on upload, and it keeps the Seurat path on the same
# scale as the raw-counts path: averaging log-normalised values gives a
# different clone profile than averaging counts.
seurat_expression_matrix <- function(data, assay = NULL, layer = NULL) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(data)
  }
  get_layer <- function(nm) {
    tryCatch(
      Seurat::GetAssayData(data, assay = assay, layer = nm),
      error = function(e) {
        Seurat::GetAssayData(data, assay = assay, slot = nm)
      }
    )
  }

  if (!is.null(layer)) {
    return(get_layer(layer))
  }
  for (candidate in c("counts", "data")) {
    out <- tryCatch(get_layer(candidate), error = function(e) NULL)
    if (!is.null(out) && nrow(out) > 0L) {
      return(out)
    }
  }
  cli::cli_abort(c(
    "The Seurat object has neither a {.field counts} nor a {.field data} layer.",
    "i" = "Assay {.val {assay}}; supply an object that carries an expression matrix."
  ))
}


#' Prepare patient data for PERCEPTIONx prediction
#'
#' End-to-end preprocessing pipeline that takes raw single-cell expression data
#' and produces a rank-normalized subclone expression matrix and clone counts
#' table, ready for direct use with \code{predict_drugs()} and
#' \code{predict_patients()}.
#'
#' The pipeline performs:
#' \enumerate{
#'   \item Seurat clustering to define transcriptional subclones
#'   \item Cell-to-patient and cell-to-clone annotation
#'   \item Clone-level mean expression computation
#'   \item Rank normalization of clone expression
#'   \item Clone abundance table construction
#' }
#'
#' @param method Character. Dimensionality reduction method passed to
#'        \code{\link{annotate_clones}()}. One of \code{"umap"} (default) or
#'        \code{"tsne"}. UMAP is faster and preserves global structure;
#'        t-SNE emphasizes local neighborhoods.
#' @param expression_matrix Matrix. Gene expression matrix with genes as rows
#'        and cells as columns.
#' @param patient_mapping List or data frame. Patient-cell mapping in one of two formats:
#'   \describe{
#'     \item{List format}{Named list where each element is a patient ID and
#'        contains a character vector of cell IDs. Example:
#'        \code{list(Patient_1 = c("Cell_1", "Cell_2"), Patient_2 = c("Cell_3"))}}
#'     \item{Data frame format}{Metadata with cell ID and patient ID columns.
#'        Specify column names via \code{cell_col} and \code{patient_col}.
#'        Example: \code{data.frame(cell_id = c("Cell_1", "Cell_2"), patient_id = c("P1", "P1"))}}
#'   }
#'   If NULL and \code{clone_labels} carries a \code{patient} (or
#'   \code{patient_id}) column, the grouping is taken from those labels.
#'   Otherwise all cells are assigned to a single patient "patient1".
#' @param cell_col Character. Cell ID column name in patient_mapping data frame.
#'        Default = "cell_id". Only used when patient_mapping is a data frame.
#' @param patient_col Character. Patient ID column name in patient_mapping data frame.
#'        Default = "patient_id". Only used when patient_mapping is a data frame.
#' @param genes_to_use Character vector. Genes to retain in the output matrix.
#'        If NULL, all genes in the expression matrix are used.
#' @param seurat_resolution Numeric. Clustering resolution. Default = 0.8.
#' @param seurat_dims Integer. PCA dimensions for clustering. Default = 10.
#' @param seurat_nfeatures Integer. Variable features count. Default = 2000.
#' @param seurat_min_cells Integer. Minimum cells per feature. Default = 3.
#' @param seurat_min_features Integer. Minimum features per cell. Default = 200.
#'        Auto-adjusted to 10% of gene count if the expression matrix has fewer genes.
#' @param seurat_seed Integer. Random seed. Default = 42.
#'
#' @return A named list with:
#' \describe{
#'   \item{clone_expression_rnorm}{Matrix. Rank-normalized clone-level expression
#'         (genes as rows, patient_clone as columns). Ready for \code{predict_drugs()}.}
#'   \item{clone_counts}{Data frame. Clone abundance per patient. Ready for
#'         \code{predict_patients()}.}
#'   \item{cell_clone_map}{Data frame. Cell-to-clone mapping with columns
#'         cell_id, clone_id, patient, dim_1, dim_2.}
#'   \item{clone_viability_df_template}{Data frame. Template with patient and clone_id
#'         columns, ready to merge with \code{predict_drugs()} output.}
#'   \item{umap_coords}{Data frame. 2D embedding coordinates per cell
#'         (cell_id, dim_1, dim_2). Ready for \code{plot_tsne_response()}.}
#'   \item{reduction_method}{Character. The method used (\code{"umap"} or
#'         \code{"tsne"}).}
#' }
#'
#' @examples
#' \dontrun{
#'   # List format (same as Rmd)
#'   prepared <- prepare_data(
#'     expression_matrix = patient_scRNA,
#'     patient_mapping = cell_names_list,
#'     genes_to_use = GOI
#'   )
#'
#'   # Or data frame format (from metadata)
#'   metadata <- data.frame(cell_id = colnames(patient_scRNA), patient_id = patient_ids)
#'   prepared <- prepare_data(patient_scRNA, metadata)
#'
#'   # Parse patient ID from cell names (e.g., "P11_M_Barcode" -> "P11")
#'   metadata <- data.frame(Cell = c("P11_M_Barcode1", "P12_M_Barcode2"))
#'   prepared <- prepare_data(
#'     patient_scRNA, metadata,
#'     cell_col = "Cell",           # Custom column name
#'     parse_patient = TRUE,        # Parse from Cell column
#'     patient_sep = "_",           # Split by "_"
#'     patient_pos = 1              # Take first element
#'   )
#'   # Result: Patient IDs = "P11", "P12"
#'
#'   # Use directly with prediction functions
#'   clone_pred <- predict_drugs(models, prepared$clone_expression_rnorm)
#'   patient_pred <- predict_patients(clone_pred, prepared)
#'
#'   # After preprocess_sc(): the labels alone are enough, because they already
#'   # carry the patient of every cell. The unfiltered matrix can be passed as
#'   # is; cells without a label are dropped.
#'   labels <- preprocess_sc(counts, cell_patients)
#'   prepared <- prepare_data(
#'     expression_matrix = counts,
#'     clone_labels = labels$cell_labels
#'   )
#' }
#'
#' @param parse_patient Logical. If TRUE, parse patient ID from cell_col using separator.
#'   Default = FALSE. Auto-enabled if patient_sep or patient_pos is provided.
#'   Useful when cell names contain patient info (e.g., "P11_M_Barcode").
#' @param patient_sep Character. Separator to split cell_col for parsing patient ID.
#'   Default = "_". Providing this parameter auto-enables parse_patient.
#' @param patient_pos Integer. Position of patient ID after splitting.
#'   Default = 1 (first element). Providing this parameter auto-enables parse_patient.
#' @param skip_clustering Logical. If TRUE, skip the Seurat clustering step and
#'   treat every column of \code{expression_matrix} as one pre-defined clone.
#'   Use this when you already have a clone-level expression matrix (e.g. from a
#'   published study). Rank normalization and clone counts are still applied.
#'   No UMAP/t-SNE embedding is produced in this mode.
#' @param progress_cb Optional function called as \code{progress_cb(phase, i, n,
#'   drug)} with the current pipeline stage (e.g. \code{"seurat-pca"},
#'   \code{"mapping"}). Used to show live progress in the Shiny app.
#' @param data Optional \code{Seurat} object. When supplied, the expression
#'   matrix is read from \code{assay}/\code{layer}, clone labels from
#'   \code{clone_col} (or \code{Idents()}), and patient IDs from
#'   \code{patient_col}. Use this instead of \code{expression_matrix} plus
#'   \code{patient_mapping}.
#' @param clone_labels Optional data frame with \code{cell_id} and
#'   \code{clone_id} columns mapping every cell to a clone. Skips clustering.
#'   This is the \code{cell_labels} element returned by \code{preprocess_sc()},
#'   and its \code{patient} column makes \code{patient_mapping} optional. Cells
#'   of \code{expression_matrix} without a label (for example cells removed by
#'   QC) are dropped, so the unfiltered matrix can be passed directly.
#' @param assay Character. Assay to read when \code{data} is a Seurat object.
#'   Default = \code{Seurat::DefaultAssay(data)}.
#' @param layer Character. Layer (Seurat v5) or slot (Seurat v4) to read.
#'   Default = \code{NULL}: the \code{counts} layer when the object has one,
#'   otherwise \code{data}. Set it explicitly to override, e.g.
#'   \code{layer = "data"} to use normalised values.
#' @param clone_col Character. Metadata column in \code{data} holding clone
#'   labels. When \code{NULL}, \code{Seurat::Idents(data)} is used.
#'
#' @export
prepare_data <- function(
  method = c("umap", "tsne"),
  expression_matrix = NULL,
  patient_mapping = NULL,
  cell_col = "cell_id",
  patient_col = "patient_id",
  parse_patient = FALSE,
  patient_sep = "_",
  patient_pos = 1L,
  genes_to_use = NULL,
  seurat_resolution = 0.8,
  seurat_dims = 10L,
  seurat_nfeatures = 2000L,
  seurat_min_cells = 3L,
  seurat_min_features = 200L,
  seurat_seed = 42L,
  skip_clustering = FALSE,
  progress_cb = NULL,
  data = NULL,
  clone_labels = NULL,
  assay = NULL,
  layer = NULL,
  clone_col = NULL
) {
  method <- rlang::arg_match(method)
  notify_stage <- function(phase, i, n) {
    if (is.function(progress_cb)) progress_cb(phase, i, n, "")
  }
  cli::cli_h2("PERCEPTIONx Patient Data Preparation")
  cli::cli_ul("Reduction method: {.val {toupper(method)}}")

  # --- Seurat object input: unpack matrix + labels -------------------------
  # Keeps the rest of this function unchanged: a Seurat object is only a
  # structured container for an expression matrix plus cell labels, so it is
  # reduced to exactly the two inputs the pipeline already understands.
  label_source <- if (isTRUE(skip_clustering)) {
    "provided_matrix_columns"
  } else {
    "auto_cluster"
  }

  if (!is.null(data)) {
    if (!is.null(expression_matrix)) {
      cli::cli_abort(
        "Supply only one of {.arg data} and {.arg expression_matrix}."
      )
    }
    if (!inherits(data, "Seurat")) {
      cli::cli_abort(
        "{.arg data} must be a {.cls Seurat} object, not {.cls {class(data)}}."
      )
    }

    if (is.null(assay)) {
      assay <- Seurat::DefaultAssay(data)
    }
    expression_matrix <- seurat_expression_matrix(data, assay, layer)

    metadata <- data[[]]

    if (is.null(clone_col)) {
      clone_ids <- as.character(Seurat::Idents(data))
      label_source <- "seurat_idents"
      if (all(grepl("^[0-9]+$", clone_ids))) {
        cli::cli_warn(c(
          "Clone labels are Seurat cluster numbers ({.val {head(unique(clone_ids), 3)}}), not validated clones.",
          "i" = "Pass {.arg clone_col} to use a curated label column instead."
        ))
      }
    } else {
      if (!clone_col %in% colnames(metadata)) {
        cli::cli_abort(c(
          "Column {.val {clone_col}} not found in the Seurat object metadata.",
          "i" = "Available columns: {.val {colnames(metadata)}}."
        ))
      }
      clone_ids <- as.character(metadata[[clone_col]])
      label_source <- "seurat_metadata"
    }

    if (!patient_col %in% colnames(metadata)) {
      cli::cli_abort(c(
        "Column {.val {patient_col}} not found in the Seurat object metadata.",
        "i" = "Use {.arg patient_col} to name the column that holds patient IDs.",
        "i" = "Available columns: {.val {colnames(metadata)}}."
      ))
    }

    clone_labels <- data.frame(
      cell_id = colnames(expression_matrix),
      clone_id = clone_ids,
      stringsAsFactors = FALSE
    )
    # Name the patient column exactly as `patient_col` so the downstream
    # split() below finds it without a second lookup table.
    clone_labels[[patient_col]] <- as.character(metadata[[patient_col]])
    patient_mapping <- clone_labels[, c("cell_id", patient_col)]
    skip_clustering <- TRUE
  }

  if (is.null(expression_matrix)) {
    cli::cli_abort(
      "Supply either {.arg data} (a Seurat object) or {.arg expression_matrix}."
    )
  }

  # Labels from preprocess_sc() already carry the patient of every cell, so a
  # separate cell-to-patient table would just repeat that information. An
  # explicit mapping still wins when both are supplied.
  if (is.null(patient_mapping) && !is.null(clone_labels)) {
    patient_from_labels <- intersect(
      c("patient", "patient_id"),
      colnames(clone_labels)
    )
    if (length(patient_from_labels) > 0L) {
      patient_mapping <- data.frame(
        as.character(clone_labels$cell_id),
        as.character(clone_labels[[patient_from_labels[1]]]),
        stringsAsFactors = FALSE
      )
      names(patient_mapping) <- c(cell_col, patient_col)
    } else {
      cli::cli_warn(c(
        "{.arg clone_labels} has no {.field patient} column and no {.arg patient_mapping} was supplied.",
        "i" = "All cells will be treated as one patient."
      ))
    }
  }

  # --- Auto-enable parse_patient if patient_sep or patient_pos is provided ---
  if (!parse_patient) {
    # Check if user provided non-default patient_sep or patient_pos
    # Default: patient_sep="_", patient_pos=1
    # If user explicitly changed either, auto-enable parsing
    caller_args <- as.list(match.call())
    if (
      "patient_sep" %in%
        names(caller_args) ||
        "patient_pos" %in% names(caller_args)
    ) {
      parse_patient <- TRUE
      cli::cli_alert_info(
        "Auto-enabling parse_patient (patient_sep or patient_pos provided)"
      )
    }
  }

  # --- Convert patient_mapping to list format if data frame ---
  if (
    is.data.frame(patient_mapping) || data.table::is.data.table(patient_mapping)
  ) {
    # Copy, so that set-by-reference edits below never touch the caller's object.
    # as.data.table() returns the same object when the input is already a
    # data.table, so an explicit copy() is required for that case.
    patient_mapping <- data.table::copy(
      data.table::as.data.table(patient_mapping)
    )

    # Optional per-clone abundance column (count / cells / n_cells /
    # n_cells_per_clone / abundance): expand each row to that many rows so that
    # clone weights reflect TRUE cell numbers instead of equal 1/n per clone.
    # The Shiny app performs the same expansion when a mapping is uploaded
    # (inst/shiny/app/R/mod_data.R); doing it here as well makes the documented
    # behaviour available to scripted users who call prepare_data() directly.
    # The abundance column is consumed here and not passed downstream.
    count_col <- intersect(
      c("count", "cells", "n_cells", "n_cells_per_clone", "abundance"),
      tolower(colnames(patient_mapping))
    )
    if (length(count_col) > 0L) {
      cnt_col <- colnames(patient_mapping)[
        tolower(colnames(patient_mapping)) == count_col[1L]
      ][1L]
      cnt <- suppressWarnings(as.numeric(patient_mapping[[cnt_col]]))
      cnt[is.na(cnt) | cnt < 1L] <- 1L
      patient_mapping <- patient_mapping[rep(
        seq_len(nrow(patient_mapping)),
        cnt
      )]
      patient_mapping[, (cnt_col) := NULL]
      cli::cli_progress_step(
        "Expanded mapping by '{cnt_col}' to {nrow(patient_mapping)} rows (true clone abundance)"
      )
    }

    if (!cell_col %chin% colnames(patient_mapping)) {
      cli::cli_abort(
        "patient_mapping data frame must have a column named '{cell_col}'"
      )
    }

    # Parse patient ID from cell_col if requested
    if (parse_patient) {
      cli::cli_progress_step(
        "Parsing patient IDs from '{cell_col}' column using separator '{patient_sep}'"
      )
      cell_names <- patient_mapping[[cell_col]]
      parsed_parts <- strsplit(
        as.character(cell_names),
        patient_sep,
        fixed = TRUE
      )

      # Extract patient ID at specified position
      patient_ids_parsed <- sapply(parsed_parts, function(x) {
        if (length(x) >= patient_pos) {
          return(x[patient_pos])
        } else {
          return(NA_character_)
        }
      })

      # Check for parsing failures
      if (anyNA(patient_ids_parsed)) {
        cli::cli_warn(
          "{sum(is.na(patient_ids_parsed))} cells could not be parsed for patient ID"
        )
      }

      # Add parsed patient IDs to the data frame
      patient_mapping[, patient_id_parsed := patient_ids_parsed]
      sample_cell_names <- split(
        patient_mapping[[cell_col]],
        patient_mapping[["patient_id_parsed"]]
      )
      cli::cli_progress_step(
        "Parsed {length(unique(patient_ids_parsed))} unique patient IDs"
      )
    } else {
      # Use existing patient_col
      if (!patient_col %chin% colnames(patient_mapping)) {
        cli::cli_abort(
          c(
            "patient_mapping data frame must have a column named '{patient_col}'",
            "i" = "Use parse_patient = TRUE to extract patient ID from cell names"
          )
        )
      }
      sample_cell_names <- split(
        patient_mapping[[cell_col]],
        patient_mapping[[patient_col]]
      )
    }
    cli::cli_progress_step(
      "Converted metadata data frame to list format: {length(sample_cell_names)} patients"
    )
  } else {
    sample_cell_names <- patient_mapping
  }

  # --- Step 1: Define clones ------------------------------------------------
  if (!is.null(clone_labels)) {
    notify_stage("clones", 1L, 5L)
    cli::cli_progress_step("[1/5] Using supplied cell-to-clone labels...")

    if (is.null(data)) {
      label_source <- attr(clone_labels, "label_source")
      if (is.null(label_source)) {
        label_source <- "user_labels"
      }
    }

    if (!all(c("cell_id", "clone_id") %in% colnames(clone_labels))) {
      cli::cli_abort(c(
        "{.arg clone_labels} must have columns {.val cell_id} and {.val clone_id}.",
        "x" = "Got columns: {.val {colnames(clone_labels)}}."
      ))
    }

    cell_clone_map <- data.table::data.table(
      cell_id = as.character(clone_labels$cell_id),
      clone_id = as.character(clone_labels$clone_id)
    )

    # Carry precomputed 2D coordinates (e.g. from preprocess_sc()) through to
    # umap_coords, so the app can still draw an embedding in clone-label mode.
    if (all(c("dim_1", "dim_2") %in% colnames(clone_labels))) {
      cell_clone_map$dim_1 <- as.numeric(clone_labels$dim_1)
      cell_clone_map$dim_2 <- as.numeric(clone_labels$dim_2)
    }

    labelled <- cell_clone_map$cell_id %chin% colnames(expression_matrix)
    if (!all(labelled)) {
      cli::cli_warn(
        "{sum(!labelled)} labelled cell{?s} are absent from the expression matrix and will be dropped."
      )
      cell_clone_map <- cell_clone_map[labelled, ]
    }
    unlabelled <- setdiff(colnames(expression_matrix), cell_clone_map$cell_id)
    if (length(unlabelled) > 0L) {
      cli::cli_alert_info(
        "Dropping {length(unlabelled)} cell{?s} without a clone label (e.g. removed by QC)."
      )
    }

    cli::cli_progress_step(
      "Using {nrow(cell_clone_map)} labelled cells across {length(unique(cell_clone_map$clone_id))} clones"
    )
  } else if (isTRUE(skip_clustering)) {
    notify_stage("clones", 1L, 5L)
    cli::cli_progress_step(
      "[1/5] Skipping Seurat clustering -- each column is treated as one clone..."
    )
    cn <- colnames(expression_matrix)
    if (is.null(cn)) {
      cli::cli_abort(
        "expression_matrix must have column names when skip_clustering = TRUE (each column is expected to be one clone)."
      )
    }
    # Clone-level input with a cell-count-resolved mapping: if the patient
    # mapping repeats clone ids (one row per real cell, e.g. "Kydar01_c1"
    # appearing 500 times), build the cell_clone_map from the mapping rows so
    # clone abundances reflect TRUE cell counts (matching the paper's figures)
    # instead of equal weights per clone. If the mapping is 1:1 with the matrix
    # columns (legacy equal-weight case) the behavior is unchanged.
    sample_cell_ids <- if (is.list(sample_cell_names)) {
      unlist(sample_cell_names, use.names = FALSE)
    } else if (is.character(sample_cell_names)) {
      sample_cell_names
    } else {
      NULL
    }

    if (
      !is.null(sample_cell_ids) &&
        anyDuplicated(sample_cell_ids) > 0L &&
        all(cn %in% sample_cell_ids)
    ) {
      cell_clone_map <- data.table::data.table(
        cell_id = sample_cell_ids,
        clone_id = sample_cell_ids
      )
      cell_id <- NULL # ease check NOTE
      cell_clone_map <- cell_clone_map[cell_id %in% cn]
      cli::cli_progress_step(
        "Using cell-count-resolved mapping: {nrow(cell_clone_map)} cells across {length(unique(cell_clone_map$cell_id))} clones"
      )
    } else {
      # Keep the full column name as the unique clone identity (e.g. "Kydar01_c1").
      # This is universal: clone names are whatever the user's data uses. Patient
      # separation is guaranteed downstream via "patient@@clone" keys.
      cell_clone_map <- data.table::data.table(
        cell_id = cn,
        clone_id = cn
      )
      cli::cli_progress_step("Using {length(cn)} provided clone columns")
    }
  } else {
    notify_stage("seurat", 1L, 5L)
    cli::cli_progress_step("[1/5] Clustering cells via Seurat...")

    # Auto-adjust min_features and nfeatures if gene count is small
    n_genes <- nrow(expression_matrix)
    if (n_genes < seurat_min_features) {
      seurat_min_features <- max(1L, floor(n_genes * 0.1))
      cli::cli_progress_step(
        "Adjusting min_features to {seurat_min_features} (only {n_genes} genes in expression matrix)"
      )
    }
    if (n_genes < seurat_nfeatures) {
      seurat_nfeatures <- max(10L, n_genes)
      cli::cli_progress_step(
        "Adjusting nfeatures to {seurat_nfeatures} (only {n_genes} genes available)"
      )
    }

    cell_clone_map <- annotate_clones(
      method = method,
      expression_matrix = expression_matrix,
      min_cells = seurat_min_cells,
      min_features = seurat_min_features,
      nfeatures = seurat_nfeatures,
      dims = seurat_dims,
      resolution = seurat_resolution,
      seed = seurat_seed,
      # Seurat step callbacks (normalizing/pca/clustering/...) are forwarded
      # through the same channel so the UI can show live stage text.
      progress_cb = function(step) notify_stage(paste0("seurat-", step), 1L, 5L)
    )
    cli::cli_progress_step(
      "Found {length(unique(cell_clone_map$clone_id))} clones across {nrow(cell_clone_map)} cells"
    )
  }

  # --- Step 2: Build patient IDs ---
  notify_stage("mapping", 2L, 5L)
  cli::cli_progress_step("[2/5] Mapping cells to patients...")
  if (is.null(sample_cell_names)) {
    patient_ids <- rep("patient1", ncol(expression_matrix))
    names(patient_ids) <- colnames(expression_matrix)
  } else {
    patient_ids <- rep(NA_character_, ncol(expression_matrix))
    names(patient_ids) <- colnames(expression_matrix)
    for (pat_name in names(sample_cell_names)) {
      pat_cells <- sample_cell_names[[pat_name]]
      matched <- intersect(pat_cells, colnames(expression_matrix))
      patient_ids[matched] <- pat_name
    }
    # Remove cells with no patient assignment
    unassigned <- is.na(patient_ids)
    if (any(unassigned)) {
      cli::cli_warn(
        "{sum(unassigned)} cells could not be assigned to any patient. Removing them."
      )
      keep_cells <- names(patient_ids)[!unassigned]
      expression_matrix <- expression_matrix[, keep_cells]
      patient_ids <- patient_ids[!unassigned]
      cell_clone_map <- cell_clone_map[cell_id %chin% keep_cells, ]
    }
  }
  cli::cli_progress_step("Found {length(unique(patient_ids))} patients")

  # Add patient column to cell_clone_map; patient_ids is named by cell id, so the
  # lookup needs no match().
  cell_clone_map[, patient := unname(patient_ids[cell_id])]

  # In clone-level mode, refine the clone label to the per-patient category
  # (e.g. "Kydar01_c1" -> "c1") using the KNOWN patient id, so clone keys
  # become "Kydar01@@c1" instead of "Kydar01@@Kydar01_c1". This only strips
  # when the cell name actually starts with its patient id; any other naming
  # scheme keeps the full name, so the behavior is universal.
  if (isTRUE(skip_clustering)) {
    clone_id <- patient <- cell_id <- NULL # ease checking NOTE
    cell_clone_map[,
      clone_id := {
        has_prefix <- !is.na(patient) & startsWith(cell_id, patient)
        stripped <- sub(
          "^[\\._-]",
          "",
          substring(cell_id, nchar(patient) + 1L)
        )
        data.table::fifelse(has_prefix & nzchar(stripped), stripped, cell_id)
      }
    ]
  }

  # --- Step 3: Compute clone-level mean expression ---
  cli::cli_progress_step("[3/5] Computing clone-level mean expression...")
  clone_expr_list <- clone_mean_expression(
    expression_matrix = expression_matrix,
    cell_clone_map = cell_clone_map[, c("cell_id", "clone_id")],
    patient_ids = patient_ids
  )

  if (length(clone_expr_list) == 0L) {
    cli::cli_abort(
      "No clone expression matrices were produced. Check cell-clone mapping."
    )
  }

  # Merge all patients into one matrix
  clone_expr_merged <- do.call(cbind, clone_expr_list)
  cli::cli_progress_step(
    "Produced {ncol(clone_expr_merged)} clone columns across {length(clone_expr_list)} patients"
  )

  # --- Step 4: Filter genes and rank normalize ---
  cli::cli_progress_step("[4/5] Rank-normalizing clone expression...")
  if (!is.null(genes_to_use)) {
    available_genes <- intersect(genes_to_use, rownames(clone_expr_merged))
    if (length(available_genes) == 0L) {
      rlang::abort(
        "None of the provided genes_to_use were found in the expression matrix."
      )
    }
    if (length(available_genes) < length(genes_to_use)) {
      cli::cli_warn(
        "Only {length(available_genes)} / {length(genes_to_use)} genes_to_use found in expression matrix."
      )
    }
    clone_expr_merged <- clone_expr_merged[available_genes, , drop = FALSE]
  }

  clone_expression_rnorm <- rank_normalization_mat(clone_expr_merged)

  # --- Step 5: Build clone counts table ---
  notify_stage("clone-counts", 5L, 5L)
  cli::cli_progress_step("[5/5] Building clone abundance table...")
  clone_counts <- build_clone_counts(
    cell_clone_map = cell_clone_map,
    patient_ids = cell_clone_map$patient
  )

  # Build clone_viability_df template from cell_clone_map
  # Row order matches clone_expression_rnorm column order
  clone_viability_template <- parse_clone_keys(colnames(clone_expression_rnorm))

  cli::cli_progress_done()
  cli::cli_h2("Preparation complete")
  cli::cli_ul("Genes: {.val {nrow(clone_expression_rnorm)}}")
  cli::cli_ul("Clones: {.val {ncol(clone_expression_rnorm)}}")
  cli::cli_ul("Patients: {.val {nrow(clone_counts)}}")

  return(list(
    clone_expression_rnorm = clone_expression_rnorm,
    clone_counts = clone_counts,
    cell_clone_map = cell_clone_map,
    clone_viability_template = clone_viability_template,
    umap_coords = if (all(c("dim_1", "dim_2") %in% names(cell_clone_map))) {
      cell_clone_map[, c("cell_id", "dim_1", "dim_2")]
    } else {
      NULL
    },
    reduction_method = if (isTRUE(skip_clustering)) "none" else method,
    label_source = label_source
  ))
}
