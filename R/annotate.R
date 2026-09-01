# encoding: UTF-8

#' PERCEPTIONx Patient Data Annotation Functions
#'
#' Functions for annotating single-cell data with clone and patient information,
#' and preparing patient data for the prediction pipeline.
#'
#' @name annotate_perception
#' @keywords internal
NULL


# Internal: shared Seurat clustering + 2D embedding pipeline.
#
# Runs NormalizeData -> FindVariableFeatures -> ScaleData -> PCA -> neighbors ->
# clusters -> UMAP/t-SNE, with PCA dimensions capped to what the (possibly
# small) dataset can support. Used by both annotate_clones() and
# plot_seurat_clustering() so they stay in sync and share the same safeguards.
run_seurat_pipeline <- function(expression_matrix, method,
                                min_cells, min_features, nfeatures,
                                dims, resolution, seed, progress_cb = NULL,
                                cluster_algorithm = 1,
                                variable_selection = c("vst", "dispersion", "mvp")) {
  notify <- function(phase) {
    if (is.function(progress_cb)) progress_cb(phase)
  }
  variable_selection <- match.arg(variable_selection)
  set.seed(seed)

  notify("normalizing")
  so <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                   project = "PERCEPTIONx",
                                   min.cells = min_cells,
                                   min.features = min_features)
  so <- Seurat::NormalizeData(so, normalization.method = "LogNormalize",
                              scale.factor = 10000)
  notify("variable-features")
  # vst is the most accurate but by far the slowest variable-feature method
  # on large matrices; dispersion/mvp trade a little feature quality for a
  # big speedup. The Seurat pipeline below works with any of them.
  so <- Seurat::FindVariableFeatures(so, selection.method = variable_selection,
                                     nfeatures = nfeatures)
  notify("scaling")
  so <- Seurat::ScaleData(so, features = Seurat::VariableFeatures(object = so))

  # PCA requires npcs < min(nrow, ncol) of the scaled matrix -- cap dims to the
  # number actually available so small datasets do not error out.
  max_pcs <- min(ncol(so), length(Seurat::VariableFeatures(object = so))) - 1
  actual_dims <- min(dims, max_pcs)
  if (actual_dims < 1) actual_dims <- 1

  if (actual_dims < dims) {
    message("  Adjusting PCA dims from ", dims, " to ", actual_dims,
            " due to limited cells/features")
  }

  notify("pca")
  so <- Seurat::RunPCA(so, features = Seurat::VariableFeatures(object = so),
                       npcs = actual_dims)
  notify("neighbors")
  so <- Seurat::FindNeighbors(so, dims = seq_len(actual_dims))
  # verbose = FALSE: FindClusters prints "Modularity Optimizer" progress to
  # stdout, which pollutes worker logs / can leave stray files behind.
  # algorithm = 1 is the default Louvain (no extra dependency). Leiden (4)
  # is slightly faster on VERY large graphs but pulls in the 'leidenbase'
  # package — not worth the dependency at PERCEPTIONx's typical cell counts.
  notify("clustering")
  so <- Seurat::FindClusters(so, resolution = resolution,
                             algorithm = cluster_algorithm, verbose = FALSE)

  notify(if (method == "umap") "umap" else "tsne")
  if (method == "umap") {
    # RunUMAP's uwot call is single-threaded and uses heavy defaults
    # (n_neighbors = 30, n_epochs = 500, n_trees = 50) -- on small datasets
    # the Annoy index build + calibration alone cost ~1-5 s of fixed
    # overhead (measured: 5.0 s on 1500 cells, ~1.2 s on 400). The UMAP here
    # is ONLY the 2D layout for the spatial plots; clone detection comes from
    # FindClusters on the SNN graph, so we call uwot directly with reduced,
    # still-valid settings (15/200/10 -> ~2.5x faster, measured). The result
    # is attached as a normal "umap" reduction, so Embeddings(so, "umap")
    # and everything downstream work unchanged.
    pca_emb <- Seurat::Embeddings(so, "pca")[, seq_len(actual_dims), drop = FALSE]
    n_threads <- max(1L, min(4L, as.integer(parallel::detectCores()), na.rm = TRUE))
    # uwot requires n_neighbors < number of points; cap it for tiny datasets
    # so a small test matrix (<= 15 cells) does not crash with an Annoy
    # index error.
    n_neighbors <- min(15, max(2L, nrow(pca_emb) - 1L))
    coords <- uwot::umap(
      X = pca_emb, n_neighbors = n_neighbors, n_components = 2, n_epochs = 200,
      n_trees = 10, metric = "cosine", learning_rate = 1, min_dist = 0.3,
      spread = 1, set_op_mix_ratio = 1, local_connectivity = 1,
      repulsion_strength = 1, negative_sample_rate = 5,
      n_threads = n_threads, verbose = FALSE)
    rownames(coords) <- colnames(so)
    colnames(coords) <- c("UMAP_1", "UMAP_2")
    so[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = coords, key = "UMAP_", assay = Seurat::DefaultAssay(so))
  } else {
    so <- Seurat::RunTSNE(so, dims = seq_len(actual_dims), check_duplicates = FALSE)
  }

  so
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
annotate_clones <- function(method = c("umap", "tsne"),
                             expression_matrix,
                             min_cells = 3,
                             min_features = 200,
                             nfeatures = 2000,
                             dims = 10,
                             resolution = 0.8,
                             seed = 42,
                             progress_cb = NULL,
                             cluster_algorithm = 1,
                             variable_selection = c("vst", "dispersion", "mvp")) {

  method <- match.arg(method)
  variable_selection <- match.arg(variable_selection)

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for clone annotation. ",
         "Install with: install.packages('Seurat')")
  }

  so <- run_seurat_pipeline(expression_matrix, method,
                            min_cells, min_features, nfeatures,
                            dims, resolution, seed, progress_cb,
                            cluster_algorithm, variable_selection)

  cluster_ids <- Seurat::Idents(so)
  emb <- Seurat::Embeddings(so, reduction = method)

  data.frame(
    cell_id = names(cluster_ids),
    clone_id = as.character(cluster_ids),
    dim_1 = emb[, 1],
    dim_2 = emb[, 2],
    stringsAsFactors = FALSE
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
    stop("cell_clone_map must have columns 'cell_id' and 'clone_id'.")
  }

  if (length(patient_ids) != nrow(cell_clone_map)) {
    stop("Length of patient_ids must match number of rows in cell_clone_map.")
  }

  map_df <- data.frame(
    patient = patient_ids,
    clone_id = cell_clone_map$clone_id,
    stringsAsFactors = FALSE
  )

  all_clones <- unique(cell_clone_map$clone_id)
  all_patients <- unique(patient_ids)

  result <- data.frame(patients = all_patients, stringsAsFactors = FALSE)
  for (cl in all_clones) {
    result[[cl]] <- sapply(all_patients, function(pat) {
      sum(map_df$patient == pat & map_df$clone_id == cl)
    })
  }

  result
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
#'   If NULL, all cells are assigned to a single patient "patient1".
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
#'
#' @export
prepare_data <- function(method = c("umap", "tsne"),
                          expression_matrix,
                          patient_mapping = NULL,
                          cell_col = "cell_id",
                          patient_col = "patient_id",
                          parse_patient = FALSE,
                          patient_sep = "_",
                          patient_pos = 1,
                          genes_to_use = NULL,
                          seurat_resolution = 0.8,
                          seurat_dims = 10,
                          seurat_nfeatures = 2000,
                          seurat_min_cells = 3,
                          seurat_min_features = 200,
                          seurat_seed = 42,
                          skip_clustering = FALSE,
                          progress_cb = NULL) {

  method <- match.arg(method)
  notify_stage <- function(phase, i, n) {
    if (is.function(progress_cb)) progress_cb(phase, i, n, "")
  }
  message("=== PERCEPTIONx Patient Data Preparation ===")
  message("  Reduction method: ", toupper(method))

  # --- Auto-enable parse_patient if patient_sep or patient_pos is provided ---
  if (!parse_patient) {
    # Check if user provided non-default patient_sep or patient_pos
    # Default: patient_sep="_", patient_pos=1
    # If user explicitly changed either, auto-enable parsing
    caller_args <- as.list(match.call())
    if ("patient_sep" %in% names(caller_args) || "patient_pos" %in% names(caller_args)) {
      parse_patient <- TRUE
      message("  Auto-enabling parse_patient (patient_sep or patient_pos provided)")
    }
  }

  # --- Convert patient_mapping to list format if data frame ---
  if (is.data.frame(patient_mapping)) {
    if (!cell_col %in% colnames(patient_mapping)) {
      stop("patient_mapping data frame must have a column named '", cell_col, "'")
    }

    # Parse patient ID from cell_col if requested
    if (parse_patient) {
      message("  Parsing patient IDs from '", cell_col, "' column using separator '", patient_sep, "'")
      cell_names <- patient_mapping[[cell_col]]
      parsed_parts <- strsplit(as.character(cell_names), patient_sep, fixed = TRUE)

      # Extract patient ID at specified position
      patient_ids_parsed <- sapply(parsed_parts, function(x) {
        if (length(x) >= patient_pos) {
          return(x[patient_pos])
        } else {
          return(NA_character_)
        }
      })

      # Check for parsing failures
      if (any(is.na(patient_ids_parsed))) {
        warning(sum(is.na(patient_ids_parsed)), " cells could not be parsed for patient ID")
      }

      # Add parsed patient IDs to the data frame
      patient_mapping$patient_id_parsed <- patient_ids_parsed
      sample_cell_names <- split(patient_mapping[[cell_col]], patient_mapping$patient_id_parsed)
      message("  Parsed ", length(unique(patient_ids_parsed)), " unique patient IDs")
    } else {
      # Use existing patient_col
      if (!patient_col %in% colnames(patient_mapping)) {
        stop("patient_mapping data frame must have a column named '", patient_col,
             "'\n  Or use parse_patient=TRUE to extract patient ID from cell names")
      }
      sample_cell_names <- split(patient_mapping[[cell_col]], patient_mapping[[patient_col]])
    }
    message("  Converted metadata data frame to list format: ", length(sample_cell_names), " patients")
  } else {
    sample_cell_names <- patient_mapping
  }

  # --- Step 1: Define clones (Seurat clustering OR clone-level input) ---
  if (isTRUE(skip_clustering)) {
    notify_stage("clones", 1L, 5L)
    message("[1/5] Skipping Seurat clustering -- each column is treated as one clone...")
    cn <- colnames(expression_matrix)
    if (is.null(cn)) {
      stop("expression_matrix must have column names when skip_clustering = TRUE ",
           "(each column is expected to be one clone).")
    }
    # Clone-level input with a cell-count-resolved mapping: if the patient
    # mapping repeats clone ids (one row per real cell, e.g. "Kydar01_c1"
    # appearing 500 times), build the cell_clone_map from the mapping rows so
    # clone abundances reflect TRUE cell counts (matching the paper's figures)
    # instead of equal weights per clone. If the mapping is 1:1 with the matrix
    # columns (legacy equal-weight case) the behavior is unchanged.
    sample_cell_ids <- NULL
    if (is.list(sample_cell_names)) {
      sample_cell_ids <- unlist(sample_cell_names, use.names = FALSE)
    } else if (is.character(sample_cell_names)) {
      sample_cell_ids <- sample_cell_names
    }
    if (!is.null(sample_cell_ids) && anyDuplicated(sample_cell_ids) > 0 &&
        all(cn %in% sample_cell_ids)) {
      cell_clone_map <- data.frame(cell_id = sample_cell_ids,
                                   clone_id = sample_cell_ids,
                                   stringsAsFactors = FALSE)
      cell_clone_map <- cell_clone_map[cell_clone_map$cell_id %in% cn, , drop = FALSE]
      message("  Using cell-count-resolved mapping: ", nrow(cell_clone_map),
              " cells across ", length(unique(cell_clone_map$cell_id)), " clones.")
    } else {
      # Keep the full column name as the unique clone identity (e.g. "Kydar01_c1").
      # This is universal: clone names are whatever the user's data uses. Patient
      # separation is guaranteed downstream via "patient@@clone" keys.
      cell_clone_map <- data.frame(cell_id = cn, clone_id = cn,
                                   stringsAsFactors = FALSE)
      message("  Using ", length(cn), " provided clone columns.")
    }
  } else {
    notify_stage("seurat", 1L, 5L)
    message("[1/5] Clustering cells via Seurat...")

    # Auto-adjust min_features and nfeatures if gene count is small
    n_genes <- nrow(expression_matrix)
    if (n_genes < seurat_min_features) {
      seurat_min_features <- max(1, floor(n_genes * 0.1))
      message("  Adjusting min_features to ", seurat_min_features,
              " (only ", n_genes, " genes in expression matrix)")
    }
    if (n_genes < seurat_nfeatures) {
      seurat_nfeatures <- max(10, n_genes)
      message("  Adjusting nfeatures to ", seurat_nfeatures,
              " (only ", n_genes, " genes available)")
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
    message("  Found ", length(unique(cell_clone_map$clone_id)), " clones across ",
            nrow(cell_clone_map), " cells.")
  }

  # --- Step 2: Build patient IDs ---
  notify_stage("mapping", 2L, 5L)
  message("[2/5] Mapping cells to patients...")
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
      warning(sum(unassigned), " cells could not be assigned to any patient. Removing them.")
      keep_cells <- names(patient_ids)[!unassigned]
      expression_matrix <- expression_matrix[, keep_cells, drop = FALSE]
      patient_ids <- patient_ids[!unassigned]
      cell_clone_map <- cell_clone_map[cell_clone_map$cell_id %in% keep_cells, ]
    }
  }
  message("  Found ", length(unique(patient_ids)), " patients.")

  # Add patient column to cell_clone_map
  cell_clone_map$patient <- patient_ids[match(cell_clone_map$cell_id,
                                                names(patient_ids))]

  # In clone-level mode, refine the clone label to the per-patient category
  # (e.g. "Kydar01_c1" -> "c1") using the KNOWN patient id, so clone keys
  # become "Kydar01@@c1" instead of "Kydar01@@Kydar01_c1". This only strips
  # when the cell name actually starts with its patient id; any other naming
  # scheme keeps the full name, so the behavior is universal.
  if (isTRUE(skip_clustering)) {
    cell_clone_map$clone_id <- vapply(seq_len(nrow(cell_clone_map)), function(i) {
      cell <- cell_clone_map$cell_id[i]
      pat  <- cell_clone_map$patient[i]
      if (is.na(pat) || !startsWith(cell, pat)) return(cell)
      rest <- sub("^[\\._-]", "", substring(cell, nchar(pat) + 1))
      if (nchar(rest) == 0) cell else rest
    }, character(1))
  }

  # --- Step 3: Compute clone-level mean expression ---
  message("[3/5] Computing clone-level mean expression...")
  clone_expr_list <- clone_mean_expression(
    expression_matrix = expression_matrix,
    cell_clone_map = cell_clone_map[, c("cell_id", "clone_id")],
    patient_ids = patient_ids
  )

  if (length(clone_expr_list) == 0) {
    stop("No clone expression matrices were produced. Check cell-clone mapping.")
  }

  # Merge all patients into one matrix
  clone_expr_merged <- do.call(cbind, clone_expr_list)
  message("  Produced ", ncol(clone_expr_merged), " clone columns across ",
          length(clone_expr_list), " patients.")

  # --- Step 4: Filter genes and rank normalize ---
  message("[4/5] Rank-normalizing clone expression...")
  if (!is.null(genes_to_use)) {
    available_genes <- intersect(genes_to_use, rownames(clone_expr_merged))
    if (length(available_genes) == 0) {
      stop("None of the provided genes_to_use were found in the expression matrix.")
    }
    if (length(available_genes) < length(genes_to_use)) {
      warning("Only ", length(available_genes), " / ", length(genes_to_use),
              " genes_to_use found in expression matrix.")
    }
    clone_expr_merged <- clone_expr_merged[available_genes, , drop = FALSE]
  }

  clone_expression_rnorm <- rank_normalization_mat(clone_expr_merged)

  # --- Step 5: Build clone counts table ---
  notify_stage("clone-counts", 5L, 5L)
  message("[5/5] Building clone abundance table...")
  clone_counts <- build_clone_counts(
    cell_clone_map = cell_clone_map,
    patient_ids = cell_clone_map$patient
  )

  # Build clone_viability_df template from cell_clone_map
  # Row order matches clone_expression_rnorm column order
  clone_viability_template <- parse_clone_keys(colnames(clone_expression_rnorm))

  message("\n=== Preparation complete ===")
  message("  Genes: ", nrow(clone_expression_rnorm))
  message("  Clones: ", ncol(clone_expression_rnorm))
  message("  Patients: ", nrow(clone_counts))

  return(list(
    clone_expression_rnorm = clone_expression_rnorm,
    clone_counts = clone_counts,
    cell_clone_map = cell_clone_map,
    clone_viability_template = clone_viability_template,
    umap_coords = if (isTRUE(skip_clustering)) NULL else
      cell_clone_map[, c("cell_id", "dim_1", "dim_2")],
    reduction_method = if (isTRUE(skip_clustering)) "none" else method
  ))
}
