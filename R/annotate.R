#' PERCEPTION Patient Data Annotation Functions
#'
#' Functions for annotating single-cell data with clone and patient information,
#' and preparing patient data for the prediction pipeline.


#' Annotate cells with clone IDs via Seurat clustering
#'
#' Performs Seurat clustering on a single-cell expression matrix and returns
#' a mapping of each cell to its cluster (clone) ID. This matches the original
#' PERCEPTION pipeline where Seurat clusters define transcriptional subclones.
#'
#' @param expression_matrix Matrix. Gene expression matrix with genes as rows
#'        and cells as columns. Raw counts or normalized values are both accepted.
#' @param min_cells Integer. Minimum cells per feature. Default = 3.
#' @param min_features Integer. Minimum features per cell. Default = 200.
#' @param nfeatures Integer. Number of variable features. Default = 2000.
#' @param dims Integer. Number of PCA dimensions for clustering. Default = 10.
#' @param resolution Numeric. Clustering resolution. Default = 0.8.
#' @param seed Integer. Random seed for reproducibility. Default = 42.
#'
#' @return A data frame with columns: \code{cell_id} and \code{clone_id}.
#'
#' @examples
#' \dontrun{
#'   cell_clone_map <- annotate_clones(patient_expression)
#' }
#'
#' @export
annotate_clones <- function(expression_matrix,
                             min_cells = 3,
                             min_features = 200,
                             nfeatures = 2000,
                             dims = 10,
                             resolution = 0.8,
                             seed = 42) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for clone annotation. ",
         "Install with: install.packages('Seurat')")
  }

  set.seed(seed)

  so <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                   project = "PERCEPTION",
                                   min.cells = min_cells,
                                   min.features = min_features)
  so <- Seurat::NormalizeData(so, normalization.method = "LogNormalize",
                              scale.factor = 10000)
  so <- Seurat::FindVariableFeatures(so, selection.method = "vst",
                                     nfeatures = nfeatures)
  so <- Seurat::ScaleData(so)

  # Adjust PCA dimensions based on available cells/features
  # PCA requires npcs < min(nrow, ncol) of the scaled matrix
  max_pcs <- min(ncol(so), length(Seurat::VariableFeatures(object = so))) - 1
  actual_dims <- min(dims, max_pcs)
  if (actual_dims < 1) actual_dims <- 1

  if (actual_dims < dims) {
    message("  Adjusting PCA dims from ", dims, " to ", actual_dims,
            " due to limited cells/features")
  }

  so <- Seurat::RunPCA(so, features = Seurat::VariableFeatures(object = so),
                       npcs = actual_dims)
  so <- Seurat::FindNeighbors(so, dims = 1:actual_dims)
  so <- Seurat::FindClusters(so, resolution = resolution)

  cluster_ids <- Seurat::Idents(so)

  data.frame(
    cell_id = names(cluster_ids),
    clone_id = as.character(cluster_ids),
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


#' Prepare patient data for PERCEPTION prediction
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
#' @param seurat_seed Integer. Random seed. Default = 42.
#'
#' @return A named list with:
#' \describe{
#'   \item{clone_expression_rnorm}{Matrix. Rank-normalized clone-level expression
#'         (genes as rows, patient_clone as columns). Ready for \code{predict_drugs()}.}
#'   \item{clone_counts}{Data frame. Clone abundance per patient. Ready for
#'         \code{predict_patients()}.}
#'   \item{cell_clone_map}{Data frame. Cell-to-clone mapping with columns
#'         cell_id, clone_id, patient.}
#'   \item{clone_killing_df_template}{Data frame. Template with patient and clone_id
#'         columns, ready to merge with \code{predict_drugs()} output.}
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
#'
#' @export
prepare_data <- function(expression_matrix,
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
                          seurat_seed = 42) {

  message("=== PERCEPTION Patient Data Preparation ===")

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

  # --- Step 1: Seurat clustering to define subclones ---
  message("[1/5] Clustering cells via Seurat...")
  cell_clone_map <- annotate_clones(
    expression_matrix = expression_matrix,
    resolution = seurat_resolution,
    dims = seurat_dims,
    nfeatures = seurat_nfeatures,
    seed = seurat_seed
  )
  message("  Found ", length(unique(cell_clone_map$clone_id)), " clones across ",
          nrow(cell_clone_map), " cells.")

  # --- Step 2: Build patient IDs ---
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
  message("[5/5] Building clone abundance table...")
  clone_counts <- build_clone_counts(
    cell_clone_map = cell_clone_map,
    patient_ids = cell_clone_map$patient
  )

  # Build clone_killing_df template from cell_clone_map
  # Row order matches clone_expression_rnorm column order
  clone_col_names <- colnames(clone_expression_rnorm)
  # Parse patient_clone format using @@ separator
  template_patients <- sapply(strsplit(clone_col_names, "@@"), `[`, 1)
  template_clones <- sapply(strsplit(clone_col_names, "@@"), `[`, 2)

  clone_killing_template <- data.frame(
    patient = template_patients,
    clone_id = template_clones,
    stringsAsFactors = FALSE
  )

  message("\n=== Preparation complete ===")
  message("  Genes: ", nrow(clone_expression_rnorm))
  message("  Clones: ", ncol(clone_expression_rnorm))
  message("  Patients: ", nrow(clone_counts))

  return(list(
    clone_expression_rnorm = clone_expression_rnorm,
    clone_counts = clone_counts,
    cell_clone_map = cell_clone_map,
    clone_killing_template = clone_killing_template
  ))
}
