#' PERCEPTION Statistical Utilities
#'
#' @name stats_perception
#' @keywords internal
#' @importFrom stats na.omit phyper p.adjust complete.cases cor pt pnorm
NULL

#' PERCEPTION Statistical Utilities
#'
#' @name stats_perception
#' @keywords internal
#' @importFrom stats na.omit phyper p.adjust complete.cases cor pt pnorm
NULL

#' Hypergeometric test for gene list overlap
#'
#' Tests whether the overlap between two gene lists is statistically significant.
#' It answers the question: given a background set of genes (global), a reference
#' gene set (base_list), and a test gene set (test_list), is the overlap between
#' test_list and base_list more than expected by chance?
#'
#' @param test_list Character vector. The test gene list.
#' @param base_list Character vector. The reference gene list.
#' @param global Character vector. The background gene set.
#' @param lower.tail Logical. If FALSE (default), calculates the probability
#'        of observing \emph{greater than or equal to} the observed overlap
#'        (enrichment). If TRUE, calculates the probability of observing
#'        \emph{less than or equal to} the observed overlap (depletion).
#'
#' @return Numeric p-value.
#' @export
hypergeometric_test_for_twolists <- function(test_list, base_list, global, lower.tail = FALSE) {
  base_in_global <- global[na.omit(match(base_list, global))]
  overlap <- test_list[!is.na(match(test_list, base_in_global))]
  phyper(
    q = length(overlap) - 1,
    m = length(base_in_global),
    n = length(global) - length(base_in_global),
    k = length(test_list),
    lower.tail = lower.tail
  )
}



#' FDR correction for multiple testing
#'
#' Performs False Discovery Rate (FDR) correction using the Benjamini-Hochberg method.
#' This is a wrapper around \code{\link{p.adjust}} with method = "fdr".
#'
#' @param test_list A numeric vector of p-values to be corrected.
#' @return A numeric vector of FDR-adjusted p-values (q-values).
#' @export
fdrcorr <- function(test_list) {
  p.adjust(test_list, method = "fdr")
}



#' Rank-normalize each column of a matrix
#' Rank-normalize each column of a matrix
#'
#' Converts each column to ranks divided by column length, producing values in
#' the range 0 to 1 (exclusive of 0, inclusive of 1). This is a **crucial preprocessing step** when using your
#' own expression data with PERCEPTION.
#'
#' @section Why rank normalization?:
#' PERCEPTION models are trained on DepMap expression data that has been
#' rank-normalized. The model coefficients capture the relationship between
#' **relative gene expression ranks** and drug response, not absolute expression
#' values. Therefore, any new expression data fed into the model must undergo
#' the same rank normalization to ensure compatibility.
#'
#' @section How it works:
#' For each column (cell/sample), every gene's expression value is replaced by
#' its rank within that column, divided by the total number of genes:
#' \deqn{x_{ij}^{norm} = \frac{\mathrm{rank}(x_{ij})}{n}}
#' where \eqn{n} is the number of rows (genes) and ties are resolved by
#' averaging. This transforms each column into a uniform distribution over
#' \eqn{(0,1]}, making the data robust to batch effects, library size
#' differences, and outliers.
#' Converts each column to ranks divided by column length, producing values in
#' the range 0 to 1 (exclusive of 0, inclusive of 1). This is a **crucial preprocessing step** when using your
#' own expression data with PERCEPTION.
#'
#' @section Why rank normalization?:
#' PERCEPTION models are trained on DepMap expression data that has been
#' rank-normalized. The model coefficients capture the relationship between
#' **relative gene expression ranks** and drug response, not absolute expression
#' values. Therefore, any new expression data fed into the model must undergo
#' the same rank normalization to ensure compatibility.
#'
#' @section How it works:
#' For each column (cell/sample), every gene's expression value is replaced by
#' its rank within that column, divided by the total number of genes:
#' \deqn{x_{ij}^{norm} = \frac{\mathrm{rank}(x_{ij})}{n}}
#' where \eqn{n} is the number of rows (genes) and ties are resolved by
#' averaging. This transforms each column into a uniform distribution over
#' \eqn{(0,1]}, making the data robust to batch effects, library size
#' differences, and outliers.
#'
#' @section Important:
#' If you provide your own expression data that has **not** been
#' rank-normalized, predictions will be unreliable. Always run
#' \code{rank_normalization_mat()} on your raw expression matrix before passing
#' it to \code{\link{predict_drugs}()}.
#'
#' @param mat A numeric matrix. Rows = genes, columns = cells/samples.
#' @return A rank-normalized matrix of the same dimensions, with values
#'   between 0 and 1.
#' @section Important:
#' If you provide your own expression data that has **not** been
#' rank-normalized, predictions will be unreliable. Always run
#' \code{rank_normalization_mat()} on your raw expression matrix before passing
#' it to \code{\link{predict_drugs}()}.
#'
#' @param mat A numeric matrix. Rows = genes, columns = cells/samples.
#' @return A rank-normalized matrix of the same dimensions, with values
#'   between 0 and 1.
#' @export
#' @examples
#' # Raw expression matrix (genes x cells)
#' raw_expr <- matrix(c(10, 200, 50, 30, 300, 80), nrow = 3, ncol = 2)
#' rownames(raw_expr) <- c("GENE_A", "GENE_B", "GENE_C")
#'
#' # MUST normalize before prediction
#' norm_expr <- rank_normalization_mat(raw_expr)
#'
#' # Then feed into predict_drugs()
#' # predict_drugs(models, norm_expr)
#' @examples
#' # Raw expression matrix (genes x cells)
#' raw_expr <- matrix(c(10, 200, 50, 30, 300, 80), nrow = 3, ncol = 2)
#' rownames(raw_expr) <- c("GENE_A", "GENE_B", "GENE_C")
#'
#' # MUST normalize before prediction
#' norm_expr <- rank_normalization_mat(raw_expr)
#'
#' # Then feed into predict_drugs()
#' # predict_drugs(models, norm_expr)
rank_normalization_mat <- function(mat){
  apply(
    mat,
    2,
    function(x) rank(x, ties.method = "average")/length(x))
}



#' Change range to 0-1
#' Change range to 0-1
#'
#' Scales a numeric vector to the 0-1 range using the 5th and 95th percentiles
#' as thresholds to handle outliers.
#' Scales a numeric vector to the 0-1 range using the 5th and 95th percentiles
#' as thresholds to handle outliers.
#'
#' @param x A numeric vector.
#' @return A numeric vector scaled between 0 and 1.
#' @param x A numeric vector.
#' @return A numeric vector scaled between 0 and 1.
#' @export
range01 <- function(x){
  # Chossing 95% and 5% percentile as thresholds for outliers
  substitute_of_Min <- topXPercentValue(vec=x,
                                     X_percentile=5)
  substitute_of_Max <- topXPercentValue(vec=x,
                                     X_percentile=95)
  x_scaled <- (x-substitute_of_Min)/(substitute_of_Max-substitute_of_Min)
  x_scaled[x_scaled<0] = 0
  x_scaled[x_scaled>1] = 1
range01 <- function(x){
  # Chossing 95% and 5% percentile as thresholds for outliers
  substitute_of_Min <- topXPercentValue(vec=x,
                                     X_percentile=5)
  substitute_of_Max <- topXPercentValue(vec=x,
                                     X_percentile=95)
  x_scaled <- (x-substitute_of_Min)/(substitute_of_Max-substitute_of_Min)
  x_scaled[x_scaled<0] = 0
  x_scaled[x_scaled>1] = 1
  x_scaled
}



#' Compute clone-level mean expression from single-cell data
#'
#' For each patient, groups cells by clone_id and computes the mean
#' expression of each gene within each clone. This reduces single-cell
#' resolution to subclone resolution, matching the original PERCEPTION pipeline.
#'
#' @param expression_matrix Matrix or data frame. Gene expression matrix
#'        with genes as rows and cells as columns.
#' @param cell_clone_map Data frame with columns \code{cell_id} and \code{clone_id},
#'        mapping each cell to its clone. Typically from Seurat clustering.
#' @param patient_ids Character vector. Patient ID for each cell, same length as
#'        columns of expression_matrix. If NULL, all cells are assumed from one patient.
#'
#' @return A named list of matrices, one per patient. Each matrix has genes as rows
#'         and clone IDs as columns, with mean expression values.
#' @export
clone_mean_expression <- function(expression_matrix, cell_clone_map, patient_ids = NULL) {

  if (!all(c("cell_id", "clone_id") %in% colnames(cell_clone_map))) {
    stop("cell_clone_map must have columns 'cell_id' and 'clone_id'.")
  }

  cell_ids <- colnames(expression_matrix)
  if (is.null(cell_ids)) {
    stop("expression_matrix must have column names (cell IDs).")
  }

  # Filter to cells present in both expression matrix and mapping
  common_cells <- intersect(cell_ids, as.character(cell_clone_map$cell_id))
  if (length(common_cells) == 0) {
    stop("No matching cells between expression matrix and cell_clone_map.")
  }

  cell_clone_map <- cell_clone_map[match(common_cells, cell_clone_map$cell_id), ]
  expr_subset <- expression_matrix[, common_cells, drop = FALSE]

  if (is.null(patient_ids)) {
    patient_ids <- rep("patient1", ncol(expr_subset))
  }
  names(patient_ids) <- colnames(expr_subset)

  # Build per-patient clone mean expression
  unique_patients <- unique(patient_ids)
  result <- list()

  for (pat in unique_patients) {
    pat_cells <- names(patient_ids[patient_ids == pat])
    pat_map <- cell_clone_map[cell_clone_map$cell_id %in% pat_cells, ]

    if (nrow(pat_map) == 0) next

    clone_ids <- unique(pat_map$clone_id)
    clone_expr <- do.call(cbind, lapply(clone_ids, function(cl) {
      cl_cells <- as.character(pat_map$cell_id[pat_map$clone_id == cl])
      cl_cells <- intersect(cl_cells, colnames(expr_subset))
      if (length(cl_cells) == 0) return(NULL)
      if (length(cl_cells) == 1) {
        return(expr_subset[, cl_cells])
      }
      rowMeans_if_one_row(data.frame(expr_subset[, cl_cells, drop = FALSE]))
    }))

    if (is.null(clone_expr) || length(clone_expr) == 0) next

    if (!is.matrix(clone_expr)) {
      clone_expr <- matrix(clone_expr, nrow = nrow(expr_subset))
      colnames(clone_expr) <- clone_ids[1]
      rownames(clone_expr) <- rownames(expr_subset)
    } else {
      colnames(clone_expr) <- clone_ids
    }

    # Use patient_clone naming convention for downstream compatibility
    colnames(clone_expr) <- paste(pat, colnames(clone_expr), sep = "_")
    result[[pat]] <- clone_expr
  }

  return(result)
}



#' Z-score scale killing values across patients
#'
#' Applies base R \code{scale()} to each drug column across all patients,
#' centering to mean 0 and standard deviation 1. This ensures comparability
#' of predicted viability scores between different drugs before patient-level
#' aggregation. Matches the original PERCEPTION pipeline.
#'
#' @param clone_killing_df Data frame. Must have columns 'patient' and 'clone_id',
#'        plus one or more drug columns with predicted viability values.
#'
#' @return A data frame with the same structure, but drug columns z-score scaled.
#' @export
zscore_killing <- function(clone_killing_df) {

  if (!"patient" %in% colnames(clone_killing_df)) {
    stop("clone_killing_df must have a 'patient' column.")
  }

  drug_cols <- setdiff(colnames(clone_killing_df), c("patient", "clone_id"))

  if (length(drug_cols) == 0) {
    stop("No drug columns found (expected columns other than 'patient' and 'clone_id').")
  }

  clone_killing_df[drug_cols] <- lapply(clone_killing_df[drug_cols], function(col) {
    as.numeric(scale(as.numeric(col)))
  })

  return(clone_killing_df)
}



# topXPercentValue of a vector
topXPercentValue <- function(vec, X_percentile=95){
  vec=na.omit(vec)
  len=length(vec)
  vec=sort(vec)
  vec[ceiling(len*(X_percentile/100))]
}
