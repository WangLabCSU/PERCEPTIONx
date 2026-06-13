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



#' Rank-normalize each column of a metrix
#'
#' Convert each column to ranks and divide by column length.
#'
#' @param mat A numeric matrix. Rows = genes, columns = cells.
#' @return A rank-normalized matrix.
#' @export
rank_normalization_mat <- function(mat){
  apply(
    mat,
    2,
    function(x) rank(x, ties.method = "average")/length(x))
}



#' Scale vector to 0-1 range with outlier robustness
#'
#' Uses the 5th and 95th percentiles as thresholds for outliers,
#' capping values outside the range to 0 or 1.
#'
#' @param x Numeric vector.
#' @return Numeric vector scaled to the 0-1 range.
#' @export
range01 <- function(x) {
  substitute_of_Min <- topXPercentValue(vec = x, X_percentile = 5)
  substitute_of_Max <- topXPercentValue(vec = x, X_percentile = 95)
  x_scaled <- (x - substitute_of_Min) / (substitute_of_Max - substitute_of_Min)
  x_scaled[x_scaled < 0] = 0
  x_scaled[x_scaled > 1] = 1
  x_scaled
}



#' Get value at a given percentile of a vector
#'
#' Returns the value at the X-th percentile of a sorted vector (after removing NAs).
#'
#' @param vec Numeric vector.
#' @param X_percentile Integer. Percentile to retrieve (0-100). Default = 95.
#' @return Numeric value at the specified percentile.
#' @keywords internal
#' @noRd
topXPercentValue <- function(vec, X_percentile = 95) {
  vec = na.omit(vec)
  len = length(vec)
  vec = sort(vec)
  vec[ceiling(len * (X_percentile / 100))]
}
