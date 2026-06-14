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
rank_normalization_mat <- function(mat){
  apply(
    mat,
    2,
    function(x) rank(x, ties.method = "average")/length(x))
}



#' Change range to 0-1
#'
#' Scales a numeric vector to the 0-1 range using the 5th and 95th percentiles
#' as thresholds to handle outliers.
#'
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
  x_scaled
}



# topXPercentValue of a vector
topXPercentValue<-function(vec, X_percentile=95){
  vec=na.omit(vec)
  len=length(vec)
  vec=sort(vec)
  vec[ceiling(len*(X_percentile/100))]
}
