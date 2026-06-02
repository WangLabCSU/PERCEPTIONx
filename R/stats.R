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
