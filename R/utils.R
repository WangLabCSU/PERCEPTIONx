# Opposite of standard %in% function (internal, not exported)
'%!in%' <- function(x, y) !('%in%'(x, y))


#' Compute column-wise maximum
#'
#' @param colData Matrix or data frame.
#' @return Numeric vector of column maxima.
#' @keywords internal
#' @noRd
colMax <- function(colData) {
  apply(colData, MARGIN = c(2), max)
}


#' Compute column-wise median
#'
#' @param colData Matrix or data frame.
#' @return Numeric vector of column medians.
#' @keywords internal
#' @noRd
colMedian <- function(colData) {
  apply(colData, MARGIN = c(2), median)
}


#' Compute column-wise minimum
#'
#' @param colData Matrix or data frame.
#' @return Numeric vector of column minima.
#' @keywords internal
#' @noRd
colMin <- function(colData) {
  apply(colData, MARGIN = c(2), min)
}


#' Compute row-wise maximum
#'
#' @param colData Matrix or data frame.
#' @return Numeric vector of row maxima.
#' @keywords internal
#' @noRd
rowMax <- function(colData) {
  apply(colData, MARGIN = c(1), max)
}


#' Compute row-wise minimum
#'
#' @param colData Matrix or data frame.
#' @return Numeric vector of row minima.
#' @keywords internal
#' @noRd
rowMin <- function(colData) {
  apply(colData, MARGIN = c(1), min)
}


#' Cap values at 1
#'
#' Sets all values greater than 1 to 1.
#'
#' @param x Numeric vector or matrix.
#' @return Input with values capped at 1.
#' @keywords internal
#' @noRd
capping_at1 <- function(x) {
  x[x > 1] = 1
  x
}


#' Cap values at 0
#'
#' Sets all values less than 0 to 0.
#'
#' @param x Numeric vector or matrix.
#' @return Input with values capped at 0.
#' @keywords internal
#' @noRd
capping_at0 <- function(x) {
  x[x < 0] = 0
  x
}


#' Row means handling single-column edge case
#'
#' Computes row means, returning the column as-is when the matrix has only one column.
#'
#' @param mat Matrix or data frame.
#' @return Numeric vector of row means.
#' @keywords internal
#' @noRd
rowMeans_if_one_row <- function(mat) {
  if (ncol(mat) > 1) {
    return(rowMeans(mat))
  } else {
    return(mat)
  }
}


#' Count NAs in each row
#'
#' @param df Matrix or data frame.
#' @return Integer vector of NA counts per row.
#' @keywords internal
#' @noRd
count_row_NAs <- function(df) {
  apply(df, 1, function(x) sum(is.na(x)))
}


#' Error handler that returns NA instead of stopping
#'
#' Wraps an expression in tryCatch, returning NA on error.
#' Useful for continuing batch operations when individual items may fail.
#'
#' @param x Expression to evaluate.
#' @return Result of x, or NA if an error occurs.
#' @export
err_handle <- function(x) { tryCatch(x, error = function(e) NA) }


#' Custom head showing first 5 rows and columns
#'
#' @param mat Matrix or data frame.
#' @return Subset of mat (up to 5x5).
#' @keywords internal
#' @noRd
myhead <- function(mat) {
  mat[1:min(5, nrow(mat)), 1:min(5, ncol(mat))]
}


#' Subset columns by name
#'
#' @param mat Matrix or data frame.
#' @param column_Names Character vector of column names to select.
#' @return Subset of mat with matching columns.
#' @keywords internal
#' @noRd
colSubset <- function(mat, column_Names) {
  mat[, na.omit(match(column_Names, colnames(mat)))]
}


#' Subset rows by name
#'
#' @param mat Matrix or data frame.
#' @param row_Names Character vector of row names to select.
#' @return Subset of mat with matching rows.
#' @keywords internal
#' @noRd
rowSubset <- function(mat, row_Names) {
  mat[na.omit(match(row_Names, rownames(mat))), ]
}


#' Convert factor to numeric
#'
#' Safely converts a factor vector to numeric via character conversion.
#'
#' @param x Factor or vector.
#' @return Numeric vector.
#' @keywords internal
#' @noRd
factor2numeric <- function(x) {
  as.numeric(as.character(x))
}


#' Strip non-alphanumeric characters and lowercase
#'
#' Primarily used to facilitate inconsistent naming comparison (e.g., drug names).
#'
#' @param x Character vector.
#' @return Lowercase character vector with non-alphanumeric characters removed.
#' @export
stripall2match <- function(x) {
  tolower(gsub('[^A-z0-9]', '', x))
}


#' Split strings and retrieve a specific element
#'
#' Splits each string by a delimiter and returns the element at the specified position.
#'
#' @param infunc_list Character vector. Strings to split.
#' @param infunc_split_by Character. Delimiter. Default = "_".
#' @param retreiving_onject_id Integer. Position to retrieve after splitting. Default = 1.
#' @return Character vector of retrieved elements.
#' @export
strsplit_customv0 <- function(infunc_list,
                              infunc_split_by = '_',
                              retreiving_onject_id = 1) {
  sapply(strsplit(infunc_list, split = infunc_split_by), function(x) x[retreiving_onject_id])
}
