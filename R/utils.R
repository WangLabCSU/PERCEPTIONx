# Opposite of standard %in% function
'%!in%' <- function(x,y)!('%in%'(x,y))


# RowMeans functions considering the boundary case where mat has a single
# column (one cell / sample) — rowMeans would error on a one-column matrix.
rowMeans_if_one_row <- function(mat){
  if(ncol(mat)>1){
    return(rowMeans(mat))
  } else {
    return(mat)
  }
}


# Count the number of NAs in each row of a matrix
count_row_NAs <- function(df){
  apply(df, 1, function(x) sum(is.na(x)))
}


# Handing a error when a function is run and returning NA in case of error
# instead of stopping the task.
err_handle <- function(x){ tryCatch(x, error=function(e){NA}) }


# Lowercase and strip non-alphanumeric characters (keep underscores) — used
# to make inconsistent drug / gene names comparable (e.g. "Drug-1" vs "Drug_1").
stripall2match <- function(x){
  tolower(gsub('[^A-z0-9]', '', x) )
}


#' Parse clone keys in "Patient@@Clone" format
#'
#' Splits a character vector of clone identifiers (e.g. "PAT_001@@0") into a
#' data frame with \code{patient} and \code{clone_id} columns.
#'
#' @param x Character vector of clone keys in \code{"patient@@clone_id"} format.
#' @return A data frame with columns \code{patient} and \code{clone_id}.
#'
#' @seealso [build_clone_key()] for the inverse operation.
#' @export
parse_clone_keys <- function(x) {
  parts <- strsplit(as.character(x), "@@")
  data.frame(
    patient  = sapply(parts, `[`, 1),
    clone_id = sapply(parts, `[`, 2),
    stringsAsFactors = FALSE
  )
}

#' Build a clone key from patient and clone id
#'
#' The inverse of [parse_clone_keys()]: pastes patient and clone identifiers
#' into a single \code{"patient@@clone_id"} key. The \code{"@@"} separator is
#' the canonical clone-key format used throughout PERCEPTIONx.
#'
#' @param patient Character vector of patient IDs.
#' @param clone_id Character vector of clone IDs.
#' @return Character vector of \code{"patient@@clone_id"} keys.
#'
#' @seealso [parse_clone_keys()] for splitting keys back apart.
#' @export
build_clone_key <- function(patient, clone_id) {
  paste(as.character(patient), as.character(clone_id), sep = "@@")
}

