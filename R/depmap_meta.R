#' Extract lightweight metadata from a DepMap.RDS file
#'
#' Reads a DepMap.RDS and returns ONLY the small metadata the UI needs (gene
#' list, drug list, lineages, component dimensions) — never the full multi-GB
#' object. This runs inside a background worker process; the main Shiny
#' process must never deserialize the 567 MB file (it would block the UI and
#' spike its RAM).
#'
#' @param depmap_path Path to the DepMap.RDS file.
#' @param cache_file Optional path to a sidecar cache (DepMap_meta.RDS). When
#'   present and newer than \code{depmap_path} it is read directly and
#'   returned; otherwise the full file is parsed and the cache (re)written.
#'
#' @return A list with \code{loaded}, \code{genes}, \code{drugs},
#'   \code{lineages}, \code{components}.
#' @noRd
extract_depmap_meta <- function(depmap_path, cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file) && file.exists(depmap_path)) {
    # Only trust a cache that is newer than the source DepMap.RDS itself:
    # if the RDS was replaced (re-download / different version) the cached
    # genes/drugs would be silently stale. A corrupt cache falls through to
    # a full re-extraction instead of erroring. (Both files must exist —
    # file.mtime() returns NA for a missing file and the comparison would
    # then produce if(NA) instead of a clean fall-through.)
    if (file.mtime(cache_file) >= file.mtime(depmap_path)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (!is.null(cached)) return(cached)
    }
  }
  DepMap <- readRDS(depmap_path)
  if (!is.list(DepMap)) {
    stop("The file is not a list. Expected a DepMap.RDS from the PERCEPTIONx release (or an object with the same structure).")
  }
  required_fields <- c("secondary_prism", "secondary_screen_drugAnnotation",
                       "expression_rnorm", "scRNA_complete", "scRNA_subset_rnorm",
                       "CPM_scRNA_CCLE_rnorm", "annotation_20Q4",
                       "metadata_CPM_scRNA", "expression_20Q4")
  missing_fields <- setdiff(required_fields, names(DepMap))
  if (length(missing_fields) > 0) {
    stop("The RDS is missing required components: ",
         paste(missing_fields, collapse = ", "),
         ". Re-save it with these components, or use the built-in 'Download & Load' (known-good).")
  }
  nms <- names(DepMap)
  components <- lapply(nms, function(nm) {
    obj <- DepMap[[nm]]
    if (is.matrix(obj) || is.data.frame(obj)) {
      list(nrow = nrow(obj), ncol = ncol(obj),
           cols_preview = head(colnames(obj), 12))
    } else {
      list(nrow = NA_integer_, ncol = NA_integer_, cols_preview = character(0))
    }
  })
  names(components) <- nms
  lineages <- sort(unique(as.character(DepMap$annotation_20Q4$lineage)))
  lineages <- lineages[!is.na(lineages) & nzchar(lineages)]
  meta <- list(
    loaded = TRUE,
    genes  = rownames(DepMap$expression_rnorm),
    drugs  = unique(as.character(DepMap$secondary_screen_drugAnnotation$CommonName)),
    lineages = lineages,
    components = components
  )
  if (!is.null(cache_file)) saveRDS(meta, cache_file)
  rm(DepMap); gc()
  meta
}
