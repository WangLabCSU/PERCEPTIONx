#' Launch PERCEPTIONx Shiny Dashboard
#'
#' Starts the interactive PERCEPTIONx dashboard for drug response prediction.
#' The dashboard provides a graphical interface for data loading, model training,
#' prediction, and visualization.
#'
#' @return A Shiny app object (invisibly). The function is called for its side effect
#'   of launching the app in the user's browser.
#'
#' @examples
#' \dontrun{
#' library(PERCEPTIONx)
#' run_perception_app()
#' }
#'
#' @export
run_perception_app <- function() {
  appDir <- system.file("shiny", "app", package = "PERCEPTIONx")
  if (appDir == "") {
    stop("Could not find PERCEPTIONx dashboard. Try re-installing the package with: install.packages('PERCEPTIONx', repos = NULL, type = 'source')",
         call. = FALSE)
  }

  # The dashboard hard-depends on these packages. They are in Suggests (not
  # Imports) to keep the core pipeline lean, so check them here and fail fast
  # with an actionable message instead of a cryptic library() error.
  app_deps <- c("bslib", "DT", "plotly", "waiter", "thematic")
  missing <- app_deps[!vapply(app_deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("The PERCEPTIONx dashboard requires these packages, which are not installed: ",
         paste(missing, collapse = ", "),
         "\nInstall them with: install.packages(c(",
         paste(shQuote(missing), collapse = ", "), "))",
         call. = FALSE)
  }

  # Set max upload size to 1 GB (for DepMap ~567MB + models)
  options(shiny.maxRequestSize = 1024 * 1024^2)

  shiny::runApp(appDir, display.mode = "normal")
}
