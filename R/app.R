#' Launch PERCEPTION Shiny Dashboard
#'
#' Starts the interactive PERCEPTION dashboard for drug response prediction.
#' The dashboard provides a graphical interface for data loading, model training,
#' prediction, and visualization.
#'
#' @return A Shiny app object (invisibly). The function is called for its side effect
#'   of launching the app in the user's browser.
#'
#' @examples
#' \dontrun{
#' library(PERCEPTION)
#' run_perception_app()
#' }
#'
#' @export
run_perception_app <- function() {
  appDir <- system.file("shiny", "app", package = "PERCEPTION")
  if (appDir == "") {
    stop("Could not find PERCEPTION dashboard. Try re-installing the package with: install.packages('PERCEPTION', repos = NULL, type = 'source')",
         call. = FALSE)
  }

  # Set max upload size to 1 GB (for DepMap ~567MB + models)
  options(shiny.maxRequestSize = 1024 * 1024^2)

  shiny::runApp(appDir, display.mode = "normal")
}
