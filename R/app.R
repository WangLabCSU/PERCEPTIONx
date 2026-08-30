#' Launch PERCEPTIONx Shiny Dashboard
#'
#' Starts the interactive PERCEPTIONx dashboard for drug response prediction.
#' The dashboard provides a graphical interface for data loading, model training,
#' prediction, and visualization.
#'
#' Under a Shiny Server deployment the server executes this function and
#' expects a shiny app \emph{object} back (a nested \code{runApp()} would
#' fail), so the function detects Shiny Server via the \code{SHINY_PORT}
#' environment variable it sets for every R worker (the same signal the shiny
#' package itself uses internally) and returns the app object instead of
#' starting its own server. This lets the same one-liner work both locally and
#' inside \code{apps/<name>/app.R}.
#'
#' @return A Shiny app object. When run interactively (no \code{SHINY_PORT}
#'   environment variable set) the app is also launched in the user's browser.
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
  app_deps <- c("bslib", "DT", "plotly", "waiter", "thematic", "callr")
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

  # Shiny Server mode: the server runs app.R and expects a shiny app OBJECT
  # back — a nested runApp() fails there. Detect it via SHINY_PORT (Shiny
  # Server sets it for every R worker; the shiny package's own inShinyServer()
  # uses the same signal) and hand over the object; runApp() below is only for
  # local/interactive use.
  if (nzchar(Sys.getenv("SHINY_PORT", ""))) {
    return(.perception_app_object(appDir))
  }

  shiny::runApp(appDir, display.mode = "normal")
}

# Build the shiny app object from the installed app directory. Shiny Server
# deployments (apps/<name>/app.R) call run_perception_app() which routes here.
# The app.R inside inst/shiny/app ends with shinyApp(ui, server), so evaluating
# its top-level expressions yields the app object. The working directory is
# switched to the app dir during evaluation so its relative source("R/*.R")
# module loads resolve.
.perception_app_object <- function(appDir) {
  app_file <- file.path(appDir, "app.R")
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(appDir)
  env <- new.env(parent = globalenv())
  exprs <- parse(app_file)
  obj <- NULL
  for (e in exprs) obj <- eval(e, envir = env)
  if (!inherits(obj, "shiny.appobj")) {
    stop("Failed to build the PERCEPTIONx app object (app.R did not yield a shinyApp).",
         call. = FALSE)
  }
  obj
}
