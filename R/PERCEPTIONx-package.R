# nocov start

#' @title Personalized Drug Response Prediction from Single-Cell Transcriptomics
#'
#' @description An R implementation of the PERCEPTION approach (Sinha et al., 2024)
#' for predicting personalized drug response from single-cell tumor
#' transcriptomics. It trains drug-response models on large-scale cell-line
#' screens, predicts sensitivity at clone and patient level, and ships an
#' interactive visualization suite plus a Shiny web application.
#'
#' @section License:
#' MIT + file LICENSE
#'
#' @docType package
#' @name PERCEPTIONx-package
#' @aliases PERCEPTIONx
#' @keywords internal
#'
"_PACKAGE"


.onAttach <- function(libname, pkgname) {
  if (startup_message_allowed()) {
    cli::cli_process_done(.envir = .time_record)
  }
}

.onLoad <- function(libname, pkgname) {
  if (startup_message_allowed()) {
    pkg_version <- utils::packageVersion(pkgname)
    .time_record$pkgname <- pkgname
    .time_record$pkg_version <- pkg_version
    cli::cli_progress_step(
      msg = "Loading {.pkg {pkgname}} v{pkg_version}",
      msg_done = "{.pkg {pkgname}} v{pkg_version} loaded",
      spinner = TRUE,
      .envir = .time_record
    )
  }

  invisible()
}

.time_record <- new.env()


startup_message_allowed <- function() {
  allowed <- FALSE

  withRestarts(
    {
      signalCondition(structure(
        list(message = ".__startup_probe__."),
        class = c(
          "packageStartupProbe",
          "packageStartupMessage",
          "condition"
        )
      ))
      allowed <- TRUE
    },
    muffleMessage = function() NULL
  )

  allowed
}

## usethis namespace: start
#' @importFrom data.table .N
#' @importFrom data.table %chin%
#' @importFrom data.table :=
#' @importFrom data.table .SD
## usethis namespace: end
NULL

# nocov end
