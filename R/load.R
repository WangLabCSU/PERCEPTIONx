#' Get GitHub download mirrors for model files
#'
#' @keywords internal
#' @noRd
get_github_mirrors <- function() {
  c(
    # Original GitHub (default)
    "https://github.com",
    # Chinese mirrors (Verified working)
    "https://gh-proxy.com/https://github.com",
    "https://ghproxy.net/https://github.com",
    "https://moeyy.cn/gh-proxy/https://github.com",
    "https://github.akams.cn/https://github.com",
    "http://toolwa.com/github/https://github.com",
    "https://v6.gh-proxy.org/https://github.com",
    "https://gh-proxy.org/https://github.com",
    "https://ghfast.top/https://github.com",
    "https://download.githubcdn.com?url=https://github.com",
    "https://proxy.gitwarp.top/https://github.com"
  )
}



#' Download file with mirror fallback and speed monitoring
#'
#' @param urls Character vector of URLs to try
#' @param destfile Destination file path
#' @param quiet Logical, suppress progress messages. Default = FALSE.
#' @param speed_threshold Numeric, speed threshold in KB/s to warn user. Default = 10.
#' @param timeout_seconds Numeric, timeout for each download attempt in seconds. Default = 300.
#' @param retries Integer, number of retries for each mirror. Default = 0.
#'
#' @return TRUE if successful, FALSE otherwise
#' @keywords internal
#' @noRd
download_with_mirrors <- function(urls, destfile, quiet = FALSE,
                                  speed_threshold = 10,
                                  timeout_seconds = 30,
                                  retries = 0) {
  for (i in seq_along(urls)) {
    if (!quiet) message("Trying mirror ", i, "/", length(urls))

    for (attempt in 0:retries) {
      if (attempt > 0 && !quiet) message("  Retry attempt ", attempt)

      start_time <- Sys.time()
      success <- FALSE
      timeout_occurred <- FALSE

      tryCatch({
        suppressWarnings({
          download.file(urls[i], destfile, mode = "wb", quiet = quiet)
        })
        success <- TRUE
      }, error = function(e) {
        if (grepl("Timeout", e$message, ignore.case = TRUE)) {
          timeout_occurred <<- TRUE
          if (!quiet && attempt == retries) {
            message("  Mirror ", i, " timed out after ", timeout_seconds, " seconds")
          }
        } else {
          if (!quiet && attempt == retries) {
            message("  Mirror ", i, " failed: ", e$message)
          }
        }
        if (file.exists(destfile)) file.remove(destfile)
      })

      if (success && file.exists(destfile) && file.size(destfile) > 0) {
        end_time <- Sys.time()
        elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
        file_size_mb <- file.size(destfile) / 1024 / 1024
        speed_kb_s <- (file.size(destfile) / 1024) / elapsed

        if (!quiet) {
          message(sprintf("✓ Downloaded %.2f MB in %.1f sec (%.0f KB/s)",
                          file_size_mb, elapsed, speed_kb_s))
        }
        return(TRUE)
      }

      if (timeout_occurred && attempt < retries) {
        if (!quiet) message("  Retrying mirror ", i, " (timeout)...")
        next
      }
    }
  }
  return(FALSE)
}



#' Get current download mirrors for PERCEPTION
#'
#' @return Character vector of mirror URLs.
#' @export
#'
#' @examples
#' get_perception_mirrors()
get_perception_mirrors <- function() {
  getOption("PERCEPTION.download_mirrors", get_default_perception_mirrors())
}

#' Get default download mirrors
#'
#' @keywords internal
#' @noRd
get_default_perception_mirrors <- function() {
  c(
    "https://github.com",
    "https://gh-proxy.com/https://github.com",
    "https://ghproxy.net/https://github.com",
    "https://moeyy.cn/gh-proxy/https://github.com",
    "https://github.akams.cn/https://github.com",
    "http://toolwa.com/github/https://github.com",
    "https://v6.gh-proxy.org/https://github.com",
    "https://gh-proxy.org/https://github.com",
    "https://ghfast.top/https://github.com",
    "https://download.githubcdn.com?url=https://github.com",
    "https://proxy.gitwarp.top/https://github.com"
  )
}



#' Add custom download mirror
#'
#' @param url Character string. The mirror URL to add.
#' @param position Character. Where to add: "first", "last", or "before_github".
#'        Default: "first".
#'
#' @return Invisibly returns the updated mirror list.
#' @export
#'
#' @examples
#' \dontrun{
#' add_perception_mirror("https://my-mirror.com/https://github.com")
#' }
add_perception_mirror <- function(url, position = c("first", "last", "before_github")) {
  position <- match.arg(position)

  if (!grepl("^https?://", url)) {
    stop("Invalid URL. Must start with http:// or https://")
  }

  current_mirrors <- getOption("PERCEPTION.download_mirrors", get_default_perception_mirrors())
  url <- sub("/$", "", url)

  if (url %in% current_mirrors) {
    message("Mirror already exists")
    return(invisible(current_mirrors))
  }

  new_mirrors <- switch(position,
                        "first" = c(url, current_mirrors),
                        "last" = c(current_mirrors, url),
                        "before_github" = {
                          github_idx <- which(current_mirrors == "https://github.com")
                          if (length(github_idx) == 0) {
                            c(current_mirrors, url)
                          } else {
                            c(current_mirrors[1:(github_idx - 1)], url, current_mirrors[github_idx:length(current_mirrors)])
                          }
                        }
  )

  options(PERCEPTION.download_mirrors = new_mirrors)
  message("Added mirror: ", url)
  invisible(new_mirrors)
}



#' List current download mirrors
#'
#' @export
list_perception_mirrors <- function() {
  mirrors <- get_perception_mirrors()
  for (i in seq_along(mirrors)) {
    cat(sprintf("%d. %s\n", i, mirrors[i]))
  }
  invisible(mirrors)
}



#' Reset mirrors to default
#'
#' @export
reset_perception_mirrors <- function() {
  options(PERCEPTION.download_mirrors = NULL)
  message("Mirrors reset to default")
  invisible(get_default_perception_mirrors())
}



#' Load pre-built model of provided drugs
#'
#' Downloads model files from GitHub Release if not cached locally.
#' Supports automatic mirror fallback for users in different regions.
#'
#' @param ... One or more drug names (e.g., "erlotinib", "gefitinib").
#' @param dest Directory to save downloaded models. Default = "./models".
#' @param read Whether to read and return the downloaded model(s) as a named list.
#'        Default = FALSE (download only).
#' @param timeout_seconds Numeric, timeout for each download attempt in seconds. Default = 30.
#' @param retries Integer, number of retries for each mirror. Default = 0.
#'
#' @return If \code{read = TRUE}, a named list of model objects (names = drug names).
#'         If \code{read = FALSE}, invisibly returns NULL (files are saved to disk only).
#'         The returned list can be directly passed to \code{predict_drugs()},
#'         \code{compare_performance()}, or \code{get_significant_models()}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Download and load models as a named list
#' models <- load_model("erlotinib", "gefitinib", read = TRUE)
#' # models$erlotinib, models$gefitinib
#'
#' # Use directly with predict_drugs
#' pred <- predict_drugs(models, expr_rnorm)
#'
#' # Download without loading (for later use)
#' load_model("erlotinib", dest = "./models")
#' }
load_model <- function(..., dest = "./models", read = FALSE, timeout_seconds = 30, retries = 0) {
  drugs <- tolower(c(...))
  if (length(drugs) == 0) stop("Drug list is empty.")

  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE)
  }

  mirrors <- get_perception_mirrors()
  result <- list()

  for (drug in drugs) {
    file_path <- file.path(dest, paste0(drug, ".RDS"))

    if (!file.exists(file_path)) {
      urls <- paste0(mirrors, "/SunPast/PERCEPTION/releases/download/models-v1/", drug, ".RDS")
      message("Downloading model for: ", drug)

      if(!download_with_mirrors(urls, file_path, quiet = FALSE,
                                timeout_seconds = timeout_seconds,
                                retries = retries)){
        stop("Failed to download model for: ", drug)
      }

      message("-> Successfully downloaded: ", drug, ".RDS")

    } else {
      message("Found cached model: ", drug)
    }

    if (read) {
      result[[drug]] <- readRDS(file_path)
      message("-> Successfully loaded: ", drug)
    }
  }

  if (read && length(result) > 0) {
    return(result)
  }

  return(invisible(NULL))
}



#' Download filtered DepMap data
#'
#' Downloads the filtered DepMap RDS file for training models from GitHub Release.
#' The file contains bulk expression, scRNA expression, drug response annotations,
#' and cell line metadata required for PERCEPTION model training.
#' Supports automatic mirror fallback for users in different regions.
#'
#' @param dest Directory to save the downloaded file. Default = ".".
#' @param read Whether to read the data and assign to global environment as "DepMap".
#'        Default = FALSE.
#' @param timeout_seconds Numeric, timeout for each download attempt in seconds. Default = 600.
#' @param retries Integer, number of retries for each mirror. Default = 1.
#'
#' @return Invisibly returns the DepMap object if read = TRUE, otherwise NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' # Download only (no loading)
#' load_depmap()
#'
#' # Download, read, and assign to global environment
#' load_depmap(read = TRUE)
#' # Then access DepMap$expression_rnorm, DepMap$scRNA_complete, etc.
#' }
load_depmap <- function(dest = ".", read = FALSE,
                        timeout_seconds = 600, retries = 1) {
  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE)
  }

  destfile <- file.path(dest, "DepMap.RDS")

  if (!file.exists(destfile)) {
    mirrors <- get_perception_mirrors()
    urls <- paste0(mirrors, "/SunPast/PERCEPTION/releases/download/depmap/DepMap.RDS")

    message("Downloading DepMap.RDS (~567 MB). This may take several minutes...")

    if (!download_with_mirrors(urls, destfile, quiet = FALSE,
                               timeout_seconds = timeout_seconds,
                               retries = retries)) {
      message("\nAutomatic download failed. Please try manual download:")
      message("  1. Download from: https://github.com/SunPast/PERCEPTION/releases/tag/depmap")
      message("  2. Save the file to: ", destfile)
      message("  3. Then run load_depmap(read = TRUE) again")
      stop("Manual download required", call. = FALSE)
    }
    message("Successfully downloaded DepMap.RDS")
  } else {
    message("File already exists: ", destfile)
  }

  if (read) {
    message("Reading DepMap.RDS...")
    DepMap <- readRDS(destfile)
    assign("DepMap", DepMap, envir = .GlobalEnv)
    message("Assigned 'DepMap' to global environment.")
    return(invisible(DepMap))
  }

  return(invisible(NULL))
}
