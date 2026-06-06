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
#' @param dest Directory to save downloaded models. Default = ".".
#' @param read Whether to read and load the downloaded model(s) into global environment.
#' @param timeout_seconds Numeric, timeout for each download attempt in seconds. Default = 120.
#' @param retries Integer, number of retries for each mirror. Default = 0.
#'
#' @return Invisibly returns a list of model objects. Creates individual variables in global environment.
#' @export
#'
#' @examples
#' \dontrun{
#' # Download and load a single model
#' load_model("erlotinib", read = TRUE)
#'
#' # Download multiple models
#' load_model("erlotinib", "gefitinib", "osimertinib", read = TRUE)
#'
#' # Download without loading (for later use)
#' load_model("erlotinib", dest = "./models")
#' }
load_model <- function(..., dest = "./models", read = FALSE, timeout_seconds = 30, retries = 0) {
  drugs <- tolower(c(...))
  if (length(drugs) == 0) stop("Drug list is empty.")

  available_drugs <- c(...)

  invalid_drugs <- drugs[!drugs %in% available_drugs]
  if (length(invalid_drugs) > 0) {
    stop("Some input drugs do not exist: ", paste(invalid_drugs, collapse = ", "))
  }

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

      # timeout & retries
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
    list2env(result, envir = .GlobalEnv)
    return(invisible(result))
  }

  return(invisible(NULL))
}



#' Download DepMapv12.RDS
#'
#' Downloads the required DepMap RDS file for training models from Zenodo.
#' The file contains bulk expression, scRNA expression, drug response annotations,
#' and cell line metadata required for PERCEPTION model training.
#'
#' @param dest Directory to save the downloaded file. Default = ".".
#' @param read Whether to read the data and assign to global environment as "DepMap".
#'        Default = FALSE.
#' @param speed_threshold Numeric, speed threshold in KB/s to suggest manual download.
#'        Default 100.
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
load_depmap <- function(dest = ".", read = FALSE, speed_threshold = 10,
                        timeout_seconds = 600, retries = 1) {
  # Create destination directory if needed
  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE)
  }

  destfile <- file.path(dest, "DepMapv12.RDS")

  if (!file.exists(destfile)) {
    message("Downloading 883.6 MB file. This may take several minutes...")
    message("Recommend you download directly from: https://zenodo.org/record/7860559/files/DepMapv12.RDS")

    # Use official Zenodo URL
    url <- "https://zenodo.org/record/7860559/files/DepMapv12.RDS"

    # Pass speed_threshold to download function
    if (!download_with_mirrors(url, destfile, quiet = FALSE, speed_threshold = speed_threshold,
                               timeout_seconds = timeout_seconds,
                               retries = retries)) {
      # If download fails, give detailed manual download instructions
      message("\n❌ Automatic download failed. Please try manual download:")
      message("  1. Download from: https://zenodo.org/record/7860559/files/DepMapv12.RDS")
      message("  2. Save the file to: ", destfile)
      message("  3. Then run load_depmap(read = TRUE) again")
      stop("Manual download required", call. = FALSE)
    }
    message("Successfully downloaded DepMapv12.RDS")
  } else {
    message("File already exists: ", destfile)
  }

  if (read) {
    message("Reading DepMapv12.RDS...")
    DepMap <- readRDS(destfile)
    assign("DepMap", DepMap, envir = .GlobalEnv)
    message("Assigned 'DepMap' to global environment.")
    return(invisible(DepMap))
  }

  return(invisible(NULL))
}
