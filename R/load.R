#' PERCEPTIONx Data Loading
#'
#' @name load_perception
#' @keywords internal
#' @importFrom utils download.file
NULL

# Package-level environment for storing DepMap data
# (avoids assigning to .GlobalEnv which triggers R CMD check NOTE)
.depmap_env <- new.env(parent = emptyenv())

# Size of the official DepMap.RDS release asset (bytes). Used to reject
# truncated downloads / partial cache files so they are never cached or
# parsed. Bump this if the release asset is ever regenerated with a
# different size.
.depmap_release_size <- 594589700

#' Get the DepMap dataset
#'
#' Retrieves the DepMap dataset from the package-level cache.
#' This is the recommended way to access DepMap data after loading
#' with \code{load_depmap(read = TRUE)}.
#'
#' @return The DepMap list object. Errors if the data has not been loaded yet.
#' @export
#'
#' @examples
#' \dontrun{
#' load_depmap(read = TRUE)
#' depmap <- get_depmap()
#' }
get_depmap <- function() {
  if (exists("DepMap", envir = .depmap_env)) {
    return(get("DepMap", envir = .depmap_env))
  }
  # Fallback: check .GlobalEnv for backward compatibility
  if (exists("DepMap", envir = .GlobalEnv)) {
    return(get("DepMap", envir = .GlobalEnv))
  }
  stop("DepMap data not loaded. Run load_depmap(read = TRUE) first.")
}

# DepMap is accessed via get_depmap() - declare to suppress R CMD check note
utils::globalVariables("DepMap")

#' Download file with mirror fallback and speed monitoring
#'
#' @param urls Character vector of URLs to try
#' @param destfile Destination file path
#' @param quiet Logical, suppress progress messages. Default = FALSE.
#' @param speed_threshold Numeric, speed threshold in KB/s to warn user. Default = 10.
#' @param timeout_seconds Numeric, timeout for each download attempt in seconds. Default = 300.
#' @param retries Integer, number of retries for each mirror. Default = 0.
#' @param expected_size Optional expected file size in bytes. When provided, a
#'        finished transfer that is SMALLER than this is treated as an
#'        interrupted/partial download: the file is discarded and the next
#'        mirror/retry is attempted (a partial file must never be cached).
#'
#' @return TRUE if successful, FALSE otherwise
#' @keywords internal
#' @noRd
download_with_mirrors <- function(urls, destfile, quiet = FALSE,
                                  speed_threshold = 10,
                                  timeout_seconds = 30,
                                  retries = 0,
                                  expected_size = NULL) {
  # Save and restore R timeout setting
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = timeout_seconds)

  # Prefer the system curl binary via R's official download.file(method = "curl")
  # (R delegates to the curl binary internally; no shell hacks). R's built-in
  # libcurl backend can be pathologically slow on some systems (measured
  # ~3 KB/s on the production server vs ~3 MB/s for the same URL via curl).
  # Fall back to libcurl when no curl binary is on PATH.
  download_method <- if (nzchar(Sys.which("curl"))) "curl" else "libcurl"

  # Extra flags for the curl method. R's curl method runs `curl url -o file`
  # WITHOUT these by default, which is what broke DepMap downloads:
  #   -L    follow the GitHub release 302 redirect to
  #         objects.githubusercontent.com (without it curl saves the 302 page)
  #   --fail  turn HTTP 4xx/5xx into a nonzero exit (a proxy returning an HTML
  #         error page would otherwise "succeed")
  #   --ssl-no-revoke  Windows only: curl.exe uses the schannel TLS backend,
  #         which by default contacts the certificate revocation servers —
  #         unreachable behind GFW, so EVERY handshake failed with
  #         "curl: (35) schannel ... revocation". This flag skips that check.
  #         Linux curl uses OpenSSL and simply ignores the flag.
  #   -C -  resume an interrupted transfer on retries, so a 567 MB file is not
  #         re-downloaded from zero every time. The first attempt is always a
  #         fresh start (curl -o overwrites), so a stale partial file can never
  #         be silently appended to.
  # (libcurl ignores `extra`; it follows redirects, checks HTTP status and has
  # no schannel CRL issue, so the fallback path needs none of this.)
  extra_base <- "-L --fail --ssl-no-revoke"

  for (i in seq_along(urls)) {
    if (!quiet) message("Trying mirror ", i, "/", length(urls))

    for (attempt in 0:retries) {
      if (attempt > 0 && !quiet) message("  Retry attempt ", attempt)

      start_time <- Sys.time()
      success <- FALSE
      timeout_occurred <- FALSE
      incomplete <- FALSE

      tryCatch({
        suppressWarnings({
          curl_extra <- if (attempt > 0) paste(extra_base, "-C -") else extra_base
          download.file(urls[i], destfile, mode = "wb", quiet = quiet,
                        method = download_method, extra = curl_extra)
        })
        success <- TRUE
      }, error = function(e) {
        if (grepl("Timeout|timed out|reached elapsed time", e$message, ignore.case = TRUE)) {
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

      # A "successful" transfer can still be truncated by a flaky connection
      # or a misbehaving proxy: verify the size against the known release
      # asset when the caller provides one. A partial file must NEVER be
      # kept — otherwise it would poison the cache and every later parse.
      if (success && file.exists(destfile)) {
        sz <- file.size(destfile)
        if (!is.null(expected_size) && is.finite(sz) && sz < expected_size) {
          incomplete <- TRUE
          if (!quiet) message(sprintf(
            "  Incomplete download: %d of %d bytes (%.1f%%), discarding and retrying",
            sz, expected_size, 100 * sz / expected_size))
          file.remove(destfile)
          success <- FALSE
        }
      }

      if (success && file.exists(destfile) && file.size(destfile) > 0) {
        end_time <- Sys.time()
        elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
        file_size_mb <- file.size(destfile) / 1024 / 1024
        speed_kb_s <- (file.size(destfile) / 1024) / elapsed

        if (!quiet) {
          message(sprintf("Downloaded %.2f MB in %.1f sec (%.0f KB/s)",
                          file_size_mb, elapsed, speed_kb_s))
        }
        if (is.finite(speed_kb_s) && speed_kb_s < speed_threshold) {
          message(sprintf("  Warning: mirror %d is slow (%.0f KB/s < %.0f KB/s threshold). Consider a different mirror.",
                          i, speed_kb_s, speed_threshold))
        }
        return(TRUE)
      }

      if ((timeout_occurred || incomplete) && attempt < retries) {
        if (!quiet) message("  Retrying mirror ", i, "...")
        next
      }
    }
  }
  return(FALSE)
}



#' Get current download mirrors
#'
#' Returns the current mirror list, including any user-added mirrors.
#' User-added mirrors are tried first by default.
#'
#' @return Character vector of mirror URLs.
#' @export
#'
#' @examples
#' get_mirrors()
get_mirrors <- function() {
  getOption("PERCEPTIONx.download_mirrors", get_default_perception_mirrors())
}

#' Get default download mirrors (excluding GitHub primary)
#'
#' @keywords internal
#' @noRd
get_default_perception_mirrors <- function() {
  c(
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



#' Add custom download mirrors
#'
#' @param urls Character vector of mirror URLs to add.
#' @param position Character. Where to add: "first" or "last". Default: "first".
#'
#' @return Invisibly returns the updated mirror list.
#' @export
#'
#' @examples
#' \dontrun{
#' add_mirrors("https://my-mirror.com/https://github.com")
#' add_mirrors(c("https://mirror1.com/https://github.com",
#'               "https://mirror2.com/https://github.com"))
#' }
add_mirrors <- function(urls, position = c("first", "last")) {
  position <- match.arg(position)

  if (!all(grepl("^https?://", urls))) {
    stop("All URLs must start with http:// or https://")
  }

  current_mirrors <- getOption("PERCEPTIONx.download_mirrors", get_default_perception_mirrors())
  urls <- sub("/$", "", urls)

  # Remove duplicates (both against current and within input)
  urls <- unique(urls)
  urls <- urls[!urls %in% current_mirrors]

  if (length(urls) == 0) {
    message("All mirrors already exist")
    return(invisible(current_mirrors))
  }

  new_mirrors <- switch(position,
                        "first" = c(urls, current_mirrors),
                        "last" = c(current_mirrors, urls))

  options(PERCEPTIONx.download_mirrors = new_mirrors)
  message("Added ", length(urls), " mirror(s)")
  invisible(new_mirrors)
}



#' List current download mirrors
#'
#' @export
list_mirrors <- function() {
  mirrors <- get_mirrors()
  for (i in seq_along(mirrors)) {
    cat(sprintf("%d. %s\n", i, mirrors[i]))
  }
  invisible(mirrors)
}



#' Reset mirrors to default
#'
#' @export
reset_mirrors <- function() {
  options(PERCEPTIONx.download_mirrors = NULL)
  message("Mirrors reset to default")
  invisible(get_default_perception_mirrors())
}



# The 44 pre-trained model drugs shipped with the package (used by
# load_model(all = TRUE) to pre-cache the full set on deployment servers).
# Keep in sync with the drug list in inst/shiny/app/R/mod_data.R.
perception_all_drugs <- c(
  "abemaciclib", "afatinib", "axitinib", "azacitidine", "cladribine",
  "clofarabine", "cobimetinib", "dabrafenib", "dasatinib", "daunorubicin",
  "decitabine", "docetaxel", "doxorubicin", "epirubicin", "erlotinib",
  "etoposide", "gefitinib", "gemcitabine", "homoharringtonine", "ibrutinib",
  "icotinib", "ixabepilone", "lapatinib", "lenvatinib", "midostaurin",
  "niraparib", "osimertinib", "paclitaxel", "palbociclib", "ponatinib",
  "romidepsin", "sunitinib", "temsirolimus", "teniposide", "thioguanine",
  "topotecan", "trametinib", "vandetanib", "vemurafenib", "vinblastine",
  "vincristine", "vindesine", "vinflunine", "vinorelbine"
)

#' Load pre-built model of provided drugs
#'
#' Downloads model files from GitHub Release if not cached locally.
#' Supports automatic mirror fallback for users in different regions.
#'
#' @param ... One or more drug names (e.g., "erlotinib", "gefitinib").
#' @param dest Directory to save downloaded models. Default: a cache-root
#'        derived path (see \code{options(PERCEPTIONx.cache_root)}), or
#'        "./models" when no cache root is set.
#' @param read Whether to read and return the downloaded model(s) as a named list.
#'        Default = FALSE (download only).
#' @param mirror Logical. If FALSE (default), download from GitHub directly.
#'        If TRUE, use mirror sites from \code{get_mirrors()}.
#' @param mirror_url Character. A specific mirror URL to use (e.g.,
#'        \code{"https://gh-proxy.com/https://github.com"}). Overrides \code{mirror}.
#' @param timeout_seconds Numeric, timeout for each download attempt in seconds. Default = 30.
#' @param retries Integer, number of retries for each mirror. Default = 0.
#' @param force Logical. If TRUE, delete an existing cached model file before
#'        downloading it again (use after an interrupted/failed download left a
#'        corrupt file behind). Default = FALSE (keep and reuse the cache).
#' @param all Logical. If TRUE, download all 44 pre-trained models (useful to
#'        pre-cache the full model set on a deployment server). Overrides any
#'        drug names given in \code{...}. Default = FALSE.
#'
#' @return If \code{read = TRUE}, a named list of model objects (names = drug names).
#'         If \code{read = FALSE}, invisibly returns NULL (files are saved to disk only).
#'         The returned list can be directly passed to \code{predict_drugs()},
#'         \code{compare_performance()}, or \code{get_significant_models()}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Download and load models from GitHub directly
#' models <- load_model("erlotinib", "gefitinib", read = TRUE)
#'
#' # Use mirror sites
#' models <- load_model("erlotinib", "gefitinib", read = TRUE, mirror = TRUE)
#'
#' # Use a specific mirror
#' models <- load_model("erlotinib", read = TRUE,
#'                      mirror_url = "https://gh-proxy.com/https://github.com")
#' }
load_model <- function(..., dest = perception_default_model_dir(), read = FALSE, mirror = FALSE,
                       mirror_url = NULL, timeout_seconds = 30, retries = 0, force = FALSE,
                       all = FALSE) {
  drugs <- tolower(c(...))
  if (all) {
    if (length(drugs) > 0) {
      warning("all = TRUE overrides the drug names given in ...; downloading all 44 models.",
              call. = FALSE)
    }
    drugs <- perception_all_drugs
  }
  if (length(drugs) == 0) {
    stop("Drug list is empty. Pass drug names or set all = TRUE to download all 44 pre-trained models.",
         call. = FALSE)
  }

  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE)
  }

  base_urls <- if (!is.null(mirror_url)) {
    mirror_url
  } else if (mirror) {
    get_mirrors()
  } else {
    "https://github.com"
  }

  result <- list()

  for (drug in drugs) {
    file_path <- file.path(dest, paste0(drug, ".RDS"))

    # force: drop a stale/corrupt cached file so the download always restarts
    if (force && file.exists(file_path)) {
      file.remove(file_path)
      message("force = TRUE: removed cached ", drug, ".RDS, re-downloading...")
    }

    if (!file.exists(file_path)) {
      urls <- paste0(base_urls, "/WangLabCSU/PERCEPTIONx/releases/download/models-v1/", drug, ".RDS")
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
#' and cell line metadata required for PERCEPTIONx model training.
#' Supports automatic mirror fallback for users in different regions.
#'
#' @param dest Directory to save the downloaded file. Default: a cache-root
#'        derived path (see \code{options(PERCEPTIONx.cache_root)}), or "." when
#'        no cache root is set.
#' @param read Whether to read the data and assign to global environment as "DepMap".
#'        Default = FALSE.
#' @param mirror Logical. If FALSE (default), download from GitHub directly.
#'        If TRUE, use mirror sites from \code{get_mirrors()}.
#' @param mirror_url Character. A specific mirror URL to use (e.g.,
#'        \code{"https://gh-proxy.com/https://github.com"}). Overrides \code{mirror}.
#' @param timeout_seconds Numeric, timeout for each download attempt in seconds. Default = 300.
#' @param retries Integer, number of retries for each mirror. Default = 1.
#' @param force Logical. If TRUE, delete an existing cached file before
#'        downloading again (use after an interrupted/failed download left a
#'        corrupt file behind). Default = FALSE (keep and reuse the cache).
#'
#' @return Invisibly returns the DepMap object if read = TRUE, otherwise NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' # Download from GitHub directly
#' load_depmap(read = TRUE)
#'
#' # Use mirror sites
#' load_depmap(read = TRUE, mirror = TRUE)
#'
#' # Force a fresh download (removes any cached copy first)
#' load_depmap(read = TRUE, mirror = TRUE, force = TRUE)
#'
#' # Use a specific mirror
#' load_depmap(read = TRUE, mirror_url = "https://gh-proxy.com/https://github.com")
#' }
load_depmap <- function(dest = perception_default_depmap_dir(), read = FALSE, mirror = FALSE, mirror_url = NULL,
                        timeout_seconds = 300, retries = 1, force = FALSE) {
  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE)
  }

  destfile <- file.path(dest, "DepMap.RDS")

  # force: drop a stale/corrupt cached file so the download always restarts
  if (force && file.exists(destfile)) {
    file.remove(destfile)
    message("force = TRUE: removed cached DepMap.RDS, re-downloading...")
  }

  # A partial file from an interrupted download must never be treated as a
  # cache hit (it would fail at readRDS). Re-download when the on-disk file is
  # smaller than the official release asset.
  if (!file.exists(destfile) ||
      is.finite(file.size(destfile)) && file.size(destfile) < .depmap_release_size) {
    if (file.exists(destfile)) {
      message(sprintf("Existing DepMap.RDS is incomplete (%d bytes < expected %d), re-downloading...",
                      file.size(destfile), .depmap_release_size))
      file.remove(destfile)
    }
    base_urls <- if (!is.null(mirror_url)) {
      mirror_url
    } else if (mirror) {
      get_mirrors()
    } else {
      "https://github.com"
    }
    urls <- paste0(base_urls, "/WangLabCSU/PERCEPTIONx/releases/download/depmap/DepMap.RDS")

    message("Downloading DepMap.RDS (~567 MB). This may take several minutes...")

    if (!download_with_mirrors(urls, destfile, quiet = FALSE,
                               timeout_seconds = timeout_seconds,
                               retries = retries,
                               expected_size = .depmap_release_size)) {
      message("\nAutomatic download failed. Please try manual download:")
      message("  1. Download from: https://github.com/WangLabCSU/PERCEPTIONx/releases/tag/depmap")
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
    assign("DepMap", DepMap, envir = .depmap_env)
    # Also assign to .GlobalEnv for backward compatibility with existing scripts
    # that reference DepMap directly
    do.call("assign", list("DepMap", DepMap, envir = .GlobalEnv))
    message("Assigned 'DepMap' to package cache and global environment.")
    message("Use get_depmap() or DepMap to access the data.")
    return(invisible(DepMap))
  }

  return(invisible(NULL))
}
