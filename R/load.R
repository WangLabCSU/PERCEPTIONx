#' Load pre-built model of provided drugs
#'
#' Downloads model files from GitHub Release if not cached locally.
#' Default destination is present working directory.
#' And (default cannot) can automatically load the downloaded model(s) after download completion.
#' @param ... One or more drug names (e.g., "erlotinib", "gefitinib").
#' @param dest Directory to downloaded models. Default = ".".
#' @param read Whether read the downloaded model or not.
#'
#' @return Invisibly returns a list of model objects. Creates individual variables in global environment.
#' @export
load_model <- function(..., dest = ".", read = FALSE) {
  drugs <- tolower(c(...))
  if (length(drugs) == 0) stop("Drug list is empty.")

  available_drugs <- c(
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

  invalid_drugs <- drugs[!drugs %in% available_drugs]
  if (length(invalid_drugs) > 0) {
    stop("Some input drugs do not exist: ", paste(invalid_drugs, collapse = ", "))
  }

  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE)
    cat(paste0("'", dest, "' does not exist, created just now."))
  }

  base_url <- "https://github.com/SunPast/PERCEPTION/releases/download/models-v1/"
  result <- list()

  for (drug in drugs) {
    file_path <- file.path(dest, paste0(drug, ".RDS"))

    if (!file.exists(file_path)) {
      url <- paste0(base_url, drug, ".RDS")
      message("Downloading model for: ", drug)
      download.file(url, file_path, mode = "wb")
      cat(paste0("Successfully downloaded: '", drug, ".RDS'!\n"))
    } else {
      cat(paste0("Found cached model at: ", "'", file_path, "'\n"))
    }

    if (read) {
      result[[drug]] <- readRDS(file_path)
      cat(paste0("Successfully loaded ", "'", drug, ".RDS' as ", "'", drug, "'!\n"))
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
#' Downloads the required DepMap RDS file for training models from Zenodo,
#' which contains: bulk expression, scRNA expression, drug response annotations, etc.
#' Default destination is present working directory.
#' And (default cannot) can automatically load it after download completion, variable named `DepMap`.
#'
#' @param dest The path for storing the DepMap RDS. Default = ".".
#' @param read Whether read the downloaded DepMap or not.
#'
#' @return DepMap RDS file.
#' @export
load_depmap <- function(dest = ".", read = FALSE) {
  url <- "https://zenodo.org/record/7860559/files/DepMapv12.RDS"
  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE)
  }
  destfile <- file.path(dest, "DepMap.RDS")
  if (!file.exists(destfile)) {
    message("Downloading 883 MB file. This may take several minutes...")
    download.file(url, destfile, mode = "wb")
    cat("Successfully download 'DepMap.RDS'!")
  } else {
    message("The file 'DepMap.RDS' already exists.")
  }
  if (read == TRUE){
    cat("Reading DepMap.RDS...")
    DepMap <- readRDS(file.path(dest, "DepMap.RDS"))
    cat("Successfully read 'DepMap.RDS' as variable 'DepMap'!")
    list2env(DepMap, envir = .GlobalEnv)
    return(invisible(DepMap))
  }
}
