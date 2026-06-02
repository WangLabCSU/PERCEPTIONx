#' Load pre-built model of provided drugs
#'
#' There are 44 drugs and each of them corresponds to a pre-built model in the path data/model.
#' Users can input one or more drug names in lowercase to load drug models they want.
#'
#' @param ... Character vector. The drugs provided by users.
#' @return One or more RDS file(s).
#' @export
load_model <- function(...){
  drug <- c(...)
  exist_model <- list.files(system.file("extdata", "models", package="PERCEPTION"), pattern = "\\.RDS$")
  available_drugs <- gsub("\\.RDS$", "", exist_model)
  drug_amount <- length(drug)

  if (drug_amount == 0) stop("Drug list is empty.")

  if (!all(drug %in% available_drugs)) stop("Some input drugs do not exist.")
  result <- list()

  for (single_drug in drug) {
    result[[single_drug]] <- readRDS(paste0(system.file("extdata", "models", package="PERCEPTION"), "/", single_drug, ".RDS"))
    cat("Successfully loaded:", single_drug, "\n")
  }
  list2env(result, envir = .GlobalEnv)
  return(invisible(result))
}

#' Download DepMapv12.RDS
#'
#' Downloads the required DepMap RDS file for traing models from Zenodo,
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
