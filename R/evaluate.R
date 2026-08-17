#' PERCEPTIONx Model Evaluation Functions
#'
#' Functions for evaluating and comparing performance of trained drug response models.
#'
#' @name evaluate_perception
#' @keywords internal
#' @importFrom stats reformulate aggregate median
NULL


# Normalize a performance object to a c(pvalue, correlation) numeric vector.
#
# Accepts every shape produced across the package:
#   * named vector c(p.value = ..., estimate.cor = ...)  (train_models())
#   * data.frame with columns estimate.cor / p.value     (Shiny demo models)
#   * positional c(pvalue, correlation)                  (fallback)
# Returns c(pvalue = NA, correlation = NA) when nothing usable is found.
perf_to_vector <- function(perf) {
  if (is.null(perf)) return(c(pvalue = NA, correlation = NA))

  if (is.data.frame(perf)) {
    cor <- if ("estimate.cor" %in% names(perf)) perf$estimate.cor[1] else
           if ("correlation" %in% names(perf))  perf$correlation[1]  else NA
    p   <- if ("p.value" %in% names(perf)) perf$p.value[1] else NA
    return(c(pvalue = unname(p), correlation = unname(cor)))
  }

  if (is.numeric(perf) && !is.null(names(perf))) {
    cor <- if ("estimate.cor" %in% names(perf)) perf["estimate.cor"] else
           if ("correlation" %in% names(perf))  perf["correlation"]  else NA
    p   <- if ("p.value" %in% names(perf)) perf["p.value"] else NA
    return(c(pvalue = unname(p), correlation = unname(cor)))
  }

  if (is.numeric(perf) && length(perf) >= 2) {
    return(c(pvalue = as.numeric(perf[1]), correlation = as.numeric(perf[2])))
  }

  c(pvalue = NA, correlation = NA)
}


#' Compare performance of multiple trained models
#'
#' Computes and summarizes performance metrics across multiple drug models.
#' Returns performance in cross-validation, bulk test, pseudo-bulk, and scRNA datasets.
#'
#' @param model_list A named list of trained model objects from train_models().
#' @param threshold Numeric. Minimum correlation threshold for "passing" models. Default = 0.3.
#' @param verbose Logical. Whether to print summary statistics. Default = TRUE.
#'
#' @return A list containing:
#'   \describe{
#'     \item{perf_cv}{Cross-validation performance (|correlation| = sqrt(R²); p-value not available) for each drug}
#'     \item{perf_bulk}{Bulk test set performance for each drug}
#'     \item{perf_pseudo_bulk}{Pseudo-bulk performance for each drug}
#'     \item{perf_scRNA}{Single-cell RNA performance for each drug}
#'     \item{summary}{Summary statistics of models passing threshold}
#'   }
#'
#' @examples
#' \dontrun{
#'   models <- train_models(drug_list = c("abemaciclib", "erlotinib"))
#'   perf <- compare_performance(models)
#' }
#'
#' @export
compare_performance <- function(model_list, threshold = 0.3, verbose = TRUE) {

  # Validate input
  if (!is.list(model_list) || length(model_list) == 0) {
    stop("model_list must be a non-empty list of trained models.")
  }

  # Extract performance metrics from each model
  # R package uses named elements, not positional indices like original code

  perf_cv <- do.call(rbind, lapply(model_list, function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) return(c(pvalue = NA, correlation = NA))
    # CV stores sqrt(R²) = |correlation| only (no p-value available).
    c(pvalue = NA, correlation = as.numeric(x$model_performance_during_cv)[1])
  }))

  perf_bulk <- do.call(rbind, lapply(model_list, function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) return(c(pvalue = NA, correlation = NA))
    perf_to_vector(x$performance_in_bulk)
  }))

  perf_pseudo_bulk <- do.call(rbind, lapply(model_list, function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) return(c(pvalue = NA, correlation = NA))
    perf_to_vector(x$performance_in_pseudo_bulk)
  }))

  perf_scRNA <- do.call(rbind, lapply(model_list, function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) return(c(pvalue = NA, correlation = NA))
    perf_to_vector(x$performance_in_scRNA)
  }))

  # Set column names and row names
  colnames(perf_cv) <- c("pvalue", "correlation")
  colnames(perf_bulk) <- c("pvalue", "correlation")
  colnames(perf_pseudo_bulk) <- c("pvalue", "correlation")
  colnames(perf_scRNA) <- c("pvalue", "correlation")
  rownames(perf_cv) <- names(model_list)
  rownames(perf_bulk) <- names(model_list)
  rownames(perf_pseudo_bulk) <- names(model_list)
  rownames(perf_scRNA) <- names(model_list)

  # Compile results
  performance_results <- list(
    perf_cv = perf_cv,
    perf_bulk = perf_bulk,
    perf_pseudo_bulk = perf_pseudo_bulk,
    perf_scRNA = perf_scRNA
  )

  # Print summary statistics if verbose
  if (verbose) {
    message("\n=== Model Performance Summary ===")
    message("Threshold: correlation >", threshold)

    # Count models passing threshold in each dataset
    passing_counts <- sapply(performance_results, function(x) {
      sum(x[, "correlation"] > threshold, na.rm = TRUE)
    })
    message("\nModels passing threshold:")
    message("  CV (|cor|):    ", passing_counts["perf_cv"], " / ", length(model_list))
    message("  Bulk test:     ", passing_counts["perf_bulk"], " / ", length(model_list))
    message("  Pseudo-bulk:   ", passing_counts["perf_pseudo_bulk"], " / ", length(model_list))
    message("  scRNA:         ", passing_counts["perf_scRNA"], " / ", length(model_list))

    # Mean correlations
    mean_corrs <- sapply(performance_results, function(x) {
      mean(x[, "correlation"], na.rm = TRUE)
    })
    message("\nMean correlations:")
    message("  CV (|cor|):    ", round(mean_corrs["perf_cv"], 3))
    message("  Bulk test:     ", round(mean_corrs["perf_bulk"], 3))
    message("  Pseudo-bulk:   ", round(mean_corrs["perf_pseudo_bulk"], 3))
    message("  scRNA:         ", round(mean_corrs["perf_scRNA"], 3))

    # Drugs passing threshold in scRNA (most important metric)
    passing_drugs <- names(model_list)[which(perf_scRNA[, "correlation"] > threshold)]
    if (length(passing_drugs) > 0) {
      message("\nDrugs passing threshold in scRNA:")
      message("  ", paste(passing_drugs, collapse = ", "))
    }
  }

  # Return summary info
  summary_info <- list(
    threshold = threshold,
    passing_counts = passing_counts,
    mean_correlations = mean_corrs,
    passing_drugs_scRNA = names(model_list)[which(perf_scRNA[, "correlation"] > threshold)]
  )

  performance_results$summary <- summary_info

  return(performance_results)
}


#' Load performance metrics from a saved model file
#'
#' Reads only the performance portion of a saved model RDS file.
#' Useful when evaluating many models without loading full model objects (memory efficient).
#'
#' @param filepath Character. Path to the saved model RDS file.
#'                 Can be a single file path or a directory with model files.
#' @param drug_names Character vector. Optional. Specific drug names to load.
#'                   If NULL, loads all drugs in the file.
#'
#' @return A list containing performance metrics for each drug model.
#'
#' @examples
#' \dontrun{
#'   # Load performance from a single model file
#'   perf <- get_performance("models/PERCEPTIONx_models_PanCan_exPanCan_20240101_120000.RDS")
#'
#'   # Or use compare_performance on the result
#'   models <- readRDS("models/PERCEPTIONx_models_PanCan_exPanCan_20240101_120000.RDS")
#'   perf <- compare_performance(models)
#' }
#'
#' @export
get_performance <- function(filepath, drug_names = NULL) {

  # Validate filepath
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }

  # Load the model file
  model_data <- readRDS(filepath)

  # Check if it's a list of models
  if (!is.list(model_data)) {
    stop("Loaded data is not a list of models.")
  }

  # Filter by drug names if specified
  if (!is.null(drug_names)) {
    available_drugs <- intersect(drug_names, names(model_data))
    if (length(available_drugs) == 0) {
      stop("None of the specified drug_names found in the model file.")
    }
    model_data <- model_data[available_drugs]
  }

  # Extract only performance metrics (memory efficient)
  performance_only <- lapply(model_data, function(x) {
    if (is.null(x) || length(x) == 1 && is.na(x)) {
      return(list(
        model_performance_during_cv = c(NA, NA),
        performance_in_bulk = c(NA, NA),
        performance_in_pseudo_bulk = c(NA, NA),
        performance_in_scRNA = c(NA, NA)
      ))
    }
    list(
      model_performance_during_cv = x$model_performance_during_cv,
      performance_in_bulk = x$performance_in_bulk,
      performance_in_pseudo_bulk = x$performance_in_pseudo_bulk,
      performance_in_scRNA = x$performance_in_scRNA
    )
  })

  return(performance_only)
}


#' Get best performing models
#'
#' Filters model list to return only models that meet performance criteria.
#' Useful for selecting significant drug models for downstream analysis.
#'
#' @param model_list A named list of trained model objects.
#' @param min_correlation Numeric. Minimum correlation threshold. Default = 0.3.
#' @param max_pvalue Numeric. Maximum p-value threshold. Default = 0.05.
#' @param dataset Character. Which dataset to use for filtering:
#'                "scRNA" (default), "bulk", "pseudo_bulk", or "cv".
#'
#' @return A filtered list containing only models meeting the criteria.
#'
#' @examples
#' \dontrun{
#'   models <- train_models()
#'   significant_models <- get_significant_models(models, min_correlation = 0.3, max_pvalue = 0.05)
#' }
#'
#' @export
get_significant_models <- function(model_list,
                                   min_correlation = 0.3,
                                   max_pvalue = 0.05,
                                   dataset = "scRNA") {

  # Validate dataset parameter
  valid_datasets <- c("scRNA", "bulk", "pseudo_bulk", "cv")
  if (!dataset %in% valid_datasets) {
    stop("dataset must be one of: ", paste(valid_datasets, collapse = ", "))
  }

  # Get the appropriate performance element name
  perf_element <- switch(dataset,
                         "scRNA" = "performance_in_scRNA",
                         "bulk" = "performance_in_bulk",
                         "pseudo_bulk" = "performance_in_pseudo_bulk",
                         "cv" = "model_performance_during_cv")

  # Filter models. "cv" stores sqrt(R²)=|correlation| with no p-value; the other
  # datasets carry a named vector or data.frame with (p-value, correlation).
  significant_models <- model_list[vapply(model_list, function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) return(FALSE)

    if (dataset == "cv") {
      cv <- x$model_performance_during_cv
      return(is.numeric(cv) && length(cv) >= 1 && isTRUE(cv[1] > min_correlation))
    }

    perf <- perf_to_vector(x[[perf_element]])
    isTRUE(perf[["correlation"]] > min_correlation) &&
      isTRUE(perf[["pvalue"]] < max_pvalue)
  }, logical(1))]

  message("Found ", length(significant_models), " models meeting criteria (cor > ",
          min_correlation, ", p < ", max_pvalue, " in ", dataset, ")")

  return(significant_models)
}


#' Compute pseudo-bulk expression for a patient
#'
#' Calculates weighted average of clone-level expression to produce a
#' pseudo-bulk expression profile for a single patient. This is used
#' for evaluating model performance on pseudo-bulk data, not for prediction.
#'
#' @param x Integer. Index of the patient. Default = 1.
#' @param comb_viability_df Data frame. Retained for backward compatibility; no longer used (clone ids are read from the expression column names).
#' @param Clone_Counts_per_patients Data frame. Patient-clone abundance data.
#' @param clone_Level_z_expression_df Matrix. Clone-level expression data
#'        with columns named by patient-clone identifiers.
#'
#' @return Named numeric vector of pseudo-bulk expression values.
#'
#' @export
each_patient_pseudo_bulk <- function(x = 1,
                                     comb_viability_df,
                                     Clone_Counts_per_patients,
                                     clone_Level_z_expression_df) {

  patient <- Clone_Counts_per_patients$patients[x]

  # Expression columns are named "Patient@@Clone". Select this patient's clones
  # and keep clone IDs in the SAME order as the columns, so weights align with
  # columns (the previous code ordered weights by comb_viability_df instead).
  # `comb_viability_df` is retained for backward compatibility but not needed here.
  parsed <- parse_clone_keys(colnames(clone_Level_z_expression_df))
  cols <- which(parsed$patient == patient)
  col_clone_ids <- parsed$clone_id[cols]

  total_cells <- sum(unlist(Clone_Counts_per_patients[x, col_clone_ids]))
  clone_weights <- unlist(Clone_Counts_per_patients[x, col_clone_ids] / total_cells)

  sub <- clone_Level_z_expression_df[, cols, drop = FALSE]

  if (ncol(sub) > 1) {
    pseudo_bulk <- drop(sub %*% clone_weights)
  } else {
    pseudo_bulk <- sub[, 1]
  }
  pseudo_bulk
}