#' PERCEPTION Prediction Functions
#'
#' Functions for predicting drug response at cell/clone level and patient level.
#'
#' @name predict_perception
#' @keywords internal
#' @importFrom stats predict setNames
NULL


#' Predict drug response for cells or clones
#'
#' Given a trained model (or list of models) and a rank-normalized expression matrix,
#' predicts viability scores for each cell/sample across one or more drugs.
#' This function merges the former viability_from_model (single drug) and
#' killing_in_each_dataset (multi-drug) into a unified interface.
#'
#' @param model_list A named list of model objects (each with a \code{$model} element),
#'        or a single model object. From \code{train_models()} or \code{load_model()}.
#' @param expr Matrix or data frame. Rank-normalized expression matrix with
#'        genes as rows and cells/samples as columns.
#'
#' @return A matrix with cells/samples as rows and drugs as columns,
#'         containing predicted viability scores. Lower values indicate higher
#'         drug sensitivity.
#'
#' @examples
#' \dontrun{
#'   # Single drug
#'   models <- load_model("erlotinib", read = TRUE)
#'   pred <- predict_drugs(models, expr_rnorm)
#'
#'   # Multiple drugs
#'   models <- train_models(drug_list = c("abemaciclib", "erlotinib"),
#'                          cancer_type = "PanCan", exclude_cancer = "PanCan", GOI = GOI)
#'   pred <- predict_drugs(models, expr_rnorm)
#' }
#'
#' @export
predict_drugs <- function(model_list, expr) {

  # Normalize input: if single model object (not a list of models), wrap it
  if (!is.list(model_list) || is.null(names(model_list)) ||
      (!is.null(model_list$model) && !is.list(model_list[[1]]))) {
    # Single model object: wrap into a named list
    if (!is.null(model_list$model)) {
      drug_name <- attr(model_list, "drug_name")
      if (is.null(drug_name)) drug_name <- "drug1"
      model_list <- setNames(list(model_list), drug_name)
    } else {
      stop("model_list must be a named list of model objects (each with $model element), ",
           "or a single model object.")
    }
  }

  # Validate expression matrix
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("expr must be a matrix or data frame with genes as rows and samples as columns.")
  }

  # Predict for each drug
  viab_raw <- lapply(names(model_list), function(drug_name) {
    model_obj <- model_list[[drug_name]]

    # Skip NULL/NA entries
    if (is.null(model_obj) || (length(model_obj) == 1 && is.na(model_obj))) {
      warning("Skipping NULL model for: ", drug_name)
      return(NULL)
    }

    # Extract the caret model
    if (is.null(model_obj$model)) {
      warning("No $model element found for: ", drug_name)
      return(NULL)
    }

    viability_from_model_internal(drug_name, model_obj$model, expr)
  })

  names(viab_raw) <- names(model_list)

  # Remove NULL entries
  viab_raw <- viab_raw[!sapply(viab_raw, is.null)]

  if (length(viab_raw) == 0) {
    stop("All model predictions failed.")
  }

  do.call(cbind, viab_raw)
}


#' Internal: predict viability for a single drug model
#'
#' Extracts features from the model, matches them to expression matrix rows,
#' and predicts viability scores. Handles gene name format differences
#' (e.g., hyphens vs dots) via make.names() fallback.
#'
#' @param drug_name Character. Drug name (for error messages).
#' @param model A caret model object (from build_on_BULK_v2 output $model).
#' @param dataset Matrix. Expression matrix with genes as rows and
#'        cells/samples as columns. Must be rank-normalized.
#'
#' @return A named numeric vector of predicted viability scores.
#'
#' @keywords internal
viability_from_model_internal <- function(drug_name, model, dataset) {

  # Step 1: Try exact match of model features to expression matrix rows
  feature_match <- match(model$coefnames, rownames(dataset))

  # Step 2: If some features not found, try make.names() conversion
  # (DepMap gene names may use hyphens while external data uses dots after make.names)
  if (any(is.na(feature_match))) {
    dataset_rownames_made <- make.names(rownames(dataset), unique = TRUE)
    coefnames_made <- make.names(model$coefnames, unique = TRUE)
    feature_match_alt <- match(coefnames_made, dataset_rownames_made)

    # Use alternative match for previously unmatched features
    na_idx <- which(is.na(feature_match))
    for (i in na_idx) {
      if (!is.na(feature_match_alt[i])) {
        feature_match[i] <- feature_match_alt[i]
      }
    }
  }

  missing_features <- model$coefnames[is.na(feature_match)]

  if (length(missing_features) > 0) {
    stop(sprintf(
      "Model features not found in expression matrix for drug '%s': %s\nTotal missing: %d / %d features.\nThis usually means gene name formats differ between training and prediction data.",
      drug_name,
      paste(head(missing_features, 5), collapse = ", "),
      length(missing_features),
      length(model$coefnames)))
  }

  # Subset expression matrix to model features
  dataset_FOI <- data.frame(
    dataset[feature_match, ],
    row.names = rownames(dataset)[feature_match]
  )
  dataset_FOI <- data.frame(t(dataset_FOI))
  predict(model, dataset_FOI)
}


#' Predict drug response at patient level
#'
#' Aggregates clone-level drug response predictions to patient-level using
#' various weighting strategies based on clone abundance. This function merges
#' the former each_patient_killing (single drug) and each_patient_killingv2
#' (multi-drug) into a unified interface.
#'
#' @param clone_killing_matrix Data frame. Must have columns 'patient' and
#'        'clone_id', plus one or more drug killing/viability columns.
#'        Typically produced by merging \code{predict_drugs()} output with
#'        patient-clone metadata.
#' @param clone_counts Data frame. Rows are patients, first column must be
#'        'patients' with patient IDs, remaining columns are clone IDs with
#'        cell counts as values.
#' @param mode Character. Aggregation method:
#'   \describe{
#'     \item{"weighted_average"}{Weighted average of clone killing by clone abundance (default)}
#'     \item{"min"}{Minimum killing across clones (most resistant clone)}
#'     \item{"max"}{Maximum killing across clones (most sensitive clone)}
#'     \item{"weighted_max"}{Maximum of weighted killing across clones}
#'   }
#'
#' @return A data frame with patients as rows and drugs as columns,
#'         containing aggregated killing scores.
#'
#' @examples
#' \dontrun{
#'   # Step 1: Predict at clone level
#'   clone_pred <- predict_drugs(models, clone_expr_rnorm)
#'
#'   # Step 2: Add patient/clone metadata
#'   clone_killing_df <- data.frame(
#'     patient = patient_ids,
#'     clone_id = clone_ids,
#'     clone_pred
#'   )
#'
#'   # Step 3: Aggregate to patient level
#'   patient_pred <- predict_patients(clone_killing_df, clone_counts,
#'                                   mode = "weighted_average")
#' }
#'
#' @export
predict_patients <- function(clone_killing_matrix, clone_counts,
                             mode = "weighted_average") {

  # Validate inputs
  if (!"patient" %in% colnames(clone_killing_matrix)) {
    stop("clone_killing_matrix must have a 'patient' column.")
  }
  if (!"clone_id" %in% colnames(clone_killing_matrix)) {
    stop("clone_killing_matrix must have a 'clone_id' column.")
  }
  if (!"patients" %in% colnames(clone_counts)) {
    stop("clone_counts must have a 'patients' column.")
  }

  valid_modes <- c("weighted_average", "min", "max", "weighted_max")
  if (!mode %in% valid_modes) {
    stop("mode must be one of: ", paste(valid_modes, collapse = ", "))
  }

  # Identify drug columns (all columns except patient and clone_id)
  drug_cols <- setdiff(colnames(clone_killing_matrix), c("patient", "clone_id"))
  n_drugs <- length(drug_cols)

  # Aggregate for each patient
  patient_results <- do.call(rbind, lapply(1:nrow(clone_counts), function(x) {
    patient_id <- clone_counts$patients[x]

    # Get clone data for this patient
    patient_clones <- clone_killing_matrix[
      clone_killing_matrix$patient == patient_id, ]

    if (nrow(patient_clones) == 0) {
      # No clone data for this patient
      return(setNames(rep(NA, n_drugs), drug_cols))
    }

    # Compute clone weights
    clone_ids <- patient_clones$clone_id
    total_cells <- sum(clone_counts[x, clone_ids])
    if (total_cells == 0) {
      return(setNames(rep(NA, n_drugs), drug_cols))
    }
    clone_weights <- unlist(clone_counts[x, clone_ids] / total_cells)

    # Aggregate based on mode
    if (n_drugs == 1) {
      # Single drug: use comb_killing column if present, otherwise use drug column
      if ("comb_killing" %in% colnames(patient_clones)) {
        killing_values <- patient_clones$comb_killing
      } else {
        killing_values <- patient_clones[[drug_cols]]
      }

      result <- switch(mode,
        "weighted_average" = sum(killing_values * clone_weights),
        "min" = min(killing_values),
        "max" = max(killing_values),
        "weighted_max" = max(killing_values * clone_weights)
      )
      return(setNames(result, drug_cols))
    } else {
      # Multiple drugs: apply aggregation across all drug columns
      drug_matrix <- as.matrix(patient_clones[, drug_cols, drop = FALSE])

      result <- switch(mode,
        "weighted_average" = apply(drug_matrix, 2, function(col) sum(col * clone_weights)),
        "min" = apply(drug_matrix, 2, function(col) min(col)),
        "max" = apply(drug_matrix, 2, function(col) max(col)),
        "weighted_max" = apply(drug_matrix, 2, function(col) max(col * clone_weights))
      )
      return(result)
    }
  }))

  # Add patient names as row names
  rownames(patient_results) <- clone_counts$patients

  # Convert to data frame
  patient_results <- as.data.frame(patient_results)

  patient_results
}

