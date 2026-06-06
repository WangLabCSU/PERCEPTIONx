#' Predict viability from a trained PERCEPTION model
#'
#' Given a trained model and an expression matrix, this function extracts the features used by the model (via coefnames),
#' matches them to the expression matrix rows, and predicts viability scores.
#' This is the robust version from step0B that uses model$coefnames for feature matching.
#'
#' @param infunc_DrugName Character. Drug name (used for naming/reference only).
#' @param infunc_model A caret model object (from build_on_BULK_v2 output $model).
#' @param infunc_dataset Matrix. Expression matrix with genes as rows and
#'        cells/samples as columns. Must be rank-normalized.
#'
#' @return A named numeric vector of predicted viability scores, one per column
#'         in infunc_dataset.
#'
#' @examples
#' \dontrun{
#'   load_depmap(read = TRUE)
#'   models <- train_perception_models("abemaciclib", "PanCan", "PanCan", GOI)
#'   viab <- viability_from_model(
#'     infunc_DrugName = "abemaciclib",
#'     infunc_model = models$abemaciclib$model,
#'     infunc_dataset = DepMap$CPM_scRNA_CCLE_rnorm
#'   )
#' }
#'
#' @export
viability_from_model <- function(infunc_DrugName,
                                 infunc_model,
                                 infunc_dataset) {

  # Step 1: Try exact match of model features to expression matrix rows
  feature_match <- match(infunc_model$coefnames, rownames(infunc_dataset))
  
  # Step 2: If some features not found, try make.names() conversion
  # (DepMap gene names may use hyphens while external data uses dots after make.names)
  if (any(is.na(feature_match))) {
    dataset_rownames_made <- make.names(rownames(infunc_dataset), unique = TRUE)
    coefnames_made <- make.names(infunc_model$coefnames, unique = TRUE)
    feature_match_alt <- match(coefnames_made, dataset_rownames_made)
    
    # Use alternative match for previously unmatched features
    na_idx <- which(is.na(feature_match))
    for (i in na_idx) {
      if (!is.na(feature_match_alt[i])) {
        feature_match[i] <- feature_match_alt[i]
      }
    }
  }
  
  missing_features <- infunc_model$coefnames[is.na(feature_match)]
  
  if (length(missing_features) > 0) {
    stop(sprintf(
      "Model features not found in expression matrix for drug '%s': %s\nTotal missing: %d / %d features.\nThis usually means gene name formats differ between training and prediction data.",
      infunc_DrugName,
      paste(head(missing_features, 5), collapse = ", "),
      length(missing_features),
      length(infunc_model$coefnames)))
  }

  # Subset expression matrix to model features
  infunc_dataset_FOI <- data.frame(
    infunc_dataset[feature_match, ],
    row.names = rownames(infunc_dataset)[feature_match]
  )
  infunc_dataset_FOI <- data.frame(t(infunc_dataset_FOI))
  Viability_score <- predict(infunc_model, infunc_dataset_FOI)
  Viability_score
}


#' Predict killing for all drugs across a dataset
#'
#' Applies viability_from_model for each drug in a model list and returns
#' a combined matrix of viability predictions.
#'
#' @param infunc_scRNAseq_dataset_rnorm Matrix. Rank-normalized expression matrix
#'        (genes as rows, cells/samples as columns).
#' @param infunc_GOI Character vector. Genes of Interest (currently unused,
#'        kept for compatibility).
#' @param infunc_model_list Named list of model objects (each with $model element).
#'
#' @return A matrix with viability predictions for each drug (columns) and
#'         each cell/sample (rows).
#'
#' @examples
#' \dontrun{
#'   models <- train_perception_models(c("abemaciclib", "erlotinib"),
#'                                     "PanCan", "PanCan", GOI)
#'   killing_mat <- killing_in_each_dataset(
#'     infunc_scRNAseq_dataset_rnorm = DepMap$scRNA_subset_rnorm,
#'     infunc_model_list = models
#'   )
#' }
#'
#' @export
killing_in_each_dataset <- function(infunc_scRNAseq_dataset_rnorm,
                                    infunc_GOI = NULL,
                                    infunc_model_list) {

  viab_raw <- lapply(1:length(infunc_model_list), function(x)
    viability_from_model(infunc_DrugName = names(infunc_model_list[x]),
                         infunc_model = infunc_model_list[[x]]$model,
                         infunc_dataset = infunc_scRNAseq_dataset_rnorm))
  names(viab_raw) <- names(infunc_model_list)
  do.call(cbind, viab_raw)
}


#' Compute per-patient killing from clone-level predictions
#'
#' Aggregates clone-level killing predictions to patient-level using various
#' weighting strategies. Each clone's contribution is weighted by its abundance
#' in the patient.
#'
#' @param x Integer. Index of the patient in Clone_Counts_per_patients. Default = 1.
#' @param mode Character. Aggregation method:
#'   \describe{
#'     \item{"weighted_average"}{Weighted average of clone killing by clone abundance (default)}
#'     \item{"min"}{Minimum killing across clones}
#'     \item{"max"}{Maximum killing across clones}
#'     \item{"weighted_max"}{Maximum of weighted killing across clones}
#'   }
#' @param clone_killing_matrix Data frame. Must have columns 'patient' and
#'        'clone_id', plus one or more drug killing columns. Typically
#'        the output of killing_in_each_dataset merged with clone metadata.
#' @param Clone_Counts_per_patients Data frame. Rows are patients, columns
#'        include clone IDs with cell counts as values, plus a 'patients' column.
#'
#' @return Numeric. The aggregated killing score for the patient.
#'
#' @export
each_patient_killing <- function(x = 1,
                                 mode = "weighted_average",
                                 clone_killing_matrix,
                                 Clone_Counts_per_patients) {

  total_cells <- sum(Clone_Counts_per_patients[x, clone_killing_matrix[
    clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], ]$clone_id])
  clone_weights <- unlist(Clone_Counts_per_patients[x, clone_killing_matrix[
    clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], ]$clone_id] / total_cells)

  if (mode == "weighted_average") {
    weighted_killing <- sum(clone_killing_matrix[
      clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], ]$comb_killing * clone_weights)
  } else if (mode == "min") {
    weighted_killing <- min(clone_killing_matrix[
      clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], ]$comb_killing)
  } else if (mode == "max") {
    weighted_killing <- max(clone_killing_matrix[
      clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], ]$comb_killing)
  } else if (mode == "weighted_max") {
    weighted_killing <- max(clone_killing_matrix[
      clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], ]$comb_killing * clone_weights)
  }
  weighted_killing
}


#' Compute per-patient killing from clone-level predictions (multi-drug version)
#'
#' Similar to each_patient_killing but returns a vector of killing scores
#' for multiple drugs (columns in clone_killing_matrix).
#'
#' @param x Integer. Index of the patient. Default = 1.
#' @param mode Character. Aggregation method: "weighted_average", "min", or "max".
#'        Default = "weighted_average".
#' @param clone_killing_matrix Data frame. Must have columns 'patient' and
#'        'clone_id', plus drug killing columns.
#' @param Clone_Counts_per_patients Data frame. Patient-clone abundance data.
#'
#' @return Named numeric vector of killing scores, one per drug column.
#'
#' @export
each_patient_killingv2 <- function(x = 1,
                                   mode = "weighted_average",
                                   clone_killing_matrix,
                                   Clone_Counts_per_patients) {

  total_cells <- sum(Clone_Counts_per_patients[x, clone_killing_matrix[
    clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], ]$clone_id])
  clone_weights <- unlist(Clone_Counts_per_patients[x, clone_killing_matrix[
    clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], ]$clone_id] / total_cells)

  if (mode == "weighted_average") {
    killing <- apply(clone_killing_matrix[
      clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], -(1:2)], 2,
      function(x) sum(x * clone_weights))
  } else if (mode == "min") {
    killing <- apply(clone_killing_matrix[
      clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], -(1:2)], 2,
      function(x) min(x))
  } else if (mode == "max") {
    killing <- apply(clone_killing_matrix[
      clone_killing_matrix$patient == Clone_Counts_per_patients$patients[x], -(1:2)], 2,
      function(x) max(x))
  }
  killing
}


#' Compute clone abundance weights for each patient
#'
#' Calculates the proportion of each clone's cell count relative to total
#' cells for a given patient.
#'
#' @param x Integer. Index of the patient. Default = 1.
#' @param mode Character. Aggregation mode (currently unused, kept for compatibility).
#'        Default = "weighted_average".
#' @param comb_killing_df Data frame. Must have columns 'patient' and 'clone_id'.
#' @param Clone_Counts_per_patients Data frame. Patient-clone abundance data.
#'
#' @return Named numeric vector of clone weights (c1, c2, c3 format).
#'
#' @export
each_patient_clone_weights <- function(x = 1,
                                       mode = "weighted_average",
                                       comb_killing_df,
                                       Clone_Counts_per_patients) {

  # Identify clone columns (all columns except 'patients')
  clone_cols <- setdiff(colnames(Clone_Counts_per_patients), "patients")

  # Get clone IDs for this patient from comb_killing_df
  patient_clones <- comb_killing_df[
    comb_killing_df$patient == Clone_Counts_per_patients$patients[x], ]$clone_id

  # Initialize all clone weights to 0
  clone_weights <- rep(0, length(clone_cols))
  names(clone_weights) <- clone_cols

  # Fill in weights for clones that exist for this patient
  if (length(patient_clones) > 0) {
    total_cells <- sum(Clone_Counts_per_patients[x, patient_clones])
    if (total_cells > 0) {
      existing_weights <- unlist(Clone_Counts_per_patients[x, patient_clones] / total_cells)
      clone_weights[match(names(existing_weights), names(clone_weights))] <- existing_weights
    }
  }

  clone_weights
}


#' Compute pseudo-bulk expression for a patient
#'
#' Calculates weighted average of clone-level expression to produce a
#' pseudo-bulk expression profile for a single patient.
#'
#' @param x Integer. Index of the patient. Default = 1.
#' @param comb_killing_df Data frame. Must have columns 'patient' and 'clone_id'.
#' @param Clone_Counts_per_patients Data frame. Patient-clone abundance data.
#' @param clone_Level_z_expression_df Matrix. Clone-level expression data
#'        with columns named by patient-clone identifiers.
#'
#' @return Named numeric vector of pseudo-bulk expression values.
#'
#' @export
each_patient_pseudo_bulk <- function(x = 1,
                                     comb_killing_df,
                                     Clone_Counts_per_patients,
                                     clone_Level_z_expression_df) {

  total_cells <- sum(Clone_Counts_per_patients[x, comb_killing_df[
    comb_killing_df$patient == Clone_Counts_per_patients$patients[x], ]$clone_id])
  clone_weights <- unlist(Clone_Counts_per_patients[x, comb_killing_df[
    comb_killing_df$patient == Clone_Counts_per_patients$patients[x], ]$clone_id] / total_cells)
  clone_expression <- data.frame(
    clone_Level_z_expression_df[, grep(Clone_Counts_per_patients$patients[x],
                                       colnames(clone_Level_z_expression_df))])

  if (ncol(clone_expression) > 1) {
    pseudo_bulk <- rowMeans(clone_Level_z_expression_df[
      , grep(Clone_Counts_per_patients$patients[x],
             colnames(clone_Level_z_expression_df))] * clone_weights)
  } else {
    pseudo_bulk <- clone_Level_z_expression_df[
      , grep(Clone_Counts_per_patients$patients[x],
             colnames(clone_Level_z_expression_df))]
  }
  pseudo_bulk
}
