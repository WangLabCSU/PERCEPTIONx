# Models are fitted with caret::train() (class "train"), so predict() dispatches
# on predict.train. Importing the method loads caret's namespace and registers
# its S3 methods; without it, predicting from a *pre-trained* model (loaded via
# load_model(), i.e. caret never called in the session) fails with
# "no applicable method for 'predict' applied to an object of class 'train'"

#' PERCEPTIONx Prediction Functions
#'
#' Functions for predicting drug response at cell/clone level and patient level.
#'
#' @name predict_perception
#' @keywords internal
NULL


#' Predict drug response for cells or clones
#'
#' Given a trained model (or list of models) and a rank-normalized expression matrix,
#' predicts viability scores for each cell/sample across one or more drugs.
#' This function merges the former viability_from_model (single drug) and
#' viability_in_each_dataset (multi-drug) into a unified interface.
#'
#' @param model_list A named list of model objects (each with a \code{$model} element),
#'        or a single model object. From \code{train_models()} or \code{load_model()}.
#' @param expr Matrix, data frame, or data.table. Rank-normalized expression with
#'        genes as rows and cells/samples as columns. A data.table carries no row
#'        names, so its first column must hold the gene identifiers and the
#'        remaining columns the cells/samples.
#'
#' @return A matrix with cells/samples as rows and drugs as columns,
#'         containing predicted viability scores. Lower values indicate higher
#'         drug sensitivity.
#'
#' @details If a small fraction (<= 50\%) of a model's features are missing from
#'   \code{expr} (e.g. genes filtered out during scRNA QC), they are imputed
#'   with the neutral rank value 0.5 and a warning is issued. If more than half
#'   of the features are missing, prediction stops with an error.
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
  # Normalize input: wrap a single model object (one with a $model element)
  # into a named list so the rest of the function treats both cases uniformly.
  if (!is.null(model_list$model)) {
    drug_name <- attr(model_list, "drug_name")
    if (is.null(drug_name)) {
      drug_name <- "drug1"
    }
    model_list <- stats::setNames(list(model_list), drug_name)
  }

  if (!is.list(model_list) || is.null(names(model_list))) {
    cli::cli_abort(
      "model_list must be a named list of model objects (each with a {.field model} element), or a single model object."
    )
  }

  # Validate expression matrix
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    cli::cli_abort(
      "expr must be a matrix, data frame, or data.table with genes as rows and samples as columns."
    )
  }

  # Predict for each drug
  viab_raw <- lapply(names(model_list), function(drug_name) {
    model_obj <- model_list[[drug_name]]

    # Skip NULL/NA entries
    if (is.null(model_obj) || (length(model_obj) == 1L && is.na(model_obj))) {
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

  if (length(viab_raw) == 0L) {
    cli::cli_abort("All model predictions failed.")
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
  if (is.null(model$coefnames) || length(model$coefnames) == 0L) {
    cli::cli_abort(
      "Model for drug {.val {drug_name}} has no feature names ({.field coefnames} is NULL/empty). Cannot subset the expression matrix for prediction."
    )
  }

  # data.table has no row names, so feature identifiers must come from the first
  # column; matrices and data frames keep using rownames().
  if (data.table::is.data.table(dataset)) {
    gene_ids <- as.character(dataset[[1L]])
    dataset <- as.matrix(dataset[, -1L])
    rownames(dataset) <- gene_ids
  } else {
    dataset <- as.matrix(dataset)
  }

  # Step 1: Try exact match of model features to expression matrix rows
  feature_match <- match(model$coefnames, rownames(dataset))

  # Step 2: If some features not found, try make.names() conversion
  # (DepMap gene names may use hyphens while external data uses dots after make.names)
  if (anyNA(feature_match)) {
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

  if (length(missing_features) > 0L) {
    frac_missing <- length(missing_features) / length(model$coefnames)
    if (frac_missing > 0.5) {
      cli::cli_abort(
        c(
          "Model features not found in expression matrix for drug {.val {drug_name}}: {length(missing_features)} / {length(model$coefnames)} features missing ({round(100 * frac_missing)}%).",
          "i" = "This usually means gene name formats differ between training and prediction data, or the wrong expression matrix was provided."
        )
      )
    }
    # Impute the remaining features with the neutral rank-normalized value (0.5,
    # the median rank). This keeps predictions unbiased when a handful of genes
    # were filtered out during scRNA QC / gene-set coverage differences.
    warning(sprintf(
      "Imputing %d / %d missing model features for drug '%s' with the median rank value (0.5): %s",
      length(missing_features),
      length(model$coefnames),
      drug_name,
      toString(utils::head(missing_features, 10L))
    ))
    impute_mat <- matrix(
      0.5,
      nrow = length(missing_features),
      ncol = ncol(dataset),
      dimnames = list(missing_features, colnames(dataset))
    )
    dataset <- rbind(as.matrix(dataset), impute_mat)
    feature_match <- match(model$coefnames, rownames(dataset))
  }

  # Subset expression matrix to model features and transpose once into the
  # samples x features frame caret expects. Transposing the matrix directly keeps
  # the sample (row) names intact, unlike data.frame(t(data.frame(...))), which
  # runs make.names() and turns clone keys like "Kydar01@@c1" into "..".
  dataset_FOI <- as.data.frame(
    t(dataset[feature_match, , drop = FALSE]),
    row.names = colnames(dataset)
  )
  stats::predict(model, dataset_FOI)
}


#' Predict drug response at patient level
#'
#' Aggregates clone-level drug response predictions to patient-level using
#' various weighting strategies based on clone abundance. This function merges
#' the former each_patient_viability (single drug) and each_patient_viabilityv2
#' (multi-drug) into a unified interface.
#'
#' @param clone_pred Matrix from \code{predict_drugs()}, with clones as rows
#'        and drugs as columns. Row names should match clone column names from
#'        \code{prepare_data()}.
#' @param prepared_data List from \code{prepare_data()}, containing
#'        \code{$clone_viability_template} and \code{$clone_counts}. Alternatively,
#'        you can pass \code{clone_counts} directly as a data frame or data.table
#'        (legacy mode).
#' @param clone_counts Optional. Data frame or data.table with patients as rows
#'        and clone IDs as columns. Only needed if \code{prepared_data} is not a
#'        list from \code{prepare_data()}.
#' @param mode Character. Aggregation method:
#'   \describe{
#'     \item{"weighted_max"}{Maximum of weighted viability across clones. Default.}
#'     \item{"max"}{Maximum viability across clones (most resistant clone)}
#'     \item{"weighted_average"}{Weighted average of clone viability by clone abundance}
#'     \item{"min"}{Minimum viability across clones (most sensitive clone)}
#'     \item{"average"}{Average of clone viability}
#'   }
#' @param zscore Logical. Whether to z-score scale drug columns across patients
#'        before aggregation. Default = TRUE. Matches the original PERCEPTION pipeline.
#'        An already-standardized \code{comb_viability} column (the output of the
#'        IDA combination step) is never scaled again.
#'
#' @return A data frame with patients as rows and drugs as columns,
#'         containing aggregated viability scores.
#'
#' @examples
#' \dontrun{
#'   # Simple workflow: pass prepare_data() output directly
#'   prepared <- prepare_data(patient_scRNA)
#'   clone_pred <- predict_drugs(models, prepared$clone_expression)
#'   patient_pred <- predict_patients(clone_pred, prepared)
#'
#'   # Legacy workflow: manually build clone_viability_matrix
#'   clone_viability_df <- data.frame(
#'     patient = patient_ids,
#'     clone_id = clone_ids,
#'     clone_pred
#'   )
#'   patient_pred <- predict_patients(clone_viability_df, clone_counts)
#' }
#'
#' @export
predict_patients <- function(
  clone_pred,
  prepared_data,
  clone_counts = NULL,
  mode = "weighted_max",
  zscore = TRUE
) {
  # Determine input mode: simple (prepared_data is a list) vs legacy (prepared_data is clone_counts)
  # Accept both clone_viability_template and clone_viability_df_template for backward compatibility
  # Also accept any list with clone_counts + any template-like element
  template_names <- c(
    "clone_viability_template",
    "clone_viability_df_template",
    "template",
    "viability_template"
  )
  has_counts <- "clone_counts" %in% names(prepared_data)

  if (is.list(prepared_data) && !is.data.frame(prepared_data) && has_counts) {
    # Simple mode: extract from prepare_data() output
    template_name <- intersect(names(prepared_data), template_names)[1L]

    if (is.na(template_name)) {
      # No template found but has counts - try to build template from clone_pred rownames
      cli::cli_alert_info(
        "No template found, building from {.field clone_pred} rownames..."
      )
      clone_col_names <- rownames(clone_pred)
      if (length(clone_col_names) > 0L && all(grepl("@@", clone_col_names))) {
        clone_template <- parse_clone_keys(clone_col_names)
      } else {
        cli::cli_abort(
          c(
            "Cannot build template: clone_pred rownames do not follow the {.val Patient@@Clone} format.",
            "x" = "clone_pred rownames: {toString(head(rownames(clone_pred), 5))}"
          )
        )
      }
    } else {
      clone_template <- prepared_data[[template_name]]
    }

    clone_counts <- prepared_data$clone_counts

    # Merge clone_pred matrix with template to build clone_viability_matrix
    clone_viability_matrix <- data.table::as.data.table(clone_pred)
    clone_viability_matrix[, patient := clone_template$patient]
    clone_viability_matrix[, clone_id := clone_template$clone_id]
  } else if (is.data.frame(prepared_data)) {
    # Legacy mode: prepared_data is actually clone_counts, clone_pred is clone_viability_matrix
    # copy(): the `:=` edits below would otherwise modify a data.table input in
    # place, i.e. the caller's object.
    clone_viability_matrix <- data.table::copy(clone_pred)
    clone_counts <- prepared_data
  } else {
    cli::cli_abort(
      c(
        "prepared_data must be either:",
        "*" = "A list from prepare_data() (with clone_counts)",
        "*" = "A data frame of clone_counts (legacy mode, requires clone_pred to be clone_viability_matrix)",
        "i" = "class(prepared_data): {class(prepared_data)[1]}",
        "i" = "is.list: {is.list(prepared_data)}",
        "i" = "is.data.frame: {is.data.frame(prepared_data)}",
        "i" = "names(prepared_data): {paste(names(prepared_data), collapse = ', ')}",
        "i" = "has clone_counts: {has_counts}"
      )
    )
  }

  # Normalise tabular inputs to data.table. as.data.table() copies, so the
  # caller's objects are never touched by the in-place operations below.
  clone_viability_matrix <- data.table::as.data.table(clone_viability_matrix)
  clone_counts <- data.table::as.data.table(clone_counts)

  # Column names used through data.table's non-standard evaluation.
  patient <- clone_id <- clone_key <- abund_col <- n_cells <- NULL

  # Validate inputs
  if (!"patient" %chin% colnames(clone_viability_matrix)) {
    cli::cli_abort(
      "clone_viability_matrix must have a {.field patient} column."
    )
  } else if (!"clone_id" %chin% colnames(clone_viability_matrix)) {
    cli::cli_abort(
      "clone_viability_matrix must have a {.field clone_id} column."
    )
  } else if (!"patients" %chin% colnames(clone_counts)) {
    cli::cli_abort("clone_counts must have a {.field patients} column.")
  }

  valid_modes <- c("weighted_average", "min", "max", "weighted_max", "average")
  if (!mode %chin% valid_modes) {
    cli::cli_abort(
      "{.arg mode} must be one of: {.val {valid_modes}}"
    )
  }

  # Identify drug columns (all columns except patient and clone_id)
  drug_cols <- setdiff(
    colnames(clone_viability_matrix),
    c("patient", "clone_id")
  )
  n_drugs <- length(drug_cols)
  if (n_drugs == 0L) {
    cli::cli_abort(
      "No drug columns found. {.arg clone_viability_matrix} must contain columns besides {.field patient} and {.field clone_id}."
    )
  }

  # Z-score scale drug columns across patients before aggregation.
  # A pre-combined score (column name "comb_viability", produced by the IDA
  # step: per-drug z-score -> per-clone minimum) is ALREADY standardized, so
  # scaling it again would silently change the weighted maximum
  # (max(z_i * w_i) != max(x_i * w_i) for unequal weights) and therefore the
  # patient ranking. The original PERCEPTION pipeline only standardizes the raw
  # per-drug predictions, so pre-combined columns are left as-is.
  if (zscore) {
    scale_cols <- setdiff(drug_cols, "comb_viability")
    if (length(scale_cols) > 0L) {
      # Same transformation as zscore_viability(), applied in place.
      clone_viability_matrix[,
        (scale_cols) := lapply(.SD, function(col) {
          as.numeric(scale(as.numeric(col)))
        }),
        .SDcols = scale_cols
      ]
    }
  }

  # --- Clone weights in long form: one row per (patient, abundance column) ---
  patient_order <- as.character(clone_counts$patients)
  abundance <- data.table::melt(
    clone_counts,
    id.vars = "patients",
    variable.name = "abund_col",
    value.name = "n_cells",
    variable.factor = FALSE
  )
  data.table::setnames(abundance, "patients", "patient")
  abundance[, patient := as.character(patient)]

  # Abundance columns are keyed "patient@@clone" (see build_clone_counts); fall
  # back to the bare clone id so objects prepared by earlier versions keep
  # working.
  clone_viability_matrix[, patient := as.character(patient)]
  clone_viability_matrix[,
    clone_key := paste(patient, as.character(clone_id), sep = "@@")
  ]
  clone_viability_matrix[,
    abund_col := data.table::fifelse(
      clone_key %chin% names(clone_counts),
      clone_key,
      as.character(clone_id)
    )
  ]

  # Attach the cell count that weights each clone, then aggregate per patient.
  clone_viability_matrix <- abundance[
    clone_viability_matrix,
    on = c("patient", "abund_col")
  ]

  n_unmatched <- clone_viability_matrix[is.na(n_cells), .N]
  if (n_unmatched > 0L) {
    cli::cli_warn(
      "{n_unmatched} clone{?s} have no abundance column in {.arg clone_counts}; affected patient scores cannot be weighted and may be {.val NA}."
    )
  }

  aggregate_cols <- switch(
    mode,
    "weighted_average" = function(values, weights) sum(values * weights),
    "weighted_max" = function(values, weights) max(values * weights),
    "min" = function(values, weights) min(values),
    "max" = function(values, weights) max(values),
    "average" = function(values, weights) mean(values)
  )

  patient_results <- clone_viability_matrix[,
    {
      total_cells <- sum(n_cells, na.rm = TRUE)
      if (!is.finite(total_cells) || total_cells <= 0) {
        # No usable abundance for this patient
        as.list(stats::setNames(rep(NA_real_, n_drugs), drug_cols))
      } else {
        weights <- n_cells / total_cells
        as.list(stats::setNames(
          vapply(
            drug_cols,
            function(col) aggregate_cols(.SD[[col]], weights),
            numeric(1L)
          ),
          drug_cols
        ))
      }
    },
    by = "patient",
    .SDcols = drug_cols
  ]

  # Every patient listed in clone_counts gets a row, in its original order (NA
  # when no clone-level prediction exists for it).
  patient_results <- patient_results[
    match(patient_order, patient_results$patient)
  ]

  # Documented output shape: patients as row names, drugs as columns.
  patient_results <- as.data.frame(patient_results[, .SD, .SDcols = drug_cols])
  rownames(patient_results) <- patient_order

  patient_results
}


#' Combine multi-drug predictions into a single combination score
#'
#' Implements the IDA (Independent Drug Action) combination principle used in
#' PERCEPTION: each clone-level drug viability is first z-score scaled across
#' all clones, then the per-clone combined score is derived from the
#' drug-specific z-scores. The default rule \code{min} means "the most effective
#' drug in the combination drives the clone-level combined score" — the exact
#' rule hard-coded in the Shiny app's internal \code{combo_clone_frame()}.
#'
#' This function is the missing glue between [predict_drugs()] and downstream
#' patient-level aggregation / visualisation. Its output feeds directly into
#' [predict_patients()] (as the app's [predict_patients()] already skips
#' re-standardising the \code{comb_viability} column) and [plot_clone_viability()].
#'
#' @param pred Matrix, data.frame, or data.table of clone-level drug viability
#'   predictions. Typically the output of [predict_drugs()]. Rows = clones,
#'   columns = drugs. Row names must be canonical PERCEPTION clone keys
#'   (\code{"patient@@clone_id"}). Alternatively, pass a \code{data.frame} /
#'   \code{data.table} with an explicit \code{patient} column and a
#'   \code{clone_id} column.
#' @param patient Optional character vector of patient IDs, length equal to
#'   \code{nrow(pred)}. Required only when \code{pred} row names are NOT
#'   \code{"patient@@clone_id"} keys.
#' @param clone_id Optional character vector of clone IDs, length equal to
#'   \code{nrow(pred)}. Required only when passing \code{patient} explicitly.
#' @param combine_rule Character scalar, per-clone aggregation after z-scoring.
#'   Must be one of \code{"min"} (IDA default), \code{"mean"}, \code{"max"}, or
#'   \code{"median"}.
#' @param na.rm Logical, passed through to the aggregator. Default \code{TRUE}.
#' @param include_single_z Logical. If \code{TRUE}, per-drug z-scored columns
#'   are retained in the output alongside \code{comb_viability}. Default
#'   \code{FALSE}.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{patient}{Patient ID.}
#'     \item{clone_id}{Clone ID.}
#'     \item{comb_viability}{Combined viability score (already z-scored across
#'       clones, so [predict_patients()] will not standardise it again).}
#'     \item{<drug_z>}{(optional, when \code{include_single_z = TRUE}) one
#'       additional column per input drug, named like \code{<drug>.z},
#'       containing the per-drug z-scores before combination.}
#'   }
#'
#' @section Pipeline:
#' \preformatted{
#'   pred     <- predict_drugs(models, rank_expr)        # clones x drugs matrix
#'   clone_df <- combine_drugs(pred)                    # patient / clone_id / comb_viability
#'   pat      <- predict_patients(clone_df, prepared_data)
#'   plot_clone_viability(clone_df, viability_var = "comb_viability")
#' }
#'
#' The patient-level aggregated score will be identical to the one produced by
#' the Shiny app for the same drugs / clones.
#'
#' @section Caveats:
#' The default \code{min} rule amplifies noise when individual drug models are
#' very poor — the "least bad" model can still be quite bad, but \code{min}
#' picks it as the best combination driver. For exploratory analyses you can
#' pass \code{combine_rule = "mean"} or \code{combine_rule = "max"} and compare.
#'
#' When only one drug is present in \code{pred}, the combined score is simply
#' that drug's z-score.
#'
#' @seealso [predict_drugs()], [predict_patients()], [plot_clone_viability()],
#'   [zscore_viability()]
#' @export
combine_drugs <- function(
  pred,
  patient = NULL,
  clone_id = NULL,
  combine_rule = c("min", "mean", "max", "median"),
  na.rm = TRUE,
  include_single_z = FALSE
) {
  combine_rule <- rlang::arg_match(combine_rule)

  # --- 1. Normalise input to a numeric matrix with row names ----------------
  if (is.data.frame(pred)) {
    # If user passes a data.frame/data.table with patient + clone_id columns,
    # stash them and keep every remaining column as a drug.
    pred_rn <- rownames(pred)
    pred_dt <- data.table::as.data.table(pred)
    df_cols <- colnames(pred_dt)
    drop_cols <- intersect(df_cols, c("patient", "clone_id", "comb_viability"))
    drug_cols <- setdiff(df_cols, drop_cols)
    if (length(drug_cols) == 0L) {
      cli::cli_abort(
        "No drug columns found in pred. A data.frame input must contain drug columns in addition to {.field patient} / {.field clone_id}."
      )
    }
    pred_mat <- as.matrix(pred_dt[, .SD, .SDcols = drug_cols])
    # data.table has no row names, so restore the clone keys carried by `pred`.
    rownames(pred_mat) <- pred_rn
    if (is.null(patient) && "patient" %chin% df_cols) {
      patient <- as.character(pred_dt$patient)
    }
    if (is.null(clone_id) && "clone_id" %chin% df_cols) {
      clone_id <- as.character(pred_dt$clone_id)
    }
  } else if (!is.matrix(pred)) {
    cli::cli_abort(
      "pred must be a matrix, data.frame, or data.table::data.table (from {.fn predict_drugs})."
    )
  } else {
    pred_mat <- pred
  }
  pred_mat <- as.matrix(pred_mat)
  storage.mode(pred_mat) <- "double"

  if (ncol(pred_mat) < 1L) {
    cli::cli_abort("pred has zero drug columns.")
  }

  # --- 2. Resolve patient / clone_id ---------------------------------------
  if (is.null(patient) || is.null(clone_id)) {
    rn <- rownames(pred_mat)
    if (is.null(rn) || !all(grepl("@@", rn))) {
      cli::cli_abort(
        c(
          "pred row names must be canonical PERCEPTION clone keys ({.val patient@@clone_id}),",
          "i" = "or supply both {.arg patient} and {.arg clone_id}."
        )
      )
    }
    keys <- parse_clone_keys(rn)
    if (is.null(patient)) {
      patient <- keys$patient
    }
    if (is.null(clone_id)) clone_id <- keys$clone_id
  }
  if (length(patient) != nrow(pred_mat) || length(clone_id) != nrow(pred_mat)) {
    cli::cli_abort(
      "{.arg patient} and {.arg clone_id} lengths ({length(patient)}/{length(clone_id)}) must equal {.code nrow(pred)} ({nrow(pred_mat)})."
    )
  }

  # --- 3. Z-score each drug column across all clones ------------------------
  z_mat <- pred_mat
  for (j in seq_len(ncol(pred_mat))) {
    col <- pred_mat[, j]
    s <- stats::sd(col, na.rm = TRUE)
    z_mat[, j] <- if (is.na(s) || s == 0L) {
      rep(0L, length(col))
    } else {
      as.numeric(scale(col))
    }
  }

  # --- 4. Per-clone combination ---------------------------------------------
  comb <- switch(
    combine_rule,
    "min" = if (rlang::is_installed("matrixStats")) {
      matrixStats::rowMins(z_mat, na.rm = na.rm)
    } else {
      apply(z_mat, 1L, min, na.rm = na.rm)
    },
    "mean" = rowMeans(z_mat, na.rm = na.rm),
    "max" = if (rlang::is_installed("matrixStats")) {
      matrixStats::rowMaxs(z_mat, na.rm = na.rm)
    } else {
      apply(z_mat, 1L, max, na.rm = na.rm)
    },
    "median" = if (rlang::is_installed("matrixStats")) {
      matrixStats::rowMedians(z_mat, na.rm = na.rm)
    } else {
      apply(z_mat, 1L, stats::median, na.rm = na.rm)
    }
  )
  # Single-row / all-NA edge case: rowMeans returns NaN, min/max/median return NA
  if (combine_rule == "mean") {
    comb[!is.finite(comb)] <- NA_real_
  }

  # --- 5. Assemble output ---------------------------------------------------
  out <- data.table::data.table(
    patient = as.character(patient),
    clone_id = as.character(clone_id),
    comb_viability = as.numeric(comb)
  )

  if (include_single_z) {
    z_dt <- data.table::as.data.table(z_mat)
    data.table::setnames(z_dt, paste0(colnames(z_mat), ".z"))
    out <- cbind(out, z_dt)
  }

  # Documented return shape is a plain data.frame.
  # ! return data.frame for compatible with Shiny app
  as.data.frame(out)
}
