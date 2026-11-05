#' PERCEPTIONx Model Training
#'
#' @name train_perception
#' @keywords internal
NULL

#' Get specific drug response data for cell-lines
#'
#' Extracts the drug response (AUC) for a given drug from the DepMap database.
#' The function handles multiple screening batches by prioritizing MTS over HTS
#' and selecting the batch with the fewest missing values.
#'
#' @param infunc_drugName Character string. Name of the drug (e.g., "erlotinib").
#'
#' @return A named numeric vector of AUC values for all cell lines.
#'
#' @export
get_response_matrix <- function(infunc_drugName) {
  if (!exists("DepMap")) {
    cli::cli_abort(
      "DepMap data not loaded. Please run load_depmap(read = TRUE) first."
    )
  }

  # Strip BOTH sides: stripall2match() removes non-alphanumerics, so comparing
  # the stripped annotation names against a raw drug name silently failed for
  # every drug containing a hyphen/space (e.g. "5-fluorouracil" -> "5fluorouracil"
  # != "5-fluorouracil"). That yielded an empty row index, and colMeans() of the
  # resulting 0-row matrix returned an all-NaN vector, so training reported
  # "feature ranking failed" for a drug whose response data was present.
  infunc_drugName_id <- which(
    stripall2match(DepMap$secondary_screen_drugAnnotation$CommonName) ==
      stripall2match(infunc_drugName)
  )
  infunc_response <- DepMap$secondary_prism[infunc_drugName_id, ]

  if (is.matrix(infunc_response)) {
    infunc_drugName_screenID <- DepMap$secondary_screen_drugAnnotation$Screen_id[
      infunc_drugName_id
    ]
    infunc_drugName_screenID_trimmed <- substring(
      infunc_drugName_screenID,
      0L,
      3L
    )

    sum_hts <- sum(infunc_drugName_screenID_trimmed == 'HTS')
    sum_mts <- sum(infunc_drugName_screenID_trimmed == 'MTS')

    if (sum_hts == 1L && sum_mts == 1L) {
      infunc_response <- infunc_response[
        infunc_drugName_screenID_trimmed == 'MTS',
      ]
    } else if (sum_hts == 1L && sum_mts > 1L) {
      infunc_response <- infunc_response[
        infunc_drugName_screenID_trimmed == 'MTS',
      ]
      infunc_response <- infunc_response[
        count_row_NAs(infunc_response) == min(count_row_NAs(infunc_response)),
      ]
      # Fall through to the colMeans() below: when several MTS batches tie on
      # NA count they are averaged ACROSS cell lines. (rowMeans here would
      # collapse to one scalar per batch and destroy cell-line granularity.)
    } else if (sum_hts == 0L && sum_mts > 1L) {
      infunc_response <- infunc_response[
        count_row_NAs(infunc_response) == min(count_row_NAs(infunc_response)),
      ]
    }
    if (is.matrix(infunc_response)) {
      infunc_response <- colMeans(infunc_response, na.rm = TRUE)
    }
  }
  return(infunc_response)
}


#' Determine training and test cell-lines for a given drug
#'
#' Identifies which cell lines should be used for training and which should be
#' excluded (reserved for testing) based on cancer type and single-cell data availability.
#'
#' @param infunc_cancerType Character. Cancer type for training. Default = "PanCan".
#' @param infunc_drugName Character. Name of the drug.
#' @param exclude_cancer Character. Cancer type to exclude from training. Default = "PanCan".
#' @param infunc_response Named numeric vector. Drug response data.
#' @param force_add_cellLines Logical. Whether to force add additional cell lines. Default = TRUE.
#' @param force_add_cellLines_list Character vector. Cell line IDs to force add. Default = NA.
#'
#' @return A list of length 2:
#' \itemize{
#'   \item \[1\] common_cellLines: Cell lines to use for training
#'   \item \[2\] cellLines2remove: Cell lines excluded (reserved for testing)
#' }
#'
#' @export
get_cellLine_list <- function(
  infunc_cancerType = "PanCan",
  infunc_drugName,
  exclude_cancer = "PanCan",
  infunc_response,
  force_add_cellLines = TRUE,
  force_add_cellLines_list = NA
) {
  if (infunc_cancerType == "PanCan") {
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    CellLinee2exclude <- DepMap$annotation_20Q4$DepMap_ID[
      grep(exclude_cancer, DepMap$annotation_20Q4$lineage)
    ]

    if (exclude_cancer == "PanCan") {
      CellLinee2exclude <- intersect(
        CellLines_with_singleCellExp,
        DepMap$annotation_20Q4$DepMap_ID
      )
    }

    cellLines2remove <- intersect(
      CellLinee2exclude,
      CellLines_with_singleCellExp
    )

    common_cellLines <- setdiff(
      intersect(CellLines_with_drugResponse, CellLines_with_bulkExp),
      cellLines2remove
    )
    if (force_add_cellLines) {
      common_cellLines <- c(
        common_cellLines,
        intersect(CellLines_with_bulkExp, force_add_cellLines_list)
      )
    }
  } else {
    CellLines_cancer <- DepMap$annotation_20Q4$DepMap_ID[
      grep(infunc_cancerType, DepMap$annotation_20Q4$lineage)
    ]
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    cellLines2remove <- intersect(
      CellLines_cancer,
      CellLines_with_singleCellExp
    )
    common_cellLines <- setdiff(
      Reduce(
        intersect,
        list(
          CellLines_with_drugResponse,
          CellLines_with_bulkExp,
          CellLines_cancer
        )
      ),
      cellLines2remove
    )
  }
  list(common_cellLines, cellLines2remove)
}


#' Feature ranking for a single drug using bulk expression
#'
#' Computes Pearson correlation between each gene's expression and drug response,
#' then ranks genes by absolute correlation. This matches the original
#' feature_ranking_bulk function from step0B.
#'
#' @param infunc_drugName Character. Name of the drug.
#' @param infunc_cancerType Character. Cancer type for training. Default = "PanCan".
#' @param exclude_cancer Character. Cancer type to exclude. Default = "PanCan".
#' @param infunc_GOI Character vector. Genes of Interest to rank.
#' @param ncores Number of cores to use for parallelization. Default = 4.
#'
#' @return Matrix of ranked features with columns: p.value, estimate.cor (sorted by abs(cor) descending)
#' @export
feature_ranking_bulk <- function(
  infunc_drugName,
  infunc_cancerType = "PanCan",
  exclude_cancer = "PanCan",
  infunc_GOI,
  ncores = 4L
) {
  if (!exists("DepMap")) {
    cli::cli_abort(
      "DepMap data not loaded. Please run load_depmap(read = TRUE) first."
    )
  }

  # Get drug response (AUC only, as IC50/VFC data not included in DepMap)
  infunc_response <- get_response_matrix(infunc_drugName)

  # Determine training cell lines (same logic as original)
  if (infunc_cancerType == "PanCan") {
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    CellLinee2exclude <- DepMap$annotation_20Q4$DepMap_ID[
      grep(exclude_cancer, DepMap$annotation_20Q4$lineage)
    ]
    if (exclude_cancer == "PanCan") {
      CellLinee2exclude <- intersect(
        CellLines_with_singleCellExp,
        DepMap$annotation_20Q4$DepMap_ID
      )
    }
    cellLines2remove <- intersect(
      CellLinee2exclude,
      CellLines_with_singleCellExp
    )

    common_cellLines <- setdiff(
      intersect(CellLines_with_drugResponse, CellLines_with_bulkExp),
      cellLines2remove
    )
    cli::cli_alert_info("Common cell lines: {length(common_cellLines)}")
  } else {
    CellLines_cancer <- DepMap$annotation_20Q4$DepMap_ID[
      grep(infunc_cancerType, DepMap$annotation_20Q4$lineage)
    ]
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    common_cellLines <- setdiff(
      Reduce(
        intersect,
        list(
          CellLines_with_drugResponse,
          CellLines_with_bulkExp,
          CellLines_cancer
        )
      ),
      colnames(DepMap$scRNA_complete)
    )
    cli::cli_alert_info("Common cell lines: {length(common_cellLines)}")
  }

  # Compute correlations — vectorized across ALL genes in one pass.
  # The previous per-gene cor.test loop is O(n_genes) R-level function calls
  # (~5 minutes for 15k genes x 306 lines); a single cor() + pt() pass gives
  # bit-identical Pearson r and p-values in well under a second.
  infunc_response_matchedSubset <- unlist(infunc_response[common_cellLines])
  X <- DepMap$expression_rnorm[infunc_GOI, common_cellLines, drop = FALSE]
  y <- as.numeric(infunc_response_matchedSubset)

  # Same complete-case semantics as the old per-gene cor.test: drop response
  # NAs up front, then each gene uses its own complete cells.
  y_ok <- !is.na(y)
  X <- X[, y_ok, drop = FALSE]
  y <- y[y_ok]

  ok <- !is.na(X) # per-gene complete-cell mask
  n_eff <- rowSums(ok)
  r <- if (rlang::is_installed("WGCNA")) {
    as.numeric(WGCNA::cor(
      t(X),
      y,
      use = "pairwise.complete.obs",
      nThreads = ncores
    ))
  } else {
    as.numeric(stats::cor(t(X), y, use = "pairwise.complete.obs"))
  }
  df <- n_eff - 2L
  t_stat <- sqrt(df) * r / sqrt(1L - r^2L)
  pval <- 2L * stats::pt(-abs(t_stat), df)
  pval[!is.finite(t_stat) | is.na(r) | df < 1L] <- NA
  r[!is.finite(t_stat)] <- NA # constant rows (sd = 0) -> NA, like cor.test

  cor_profile <- cbind(p.value = pval, estimate.cor = r)
  rownames(cor_profile) <- rownames(X)
  cor_profile <- cor_profile[order(abs(cor_profile[, 2L]), decreasing = TRUE), ]

  return(cor_profile)
}


#' Parallel feature ranking for multiple drugs
#'
#' Runs feature_ranking_bulk for a list of drugs.
#'
#' NOTE: feature ranking is fully vectorized (~0.15 s per drug over 15k genes),
#' so this deliberately runs SERIALLY. The old Windows PSOCK path serialized
#' the 567 MB DepMap into every worker (clusterExport) and then still failed —
#' workers are clean R sessions without the package's internal helpers — before
#' silently falling back to serial; that cost ~10x for zero benefit.
#'
#' @param infunc_DrugsToUse Character vector. Drug names to rank features for.
#' @param id_cancerType Character. Cancer type. Default = "PanCan".
#' @param infunc_exclude_cancer Character. Cancer type to exclude. Default = "PanCan".
#' @param infunc_GOI Character vector. Genes of Interest.
#' @param ncores Integer. Kept for API compatibility; ignored (see NOTE above).
#'
#' @return A list of feature ranking results, one per drug.
#' @export
run_parallel_feature_ranking_bulk <- function(
  infunc_DrugsToUse,
  id_cancerType = "PanCan",
  infunc_exclude_cancer = "PanCan",
  infunc_GOI,
  ncores = 4L
) {
  idx <- seq_along(infunc_DrugsToUse)

  featuresRank_fromBulk <- lapply(
    idx,
    function(x) {
      err_handle(feature_ranking_bulk(
        infunc_drugName = infunc_DrugsToUse[x],
        infunc_cancerType = id_cancerType,
        infunc_GOI = infunc_GOI,
        exclude_cancer = infunc_exclude_cancer
      ))
      # matrix/NA
    }
  )

  names(featuresRank_fromBulk) <- infunc_DrugsToUse
  featuresRank_fromBulk # a list of matrices
}


#' Build PERCEPTIONx model on bulk expression data
#'
#' Core training function that builds a glmnet model using bulk expression,
#' then evaluates on pseudo-bulk, bulk test, and single-cell test data.
#' This is the direct port of build_on_BULK_v2 from step0B.
#'
#' @param infunc_drugName Character. Drug name.
#' @param infunc_cancerType Character. Cancer type for training. Default = "PanCan".
#' @param exclude_cancer Character. Cancer type to exclude. Default = "PanCan".
#' @param infunc_features Character vector. Feature gene names (ranked).
#' @param single_best Character. Name of the single best feature. Default = infunc_features\[1\].
#' @param k_features Integer. Number of top features to use. Default = 100.
#' @param infunc_alpha Numeric. Alpha for glmnet (not used directly, tuning via alpha_gradient). Default = 1.
#' @param model_type Character. Model type: "glmnet" or "rf". Default = "glmnet".
#' @param num_folds Integer. Number of CV folds. Default = 3.
#' @param num_tree Integer. Number of trees for random forest. Default = 500.
#' @param seed Integer. Random seed for training. Default = 1.
#' @param alpha_gradient Numeric. Step size for alpha grid in glmnet tuning. Default = 0.05.
#' @param lambda_gradient Integer. Number of lambda values in glmnet tuning grid. Default = 20.
#' @param lambda_range Numeric vector of length 2. Min and max lambda for tuning grid. Default = c(0.0001, 1).
#' @param cv_method Character. Cross-validation method for caret::trainControl. Default = "cv".
#'
#' @return A list containing: model, single_best, model_performance_during_cv,
#'         performance_in_bulk, performance_in_pseudo_bulk, performance_in_scRNA,
#'         predVSgroundTruth.
#' @export
build_on_BULK_v2 <- function(
  infunc_drugName,
  infunc_cancerType = "PanCan",
  exclude_cancer = "PanCan",
  infunc_features,
  single_best = infunc_features[1L],
  k_features = 100L,
  infunc_alpha = 1L,
  model_type = "glmnet",
  num_folds = 3L,
  num_tree = 500L,
  seed = 1L,
  alpha_gradient = 0.05,
  lambda_gradient = 20L,
  lambda_range = c(0.0001, 1L),
  cv_method = "cv"
) {
  # Get drug response
  infunc_response <- get_response_matrix(infunc_drugName = infunc_drugName)

  # Get training cell lines
  cellLine_list <- get_cellLine_list(
    infunc_cancerType = infunc_cancerType,
    infunc_drugName = infunc_drugName,
    exclude_cancer = exclude_cancer,
    infunc_response = infunc_response,
    force_add_cellLines = TRUE
  )
  common_cellLines <- cellLine_list[[1L]]

  # Get test cell lines (excluded from training)
  cellLines2remove <- cellLine_list[[2L]]

  # Check consistency between rnorm matrices
  if (nrow(DepMap$expression_rnorm) != nrow(DepMap$scRNA_subset_rnorm)) {
    warning(
      "Mismatch in initial features (rNorm sc vs bulk); skipping this model."
    )
    return(NULL)
  }

  ############
  # Train Data: on bulk
  ##########
  Train_infunc_features_id <- utils::head(
    stats::na.omit(match(infunc_features, rownames(DepMap$expression_rnorm))),
    k_features
  )
  Train_features <- t(DepMap$expression_rnorm[
    Train_infunc_features_id,
    common_cellLines
  ])
  Train_Target_Label <- infunc_response[common_cellLines]
  Train_features <- Train_features[which(!is.na(Train_Target_Label)), ]
  Train_Target_Label <- unlist(Train_Target_Label[which(
    !is.na(Train_Target_Label)
  )])
  Train_features <- Train_features[which(Train_Target_Label < 10L^100L), ]
  Train_Target_Label <- unlist(Train_Target_Label[which(
    Train_Target_Label < 10L^100L
  )])

  ############
  # Test Data 1: on pseudo-Bulk
  ##########
  cell_ID_for_test_data1 <- Reduce(
    intersect,
    list(
      cellLines2remove,
      colnames(DepMap$scRNA_subset_rnorm),
      names(infunc_response)
    )
  )
  Test_infunc_features_id <- utils::head(
    stats::na.omit(match(infunc_features, rownames(DepMap$scRNA_subset_rnorm))),
    k_features
  )
  Test_features <- t(DepMap$scRNA_subset_rnorm[
    Test_infunc_features_id,
    cell_ID_for_test_data1
  ])
  Test_Target_Label <- unlist(infunc_response)[cell_ID_for_test_data1]
  Test_features <- Test_features[which(!is.na(Test_Target_Label)), ]
  Test_Target_Label <- unlist(Test_Target_Label[which(
    !is.na(Test_Target_Label)
  )])

  ############
  # Test Data 2: on bulk
  ##########
  cell_ID_for_test_data2 <- Reduce(
    intersect,
    list(
      cellLines2remove,
      colnames(DepMap$expression_rnorm),
      names(infunc_response)
    )
  )
  Test2_infunc_features_id <- utils::head(
    stats::na.omit(match(infunc_features, rownames(DepMap$expression_rnorm))),
    k_features
  )
  Test2_features <- t(DepMap$expression_rnorm[
    Test2_infunc_features_id,
    cell_ID_for_test_data2
  ])
  Test2_Target_Label <- unlist(infunc_response)[cell_ID_for_test_data2]
  Test2_features <- Test2_features[which(!is.na(Test2_Target_Label)), ]
  Test2_Target_Label <- unlist(Test2_Target_Label[which(
    !is.na(Test2_Target_Label)
  )])

  ############
  # Test Data 3: on sc-RNA-seq (All cancer)
  ##########
  cells_ID_for_test_data3 <- DepMap$metadata_CPM_scRNA$NAME[
    DepMap$metadata_CPM_scRNA$DepMap_ID %in%
      Reduce(intersect, list(cellLines2remove, names(infunc_response)))
  ]
  cellLine_ID_for_test_data3 <- unique(
    DepMap$metadata_CPM_scRNA$DepMap_ID[
      DepMap$metadata_CPM_scRNA$DepMap_ID %in%
        Reduce(intersect, list(cellLines2remove, names(infunc_response)))
    ]
  )
  Test3_infunc_features_id <- utils::head(
    stats::na.omit(match(
      infunc_features,
      rownames(DepMap$CPM_scRNA_CCLE_rnorm)
    )),
    k_features
  )
  Test3_features <- t(DepMap$CPM_scRNA_CCLE_rnorm[
    Test3_infunc_features_id,
    cells_ID_for_test_data3
  ])
  Test3_Target_Label <- unlist(infunc_response)[as.character(
    cellLine_ID_for_test_data3
  )]

  ##########
  # ML Model
  ##########
  set.seed(seed)
  tc <- caret::trainControl(method = cv_method, number = num_folds)

  if (model_type == 'rf') {
    tunegrid <- expand.grid(.mtry = ncol(Train_features) / seq(10L, 1L, -3L))
    cv.out <- suppressWarnings(caret::train(
      Train_Target_Label ~ .,
      data = data.frame(Train_Target_Label, Train_features),
      method = 'rf',
      trControl = tc,
      tunegrid = tunegrid,
      ntree = num_tree
    ))
  } else if (model_type == 'glmnet') {
    tunegrid <- expand.grid(
      alpha = seq(0L, 1L, alpha_gradient),
      lambda = seq(lambda_range[1L], lambda_range[2L], length = lambda_gradient)
    )
    cv.out <- suppressWarnings(caret::train(
      Train_Target_Label ~ .,
      data = data.frame(Train_Target_Label, Train_features),
      method = 'glmnet',
      trControl = tc,
      metric = "RMSE",
      tuneGrid = tunegrid,
      savePredictions = "final"
    ))
  }

  ############
  # Performance
  ##########
  # Best performance during CV. NOTE: sqrt(R²) = |correlation| — R² discards the
  # sign, so this metric is the absolute value of the CV correlation (no p-value).
  model_performance_cv <- max(sqrt(cv.out$results$Rsquared), na.rm = TRUE)

  # Test on bulk (test 2)
  Test_pred_bulk <- stats::predict(cv.out, Test2_features)
  performance_in_bulk <- unlist(stats::cor.test(
    Test_pred_bulk,
    Test2_Target_Label,
    method = "pearson"
  )[c(3L, 4L)])

  # Test on pseudo-Bulk (test 1)
  Test_pred_Msc <- predict(cv.out, newdata = Test_features)
  performance_in_pseudo_bulk <- unlist(stats::cor.test(
    Test_pred_Msc,
    Test_Target_Label,
    method = "pearson"
  )[c(3L, 4L)])

  # Test on scRNA-seq (test 3)
  Test_pred_sc <- stats::predict(cv.out, newdata = Test3_features)
  cells_mapping2_cellLine <- DepMap$metadata_CPM_scRNA$DepMap_ID[
    DepMap$metadata_CPM_scRNA$NAME %in% names(Test_pred_sc)
  ]
  Test_pred_sc_mean_byCellLine <- stats::aggregate(
    Test_pred_sc ~ cells_mapping2_cellLine,
    data.frame(Test_pred_sc, cells_mapping2_cellLine),
    function(x) mean(x, na.rm = TRUE)
  )
  performance_in_scRNA <- unlist(stats::cor.test(
    Test_pred_sc_mean_byCellLine$Test_pred_sc,
    Test3_Target_Label[
      as.character(Test_pred_sc_mean_byCellLine$cells_mapping2_cellLine)
    ],
    method = "pearson"
  )[c(3L, 4L)])

  # Predictions vs ground truth
  predVSgroundTruth <- list(
    pred_gt_scRNA = data.frame(
      Test_pred_sc = Test_pred_sc_mean_byCellLine$Test_pred_sc,
      Observed = Test3_Target_Label[as.character(
        Test_pred_sc_mean_byCellLine$cells_mapping2_cellLine
      )]
    ),
    pred_gt_bulk = data.frame(
      Test_pred_Msc = Test_pred_bulk,
      Observed = Test2_Target_Label
    ),
    pred_gt_mscRNA = data.frame(
      Test_pred_Msc = Test_pred_Msc,
      Observed = Test_Target_Label
    )
  )

  ############
  # Return
  ##########
  list(
    model = cv.out,
    single_best = single_best,
    model_performance_during_cv = model_performance_cv,
    performance_in_bulk = performance_in_bulk,
    performance_in_pseudo_bulk = performance_in_pseudo_bulk,
    performance_in_scRNA = performance_in_scRNA,
    predVSgroundTruth = predVSgroundTruth
  )
}


#' Train PERCEPTIONx models for multiple drugs
#'
#' This function runs the complete PERCEPTIONx training pipeline for a list of drugs.
#' It performs feature ranking, model training with hyperparameter tuning across
#' different k_features values, and selects the best model per drug based on
#' single-cell test performance.
#'
#' @param drug_list Character vector. Names of drugs to train models for.
#'        If NULL, uses the 44 FDA-approved drugs from the paper.
#' @param cancer_type Character. Cancer type for training. Default = "PanCan".
#' @param exclude_cancer Character. Cancer type to exclude from training. Default = "PanCan".
#' @param GOI Character vector. Genes of Interest to use as features.
#' @param k_features_values Numeric vector. Feature counts to test during tuning.
#'        If NULL, automatically calculated from expression_rnorm dimensions.
#' @param ncores Integer. Number of CPU cores for parallel feature ranking. Default = 4.
#' @param output_dir Character. Directory to save trained model RDS file. Default = "./models".
#' @param model_type Character. Model type: "glmnet" or "rf". Default = "glmnet".
#' @param num_folds Integer. Number of CV folds. Default = 3.
#' @param num_tree Integer. Number of trees for random forest. Default = 500.
#' @param seed Integer. Random seed for training. Default = 1.
#' @param alpha_gradient Numeric. Step size for alpha grid in glmnet tuning. Default = 0.05.
#' @param lambda_gradient Integer. Number of lambda values in glmnet tuning grid. Default = 20.
#' @param lambda_range Numeric vector of length 2. Min and max lambda for tuning grid. Default = c(0.0001, 1).
#' @param cv_method Character. Cross-validation method. Default = "cv".
#' @param progress_cb Function. Optional progress callback called as
#'        \code{progress_cb(phase, i, n, drug)} where phase is "rank", "train",
#'        or "done". Useful for progress bars in interactive sessions.
#'
#' @return A named list of trained model objects, one per drug. Also saved as a single RDS file.
#'
#' @examples
#' \dontrun{
#'   load_depmap(read = TRUE)
#'   available_genes <- intersect(rownames(DepMap$expression_20Q4),
#'                                rownames(DepMap$scRNA_complete))
#'   set.seed(123)
#'   GOI_100 <- sample(available_genes, 100)
#'   train_models("abemaciclib", "PanCan", "PanCan", GOI_100, ncores = 1)
#' }
#'
#' @export
train_models <- function(
  drug_list = NULL,
  cancer_type = "PanCan",
  exclude_cancer = "PanCan",
  GOI = NULL,
  k_features_values = NULL,
  ncores = 4L,
  output_dir = "./models",
  model_type = "glmnet",
  num_folds = 3L,
  num_tree = 500L,
  seed = 1L,
  alpha_gradient = 0.05,
  lambda_gradient = 20L,
  lambda_range = c(0.0001, 1L),
  cv_method = "cv",
  progress_cb = NULL
) {
  # progress_cb: optional function(phase, i, n, drug) called as training
  # advances — phase is "rank" (feature ranking done) or "train" (per drug,
  # i/n with the drug name). Lets a Shiny app drive a progress bar during the
  # long synchronous training loop.

  # Hard cap on CPU cores: on a shared deployment server the UI would otherwise
  # expose ALL server cores to every user, letting one session hog the machine.
  ncores <- min(max(1L, as.integer(ncores)), 4L)

  # ============================================================================
  # 1. Check DepMap data
  # ============================================================================
  if (!exists("DepMap")) {
    cli::cli_abort(
      "DepMap data not loaded. Please run load_depmap(read = TRUE) first."
    )
  }

  # ============================================================================
  # 2. Set up GOI (Genes of Interest)
  # ============================================================================
  if (is.null(GOI)) {
    GOI <- intersect(
      rownames(DepMap$expression_20Q4),
      rownames(DepMap$scRNA_complete)
    )
    cli::cli_alert_info(
      "Using all {.val {length(GOI)}} intersecting genes as features."
    )
  } else {
    available_genes <- intersect(GOI, rownames(DepMap$expression_20Q4))
    if (length(available_genes) == 0L) {
      cli::cli_abort(
        "None of the provided GOI genes were found in DepMap expression data."
      )
    }
    if (length(available_genes) < length(GOI)) {
      cli::cli_warn(
        "Some GOI genes not found in DepMap: {toString(setdiff(GOI, available_genes))}"
      )
    }
    GOI <- available_genes
    cli::cli_alert_info(
      "Using user-specified {length(GOI)} genes as features."
    )
  }

  # ============================================================================
  # 3. Determine drug list
  # ============================================================================
  if (is.null(drug_list)) {
    drug_list <- perception_all_drugs
  }

  # ============================================================================
  # 4. Set k_features candidates
  # ============================================================================
  if (is.null(k_features_values)) {
    n_genes_total <- nrow(DepMap$expression_rnorm)
    k_features_values <- n_genes_total * seq(0.0005, 0.02, length = 5L)
    k_features_values <- round(k_features_values)
    cli::cli_alert_info(
      "Testing feature counts: {.val {toString(k_features_values)}}"
    )
  }

  # ============================================================================
  # 5. Create output directory if needed
  # ============================================================================
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  training_start <- Sys.time()

  # ============================================================================
  # 6. Step 1: Parallel feature ranking for all drugs
  # ============================================================================
  cli::cli_progress_step(
    "Step 1/2: Feature ranking for {length(drug_list)} drugs..."
  )
  features_list <- run_parallel_feature_ranking_bulk(
    infunc_DrugsToUse = drug_list,
    id_cancerType = cancer_type,
    infunc_exclude_cancer = exclude_cancer,
    infunc_GOI = GOI,
    ncores = ncores
  )
  if (!is.null(progress_cb)) {
    progress_cb("rank", 0L, length(drug_list), "")
  }

  # ============================================================================
  # 7. Step 2: Train models for each drug with hyperparameter tuning
  #    (matches original Step1_Build_Perception_models.Rmd logic exactly)
  # ============================================================================
  cli::cli_progress_step(
    "Step 2/2: Training models with hyperparameter tuning..."
  )

  for_output_lung_Test_vglm <- vector("list", length(drug_list))
  Performance_by_features <- list()

  for (i in seq_along(drug_list)) {
    drug <- drug_list[i]
    cli::cli_progress_step("{i}/{length(drug_list)} Processing: {.val {drug}}")
    if (!is.null(progress_cb)) {
      progress_cb("train", i, length(drug_list), drug)
    }

    # Check if feature ranking succeeded
    features <- features_list[[i]]
    if (is.null(features) || (length(features) == 1L && is.na(features))) {
      cli::cli_warn("Skipping {drug} - feature ranking failed.")
      for_output_lung_Test_vglm[i] <- list(NULL)
      next
    }

    # Check if drug exists in DepMap response data
    # Strip both sides here too, for the same reason as get_response_matrix():
    # a raw comparison rejects every drug name containing a hyphen or space.
    drug_match <- which(
      stripall2match(DepMap$secondary_screen_drugAnnotation$CommonName) ==
        stripall2match(drug)
    )
    if (length(drug_match) == 0L) {
      cli::cli_warn("Skipping {drug} - drug not found in DepMap response data.")
      for_output_lung_Test_vglm[i] <- list(NULL)
      next
    }

    # Iterate over k_features values (exactly like original script)
    Raw_models_output <- lapply(
      seq_along(k_features_values),
      function(id_counter) {
        infunc_k_features_grid <- k_features_values[id_counter]
        tryCatch(
          build_on_BULK_v2(
            infunc_drugName = drug,
            infunc_cancerType = cancer_type,
            exclude_cancer = exclude_cancer,
            infunc_features = rownames(features),
            single_best = rownames(features)[1L],
            k_features = infunc_k_features_grid,
            model_type = model_type,
            num_folds = num_folds,
            num_tree = num_tree,
            seed = seed,
            alpha_gradient = alpha_gradient,
            lambda_gradient = lambda_gradient,
            lambda_range = lambda_range,
            cv_method = cv_method
          ),
          error = function(e) {
            cli::cli_warn(
              "build_on_BULK_v2 failed for {drug} (k = {infunc_k_features_grid}): {e$message}"
            )
            NULL
          }
        )
      }
    )

    # Check if any models were successfully built
    if (
      length(Raw_models_output) == 0L ||
        all(vapply(Raw_models_output, base::is.null, logical(1L)))
    ) {
      cli::cli_warn("Skipping {drug} - all model builds failed.")
      for_output_lung_Test_vglm[i] <- list(NULL)
      next
    }

    # Keep only successful builds (some may be NULL or the "Mismatch" string).
    valid <- vapply(
      Raw_models_output,
      function(x) {
        is.list(x) &&
          !is.null(x$performance_in_scRNA) &&
          length(x$performance_in_scRNA) >= 2L
      },
      logical(1L)
    )

    if (!any(valid)) {
      cli::cli_warn("Skipping {drug} - all model builds failed.")
      for_output_lung_Test_vglm[i] <- list(NULL)
      next
    }

    # Print and store the performance (rows: p-value, correlation)
    perf <- vapply(
      Raw_models_output[valid],
      function(x) x$performance_in_scRNA,
      numeric(2L)
    )
    if (rlang::is_installed("kableExtra")) {
      kableExtra::kable(perf)
    } else {
      print(perf)
    }
    Performance_by_features[[i]] <- perf

    # Store the Tuned model (select based on highest correlation estimate)
    best_pos <- which.max(perf[2L, ])
    for_output_lung_Test_vglm[[i]] <- Raw_models_output[valid][[best_pos]]
    cli::cli_progress_done()
  }
  # Clean up
  gc(verbose = FALSE)

  # All drugs processed — final tick so a progress bar reaches 100%.
  if (!is.null(progress_cb)) {
    progress_cb("done", length(drug_list), length(drug_list), "")
  }

  # Assign names and filter out NULL entries (failed drugs)
  names(for_output_lung_Test_vglm) <- drug_list
  for_output_lung_Test_vglm <- for_output_lung_Test_vglm[
    !vapply(for_output_lung_Test_vglm, base::is.null, logical(1L))
  ]
  Tuned_models_output <- for_output_lung_Test_vglm

  timestamp <- format(training_start, "%Y%m%d_%H%M%S")
  save_filename <- paste0(
    "PERCEPTIONx_models_",
    cancer_type,
    "_ex",
    exclude_cancer,
    "_",
    timestamp,
    ".RDS"
  )
  save_path <- file.path(output_dir, save_filename)

  # Strip formula/call environment references before serialization. The model's
  # formula/call point into the build frame -> package env -> (inside Shiny)
  # session/globals; saveRDS traverses that chain and would serialize the whole
  # 7 GB DepMap, stalling for minutes. Only the DISK copy is stripped — the
  # models returned to the caller keep call/terms intact.
  strip_model_env <- function(m) {
    if (is.list(m) && inherits(m$model, "train")) {
      m$model$call <- NULL
      if (!is.null(m$model$terms)) {
        attr(m$model$terms, ".Environment") <- baseenv()
      }
    }
    m
  }
  out_for_save <- lapply(Tuned_models_output, strip_model_env)
  saveRDS(out_for_save, save_path)
  cli::cli_alert_success("All models saved to: {.path {save_path}}")

  successful <- sum(
    !vapply(
      Tuned_models_output,
      function(x) length(x) == 1L && is.na(x),
      logical(1L)
    )
  )
  cli::cli_alert_success(
    "Training complete. Successful models: {successful} / {length(drug_list)}"
  )

  return(invisible(Tuned_models_output))
}
