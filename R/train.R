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
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }

  infun_response_matrix <- DepMap$secondary_prism
  infunc_drugName_id <- which(stripall2match(DepMap$secondary_screen_drugAnnotation$CommonName) == infunc_drugName)
  infunc_response <- DepMap$secondary_prism[infunc_drugName_id, ]

  if (is.matrix(infunc_response)) {
    infunc_drugName_screenID <- DepMap$secondary_screen_drugAnnotation$Screen_id[infunc_drugName_id]
    infunc_drugName_screenID_trimmed <- substring(infunc_drugName_screenID, 0, 3)

    if (sum(infunc_drugName_screenID_trimmed == 'HTS') == 1 &
        sum(infunc_drugName_screenID_trimmed == 'MTS') == 1) {
      infunc_response <- infunc_response[infunc_drugName_screenID_trimmed == 'MTS', ]
    } else if (sum(infunc_drugName_screenID_trimmed == 'HTS') == 1 &
               sum(infunc_drugName_screenID_trimmed == 'MTS') > 1) {
      infunc_response <- infunc_response[infunc_drugName_screenID_trimmed == 'MTS', ]
      infunc_response <- infunc_response[count_row_NAs(infunc_response) == min(count_row_NAs(infunc_response)), ]
      if (is.matrix(infunc_response)) {
        infunc_response <- rowMeans(infunc_response, na.rm = TRUE)
      }
    } else if (sum(infunc_drugName_screenID_trimmed == 'HTS') == 0 &
               sum(infunc_drugName_screenID_trimmed == 'MTS') > 1) {
      infunc_response <- infunc_response[count_row_NAs(infunc_response) == min(count_row_NAs(infunc_response)), ]
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
#'   \item Element 1: common_cellLines - Cell lines to use for training
#'   \item Element 2: cellLines2remove - Cell lines excluded (reserved for testing)
#' }
#'
#' @export
get_cellLine_list <- function(infunc_cancerType = "PanCan",
                              infunc_drugName,
                              exclude_cancer = "PanCan",
                              infunc_response,
                              force_add_cellLines = TRUE,
                              force_add_cellLines_list = NA) {
  if (infunc_cancerType == "PanCan") {
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    CellLinee2exclude <- DepMap$annotation_20Q4$DepMap_ID[
      grep(exclude_cancer, DepMap$annotation_20Q4$lineage)]

    if (exclude_cancer == "PanCan") {
      CellLinee2exclude <- intersect(CellLines_with_singleCellExp,
                                     DepMap$annotation_20Q4$DepMap_ID)
    }

    cellLines2remove <- intersect(CellLinee2exclude, CellLines_with_singleCellExp)

    common_cellLines <- setdiff(intersect(CellLines_with_drugResponse, CellLines_with_bulkExp),
                                cellLines2remove)
    if (force_add_cellLines) {
      common_cellLines <- c(common_cellLines, intersect(CellLines_with_bulkExp, force_add_cellLines_list))
    }
  } else {
    CellLines_cancer <- DepMap$annotation_20Q4$DepMap_ID[
      grep(infunc_cancerType, DepMap$annotation_20Q4$lineage)]
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    cellLines2remove <- intersect(CellLines_cancer, CellLines_with_singleCellExp)
    common_cellLines <- setdiff(Reduce(intersect, list(CellLines_with_drugResponse,
                                                       CellLines_with_bulkExp,
                                                       CellLines_cancer)),
                                cellLines2remove)
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
#'
#' @return Matrix of ranked features with columns: p.value, estimate.cor (sorted by abs(cor) descending)
#' @export
feature_ranking_bulk <- function(infunc_drugName,
                                 infunc_cancerType = "PanCan",
                                 exclude_cancer = "PanCan",
                                 infunc_GOI) {

  if (!exists("DepMap")) {
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }

  # Get drug response (AUC only, as IC50/VFC data not included in DepMap)
  infunc_response <- get_response_matrix(infunc_drugName)

  # Determine training cell lines (same logic as original)
  if (infunc_cancerType == "PanCan") {
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    CellLinee2exclude <- DepMap$annotation_20Q4$DepMap_ID[
      grep(exclude_cancer, DepMap$annotation_20Q4$lineage)]
    if (exclude_cancer == "PanCan") {
      CellLinee2exclude <- intersect(CellLines_with_singleCellExp,
                                     DepMap$annotation_20Q4$DepMap_ID)
    }
    cellLines2remove <- intersect(CellLinee2exclude, CellLines_with_singleCellExp)

    common_cellLines <- setdiff(intersect(CellLines_with_drugResponse, CellLines_with_bulkExp),
                                cellLines2remove)
    print(length(common_cellLines))
  } else {
    CellLines_cancer <- DepMap$annotation_20Q4$DepMap_ID[
      grep(infunc_cancerType, DepMap$annotation_20Q4$lineage)]
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    common_cellLines <- setdiff(Reduce(intersect, list(CellLines_with_drugResponse,
                                                       CellLines_with_bulkExp,
                                                       CellLines_cancer)),
                                colnames(DepMap$scRNA_complete))
    print(length(common_cellLines))
  }

  # Compute correlations
  infunc_response_matchedSubset <- unlist(infunc_response[common_cellLines])
  cor_profile <- t(apply(DepMap$expression_rnorm[infunc_GOI, common_cellLines], 1, function(x)
    unlist(cor.test_trimmed_v0(x, infunc_response_matchedSubset, method = "pearson"))))
  cor_profile <- cor_profile[order(abs(cor_profile[, 2]), decreasing = TRUE), ]

  return(cor_profile)
}


#' Parallel feature ranking for multiple drugs
#'
#' Runs feature_ranking_bulk in parallel for a list of drugs.
#'
#' @param infunc_DrugsToUse Character vector. Drug names to rank features for.
#' @param id_cancerType Character. Cancer type. Default = "PanCan".
#' @param infunc_exclude_cancer Character. Cancer type to exclude. Default = "PanCan".
#' @param infunc_GOI Character vector. Genes of Interest.
#' @param ncores Integer. Number of cores for parallel processing. Default = 4.
#'
#' @return A list of feature ranking results, one per drug.
#' @export
run_parallel_feature_ranking_bulk <- function(infunc_DrugsToUse,
                                              id_cancerType = "PanCan",
                                              infunc_exclude_cancer = "PanCan",
                                              infunc_GOI,
                                              ncores = 4) {

  featuresRank_fromBulk <- parallel::mclapply(1:length(infunc_DrugsToUse), function(x)
    err_handle(feature_ranking_bulk(infunc_drugName = infunc_DrugsToUse[x],
                                    infunc_cancerType = id_cancerType,
                                    infunc_GOI = infunc_GOI,
                                    exclude_cancer = infunc_exclude_cancer)),
    mc.cores = ncores)

  names(featuresRank_fromBulk) <- infunc_DrugsToUse
  featuresRank_fromBulk
}


#' Build PERCEPTION model on bulk expression data
#'
#' Core training function that builds a glmnet model using bulk expression,
#' then evaluates on pseudo-bulk, bulk test, and single-cell test data.
#' This is the direct port of build_on_BULK_v2 from step0B.
#'
#' @param infunc_drugName Character. Drug name.
#' @param infunc_cancerType Character. Cancer type for training. Default = "PanCan".
#' @param exclude_cancer Character. Cancer type to exclude. Default = "PanCan".
#' @param infunc_features Character vector. Feature gene names (ranked).
#' @param single_best Character. Name of the single best feature. Default is the first element of infunc_features.
#' @param k_features Integer. Number of top features to use. Default = 100.
#' @param infunc_alpha Numeric. Alpha for glmnet (not used directly, tuning via alpha_gradient). Default = 1.
#' @param model_type Character. Model type: "glmnet" or "rf". Default = "glmnet".
#' @param alpha_gradient Numeric. Step size for alpha grid in glmnet tuning. Default = 0.05.
#' @param lambda_gradient Integer. Number of lambda values in glmnet tuning grid. Default = 20.
#' @param lambda_range Numeric vector of length 2. Min and max lambda for tuning grid. Default = c(0.0001, 1).
#' @param cv_method Character. Cross-validation method for caret::trainControl. Default = "cv".
#'
#' @return A list containing: model, single_best, model_performance_during_cv,
#'         performance_in_bulk, performance_in_pseudo_bulk, performance_in_scRNA,
#'         predVSgroundTruth.
#' @export
build_on_BULK_v2 <- function(infunc_drugName,
                             infunc_cancerType = "PanCan",
                             exclude_cancer = "PanCan",
                             infunc_features,
                             single_best = infunc_features[1],
                             k_features = 100,
                             infunc_alpha = 1,
                             model_type = "glmnet",
                             alpha_gradient = 0.05,
                             lambda_gradient = 20,
                             lambda_range = c(0.0001, 1),
                             cv_method = "cv") {

  # Get drug response
  infunc_response <- get_response_matrix(infunc_drugName = infunc_drugName)

  # Get training cell lines
  common_cellLines <- get_cellLine_list(infunc_cancerType = infunc_cancerType,
                                        infunc_drugName = infunc_drugName,
                                        exclude_cancer = exclude_cancer,
                                        infunc_response = infunc_response,
                                        force_add_cellLines = TRUE)[[1]]

  # Get test cell lines (excluded from training)
  cellLines2remove <- get_cellLine_list(infunc_cancerType = infunc_cancerType,
                                        infunc_drugName = infunc_drugName,
                                        exclude_cancer = exclude_cancer,
                                        infunc_response = infunc_response)[[2]]

  # Check consistency between rnorm matrices
  if (nrow(DepMap$expression_rnorm) != nrow(DepMap$scRNA_subset_rnorm)) {
    return('Mismatch in initial features: rNorm sc vs bulk')
  }

  ############
  # Train Data: on bulk
  ##########
  Train_infunc_features_id <- head(na.omit(match(infunc_features,
                                                  rownames(DepMap$expression_rnorm))), k_features)
  Train_features <- t(DepMap$expression_rnorm[Train_infunc_features_id, common_cellLines])
  Train_Target_Label <- infunc_response[common_cellLines]
  Train_features <- Train_features[which(!is.na(Train_Target_Label)), ]
  Train_Target_Label <- unlist(Train_Target_Label[which(!is.na(Train_Target_Label))])
  Train_features <- Train_features[which(Train_Target_Label < 10^100), ]
  Train_Target_Label <- unlist(Train_Target_Label[which(Train_Target_Label < 10^100)])

  ############
  # Test Data 1: on pseudo-Bulk
  ##########
  cell_ID_for_test_data1 <- Reduce(intersect, list(cellLines2remove,
                                                    colnames(DepMap$scRNA_subset_rnorm),
                                                    names(infunc_response)))
  Test_infunc_features_id <- head(na.omit(match(infunc_features,
                                                 rownames(DepMap$scRNA_subset_rnorm))), k_features)
  Test_features <- t(DepMap$scRNA_subset_rnorm[Test_infunc_features_id, cell_ID_for_test_data1])
  Test_Target_Label <- unlist(infunc_response)[cell_ID_for_test_data1]
  Test_features <- Test_features[which(!is.na(Test_Target_Label)), ]
  Test_Target_Label <- unlist(Test_Target_Label[which(!is.na(Test_Target_Label))])

  ############
  # Test Data 2: on bulk
  ##########
  cell_ID_for_test_data2 <- Reduce(intersect, list(cellLines2remove,
                                                    colnames(DepMap$expression_rnorm),
                                                    names(infunc_response)))
  Test2_infunc_features_id <- head(na.omit(match(infunc_features,
                                                  rownames(DepMap$expression_rnorm))), k_features)
  Test2_features <- t(DepMap$expression_rnorm[Test2_infunc_features_id, cell_ID_for_test_data2])
  Test2_Target_Label <- unlist(infunc_response)[cell_ID_for_test_data2]
  Test2_features <- Test2_features[which(!is.na(Test2_Target_Label)), ]
  Test2_Target_Label <- unlist(Test2_Target_Label[which(!is.na(Test2_Target_Label))])

  ############
  # Test Data 3: on sc-RNA-seq (All cancer)
  ##########
  cells_ID_for_test_data3 <- DepMap$metadata_CPM_scRNA$NAME[
    DepMap$metadata_CPM_scRNA$DepMap_ID %in% Reduce(intersect, list(cellLines2remove,
                                                                     names(infunc_response)))]
  cellLine_ID_for_test_data3 <- unique(
    DepMap$metadata_CPM_scRNA$DepMap_ID[DepMap$metadata_CPM_scRNA$DepMap_ID %in%
      Reduce(intersect, list(cellLines2remove, names(infunc_response)))])
  Test3_infunc_features_id <- head(na.omit(match(infunc_features,
                                                  rownames(DepMap$CPM_scRNA_CCLE_rnorm))), k_features)
  Test3_features <- t(DepMap$CPM_scRNA_CCLE_rnorm[Test3_infunc_features_id,
                                                   cells_ID_for_test_data3])
  Test3_Target_Label <- unlist(infunc_response)[as.character(cellLine_ID_for_test_data3)]

  ##########
  # ML Model
  ##########
  set.seed(1)
  tc <- caret::trainControl(method = cv_method, number = 3)

  if (model_type == 'rf') {
    tunegrid <- expand.grid(.mtry = ncol(Train_features) / seq(10, 1, -3))
    cv.out <- suppressWarnings(caret::train(Train_Target_Label ~ .,
                           data = data.frame(Train_Target_Label, Train_features),
                           method = 'rf',
                           trControl = tc,
                           tunegrid = tunegrid,
                           ntree = 500))
  } else if (model_type == 'glmnet') {
    tunegrid <- expand.grid(alpha = seq(0, 1, alpha_gradient),
                            lambda = seq(lambda_range[1], lambda_range[2], length = lambda_gradient))
    cv.out <- suppressWarnings(caret::train(Train_Target_Label ~ .,
                           data = data.frame(Train_Target_Label, Train_features),
                           method = 'glmnet',
                           trControl = tc,
                           metric = "RMSE",
                           tuneGrid = tunegrid,
                           savePredictions = "final"))
  }

  ############
  # Performance
  ##########
  # Best performance during CV
  model_performance_cv <- max(sqrt(cv.out$results$Rsquared), na.rm = TRUE)

  # Test on bulk (test 2)
  Test_pred_bulk <- predict(cv.out, Test2_features)
  performance_in_bulk <- unlist(cor.test(Test_pred_bulk,
                                         Test2_Target_Label, method = "pearson")[c(3, 4)])

  # Test on pseudo-Bulk (test 1)
  Test_pred_Msc <- predict(cv.out, newdata = Test_features)
  performance_in_pseudo_bulk <- unlist(cor.test(Test_pred_Msc,
                                                Test_Target_Label, method = "pearson")[c(3, 4)])

  # Test on scRNA-seq (test 3)
  Test_pred_sc <- predict(cv.out, newdata = Test3_features)
  cells_mapping2_cellLine <- DepMap$metadata_CPM_scRNA$DepMap_ID[
    DepMap$metadata_CPM_scRNA$NAME %in% names(Test_pred_sc)]
  Test_pred_sc_mean_byCellLine <- aggregate(Test_pred_sc ~ cells_mapping2_cellLine,
                                            data.frame(Test_pred_sc, cells_mapping2_cellLine),
                                            function(x) mean(x, na.rm = TRUE))
  performance_in_scRNA <- unlist(cor.test(Test_pred_sc_mean_byCellLine$Test_pred_sc,
                                          Test3_Target_Label[
                                            as.character(Test_pred_sc_mean_byCellLine$cells_mapping2_cellLine)],
                                          method = "pearson")[c(3, 4)])

  # Predictions vs ground truth
  predVSgroundTruth <- list(
    pred_gt_scRNA = data.frame(Test_pred_sc = Test_pred_sc_mean_byCellLine$Test_pred_sc,
                               Observed = Test3_Target_Label[as.character(Test_pred_sc_mean_byCellLine$cells_mapping2_cellLine)]),
    pred_gt_bulk = data.frame(Test_pred_Msc = Test_pred_bulk,
                              Observed = Test2_Target_Label),
    pred_gt_mscRNA = data.frame(Test_pred_Msc = Test_pred_Msc,
                                Observed = Test_Target_Label)
  )

  ############
  # Return
  ##########
  toreturn <- list(model = cv.out,
                   single_best = single_best,
                   model_performance_during_cv = model_performance_cv,
                   performance_in_bulk = performance_in_bulk,
                   performance_in_pseudo_bulk = performance_in_pseudo_bulk,
                   performance_in_scRNA = performance_in_scRNA,
                   predVSgroundTruth = predVSgroundTruth)
  toreturn
}


#' Train PERCEPTION models for multiple drugs
#'
#' This function runs the complete PERCEPTION training pipeline for a list of drugs.
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
#' @param alpha_gradient Numeric. Step size for alpha grid in glmnet tuning. Default = 0.05.
#' @param lambda_gradient Integer. Number of lambda values in glmnet tuning grid. Default = 20.
#' @param lambda_range Numeric vector of length 2. Min and max lambda for tuning grid. Default = c(0.0001, 1).
#' @param cv_method Character. Cross-validation method. Default = "cv".
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
train_models <- function(drug_list = NULL,
                                    cancer_type = "PanCan",
                                    exclude_cancer = "PanCan",
                                    GOI = NULL,
                                    k_features_values = NULL,
                                    ncores = 4,
                                    output_dir = "./models",
                                    model_type = "glmnet",
                                    alpha_gradient = 0.05,
                                    lambda_gradient = 20,
                                    lambda_range = c(0.0001, 1),
                                    cv_method = "cv") {

  # ============================================================================
  # 1. Check DepMap data
  # ============================================================================
  if (!exists("DepMap")) {
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }

  # ============================================================================
  # 2. Set up GOI (Genes of Interest)
  # ============================================================================
  if (is.null(GOI)) {
    GOI <- intersect(rownames(DepMap$expression_20Q4),
                     rownames(DepMap$scRNA_complete))
    message("Using all ", length(GOI), " intersecting genes as features.")
  } else {
    available_genes <- intersect(GOI, rownames(DepMap$expression_20Q4))
    if (length(available_genes) == 0) {
      stop("None of the provided GOI genes were found in DepMap expression data.")
    }
    if (length(available_genes) < length(GOI)) {
      warning("Some GOI genes not found in DepMap: ",
              paste(setdiff(GOI, available_genes), collapse = ", "))
    }
    GOI <- available_genes
    message("Using user-specified ", length(GOI), " genes as features.")
  }

  # ============================================================================
  # 3. Determine drug list
  # ============================================================================
  if (is.null(drug_list)) {
    drug_list <- c(
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
  }

  # ============================================================================
  # 4. Set k_features candidates
  # ============================================================================
  if (is.null(k_features_values)) {
    n_genes_total <- nrow(DepMap$expression_rnorm)
    k_features_values <- n_genes_total * seq(0.0005, 0.02, length = 5)
    k_features_values <- round(k_features_values)
    message("Testing feature counts: ", paste(k_features_values, collapse = ", "))
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
  message("Step 1/2: Feature ranking for ", length(drug_list), " drugs...")
  features_list <- run_parallel_feature_ranking_bulk(
    infunc_DrugsToUse = drug_list,
    id_cancerType = cancer_type,
    infunc_exclude_cancer = exclude_cancer,
    infunc_GOI = GOI,
    ncores = ncores
  )

  # ============================================================================
  # 7. Step 2: Train models for each drug with hyperparameter tuning
  #    (matches original Step1_Build_Perception_models.Rmd logic exactly)
  # ============================================================================
  message("Step 2/2: Training models with hyperparameter tuning...")

  tuned_models <- list()
  Performance_by_features <- list()

  for (i in seq_along(drug_list)) {
    drug <- drug_list[i]
    message("  Processing: ", drug)

    # Check if feature ranking succeeded
    features <- features_list[[i]]
    if (is.null(features) || (length(features) == 1 && is.na(features))) {
      warning("  Skipping ", drug, " - feature ranking failed.")
      tuned_models[[i]] <- NULL
      next
    }

    # Check if drug exists in DepMap response data
    drug_match <- which(stripall2match(DepMap$secondary_screen_drugAnnotation$CommonName) == drug)
    if (length(drug_match) == 0) {
      warning("  Skipping ", drug, " - drug not found in DepMap response data.")
      tuned_models[[i]] <- NULL
      next
    }

    # Iterate over k_features values (exactly like original script)
    Raw_models_output <- list()
    id_counter <- 1
    build_failed <- FALSE

    for (infunc_k_features_grid in k_features_values) {
      Raw_models_output[[id_counter]] <- tryCatch(
        build_on_BULK_v2(
          infunc_drugName = drug,
          infunc_cancerType = cancer_type,
          exclude_cancer = exclude_cancer,
          infunc_features = rownames(features),
          single_best = rownames(features)[1],
          k_features = infunc_k_features_grid,
          model_type = model_type,
          alpha_gradient = alpha_gradient,
          lambda_gradient = lambda_gradient,
          lambda_range = lambda_range,
          cv_method = cv_method
        ),
        error = function(e) {
          warning("  build_on_BULK_v2 failed for ", drug, " (k=", infunc_k_features_grid, "): ", e$message)
          NULL
        }
      )
      id_counter <- id_counter + 1
    }

    # Check if any models were successfully built
    if (length(Raw_models_output) == 0 || all(sapply(Raw_models_output, is.null))) {
      warning("  Skipping ", drug, " - all model builds failed.")
      tuned_models[[i]] <- NULL
      next
    }

    # Print and store the performance
    print(sapply(Raw_models_output, function(x) x$performance_in_scRNA))
    Performance_by_features[[i]] <- lapply(Raw_models_output, function(x) x$performance_in_scRNA)

    # Store the Tuned model (select based on highest correlation estimate)
    tuned_models[[i]] <- err_handle(
      Raw_models_output[[which.max(sapply(Raw_models_output, function(x) x$performance_in_scRNA)[2, ])]])

    # Clean up
    rm(Raw_models_output)
    gc()
  }

  # Assign names and filter out NULL entries (failed drugs)
  names(tuned_models) <- drug_list
  tuned_models <- tuned_models[!sapply(tuned_models, is.null)]
  Tuned_models_output <- tuned_models

  # ============================================================================
  # 8. Save all models to a single RDS file
  # ============================================================================
  timestamp <- format(training_start, "%Y%m%d_%H%M%S")
  save_filename <- paste0("PERCEPTION_models_", cancer_type,
                          "_ex", exclude_cancer, "_", timestamp, ".RDS")
  save_path <- file.path(output_dir, save_filename)

  saveRDS(Tuned_models_output, save_path)
  message("\nAll models saved to: ", save_path)

  successful <- sum(!sapply(Tuned_models_output, function(x) length(x) == 1 && is.na(x)))
  message("Training complete. Successful models: ", successful, " / ", length(drug_list))

  return(invisible(Tuned_models_output))
}
