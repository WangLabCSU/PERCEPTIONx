#' Train PERCEPTION models for multiple drugs
#'
#' This function runs the complete PERCEPTION training pipeline for a list of drugs.
#' It performs feature ranking, model training with hyperparameter tuning, and saves the final models.
#'
#' @param drug_list Character vector. Names of drugs to train models for.
#'        If NULL, uses the 44 drugs that were found to be predictive in the paper.
#' @param cancer_type Character. Cancer type for training. Use "PanCan" for pan-cancer
#'        model, or specific type like "lung", "breast". Default = "PanCan".
#' @param exclude_cancer Character. Cancer type to exclude from training. These cell
#'        lines are reserved for independent testing. Default = "PanCan".
#' @param k_features_values Numeric vector. Feature counts to test during tuning.
#'        If NULL, automatically calculated as 0.05%, 0.5%, 1%, 1.5%, 2% of total genes.
#' @param ncores Integer. Number of CPU cores for parallel feature ranking.
#' @param output_dir Character. Directory to save trained model RDS files.
#'        Default = "./models".
#'
#' @return A named list of trained model objects, one per drug. List names are drug names.
#'         Models are also saved as individual RDS files in output_dir.
#'
#' @details
#' The training pipeline consists of two main steps:
#' \enumerate{
#'   \item \strong{Feature ranking}: Computes Pearson correlation between each gene's
#'         expression and drug response, then ranks genes by absolute correlation.
#'         This step is parallelized across drugs.
#'   \item \strong{Model training}: For each drug, trains elastic net models with
#'         different k_features values and selects the best one based on performance
#'         on single-cell test data.
#' }
#'
#' The function automatically handles:
#' \itemize{
#'   \item Missing drug response data (skips drugs that fail)
#'   \item Multiple feature count candidates (selects best performing)
#'   \item Saving intermediate results to disk
#' }
#'
#' @seealso
#' \code{\link{load_depmap}} for loading DepMap data,
#' \code{\link{run_parallel_feature_ranking_bulk}} for parallel feature ranking,
#' \code{\link{build_on_BULK_v2}} for single-drug model training.
#'
#' @examples
#' \dontrun{
#'   # Load DepMap data first (required)
#'   load_depmap(read = TRUE)
#'
#'   # Train models for default 44 drugs
#'   models <- train_perception_models(ncores = 8)
#'
#'   # Train models for specific drugs only
#'   models <- train_perception_models(
#'       drug_list = c("erlotinib", "gefitinib"),
#'       ncores = 4
#'   )
#'
#'   # Train lung-specific models excluding breast cancer
#'   models <- train_perception_models(
#'       cancer_type = "lung",
#'       exclude_cancer = "breast",
#'       ncores = 8
#'   )
#' }
#'
#' @export
train_perception_models <- function(drug_list = NULL,
                                    cancer_type = "PanCan",
                                    exclude_cancer = "PanCan",
                                    k_features_values = NULL,
                                    ncores = 4,
                                    output_dir = "./models") {

  # 1. Check DepMap data
  if (!exists("DepMap")) {
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }

  # 2. Determine drug list
  if (is.null(drug_list)) {
    # Use the 44 drugs that had predictive models in the paper
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

  # 3. Set k_features candidates (feature counts to test)
  if (is.null(k_features_values)) {
    n_genes <- nrow(DepMap$expression_rnorm)
    # Test 5 different feature counts: 0.05%, 0.5%, 1%, 1.5%, 2% of total genes
    k_features_values <- n_genes * c(0.0005, 0.005, 0.01, 0.015, 0.02)
    k_features_values <- round(k_features_values)
    message("Testing feature counts: ", paste(k_features_values, collapse = ", "))
  }

  # 4. Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 5. Step 1: Parallel feature ranking for all drugs
  message("Step 1/2: Feature ranking for ", length(drug_list), " drugs...")
  features_list <- run_parallel_feature_ranking_bulk(
    inpara_drugName = drug_list,
    inpara_cancerType = cancer_type,
    inpara_exclude_cancer = exclude_cancer,
    ncores = ncores
  )

  # 6. Step 2: Train models for each drug with hyperparameter tuning
  message("Step 2/2: Training models with hyperparameter tuning...")
  all_models <- list()

  for (i in seq_along(drug_list)) {
    drug <- drug_list[i]
    message("  Processing: ", drug)

    # Check if feature ranking succeeded
    features <- features_list[[drug]]
    if (is.null(features) || (length(features) == 1 && is.na(features))) {
      warning("  Skipping ", drug, " - feature ranking failed.")
      next
    }

    # Try different k_features values and select the best
    best_model <- NULL
    best_perf <- -Inf

    for (k in k_features_values) {
      message("    Testing k_features = ", k)
      model <- build_on_BULK_v2(
        infunc_drugName = drug,
        infunc_cancerType = cancer_type,
        infunc_features = rownames(features),
        k_features = k,
        exclude_cancer = exclude_cancer
      )

      # Check if model training succeeded and extract performance
      if (!is.na(model$performance_in_scRNA[2])) {
        current_perf <- model$performance_in_scRNA[2]
        message("      Performance (scRNA correlation): ", round(current_perf, 4))
        if (current_perf > best_perf) {
          best_perf <- current_perf
          best_model <- model
        }
      }
    }

    # Save the best model
    if (!is.null(best_model)) {
      all_models[[drug]] <- best_model
      model_name <- paste(drug, cancer_type, paste0("ex-", exclude_cancer), sep = "_")
      save_path <- file.path(output_dir, paste0(model_name, ".rds"))
      saveRDS(best_model, save_path)
      message("  Best model saved to: ", save_path)
      message("    Best performance: ", round(best_perf, 4))
    } else {
      warning("  Failed to train model for: ", drug)
    }
  }

  # Report summary
  successful <- length(Filter(Negate(is.null), all_models))
  message("\nTraining complete. Successful models: ", successful, " / ", length(drug_list))

  return(invisible(all_models))
}





#' Get specific drug response data for cell-lines
#'
#' Extracts the drug response (AUC) for a given drug from the DepMap database.
#' The function handles multiple screening batches by prioritizing MTS over HTS
#' and selecting the batch with the fewest missing values.
#' The output is a named numeric vector, whose names are cell line IDs (e.g., "ACH-000001")
#' and values are AUC (lower = more sensitive).
#'
#' @param infunc_drugName Character string. Name of the drug (e.g., "erlotinib").
#'
#' @return A named numeric vector of AUC values for all cell lines that have
#' response data for the specified drug.
#'
#' @examples
#' \dontrun{
#' response <- get_response_matrix("erlotinib")
#' head(response, 3)
#' # ACH-000001 ACH-000002 ACH-000003
#' # 0.234 0.456 0.678
#' }
#'
#' @export
get_response_matrix <- function(infunc_drugName){
  if (!exists("DepMap")) {
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }

  # Mode 1: All HTS
  # Mode 2: Priortize MTS
  infun_response_matrix <- DepMap$secondary_prism
  infunc_drugName_id <- which(stripall2match(DepMap$secondary_screen_drugAnnotation$CommonName)==infunc_drugName)
  infunc_response <- DepMap$secondary_prism[infunc_drugName_id,]
  # is.matrix(infunc_response) - denotes if there exists a matrix with multiple screen
  # If there are multiple measures of response;
  if(is.matrix(infunc_response)){
    # Choose screen of interest
    # Rule: Priority, MTS over HTS; MTS with least number of NAs
    infunc_drugName_screenID <- DepMap$secondary_screen_drugAnnotation$Screen_id[infunc_drugName_id]
    infunc_drugName_screenID_trimmed <- substring(infunc_drugName_screenID, 0, 3)
    # if single MTS vs single HTS
    if(sum(infunc_drugName_screenID_trimmed == 'HTS') == 1 &
       sum(infunc_drugName_screenID_trimmed == 'MTS') == 1){
      infunc_response <- infunc_response[infunc_drugName_screenID_trimmed=='MTS', ]
      # if HTS vs multiple MTS
    } else if(sum(infunc_drugName_screenID_trimmed == 'HTS') == 1 &
              sum(infunc_drugName_screenID_trimmed == 'MTS')>1){
      infunc_response <- infunc_response[infunc_drugName_screenID_trimmed == 'MTS', ]
      infunc_response[count_row_NAs(infunc_response) == min(count_row_NAs(infunc_response)),]
      # If above the above is still a tie btw MTS screens; Take a mean
      if(is.matrix(infunc_response)){ infunc_response <- rowMeans(infunc_response, na.rm = T) }
      # if multiple MTS are competing
    } else if(sum(infunc_drugName_screenID_trimmed == 'HTS') == 0 &
              sum(infunc_drugName_screenID_trimmed == 'MTS') > 1){
      infunc_response <- infunc_response[
        count_row_NAs(infunc_response) == min(count_row_NAs(infunc_response)), ]
      # If above the above is still a tie btw MTS screens; Take a mean
    }
    if(is.matrix(infunc_response)){ infunc_response <- colMeans(infunc_response, na.rm = T) }
  }
  infunc_response
}


#' Determine training and test cell-lines for a given drug
#'
#' Identifies which cell lines should be used for training and which should be
#' excluded (reserved for testing) based on cancer type and single-cell data availability.
#' Cell lines with single-cell data are excluded from
#' training to enable independent validation of sc-level predictions.
#'
#' @param infunc_cancerType Character. Cancer type for training. Use "PanCan"
#'        for pan-cancer model, or specific cancer type like "lung".
#' @param infunc_drugName Character. Name of the drug.
#' @param exclude_cancer Character. Cancer type to exclude from training.
#'        When infunc_cancerType = "PanCan", set to "PanCan" to exclude all
#'        cell lines with single-cell data, or specify a cancer name
#'        (e.g., "lung") to exclude only that cancer type.
#' @param infunc_response Named numeric vector. Drug response data from
#'        \code{\link{get_response_matrix}}. Names are cell-line IDs.
#' @param force_add_cellLines Logical. Whether to force add additional cell
#'        lines to the training set. Default = TRUE.
#' @param force_add_cellLines_list Character vector. Cell line IDs to force
#'        add when force_add_cellLines = TRUE. Default = NA.
#'
#' @return A list of length 2:
#' \itemize{
#'   \item [[1]] common_cellLines: Cell lines to use for training
#'   \item [[2]] cellLines2remove: Cell lines excluded (reserved for testing)
#' }
#'
#' @examples
#' \dontrun{
#'   # Get drug response first
#'   response <- get_response_matrix("erlotinib")
#'
#'   # Pan-cancer model: exclude all scRNA-seq cell lines
#'   cells <- get_cellLine_list(
#'       infunc_cancerType = "PanCan",
#'       infunc_drugName = "erlotinib",
#'       exclude_cancer = "PanCan",
#'       infunc_response = response
#'   )
#'   train_cells <- cells[[1]]
#'   test_cells <- cells[[2]]
#'
#'   # Pan-cancer model: exclude only lung cancer cell lines
#'   cells <- get_cellLine_list(
#'       infunc_cancerType = "PanCan",
#'       infunc_drugName = "erlotinib",
#'       exclude_cancer = "lung",
#'       infunc_response = response
#'   )
#' }
#'
#' @seealso \code{\link{get_response_matrix}} for obtaining drug response data
#'
#' @export
get_cellLine_list <- function(infunc_drugName,
                              infunc_cancerType = "PanCan",
                              exclude_cancer = "PanCan",
                              infunc_response,
                              force_add_cellLines = TRUE,
                              force_add_cellLines_list = NA){
  if(infunc_cancerType == "PanCan"){
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    # Exclude 'CellLinee2exclude' cell-lines with scRNA-seq
    CellLinee2exclude <- DepMap$annotation_20Q4$DepMap_ID[
      grep(exclude_cancer, DepMap$annotation_20Q4$lineage)]

    if(exclude_cancer == "PanCan"){
      CellLinee2exclude <- intersect(CellLines_with_singleCellExp, DepMap$annotation_20Q4$DepMap_ID)
    }

    cellLines2remove <- intersect(CellLinee2exclude, CellLines_with_singleCellExp)

    common_cellLines <- setdiff(intersect(CellLines_with_drugResponse, CellLines_with_bulkExp),
                             cellLines2remove)
    if(force_add_cellLines){
      common_cellLines <- c(common_cellLines, intersect(CellLines_with_bulkExp, force_add_cellLines_list))
    }
  } else{
    CellLines_cancer <- DepMap$annotation_20Q4$DepMap_ID[
      grep(infunc_cancerType, DepMap$annotation_20Q4$lineage)]
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)

    # Remove only target cancer with scRNA-seq
    cellLines2remove <- intersect(CellLines_cancer, CellLines_with_singleCellExp)
    common_cellLines <- setdiff(Reduce(intersect, list(CellLines_with_drugResponse,
                                                    CellLines_with_bulkExp,
                                                    CellLines_cancer)),
                             cellLines2remove)
  }
  list(common_cellLines, cellLines2remove)
}



# This function takes, drugname, cancer type, all gene features, response measure (AUC vs viability),
# which cancer type cell lines to exclude from analysis.
feature_ranking_bulk <- function(infunc_drugName,
                                 infunc_cancerType = "PanCan",
                                 exclude_cancer = "PanCan",
                                 infunc_GOI = NULL){
  # Check DepMap data
  if (!exists("DepMap")) {
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }

  if (is.null(infunc_GOI)) {
    infunc_GOI <- intersect(rownames(DepMap$expression_20Q4),
                                              rownames(DepMap$scRNA_complete))
  }

  # Choose the resp matrix
  infun_response_matrix <- DepMap$secondary_prism
  infunc_drugName_id <- which(stripall2match(DepMap$secondary_screen_drugAnnotation$CommonName)==infunc_drugName)
  infunc_response <- DepMap$secondary_prism[infunc_drugName_id, ]

  # is.matrix(infunc_response) - denotes if there exists a matrix with multiple screen
  # If there are multiple measures of response;
  if(is.matrix(infunc_response)){
    # Choose screen of interest
    # Rule: Priority, MTS over HTS; MTS with least number of NAs
    infunc_drugName_screenID <- DepMap$secondary_screen_drugAnnotation$Screen_id[infunc_drugName_id]
    infunc_drugName_screenID_trimmed <- substring(infunc_drugName_screenID, 0, 3)
    # if single MTS vs single HTS
    if(sum(infunc_drugName_screenID_trimmed == 'HTS') == 1 &
       sum(infunc_drugName_screenID_trimmed == 'MTS') == 1){
      infunc_response <- infunc_response[infunc_drugName_screenID_trimmed == 'MTS',]
      # if HTS vs multiple MTS
    } else if(sum(infunc_drugName_screenID_trimmed == 'HTS') == 1 &
              sum(infunc_drugName_screenID_trimmed == 'MTS') > 1){
      infunc_response <- infunc_response[infunc_drugName_screenID_trimmed == 'MTS',]
      infunc_response[count_row_NAs(infunc_response) == min(count_row_NAs(infunc_response)),]
      # If above the above is still a tie btw MTS screens; Take a mean
      if(is.matrix(infunc_response)){ infunc_response=rowMeans(infunc_response, na.rm = T) }
      # if multiple MTS are competing
    } else if(sum(infunc_drugName_screenID_trimmed == 'HTS') == 0 &
              sum(infunc_drugName_screenID_trimmed == 'MTS') > 1){
      infunc_response <- infunc_response[
        count_row_NAs(infunc_response) == min(count_row_NAs(infunc_response)), ]
      # If above the above is still a tie btw MTS screens; Take a mean
    }
    if(is.matrix(infunc_response)){ infunc_response <- colMeans(infunc_response, na.rm = T) }
  }
  # Contrasting to the previous, we want to remove
  # the cell lines which are present in scRNA-seq
  if(infunc_cancerType == "PanCan"){
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)
    # Exclude 'CellLinee2exclude' cell lines with scRNA-seq
    CellLinee2exclude <- DepMap$annotation_20Q4$DepMap_ID[
      grep(exclude_cancer, DepMap$annotation_20Q4$lineage)]
    if(exclude_cancer == "PanCan"){
      CellLinee2exclude <- intersect(CellLines_with_singleCellExp,
                                  DepMap$annotation_20Q4$DepMap_ID)
    }
    cellLines2remove <- intersect(CellLinee2exclude, CellLines_with_singleCellExp)

    common_cellLines <- setdiff(intersect(CellLines_with_drugResponse, CellLines_with_bulkExp),
                             cellLines2remove)
    cat("Drug:", infunc_drugName, "- Training cell lines:", length(common_cellLines), "\n")
  } else{
    CellLines_cancer=DepMap$annotation_20Q4$DepMap_ID[
      grep(infunc_cancerType, DepMap$annotation_20Q4$lineage)]
    CellLines_with_drugResponse <- names(infunc_response)
    CellLines_with_bulkExp <- colnames(DepMap$expression_rnorm)
    CellLines_with_singleCellExp <- colnames(DepMap$scRNA_complete)
    length(CellLines_with_singleCellExp)

    common_cellLines <- setdiff(Reduce(intersect, list(CellLines_with_drugResponse,
                                                    CellLines_with_bulkExp,
                                                    CellLines_cancer)),
                             colnames(DepMap$scRNA_complete))
    cat("Drug:", infunc_drugName, " - Training cell lines:", length(common_cellLines), "\n")
  }

  infunc_response_matchedSubset <- unlist(infunc_response[common_cellLines])
  # Pearson correlation is computed
  # Spearman cor can be compute using providing a ranked vector
  cor_profile <- t(apply(DepMap$expression_rnorm[infunc_GOI, common_cellLines], 1, function(x)
    unlist(cor.test_trimmed_v0(x, infunc_response_matchedSubset, method='p')) ))
  cor_profile <- cor_profile[order(abs(cor_profile[,2]), decreasing = T),]
  return(cor_profile)
}




#' Parallel feature ranking for multiple drugs
#'
#' Runs `feature_ranking_bulk` in parallel for a list of drugs.
#' Note: On Windows, parallel execution is disabled due to technical limitations.
#' Windows users will run in serial mode automatically.
#'
#' @param inpara_drugName Character vector. Names of drugs to analyze.
#' @param inpara_cancerType Character. Cancer type for training.
#'        Use "PanCan" for pan-cancer model, or specific type like "lung", "breast".
#'        Default "PanCan".
#' @param inpara_exclude_cancer Character. Cancer type to exclude from training.
#'        These cell lines are reserved for independent testing.
#'        Default "PanCan" (exclude all cell lines with scRNA-seq data).
#' @param inpara_GOI Character vector. Genes of interest to rank. If NULL,
#'        uses the intersection of genes in bulk expression and scRNA-seq matrices.
#'        Default NULL.
#' @param ncores Integer. Number of CPU cores to use for parallel computation.
#'        This parameter only affects Unix/Mac systems. On Windows, execution is
#'        always serial regardless of this value. Default 4.
#'
#' @return A named list. Each element is the feature ranking result from
#'         `feature_ranking_bulk` for one drug. List names are drug names.
#'
#' Platform-specific behavior:
#' \itemize{
#'   \item On Unix/Mac: True parallel execution using `mclapply` from the `parallel` package.
#'   \item On Windows: Serial execution using `lapply` (parallel is not supported due to
#'         limitations in exporting large data objects to worker processes).
#' }
#'
#' Individual drug failures are caught by `err_handle()` and return NA,
#' preventing the entire process from stopping.
#'
#' @seealso
#' \code{\link{feature_ranking_bulk}} for single-drug feature ranking.
#'
#' @examples
#' \dontrun{
#'   # Load DepMap data first
#'   load_depmap(read = TRUE)
#'
#'   # Define drugs to analyze
#'   drugs <- c("erlotinib", "gefitinib", "osimertinib")
#'
#'   # Run feature ranking (parallel on Unix/Mac, serial on Windows)
#'   features_list <- run_parallel_feature_ranking_bulk(
#'       inpara_drugName = drugs,
#'       ncores = 4
#'   )
#'
#'   # Access results for a specific drug
#'   head(features_list[["erlotinib"]])
#' }
#'
#' @export
run_parallel_feature_ranking_bulk <- function(inpara_drugName,
                                              inpara_cancerType = "PanCan",
                                              inpara_exclude_cancer = "PanCan",
                                              inpara_GOI = NULL,
                                              ncores = 4){
  # Check DepMap data
  if (!exists("DepMap")) {
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }
  # Check gene of interest
  if (is.null(inpara_GOI)) {
    inpara_GOI <- intersect(rownames(DepMap$expression_20Q4),
                            rownames(DepMap$scRNA_complete))
  }

  # Check operating system and choose parallel backend
  if (.Platform$OS.type == "windows") {
    # Windows: use `parLapply`
    if (.Platform$OS.type == "windows") {
      message("Windows detected: Running in serial mode (ncores forced to 1)")
      featuresRank_fromBulk <- lapply(inpara_drugName, function(drug) {
        err_handle(feature_ranking_bulk(
          infunc_drugName = drug,
          infunc_cancerType = inpara_cancerType,
          exclude_cancer = inpara_exclude_cancer,
          infunc_GOI = inpara_GOI
        ))
      })
    }else {
      # Unix/Mac: use `mclapply`
      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Package 'parallel' is required.")
      }
      featuresRank_fromBulk <- parallel::mclapply(inpara_drugName, function(drug) {
        err_handle(feature_ranking_bulk(
          infunc_drugName = drug,
          infunc_cancerType = inpara_cancerType,
          exclude_cancer = inpara_exclude_cancer,
          infunc_GOI = inpara_GOI
        ))
      }, mc.cores = ncores)
    }
  }
  names(featuresRank_fromBulk) <- inpara_drugName
  return(featuresRank_fromBulk)
}



#' Train a PERCEPTION model for a single drug
#'
#' This function trains a PERCEPTION model for a given drug using bulk expression data,
#' then tunes hyperparameters using single-cell data. It evaluates performance on
#' multiple test sets: bulk, pseudo-bulk, and single-cell.
#'
#' @param infunc_drugName Character. Name of the drug (e.g., "erlotinib").
#' @param infunc_cancerType Character. Cancer type for training. Default = "PanCan".
#' @param infunc_features Character vector. Ranked feature genes from feature_ranking_bulk.
#'        Use rownames(features) as input.
#' @param k_features Integer. Number of top features to use. Default = 100.
#' @param model_type Character. "glmnet" for elastic net, "rf" for random forest.
#'        Default = "glmnet".
#' @param exclude_cancer Character. Cancer type to exclude from training. Default = "PanCan".
#' @param alpha_gradient Numeric. Step size for alpha grid search. Default = 0.05.
#' @param lambda_gradient Integer. Number of lambda values to try. Default = 20.
#' @param cv_method Character. Cross-validation method. Default = "cv".
#'
#' @return A list containing:
#' \itemize{
#'   \item model: The trained caret model object
#'   \item model_performance_during_cv: Best R^2 from cross-validation
#'   \item performance_in_bulk: Pearson correlation on bulk test set
#'   \item performance_in_pseudo_bulk: Pearson correlation on pseudo-bulk test set
#'   \item performance_in_scRNA: Pearson correlation on single-cell test set
#'   \item predVSgroundTruth: Data frames of predictions vs observations
#' }
#'
#' @export
build_on_BULK_v2 <- function(infunc_drugName,
                             infunc_cancerType = "PanCan",
                             infunc_features,
                             k_features = 100,
                             model_type = "glmnet",
                             exclude_cancer = "PanCan",
                             alpha_gradient = 0.05,
                             lambda_gradient = 20,
                             cv_method = "cv") {

  # Check DepMap data
  if (!exists("DepMap")) {
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }

  # Step 1: Get drug response
  infunc_response <- get_response_matrix(infunc_drugName)

  # Step 2: Get training and test cell lines
  cell_lines <- get_cellLine_list(
    infunc_cancerType = infunc_cancerType,
    infunc_drugName = infunc_drugName,
    exclude_cancer = exclude_cancer,
    infunc_response = infunc_response
  )
  common_cellLines <- cell_lines[[1]]   # Training set
  cellLines2remove <- cell_lines[[2]]   # Test set (excluded from training)

  # Step 3: Check consistency between bulk and scRNA matrices
  if (nrow(DepMap$expression_rnorm) != nrow(DepMap$scRNA_subset_rnorm)) {
    stop("Mismatch in initial features: bulk vs scRNA")
  }

  #===========================================================================
  # Train Data: Bulk expression
  #===========================================================================

  # Get feature indices for top k_features genes
  feature_idx <- match(infunc_features, rownames(DepMap$expression_rnorm))
  feature_idx <- na.omit(feature_idx)[1:k_features]

  # Training features (cells as rows, genes as columns)
  Train_features <- t(DepMap$expression_rnorm[feature_idx, common_cellLines])
  Train_Target_Label <- infunc_response[common_cellLines]

  # Remove NAs and infinite values
  valid_idx <- which(!is.na(Train_Target_Label) & Train_Target_Label < 1e100)
  Train_features <- Train_features[valid_idx, ]
  Train_Target_Label <- unlist(Train_Target_Label[valid_idx])

  #===========================================================================
  # Test Data 1: Pseudo-bulk (aggregated from single-cell)
  #===========================================================================

  test1_cells <- Reduce(intersect, list(
    cellLines2remove,
    colnames(DepMap$scRNA_subset_rnorm),
    names(infunc_response)
  ))

  feature_idx_test1 <- match(infunc_features, rownames(DepMap$scRNA_subset_rnorm))
  feature_idx_test1 <- na.omit(feature_idx_test1)[1:k_features]

  Test_features <- t(DepMap$scRNA_subset_rnorm[feature_idx_test1, test1_cells])
  Test_Target_Label <- unlist(infunc_response[test1_cells])

  valid_idx1 <- which(!is.na(Test_Target_Label))
  Test_features <- Test_features[valid_idx1, ]
  Test_Target_Label <- unlist(Test_Target_Label[valid_idx1])

  #===========================================================================
  # Test Data 2: Independent bulk
  #===========================================================================

  test2_cells <- Reduce(intersect, list(
    cellLines2remove,
    colnames(DepMap$expression_rnorm),
    names(infunc_response)
  ))

  feature_idx_test2 <- match(infunc_features, rownames(DepMap$expression_rnorm))
  feature_idx_test2 <- na.omit(feature_idx_test2)[1:k_features]

  Test2_features <- t(DepMap$expression_rnorm[feature_idx_test2, test2_cells])
  Test2_Target_Label <- unlist(infunc_response[test2_cells])

  valid_idx2 <- which(!is.na(Test2_Target_Label))
  Test2_features <- Test2_features[valid_idx2, ]
  Test2_Target_Label <- unlist(Test2_Target_Label[valid_idx2])



  #===========================================================================
  # Test Data 3: Real single-cell data (cell-level)
  #===========================================================================

  # Get cell IDs for test cell lines
  test3_cell_lines <- Reduce(intersect, list(cellLines2remove, names(infunc_response)))
  test3_cells <- as.character(DepMap$metadata_CPM_scRNA$NAME)[
    DepMap$metadata_CPM_scRNA$DepMap_ID %in% test3_cell_lines
  ]
  test3_cell_lines_unique <- unique(test3_cell_lines)

  feature_idx_test3 <- match(infunc_features, rownames(DepMap$CPM_scRNA_CCLE_rnorm))
  feature_idx_test3 <- na.omit(feature_idx_test3)[1:k_features]

  Test3_features <- t(DepMap$CPM_scRNA_CCLE_rnorm[feature_idx_test3, test3_cells])
  Test3_Target_Label <- unlist(infunc_response)[as.character(test3_cell_lines_unique)]

  #===========================================================================
  # Train ML Model
  #===========================================================================

  set.seed(1)

  tc <- caret::trainControl(method = cv_method, number = 3)

  if (model_type == "glmnet") {
    tunegrid <- expand.grid(
      alpha = seq(0, 1, alpha_gradient),
      lambda = seq(0.0001, 1, length = lambda_gradient)
    )
    cv.out <- suppressWarnings(
      caret::train(
        Train_Target_Label ~ .,
        data = data.frame(Train_Target_Label, Train_features),
        method = "glmnet",
        trControl = tc,
        metric = "RMSE",
        tuneGrid = tunegrid,
        savePredictions = "final"
      )
    )
  } else if (model_type == "rf") {
    tunegrid <- expand.grid(.mtry = ncol(Train_features) / seq(10, 1, -3))
    cv.out <- suppressWarnings(
      caret::train(
        Train_Target_Label ~ .,
        data = data.frame(Train_Target_Label, Train_features),
        method = "rf",
        trControl = tc,
        tunegrid = tunegrid,
        ntree = 500
      )
    )
  } else {
    stop("model_type must be 'glmnet' or 'rf'")
  }
  #===========================================================================
  # Evaluate Performance
  #===========================================================================

  # Cross-validation performance
  model_performance_cv <- max(sqrt(cv.out$results$Rsquared), na.rm = TRUE)

  # Test on independent bulk
  Test_pred_bulk <- predict(cv.out, Test2_features)
  perf_bulk <- unlist(stats::cor.test(Test_pred_bulk, Test2_Target_Label, method = "pearson")[c(3, 4)])
  names(perf_bulk) <- c("p.value", "estimate.cor")

  # Test on pseudo-bulk
  Test_pred_pseudo <- predict(cv.out, newdata = Test_features)
  perf_pseudo <- unlist(stats::cor.test(Test_pred_pseudo, Test_Target_Label, method = "pearson")[c(3, 4)])
  names(perf_pseudo) <- c("p.value", "estimate.cor")

  # Test on single-cell (aggregate predictions by cell line)
  Test_pred_sc <- predict(cv.out, newdata = Test3_features)
  metadata_clean <- DepMap$metadata_CPM_scRNA[-1, ]
  NAME_clean <- as.character(metadata_clean$NAME)
  DepMap_ID_clean <- metadata_clean$DepMap_ID

  cell_mapping <- DepMap_ID_clean[NAME_clean %in% names(Test_pred_sc)]

  pred_by_cellline <- aggregate(Test_pred_sc ~ cell_mapping,
                                data = data.frame(Test_pred_sc, cell_mapping),
                                function(x) mean(x, na.rm = TRUE))

  perf_sc <- unlist(stats::cor.test(pred_by_cellline$Test_pred_sc,
                                    Test3_Target_Label[as.character(pred_by_cellline$cell_mapping)],
                                    method = "pearson")[c(3, 4)])
  names(perf_sc) <- c("p.value", "estimate.cor")

  #===========================================================================
  # Prepare output
  #===========================================================================

  predVSgroundTruth <- list(
    pred_gt_scRNA = data.frame(
      Test_pred_sc = pred_by_cellline$Test_pred_sc,
      Observed = Test3_Target_Label[as.character(pred_by_cellline$cell_mapping)]
    ),
    pred_gt_bulk = data.frame(
      Test_pred_bulk = Test_pred_bulk,
      Observed = Test2_Target_Label
    ),
    pred_gt_pseudo = data.frame(
      Test_pred_pseudo = Test_pred_pseudo,
      Observed = Test_Target_Label
    )
  )

  # Return results
  toreturn <- list(
    model = cv.out,
    model_performance_during_cv = model_performance_cv,
    performance_in_bulk = perf_bulk,
    performance_in_pseudo_bulk = perf_pseudo,
    performance_in_scRNA = perf_sc,
    predVSgroundTruth = predVSgroundTruth
  )

  return(toreturn)
}



#' Build PERCEPTION models in parallel for multiple drugs
#'
#' This function takes feature ranking results from `run_parallel_feature_ranking_bulk`
#' and trains PERCEPTION models for each drug in parallel using `lapply` (or can be
#' easily modified for true parallel execution).
#'
#' @param features_fromBulk A list of feature ranking results from
#'        `run_parallel_feature_ranking_bulk`. Each element is a matrix
#'        with genes as rows and columns "p.value" and "estimate.cor".
#' @param id_cancerType Character. Cancer type for training. Default = "PanCan".
#' @param infunc_drugsList Character vector. Names of drugs to build models for.
#' @param infunc_k_features Integer. Number of top features to use for each model.
#' @param infunc_model_type Character. "glmnet" or "rf". Default = "glmnet".
#' @param infunc_exclude_cancer Character. Cancer type to exclude from training.
#'
#' @return A list of trained model objects (output from `build_on_BULK_v2`),
#'         one per drug. If a drug fails, the corresponding element is NA.
#'
#' @seealso
#' \code{\link{build_on_BULK_v2}} for single-drug model training.
#'
#' @export
build_biomarker_step2_bulk <- function(features_fromBulk,
                                       id_cancerType = "PanCan",
                                       infunc_drugsList,
                                       infunc_k_features,
                                       infunc_model_type = "glmnet",
                                       infunc_exclude_cancer = "PanCan") {

  # Check DepMap data
  if (!exists("DepMap")) {
    stop("DepMap data not loaded. Please run load_depmap(read = TRUE) first.")
  }

  # Ensure drug list is character vector
  infunc_drugsList <- as.character(infunc_drugsList)

  # Build models for each drug
  model_for <- lapply(seq_along(infunc_drugsList), function(i) {
    err_handle(
      build_on_BULK_v2(
        infunc_drugName = infunc_drugsList[i],
        infunc_cancerType = id_cancerType,
        infunc_features = rownames(features_fromBulk[[i]]),
        k_features = infunc_k_features,
        model_type = infunc_model_type,
        exclude_cancer = infunc_exclude_cancer
      )
    )
  })

  # Name the list elements by drug names
  names(model_for) <- infunc_drugsList

  return(model_for)
}
