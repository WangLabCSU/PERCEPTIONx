# =============================================================================
# PERCEPTION Unit Tests: Prediction Functions
# =============================================================================
library(PERCEPTION)
library(testthat)

test_that("predict_drugs returns matrix with correct dimensions", {
  skip_if_not_installed("glmnet")
  skip_if_not(exists("DepMap"), "DepMap data not loaded")

  models <- train_models(
    drug_list = c("erlotinib", "gefitinib"),
    cancer_type = "PanCan",
    exclude_cancer = "PanCan",
    GOI = NULL,
    ncores = 1
  )

  # Get test expression
  cellLines_test <- get_cellLine_list(
    infunc_cancerType = "PanCan",
    infunc_drugName = "erlotinib",
    exclude_cancer = "PanCan",
    infunc_response = get_response_matrix("erlotinib")
  )[[2]]

  test_cells <- DepMap$metadata_CPM_scRNA$NAME[
    DepMap$metadata_CPM_scRNA$DepMap_ID %in% cellLines_test]
  expr_test <- DepMap$CPM_scRNA_CCLE_rnorm[, test_cells, drop = FALSE]

  pred <- predict_drugs(models, expr_test)
  expect_type(pred, "double")
  expect_equal(ncol(pred), 2) # 2 drugs
  expect_equal(nrow(pred), ncol(expr_test))
})

test_that("predict_patients aggregates clone predictions", {
  # Simulate clone killing data
  clone_killing_df <- data.frame(
    patient  = c("P1", "P1", "P2", "P2", "P2"),
    clone_id = c("P1_c1", "P1_c2", "P2_c1", "P2_c2", "P2_c3"),
    erlotinib = c(-1.2, -0.5, -2.1, -0.3, -1.8),
    check.names = FALSE
  )

  clone_counts_df <- data.frame(
    P1_c1 = c(200, 0),
    P1_c2 = c(300, 0),
    P2_c1 = c(0, 100),
    P2_c2 = c(0, 250),
    P2_c3 = c(0, 150),
    patients = c("P1", "P2")
  )

  result <- predict_patients(clone_killing_df, clone_counts_df, mode = "weighted_average")
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
})

test_that("predict_patients weighted_max mode works", {
  clone_killing_df <- data.frame(
    patient  = c("P1", "P1"),
    clone_id = c("P1_c1", "P1_c2"),
    erlotinib = c(-1.2, -0.5),
    check.names = FALSE
  )

  clone_counts_df <- data.frame(
    P1_c1 = 200,
    P1_c2 = 300,
    patients = "P1"
  )

  result <- predict_patients(clone_killing_df, clone_counts_df, mode = "weighted_max")
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) == 1)
})
