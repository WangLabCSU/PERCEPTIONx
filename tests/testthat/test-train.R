# =============================================================================
# PERCEPTION Unit Tests: Model Training
# =============================================================================
library(PERCEPTION)
library(testthat)

test_that("train_models returns a list with expected structure", {
  skip_if_not_installed("glmnet")
  skip_if_not(exists("DepMap"), "DepMap data not loaded")

  models <- train_models(
    drug_list = "erlotinib",
    cancer_type = "PanCan",
    exclude_cancer = "PanCan",
    GOI = NULL,
    ncores = 1
  )

  expect_type(models, "list")
  expect_true("erlotinib" %in% names(models))
})

test_that("train_models handles invalid drug gracefully", {
  skip_if_not_installed("glmnet")
  skip_if_not(exists("DepMap"), "DepMap data not loaded")

  models <- train_models(
    drug_list = "nonexistent_drug_xyz",
    cancer_type = "PanCan",
    exclude_cancer = "PanCan",
    GOI = NULL,
    ncores = 1
  )

  expect_type(models, "list")
  # Should return NA for invalid drug
  expect_true(is.na(models[["nonexistent_drug_xyz"]]))
})

test_that("feature_ranking_bulk returns ranked features", {
  skip_if_not(exists("DepMap"), "DepMap data not loaded")

  result <- feature_ranking_bulk(
    infunc_drugName = "erlotinib",
    infunc_cancerType = "PanCan",
    exclude_cancer = "PanCan",
    GOI = NULL
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
})

test_that("get_response_matrix returns drug response data", {
  skip_if_not(exists("DepMap"), "DepMap data not loaded")

  resp <- get_response_matrix("erlotinib")
  expect_type(resp, "matrix")
  expect_true(nrow(resp) > 0)
})

test_that("get_cellLine_list returns train/test split", {
  skip_if_not(exists("DepMap"), "DepMap data not loaded")

  result <- get_cellLine_list(
    infunc_cancerType = "PanCan",
    infunc_drugName = "erlotinib",
    exclude_cancer = "PanCan",
    infunc_response = get_response_matrix("erlotinib")
  )

  expect_type(result, "list")
  expect_equal(length(result), 2)
  expect_true(length(result[[1]]) > 0) # training
  expect_true(length(result[[2]]) > 0) # test
})
