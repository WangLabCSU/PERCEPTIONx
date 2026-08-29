# =============================================================================
# PERCEPTIONx Unit Tests: Model Training
# =============================================================================
library(PERCEPTIONx)
library(testthat)

test_that("train_models returns a list with expected structure", {
  skip_if_not_installed("glmnet")
  assign("DepMap", mock_depmap(), envir = .GlobalEnv)
  on.exit(rm("DepMap", envir = .GlobalEnv), add = TRUE)

  models <- train_models(
    drug_list = "erlotinib",
    cancer_type = "PanCan",
    exclude_cancer = "PanCan",
    GOI = head(rownames(DepMap$expression_rnorm), 20),
    k_features_values = c(5, 10),   # small: default would give k = 0 on 60 genes
    ncores = 1,
    output_dir = tempdir()   # do not leave large model RDS in tests/testthat
  )

  expect_type(models, "list")
  expect_true("erlotinib" %in% names(models))
})

test_that("train_models handles invalid drug gracefully", {
  skip_if_not_installed("glmnet")
  assign("DepMap", mock_depmap(), envir = .GlobalEnv)
  on.exit(rm("DepMap", envir = .GlobalEnv), add = TRUE)

  models <- train_models(
    drug_list = "nonexistent_drug_xyz",
    cancer_type = "PanCan",
    exclude_cancer = "PanCan",
    GOI = head(rownames(DepMap$expression_rnorm), 20),
    ncores = 1,
    output_dir = tempdir()
  )

  expect_type(models, "list")
  # Invalid drug is filtered out, leaving an empty list
  expect_length(models, 0)
})

test_that("feature_ranking_bulk returns ranked features", {
  assign("DepMap", mock_depmap(), envir = .GlobalEnv)
  on.exit(rm("DepMap", envir = .GlobalEnv), add = TRUE)

  result <- feature_ranking_bulk(
    infunc_drugName = "erlotinib",
    infunc_cancerType = "PanCan",
    exclude_cancer = "PanCan",
    infunc_GOI = head(rownames(DepMap$expression_rnorm), 50)
  )

  expect_true(is.matrix(result))
  expect_true(nrow(result) > 0)
})

test_that("get_response_matrix returns drug response data", {
  assign("DepMap", mock_depmap(), envir = .GlobalEnv)
  on.exit(rm("DepMap", envir = .GlobalEnv), add = TRUE)

  resp <- get_response_matrix("erlotinib")
  expect_type(resp, "double")
  expect_true(length(resp) > 0)
})

test_that("get_cellLine_list returns train/test split", {
  assign("DepMap", mock_depmap(), envir = .GlobalEnv)
  on.exit(rm("DepMap", envir = .GlobalEnv), add = TRUE)

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
