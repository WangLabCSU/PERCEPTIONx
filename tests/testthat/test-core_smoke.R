# =============================================================================
# PERCEPTIONx core-pipeline smoke test on mock DepMap data
# -----------------------------------------------------------------------------
# Runs the REAL code path: feature ranking -> train_models -> predict_drugs ->
# predict_patients against a tiny fake DepMap (see helper-mock_depmap.R).
# This is a regression safety net for CI, where the real 567 MB DepMap bundle
# is unavailable and would otherwise skip every heavy test.
# =============================================================================

test_that("core pipeline runs end-to-end on mock DepMap", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("caret")

  assign("DepMap", mock_depmap(), envir = .GlobalEnv)
  on.exit(rm("DepMap", envir = .GlobalEnv), add = TRUE)

  genes <- rownames(DepMap$expression_rnorm)

  models <- train_models(
    drug_list       = "erlotinib",
    cancer_type     = "PanCan",
    exclude_cancer  = "PanCan",
    GOI             = head(genes, 20),
    k_features_values = c(5, 10),   # small: 60 genes would otherwise give k = 0
    ncores          = 1,
    model_type      = "glmnet",
    num_folds       = 2,
    alpha_gradient  = 0.5,          # shrink the tuning grid for speed
    lambda_gradient = 5,
    output_dir      = tempdir()
  )

  expect_type(models, "list")
  expect_true("erlotinib" %in% names(models))

  # Predict clone-level viability on a mock single-cell matrix.
  expr_test <- DepMap$scRNA_subset_rnorm[, 1:2, drop = FALSE]
  pred <- predict_drugs(models, expr_test)
  expect_type(pred, "double")
  expect_equal(ncol(pred), 1)
  expect_equal(nrow(pred), 2)

  # Aggregate to patient level with the weighted_max default.
  clone_ids <- colnames(expr_test)   # clone/sample names from the mock matrix
  clone_viability_df <- data.frame(
    patient  = rep("P1", length(clone_ids)),
    clone_id = clone_ids,
    erlotinib = as.numeric(pred[, 1]),
    check.names = FALSE
  )
  clone_counts <- data.frame(patients = "P1", check.names = FALSE)
  clone_counts[[clone_ids[1]]] <- 100
  clone_counts[[clone_ids[2]]] <- 200

  patient_pred <- predict_patients(clone_viability_df, clone_counts)
  expect_s3_class(patient_pred, "data.frame")
  expect_equal(nrow(patient_pred), 1)
})
