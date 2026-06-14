# =============================================================================
# PERCEPTION Unit Tests: Plotting Functions
# =============================================================================
library(PERCEPTION)
library(testthat)
library(ggplot2)

test_that("plot_tsne_response returns ggplot object", {
  tsne_data <- data.frame(
    X = rnorm(50), Y = rnorm(50),
    killing_scaled = runif(50)
  )
  p <- plot_tsne_response(tsne_data, color_var = "killing_scaled")
  expect_s3_class(p, "ggplot")
})

test_that("plot_tsne_response errors on missing columns", {
  bad_data <- data.frame(X = 1:10, Y = 1:10)
  expect_error(plot_tsne_response(bad_data, color_var = "killing_scaled"))
})

test_that("plot_clone_distribution returns ggplot object", {
  clone_dist <- data.frame(
    patients = c("P1", "P1", "P2", "P2"),
    clones   = c("c1", "c2", "c1", "c2"),
    weights  = c(0.6, 0.4, 0.3, 0.7)
  )
  p <- plot_clone_distribution(clone_dist)
  expect_s3_class(p, "ggplot")
})

test_that("plot_clone_killing returns ggplot object", {
  clone_killing <- data.frame(
    patient     = c("P1", "P1", "P2"),
    clone_id    = c("P1_c1", "P1_c2", "P2_c1"),
    comb_killing = c(-1.5, -0.8, -2.1),
    weights     = c(0.6, 0.4, 1.0)
  )
  p <- plot_clone_killing(clone_killing, killing_var = "comb_killing",
                          weights_var = "weights")
  expect_s3_class(p, "ggplot")
})

test_that("plot_clone_killing works without weights", {
  clone_killing <- data.frame(
    patient     = c("P1", "P1"),
    clone_id    = c("P1_c1", "P1_c2"),
    comb_killing = c(-1.5, -0.8)
  )
  p <- plot_clone_killing(clone_killing, killing_var = "comb_killing")
  expect_s3_class(p, "ggplot")
})

test_that("plot_response_boxplot returns ggplot object", {
  exp_vs_pred <- data.frame(
    response = factor(c("R", "R", "NR", "NR", "R", "NR")),
    predicted_killing = c(-2.1, -1.5, -0.3, -0.1, -1.8, -0.5)
  )
  p <- plot_response_boxplot(exp_vs_pred)
  expect_s3_class(p, "ggplot")
})

test_that("plot_roc_curve returns ggplot object", {
  # Need enough points for ROC
  set.seed(123)
  n <- 50
  response <- factor(sample(c("R", "NR"), n, replace = TRUE), levels = c("R", "NR"))
  predictor <- rnorm(n, mean = ifelse(response == "R", -1.5, -0.3))
  p <- plot_roc_curve(response = response, predictor = predictor, smooth_curve = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_model_performance returns ggplot object", {
  skip_if_not(exists("DepMap"), "DepMap data not loaded")

  models <- train_models(
    drug_list = "erlotinib",
    cancer_type = "PanCan",
    exclude_cancer = "PanCan",
    GOI = NULL,
    ncores = 1
  )
  p <- plot_model_performance(models)
  expect_s3_class(p, "ggplot")
})
