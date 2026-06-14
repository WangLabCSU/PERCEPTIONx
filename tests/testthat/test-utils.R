# =============================================================================
# PERCEPTION Unit Tests: Utility Functions
# =============================================================================
library(PERCEPTION)
library(testthat)

test_that("rank_normalization_mat works on a simple matrix", {
  mat <- matrix(c(3, 1, 2, 6, 4, 5), nrow = 2, byrow = TRUE)
  result <- rank_normalization_mat(mat)
  expect_type(result, "double")
  expect_equal(dim(result), dim(mat))
  # Each column should be rank-normalized (0-1 range)
  for (j in seq_len(ncol(result))) {
    vals <- result[, j]
    expect_true(all(vals >= 0 & vals <= 1))
  }
})

test_that("rank_normalization_mat handles single-row matrix", {
  mat <- matrix(c(1, 2, 3), nrow = 1)
  result <- rank_normalization_mat(mat)
  expect_type(result, "double")
  expect_equal(length(result), 3)
})

test_that("range01 scales vector to [0, 1]", {
  x <- c(10, 20, 30, 40, 50)
  result <- range01(x)
  expect_equal(min(result), 0)
  expect_equal(max(result), 1)
})

test_that("range01 handles constant vector", {
  x <- rep(5, 10)
  result <- range01(x)
  expect_true(all(is.nan(result)))
})

test_that("fdrcorr returns correct structure", {
  pvals <- c(0.01, 0.05, 0.1, 0.5, 0.9)
  result <- fdrcorr(pvals)
  expect_type(result, "double")
  expect_equal(length(result), length(pvals))
  # FDR-adjusted p-values should be >= original
  expect_true(all(result >= pvals - 1e-10))
})

test_that("%!in% operator works correctly", {
  expect_true(3 %!in% c(1, 2, 4))
  expect_false(2 %!in% c(1, 2, 3))
})

test_that("stripall2match normalizes strings", {
  result <- stripall2match("Hello World 123")
  expect_type(result, "character")
  # Should be lowercase with no spaces
  expect_false(grepl(" ", result))
})

test_that("strsplit_customv0 splits and extracts correctly", {
  result <- strsplit_customv0("Patient_A_clone1", "_", 1)
  expect_equal(result, "Patient")
  result2 <- strsplit_customv0("Patient_A_clone1", "_", 3)
  expect_equal(result2, "clone1")
})
