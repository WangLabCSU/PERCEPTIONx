# =============================================================================
# PERCEPTIONx fast correlation test smoke tests
# =============================================================================

test_that("cor.test_trimmed_v0 pearson matches base cor.test", {
  set.seed(3)
  x <- rnorm(20)
  y <- 0.5 * x + rnorm(20)

  res <- cor.test_trimmed_v0(x, y)
  expect_s3_class(res, "htest")
  expect_true("p.value" %in% names(res))
  expect_true("estimate" %in% names(res))

  base <- cor.test(x, y, method = "pearson")
  expect_equal(res$p.value, base$p.value, tolerance = 1e-8)
  expect_equal(res$estimate, base$estimate, tolerance = 1e-8)
})

test_that("cor.test_trimmed_v0 formula interface works", {
  df <- data.frame(a = rnorm(15), b = rnorm(15))
  res <- cor.test_trimmed_v0(~ a + b, data = df)
  expect_s3_class(res, "htest")
  expect_true("p.value" %in% names(res))
})

test_that("cor.test_trimmed_v0 validates input", {
  expect_error(cor.test_trimmed_v0(1:5, 1:6), "same length")
  expect_error(cor.test_trimmed_v0(c("a", "b"), 1:2), "numeric")
  expect_error(cor.test_trimmed_v0(1:2, 1:2), "not enough finite")
})
