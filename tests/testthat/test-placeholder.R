# Placeholder for PERCEPTIONx tests
# Add test files here as needed, e.g., test-train.R, test-predict.R

test_that("package loads correctly", {
  expect_true(requireNamespace("PERCEPTIONx", quietly = TRUE))
})
