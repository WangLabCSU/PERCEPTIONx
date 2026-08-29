# =============================================================================
# PERCEPTIONx extract_depmap_meta smoke tests
# =============================================================================

test_that("extract_depmap_meta reads metadata from a mock DepMap RDS", {
  rds <- tempfile(fileext = ".RDS")
  saveRDS(mock_depmap(), rds)
  on.exit(unlink(rds), add = TRUE)

  meta <- PERCEPTIONx:::extract_depmap_meta(rds)

  expect_true(meta$loaded)
  expect_type(meta$genes, "character")
  expect_true(length(meta$genes) == 60)
  expect_equal(meta$drugs, c("erlotinib", "gefitinib"))
  expect_equal(meta$lineages, "Lung")
  expect_true(all(c("expression_rnorm", "secondary_prism") %in% names(meta$components)))
})

test_that("extract_depmap_meta uses a fresh cache file", {
  rds <- tempfile(fileext = ".RDS")
  saveRDS(mock_depmap(), rds)
  cache <- tempfile(fileext = ".RDS")
  on.exit(unlink(c(rds, cache)), add = TRUE)

  meta <- PERCEPTIONx:::extract_depmap_meta(rds, cache_file = cache)
  expect_true(meta$loaded)
  expect_true(file.exists(cache))

  # Second call reads the (newer) cache instead of the RDS.
  meta2 <- PERCEPTIONx:::extract_depmap_meta(rds, cache_file = cache)
  expect_equal(meta2, meta)
})

test_that("extract_depmap_meta errors on a non-list RDS", {
  rds <- tempfile(fileext = ".RDS")
  saveRDS(matrix(1:4, 2), rds)
  on.exit(unlink(rds), add = TRUE)

  expect_error(PERCEPTIONx:::extract_depmap_meta(rds), "not a list")
})
