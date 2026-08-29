# =============================================================================
# PERCEPTIONx load/mirror utilities smoke tests
# -----------------------------------------------------------------------------
# get_depmap() reads the mock DepMap from the global env; mirror helpers are
# pure functions. Download paths (load_model/load_depmap) are not exercised.
# =============================================================================

test_that("get_depmap returns the loaded DepMap object", {
  assign("DepMap", mock_depmap(), envir = .GlobalEnv)
  on.exit(rm("DepMap", envir = .GlobalEnv), add = TRUE)

  expect_true(is.list(get_depmap()))
  expect_true("expression_rnorm" %in% names(get_depmap()))
})

test_that("get_depmap errors when no DepMap is loaded", {
  # Guard: make sure nothing is in the global env already.
  was_loaded <- exists("DepMap", envir = .GlobalEnv)
  on.exit(if (was_loaded) assign("DepMap", get("DepMap", envir = .GlobalEnv), envir = .GlobalEnv),
          add = TRUE)
  rm("DepMap", envir = .GlobalEnv)

  expect_error(get_depmap(), "DepMap data not loaded")
})

test_that("mirror helpers return character vectors", {
  expect_type(get_mirrors(), "character")
  expect_true(length(get_mirrors()) > 0)

  mirrors <- c("https://mirror-a.com/https://github.com",
               "https://mirror-b.com/https://github.com")
  add_mirrors(mirrors)
  expect_true(all(mirrors %in% get_mirrors()))

  reset_mirrors()
  expect_false(any(mirrors %in% get_mirrors()))
})
