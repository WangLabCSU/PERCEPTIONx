# =============================================================================
# PERCEPTIONx evaluate functions smoke tests
# -----------------------------------------------------------------------------
# compare_performance() and get_significant_models() only need a list of model
# objects (never DepMap), so we feed them hand-built minimal models.
# =============================================================================

# Minimal model object shaped like build_on_BULK_v2() output.
fake_model <- function(cv = 0.5,
                       bulk = c(p.value = 0.1, estimate.cor = 0.5),
                       pseudo = c(p.value = 0.1, estimate.cor = 0.4),
                       sc = c(p.value = 0.1, estimate.cor = 0.3)) {
  list(
    model_performance_during_cv = cv,
    performance_in_bulk         = bulk,
    performance_in_pseudo_bulk  = pseudo,
    performance_in_scRNA        = sc
  )
}

test_that("compare_performance returns per-dataset metrics", {
  models <- list(
    drugA = fake_model(cv = 0.6, sc = c(p.value = 0.01, estimate.cor = 0.8)),
    drugB = fake_model(cv = 0.2, sc = c(p.value = 0.9, estimate.cor = 0.05))
  )

  perf <- compare_performance(models, verbose = FALSE)

  expect_type(perf, "list")
  expect_true(all(c("perf_cv", "perf_scRNA", "summary") %in% names(perf)))
  expect_equal(nrow(perf$perf_scRNA), 2)
  # Row names follow the model list names.
  expect_equal(rownames(perf$perf_scRNA), c("drugA", "drugB"))
})

test_that("get_significant_models keeps models above thresholds", {
  models <- list(
    good = fake_model(cv = 0.6, sc = c(p.value = 0.01, estimate.cor = 0.8)),
    bad  = fake_model(cv = 0.2, sc = c(p.value = 0.9, estimate.cor = 0.05))
  )

  sig <- get_significant_models(models, min_correlation = 0.3, max_pvalue = 0.05)

  expect_type(sig, "list")
  expect_equal(names(sig), "good")
})
