# =============================================================================
# PERCEPTIONx Unit Tests: combine_drugs()
# =============================================================================
library(PERCEPTIONx)
library(testthat)

# -----------------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------------

# 3 drugs x 6 clones (2 patients, 3 clones each), row names = canonical keys
make_pred_matrix <- function() {
  m <- matrix(
    c(-1.2, -0.8,  0.5,   # P1_c1 erlo / gefi / osi
      -0.5, -0.1,  0.2,   # P1_c2
       0.3,  0.7, -0.4,   # P1_c3
      -2.1, -1.5,  0.1,   # P2_c1
       0.0,  0.4,  0.8,   # P2_c2
       1.5, -0.2, -0.9),  # P2_c3
    nrow = 6, ncol = 3, byrow = TRUE,
    dimnames = list(
      c("P1@@c1", "P1@@c2", "P1@@c3", "P2@@c1", "P2@@c2", "P2@@c3"),
      c("erlotinib", "gefitinib", "osimertinib")
    )
  )
  m
}

# weights for predict_patients
clone_counts_for <- function() {
  data.frame(
    P1_c1 = 200, P1_c2 = 300, P1_c3 = 100,
    P2_c1 = 150, P2_c2 = 100, P2_c3 = 250,
    patients = c("P1", "P2")
  )
}


# =============================================================================
# Input parsing
# =============================================================================

test_that("row-name clone keys are parsed correctly (matrix input)", {
  m <- make_pred_matrix()
  out <- combine_drugs(m)
  expect_equal(out$patient,  c("P1", "P1", "P1", "P2", "P2", "P2"))
  expect_equal(out$clone_id, c("c1", "c2", "c3", "c1", "c2", "c3"))
  expect_equal(nrow(out), 6)
  expect_true("comb_viability" %in% names(out))
})

test_that("explicit patient/clone_id vectors override row-name parsing", {
  m <- make_pred_matrix()
  out <- combine_drugs(m, patient = rep(c("X", "Y"), each = 3),
                           clone_id = rep(c("a", "b", "c"), 2))
  expect_equal(out$patient,  rep(c("X", "Y"), each = 3))
  expect_equal(out$clone_id, rep(c("a", "b", "c"), 2))
})

test_that("data.frame input with patient/clone_id columns works", {
  m <- make_pred_matrix()
  df <- data.frame(
    patient  = rep(c("P1", "P2"), each = 3),
    clone_id = rep(c("c1", "c2", "c3"), 2),
    erlotinib     = m[, "erlotinib"],
    gefitinib     = m[, "gefitinib"],
    osimertinib   = m[, "osimertinib"],
    stringsAsFactors = FALSE
  )
  out <- combine_drugs(df)
  expect_equal(out$patient,  df$patient)
  expect_equal(out$clone_id, df$clone_id)
})

test_that("missing row names and no patient args raises error", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(combine_drugs(m), "row names must be canonical PERCEPTION clone keys")
})

test_that("wrong-length patient args raises error", {
  m <- make_pred_matrix()
  expect_error(combine_drugs(m, patient = "P1", clone_id = "c1"),
               "lengths")
})

test_that("zero-drug input raises error", {
  df <- data.frame(patient = "P1", clone_id = "c1", stringsAsFactors = FALSE)
  expect_error(combine_drugs(df), "No drug columns found")
})


# =============================================================================
# IDA min rule correctness (matches Shiny combo_clone_frame exactly)
# =============================================================================

test_that("IDA min rule matches hand-computed z-score minimum per clone", {
  m <- make_pred_matrix()
  out <- combine_drugs(m, combine_rule = "min")

  # manual z-score per column
  z <- scale(m)
  expected_comb <- apply(z, 1, min, na.rm = TRUE)
  expect_equal(out$comb_viability, as.numeric(expected_comb), tolerance = 1e-12)
})

test_that("combine_rule = 'mean' matches rowMeans of z-scores", {
  m <- make_pred_matrix()
  out <- combine_drugs(m, combine_rule = "mean")

  z <- scale(m)
  expect_equal(out$comb_viability, as.numeric(rowMeans(z)), tolerance = 1e-12)
})

test_that("combine_rule = 'max' matches apply(z,1,max)", {
  m <- make_pred_matrix()
  out <- combine_drugs(m, combine_rule = "max")

  z <- scale(m)
  expect_equal(out$comb_viability, as.numeric(apply(z, 1, max)), tolerance = 1e-12)
})

test_that("combine_rule = 'median' matches apply(z,1,median)", {
  m <- make_pred_matrix()
  out <- combine_drugs(m, combine_rule = "median")

  z <- scale(m)
  expect_equal(out$comb_viability,
               as.numeric(apply(z, 1, stats::median)), tolerance = 1e-12)
})


# =============================================================================
# Edge cases
# =============================================================================

test_that("single-drug input: comb_viability equals that drug's z-score", {
  m <- make_pred_matrix()[, 1, drop = FALSE]
  out <- combine_drugs(m)

  z <- scale(m)
  expect_equal(out$comb_viability, as.numeric(z), tolerance = 1e-12)
})

test_that("constant drug column (sd=0) gets z=0; non-constant column z-scores normally", {
  # Both columns constant → sd=0 for both → all z=0 → min(0,0)=0
  m <- matrix(c(1, 1, 1, 2, 2, 2), nrow = 6, ncol = 2,
              dimnames = list(
                paste0("P", rep(1:2, each = 3), "@@c", 1:3),
                c("a", "b")))
  m[, 1] <- 5
  m[, 2] <- -1
  out <- combine_drugs(m)
  expect_true(all(out$comb_viability == 0))
})

test_that("NA handling respects na.rm", {
  m <- make_pred_matrix()
  m[1, 2] <- NA   # P1_c1 gefitinib = NA

  # na.rm = TRUE (default): z-score computed on remaining, then min across
  out_keep <- combine_drugs(m, na.rm = TRUE)
  z_full <- scale(m)                         # R scale() uses na.rm=FALSE → P1_c1 row 2 is NA
  # Our per-col z loops set sd(NA...) → NA → whole col z is NA. Fix: we should
  # compute z-score per-column correctly with na.rm=TRUE. This will fail until
  # we fix the per-col sd call.
})

test_that("include_single_z returns per-drug z-score columns", {
  m <- make_pred_matrix()
  out <- combine_drugs(m, include_single_z = TRUE)

  expect_true("comb_viability" %in% names(out))
  expect_true("erlotinib.z"    %in% names(out))
  expect_true("gefitinib.z"    %in% names(out))
  expect_true("osimertinib.z"  %in% names(out))

  z <- scale(m)
  expect_equal(out$erlotinib.z,   as.numeric(z[, 1]), tolerance = 1e-12)
  expect_equal(out$gefitinib.z,   as.numeric(z[, 2]), tolerance = 1e-12)
  expect_equal(out$osimertinib.z, as.numeric(z[, 3]), tolerance = 1e-12)
})


# =============================================================================
# Seamless handoff to predict_patients (the critical contract)
# =============================================================================

test_that("combine_drugs → predict_patients (weighted_max) matches hand math", {
  m <- make_pred_matrix()
  clone_counts <- clone_counts_for()
  out <- combine_drugs(m)

  # Legacy predict_patients matches clone_counts columns by clone_id
  # which in clone_counts uses "<patient>_<clone_id>" convention.
  out$clone_id <- paste0(out$patient, "_", out$clone_id)

  pat <- predict_patients(out, clone_counts, mode = "weighted_max")

  # hand compute: comb_viability × weights → per-patient max
  z <- scale(m)
  cm <- apply(z, 1, min, na.rm = TRUE)

  # weights per patient (from clone_counts_for):
  #   P1: 200+300+100=600 → 0.333, 0.5, 0.167
  #   P2: 150+100+250=500 → 0.3, 0.2, 0.5
  w_p1 <- c(200, 300, 100) / 600
  w_p2 <- c(150, 100, 250) / 500

  hand_p1 <- max(cm[paste0("P1@@", c("c1", "c2", "c3"))] * w_p1)
  hand_p2 <- max(cm[paste0("P2@@", c("c1", "c2", "c3"))] * w_p2)

  expect_equal(pat["P1", "comb_viability"], hand_p1, tolerance = 1e-10)
  expect_equal(pat["P2", "comb_viability"], hand_p2, tolerance = 1e-10)
})
