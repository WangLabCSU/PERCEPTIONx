# =============================================================================
# PERCEPTIONx additional plot-function smoke tests
# -----------------------------------------------------------------------------
# Covers the remaining exported plotting functions (and two small internal
# helpers) that the original test file did not reach.
# =============================================================================

# Minimal model shaped like build_on_BULK_v2() output, with enough rows for
# the validation-ROC computation (>= 6 rows per dataset).
fake_roc_model <- function(seed = 1) {
  set.seed(seed)
  mk <- function() {
    n <- 10
    obs <- rnorm(n)
    data.frame(Observed = obs, Pred = 0.5 * obs + rnorm(n))
  }
  list(predVSgroundTruth = list(
    pred_gt_bulk   = mk(),
    pred_gt_mscRNA = mk(),
    pred_gt_scRNA  = mk()
  ))
}

test_that("plot_model_roc builds a ROC plot from model objects", {
  models <- list(drugA = fake_roc_model(1), drugB = fake_roc_model(2))
  res <- plot_model_roc(models)
  expect_true(!is.null(res))
})

test_that("plot_clone_umap returns ggplot with tooltip handling", {
  tsne_data <- data.frame(
    X = rnorm(50), Y = rnorm(50),
    clone_id = sample(paste0("c", 1:4), 50, replace = TRUE)
  )
  p <- plot_clone_umap(tsne_data)
  expect_s3_class(p, "ggplot")

  expect_error(plot_clone_umap(data.frame(X = 1:5, Y = 1:5)),
               "must contain columns")
})

test_that("plot_seurat_clustering returns embedding plot and cluster ids", {
  skip_if_not_installed("Seurat")
  set.seed(7)
  genes <- sprintf("GENE%03d", seq_len(60))
  cells <- paste0("C", seq_len(30))
  expr <- matrix(rpois(60 * 30, lambda = 5), nrow = 60,
                 dimnames = list(genes, cells))

  res <- plot_seurat_clustering("umap", expr,
                                min_cells = 3, min_features = 5,
                                nfeatures = 60, dims = 5)
  expect_type(res, "list")
  expect_true(all(c("seurat_object", "embedding_plot", "cluster_ids") %in% names(res)))
})

test_that("plot_patient_response_panel arranges all panels", {
  clone_dist <- data.frame(
    patients = c("P1", "P1", "P2", "P2"),
    clones   = c("c1", "c2", "c1", "c2"),
    weights  = c(0.6, 0.4, 0.7, 0.3),
    response = c("R", "R", "NR", "NR")
  )
  clone_viab <- data.frame(
    patient = c("P1", "P1", "P2", "P2"),
    clone_id = c("P1_c1", "P1_c2", "P2_c1", "P2_c2"),
    comb_viability = c(-1.2, -0.5, -2.1, -0.3),
    weights = c(0.6, 0.4, 0.7, 0.3),
    response = c("R", "R", "NR", "NR")
  )
  exp_pred <- data.frame(
    response = factor(c("R", "R", "NR", "NR", "R", "NR")),
    predicted_viability = c(-2.1, -1.5, -0.3, -0.1, -1.8, -0.5)
  )

  panel <- plot_patient_response_panel(clone_dist, clone_viab, exp_pred)
  expect_true(!is.null(panel))
})

test_that("plot_tsne_biomarker_viability arranges two UMAP panels", {
  tsne_data <- data.frame(
    X = rnorm(50), Y = rnorm(50),
    biomarker_scaled = runif(50),
    viability_scaled = runif(50)
  )
  res <- plot_tsne_biomarker_viability(tsne_data)
  expect_true(!is.null(res))
})

test_that("internal response tag and p-value helpers behave", {
  expect_equal(PERCEPTIONx:::resp_short_tag("Responder"), "R")
  expect_equal(PERCEPTIONx:::resp_short_tag("RESISTANT"), "NR")
  expect_equal(PERCEPTIONx:::resp_short_tag("PD"), "PD")

  expect_equal(PERCEPTIONx:::fmt_pval(0.5), "0.500")
  expect_equal(PERCEPTIONx:::fmt_pval(NA), "NA")
})

test_that("export_plot_cairo writes a file", {
  p <- ggplot2::ggplot(data.frame(x = 1:5, y = 1:5),
                       ggplot2::aes(x, y)) + ggplot2::geom_point()
  f <- tempfile(fileext = ".png")
  export_plot_cairo(f, p, format = "png")
  expect_true(file.exists(f))
  expect_true(file.info(f)$size > 0)
  unlink(f)
})
