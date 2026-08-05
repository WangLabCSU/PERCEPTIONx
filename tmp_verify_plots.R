# Full plot-function verification for PERCEPTION
suppressMessages(devtools::load_all("c:/Users/Lenovo/Desktop/PERCEPTION", quiet = TRUE))
suppressMessages(library(ggplot2))

ok <- function(name, expr) {
  res <- tryCatch({
    expr
    cat(sprintf("OK   %s\n", name))
    TRUE
  }, error = function(e) {
    cat(sprintf("FAIL %s : %s\n", name, conditionMessage(e)))
    FALSE
  })
  invisible(res)
}

# --- 1. plot_clone_distribution ---
cd <- data.frame(
  patients = c("P1","P1","P2","P2","P3","P3","P4","P4"),
  clones   = c("c1","c2","c1","c2","c1","c2","c1","c2"),
  weights  = c(.6,.4,.3,.7,.5,.5,.8,.2),
  response = c("R","R","NR","NR","R","R","NR","NR"),
  stringsAsFactors = FALSE
)
ok("DIST-4-BASIC", {
  p <- PERCEPTION::plot_clone_distribution(cd)
  stopifnot(inherits(p, "ggplot")); ggplot_build(p)
})
ok("DIST-4-FACET", {
  p <- PERCEPTION::plot_clone_distribution(cd, response_var = "response")
  stopifnot(inherits(p, "ggplot")); ggplot_build(p)
})
cd14 <- do.call(rbind, lapply(1:14, function(i)
  data.frame(patients = sprintf("PAT_%02d", i),
             clones = c("c1","c2"), weights = c(.5,.5),
             response = rep(c("R","NR"), 7)[i], stringsAsFactors = FALSE)))
ok("DIST-14-FACET", {
  p <- PERCEPTION::plot_clone_distribution(cd14, response_var = "response")
  stopifnot(inherits(p, "ggplot")); ggplot_build(p)
})

# --- 2. plot_clone_killing (lollipop) ---
ck <- data.frame(
  patient = c("P1","P1","P2","P2","P3","P3"),
  clone_id = c("P1_c1","P1_c2","P2_c1","P2_c2","P3_c1","P3_c2"),
  comb_killing = c(-1.5,-0.8,-2.1,-0.5,-1.2,0.3),
  weights = c(.6,.4,.3,.7,.5,.5),
  response = c("R","R","NR","NR","R","R"),
  stringsAsFactors = FALSE
)
ok("LOLLIPOP-6-SINGLE-ROW", {
  p <- PERCEPTION::plot_clone_killing(ck, killing_var = "comb_killing",
                                      weights_var = "weights", response_var = "response",
                                      drug = "abemaciclib")
  ggplot_build(p)
  f <- file.path(tempdir(), "lolli_6.png")
  grDevices::png(f, width = 9, height = 4, units = "in", res = 96)
  print(p); grDevices::dev.off()
  stopifnot(file.exists(f) && file.info(f)$size > 0)
})
ck20 <- do.call(rbind, lapply(1:20, function(i)
  data.frame(patient = sprintf("P%02d", i), clone_id = sprintf("P%02d_c1", i),
             comb_killing = runif(1, -2, 1), stringsAsFactors = FALSE)))
ok("LOLLIPOP-20-WRAP", {
  p <- PERCEPTION::plot_clone_killing(ck20, killing_var = "comb_killing")
  ggplot_build(p)
})
ok("LOLLIPOP-VIRIDIS-FALSE", {
  p <- PERCEPTION::plot_clone_killing(ck, killing_var = "comb_killing",
                                      viridis_scale = FALSE)
  ggplot_build(p)
})

# --- 3. plot_response_boxplot ---
epv <- data.frame(
  response = factor(c("R","R","NR","NR","R","NR","R","NR")),
  predicted_killing = c(-2.1,-1.5,-0.3,-0.1,-1.8,-0.5,-2.2,-0.2)
)
ok("BOXPLOT", {
  p <- PERCEPTION::plot_response_boxplot(epv)
  ggplot_build(p)
  f <- file.path(tempdir(), "boxplot.png")
  grDevices::png(f, width = 5, height = 5, units = "in", res = 96)
  print(p); grDevices::dev.off()
  stopifnot(file.exists(f) && file.info(f)$size > 0)
})

# --- 4. plot_roc_curve ---
set.seed(123); n <- 60
resp <- factor(sample(c("R","NR"), n, replace = TRUE), levels = c("R","NR"))
predv <- rnorm(n, mean = ifelse(resp == "R", -1.5, -0.3))
ok("ROC-UNSMOOTHED", {
  p <- PERCEPTION::plot_roc_curve(response = resp, predictor = predv, smooth_curve = FALSE)
  ggplot_build(p)
})
ok("ROC-SMOOTHED", {
  p <- PERCEPTION::plot_roc_curve(response = resp, predictor = predv, smooth_curve = TRUE)
  ggplot_build(p)
})

# --- 5. plot_tsne_response ---
ts <- data.frame(X = rnorm(300), Y = rnorm(300), killing_scaled = runif(300))
ok("TSNE-VIRIDIS", {
  p <- PERCEPTION::plot_tsne_response(ts, color_var = "killing_scaled")
  ggplot_build(p)
})
ok("TSNE-DIVERGING", {
  p <- PERCEPTION::plot_tsne_response(ts, color_var = "killing_scaled", palette = "diverging")
  ggplot_build(p)
})
ts_big <- data.frame(X = rnorm(60000), Y = rnorm(60000), killing_scaled = runif(60000))
ok("TSNE-60K-DOWNSAMPLE", {
  p <- PERCEPTION::plot_tsne_response(ts_big, color_var = "killing_scaled")
  stopifnot(nrow(ggplot_build(p)$data[[1]]) <= 50000)
})

# --- 6. plot_model_performance (fake perf list) ---
make_perf <- function(r) list(performance_in_scRNA = data.frame(estimate.cor = r),
                              performance_in_bulk = data.frame(estimate.cor = r * 0.9),
                              performance_in_pseudo_bulk = data.frame(estimate.cor = r * 0.8))
perf_list <- setNames(lapply(c(0.1, 0.3, 0.45, 0.55, 0.7), make_perf), paste0("drug", 1:5))
ok("MODEL-PERF", {
  p <- PERCEPTION::plot_model_performance(perf_list)
  ggplot_build(p)
})

# --- 7. export_plot_cairo ---
ok("EXPORT-CAIRO-PNG", {
  p <- PERCEPTION::plot_clone_distribution(cd, response_var = "response")
  f <- file.path(tempdir(), "t_export.png")
  PERCEPTION::export_plot_cairo(file = f, plot = p, format = "png", res = 96, width = 8, height = 6)
  stopifnot(file.exists(f) && file.info(f)$size > 0)
})
ok("EXPORT-CAIRO-PDF", {
  p <- PERCEPTION::plot_clone_killing(ck, killing_var = "comb_killing")
  f <- file.path(tempdir(), "t_export.pdf")
  PERCEPTION::export_plot_cairo(file = f, plot = p, format = "pdf", width = 8, height = 4)
  stopifnot(file.exists(f) && file.info(f)$size > 0)
})

cat("=== ALL PLOT CHECKS DONE ===\n")
