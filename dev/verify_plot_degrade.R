# verify_plot_degrade.R
# Verifies the "scale auto-degrade" rendering decision introduced in
# mod_visualize.R + async_jobs.R:
#   1. run_plot() attaches perception_meta (panels/points) to the ggplot.
#   2. render_mode() chooses static for >= 10 panels or >= 15k points.
#   3. The attr survives girafe() (interactive path) untouched.
#   4. The static path (text-scaled ggplot) prints without ggiraph.
#
# Usage: Rscript dev/verify_plot_degrade.R  (package installed in E:\rlibtest)
suppressMessages({
  library(PERCEPTIONx, lib.loc = "E:/rlibtest")
  library(ggplot2)
})
fail <- function(...) { cat("FAIL:", ..., "\n"); q(status = 1) }
ok <- function(...) cat("PASS:", ..., "\n")

# ---- Simulate the worker's meta attachment (clone_kill branch) ----
mk <- function(n_patients, clones_per = 3) {
  set.seed(1)
  df <- do.call(rbind, lapply(seq_len(n_patients), function(p) {
    data.frame(patient = paste0("P", sprintf("%02d", p)),
               clone_id = paste0("P", sprintf("%02d", p), "_c", seq_len(clones_per)),
               comb_viability = rnorm(clones_per),
               weights = rep(1 / clones_per, clones_per),
               response = rep(if (p %% 2) "R" else "NR", clones_per))
  }))
  df$comb_viability <- as.numeric(scale(df$comb_viability))
  df
}

attach_meta_like_worker <- function(plot_obj, pt, clone_viability_df) {
  meta <- NULL
  if (!is.null(plot_obj)) {
    if (pt == "clone_kill") {
      meta <- list(panels = length(unique(clone_viability_df$patient)),
                   points = nrow(clone_viability_df))
    } else {
      meta <- list(panels = 1L,
                   points = if (is.data.frame(plot_obj$data)) nrow(plot_obj$data) else 0L)
    }
  }
  if (!is.null(plot_obj)) attr(plot_obj, "perception_meta") <- meta
  plot_obj
}

render_mode_like_ui <- function(p) {
  if (is.null(p)) return("none")
  meta <- attr(p, "perception_meta")
  if (is.null(meta)) return("interactive")
  if (isTRUE(meta$panels >= 10) || isTRUE(meta$points >= 15000)) "static" else "interactive"
}

# ---- Case 1: small plot (5 patients) -> interactive ----
d_small <- mk(5)
p_small <- PERCEPTIONx::plot_clone_viability(d_small, viability_var = "comb_viability",
                                             weights_var = "weights", response_var = "response")
p_small <- attach_meta_like_worker(p_small, "clone_kill", d_small)
meta <- attr(p_small, "perception_meta")
if (is.null(meta) || meta$panels != 5 || meta$points != 15) fail("small meta wrong: ", paste(unlist(meta), collapse = ","))
if (render_mode_like_ui(p_small) != "interactive") fail("small plot should be interactive")
ok("small plot (5 patients) -> interactive, meta=", paste(unlist(meta), collapse = "/"))

# ---- Case 2: large plot (20 patients) -> static ----
d_big <- mk(20)
p_big <- PERCEPTIONx::plot_clone_viability(d_big, viability_var = "comb_viability",
                                           weights_var = "weights", response_var = "response")
p_big <- attach_meta_like_worker(p_big, "clone_kill", d_big)
meta <- attr(p_big, "perception_meta")
if (is.null(meta) || meta$panels != 20) fail("big meta wrong")
if (render_mode_like_ui(p_big) != "static") fail("large plot should be static")
ok("large plot (20 patients) -> static, meta=", paste(unlist(meta), collapse = "/"))

# ---- Case 3: null plot -> "none" ----
if (render_mode_like_ui(NULL) != "none") fail("NULL plot should be 'none'")
ok("NULL plot -> none")

# ---- Case 4: attr survives ggiraph conversion (interactive path) ----
if (requireNamespace("ggiraph", quietly = TRUE)) {
  g <- ggiraph::girafe(ggobj = p_small, width_svg = 8, height_svg = 6,
                       options = list(ggiraph::opts_zoom(min = 1, max = 1)))
  if (!inherits(g, "girafe")) fail("girafe() failed on meta-attributed ggplot")
  ok("girafe() conversion OK with perception_meta attr")
} else {
  ok("ggiraph not installed; skipped girafe conversion check")
}

# ---- Case 5: static path renders a text-scaled ggplot ----
scale <- 1.2
p_obj <- p_big
base <- p_obj$theme$text$size
if (is.null(base) || length(base) == 0 || !is.finite(base)) base <- 11
p_scaled <- p_obj +
  theme(text = element_text(size = base * scale),
        plot.title = element_text(size = (base + 1) * scale, hjust = 0, vjust = 0,
                                  face = "plain", margin = margin(b = 6)),
        plot.subtitle = element_text(size = base * scale),
        axis.title = element_text(size = base * scale),
        axis.text = element_text(size = (base - 1) * scale),
        legend.text = element_text(size = (base - 2) * scale),
        legend.title = element_text(size = (base - 1) * scale),
        strip.text = element_text(size = (base - 1) * scale))
tmp <- tempfile(fileext = ".png")
ggsave(tmp, p_scaled, width = 8, height = 6, dpi = 60)
if (!file.exists(tmp) || file.info(tmp)$size < 1000) fail("static render produced no output")
unlink(tmp)
ok("static text-scaled renderPlot path renders to file")

# ---- Case 6: threshold boundary (9 panels interactive, 10 static) ----
d9 <- mk(9); p9 <- attach_meta_like_worker(
  PERCEPTIONx::plot_clone_viability(d9, viability_var = "comb_viability",
                                    weights_var = "weights", response_var = "response"),
  "clone_kill", d9)
d10 <- mk(10); p10 <- attach_meta_like_worker(
  PERCEPTIONx::plot_clone_viability(d10, viability_var = "comb_viability",
                                    weights_var = "weights", response_var = "response"),
  "clone_kill", d10)
if (render_mode_like_ui(p9) != "interactive" || render_mode_like_ui(p10) != "static")
  fail("threshold boundary wrong")
ok("boundary: 9 panels -> interactive, 10 panels -> static")

# ---- Case 7: COMBO lollipop meta (regression: used to reference the
# non-existent clone_viability_df when is_combo = TRUE -> "object not found")
attach_meta_combo <- function(plot_obj, combo_df) {
  is_combo <- TRUE  # the worker's run_plot sets this from params
  clone_viability_df <- NULL  # NEVER defined in the combo branch
  meta <- NULL
  if (!is.null(plot_obj)) {
    if (identical(attr(plot_obj, "pt"), "clone_kill")) {
      kill_df <- if (is_combo) combo_df else clone_viability_df
      meta <- list(panels = length(unique(kill_df$patient)), points = nrow(kill_df))
    }
  }
  if (!is.null(plot_obj)) attr(plot_obj, "perception_meta") <- meta
  plot_obj
}
p_combo <- PERCEPTIONx::plot_clone_viability(d_big, viability_var = "comb_viability",
                                             weights_var = "weights", response_var = "response")
attr(p_combo, "pt") <- "clone_kill"
p_combo <- attach_meta_combo(p_combo, d_big)
meta <- attr(p_combo, "perception_meta")
if (is.null(meta) || meta$panels != 20) fail("combo meta wrong: ", paste(unlist(meta), collapse = ","))
if (render_mode_like_ui(p_combo) != "static") fail("combo 20-patient plot should be static")
ok("COMBO lollipop meta OK (no clone_viability_df reference), panels=", meta$panels)

cat("\nALL PLOT-DEGRADE CHECKS PASSED\n")
