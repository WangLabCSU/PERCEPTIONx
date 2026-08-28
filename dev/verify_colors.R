# verify_colors.R — verify the colour scales:
#   1. Gene expression (0-1) uses a grey -> red sequential ramp (0 = grey,
#      no blue at all).
#   2. Viability / lollipop use DATA-DRIVEN diverging limits — values outside
#      any fixed window (e.g. -2.2) still get a real colour, never grey/NA.
# Usage: Rscript dev/verify_colors.R
suppressMessages({ devtools::load_all(".", quiet = TRUE) })
suppressMessages({ library(ggplot2) })

set.seed(1)
xy <- data.frame(X = rnorm(50), Y = rnorm(50))

# 1. Gene Expression: light grey -> red, low = grey, high = red, no blue/white mid.
p_gene <- PERCEPTIONx::plot_tsne_response(
  data.frame(X = xy$X, Y = xy$Y, expression = runif(50)),
  color_var = "expression",
  colors = c("#e0e0e0", "#c13232"))
sc <- ggplot_build(p_gene)$plot$scales$get_scales("colour")
# gradientn stores the ramp as a palette() function over RESCALED positions
# (0 = low limit, 1 = high limit); sample those positions.
hex_at <- function(sc) {
  rgb(t(grDevices::col2rgb(sc$palette(c(0, 0.5, 1)))), maxColorValue = 255)
}
cols_g <- tolower(hex_at(sc))
cat("gene expr: colors at low/mid/high =", cols_g, "\n")
# endpoints are exact (light grey -> red); midpoint is interpolated in Lab
# space, so only require it stays a red-ish tint, never blue/white.
stopifnot(identical(cols_g[c(1, 3)], c("#e0e0e0", "#c13232")))
mid_rgb <- grDevices::col2rgb(cols_g[2])
stopifnot(mid_rgb[1] >= mid_rgb[3])  # red channel >= blue channel

# 1b. Constant-0 expression must map to the GREY end (not red).
p_zero <- PERCEPTIONx::plot_tsne_response(
  data.frame(X = xy$X, Y = xy$Y, expression = rep(0, 50)),
  color_var = "expression",
  colors = c("#e0e0e0", "#c13232"),
  limits = c(0, 1))
b0 <- ggplot_build(p_zero)
cat("gene expr all-0: point colour =", unique(b0$data[[1]]$colour), "\n")
stopifnot(identical(unique(b0$data[[1]]$colour), "#e0e0e0"))

# 2. Drug Viability: diverging, DATA-DRIVEN limits (not pinned to -2/2).
set.seed(2)
dat_viab <- data.frame(X = xy$X, Y = xy$Y,
                       viability_scaled = c(rnorm(49), -2.2))
p_viab <- PERCEPTIONx::plot_tsne_response(dat_viab, color_var = "viability_scaled",
                                          palette = "diverging", midpoint = 0)
sc2 <- ggplot_build(p_viab)$plot$scales$get_scales("colour")
stopifnot(is.null(sc2$limits))  # data-driven
b2 <- ggplot_build(p_viab)
stopifnot(!anyNA(b2$data[[1]]$colour))  # -2.2 keeps a real colour
cat("viability: data-driven limits OK; -2.2 maps to",
    b2$data[[1]]$colour[nrow(b2$data[[1]])], "\n")

# 3. Lollipop: diverging, DATA-DRIVEN limits; -2.2 keeps a colour.
kill <- data.frame(
  patient = rep(paste0("P", 1:5), each = 3),
  clone_id = paste0("P", rep(1:5, each = 3), "_c", 1:3),
  comb_viability = c(rnorm(14), -2.2), weights = 1/3,
  stringsAsFactors = FALSE)
p_kill <- PERCEPTIONx::plot_clone_viability(kill, viability_var = "comb_viability",
                                            weights_var = "weights")
sc3 <- ggplot_build(p_kill)$plot$scales$get_scales("colour")
stopifnot(is.null(sc3$limits))
b3 <- ggplot_build(p_kill)
# layer order: hline(1), segment(2), points(3) — points carry the mapping
stopifnot(!anyNA(b3$data[[3]]$colour))
cat("lollipop: data-driven limits OK; -2.2 maps to",
    b3$data[[3]]$colour[nrow(b3$data[[3]])], "\n")

# 4. viridis_scale = TRUE still available (no regression)
p_v <- PERCEPTIONx::plot_clone_viability(kill, viability_var = "comb_viability",
                                         weights_var = "weights", viridis_scale = TRUE)
sc4 <- ggplot_build(p_v)$plot$scales$get_scales("colour")
stopifnot(!is.null(sc4$palette) || inherits(sc4, "ScaleContinuous"))
cat("viridis_scale=TRUE still works\n")

cat("\nCOLOR VERIFY PASSED\n")
