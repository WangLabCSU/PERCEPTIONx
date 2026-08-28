# get_site_images.R — one-off: download the tutorial screenshots from the
# external img.remit.ee hotlinks into man/figures/ so both GitHub and the
# pkgdown site render them reliably.
md <- readLines("vignettes/shiny_app.Rmd", warn = FALSE)
pat <- "https://img\\.remit\\.ee/api/file/[A-Za-z0-9_-]+\\.png"
urls <- unique(unlist(regmatches(md, gregexpr(pat, md, perl = TRUE))))
stopifnot(length(urls) >= 9)
names <- c(
  "shiny-home.png", "shiny-data-clustering.png", "shiny-predict-heatmap.png",
  "shiny-clone-distribution.png", "shiny-lollipop.png", "shiny-roc.png",
  "shiny-boxplot.png", "shiny-model-performance.png", "shiny-umap-gene.png",
  "shiny-umap-viability.png")
names <- names[seq_along(urls)]
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
for (i in seq_along(urls)) {
  dest <- file.path("man/figures", names[i])
  ok <- tryCatch({
    download.file(urls[i], dest, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  sz <- if (file.exists(dest)) file.size(dest) else 0
  cat(sprintf("%-34s %8d bytes  %s\n", names[i], sz, if (ok && sz > 10000) "OK" else "FAIL"))
}
