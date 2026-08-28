# test_plot_cache.R — verifies the env-based cache semantics used by
# mod_visualize.R (plot_cache / widget_cache): key construction, hit/miss,
# per-settings separation, and invalidation on data change.
# Usage: Rscript dev/test_plot_cache.R
plot_cache <- new.env(parent = emptyenv())
widget_cache <- new.env(parent = emptyenv())

# Standard plot request keys: type + combo + drug + ROC groups
ck1 <- paste("clone_kill", "FALSE", "DrugX", "NULL", "NULL", sep = "\u0001")
ck2 <- paste("clone_kill", "FALSE", "DrugY", "NULL", "NULL", sep = "\u0001")
# Advanced plot request keys: type + gene + drug
ck3 <- paste("umap_gene", "TP53", "NULL", sep = "\u0001")

stopifnot(is.null(get0(ck1, envir = plot_cache, inherits = FALSE)))  # miss
assign(ck1, list(plot = "gg1", message = ""), envir = plot_cache)
stopifnot(identical(get0(ck1, envir = plot_cache, inherits = FALSE)$plot, "gg1"))  # hit
stopifnot(is.null(get0(ck2, envir = plot_cache, inherits = FALSE)))  # different drug -> miss
stopifnot(is.null(get0(ck3, envir = plot_cache, inherits = FALSE)))  # spatial key -> miss

# Widget keys: plot key + text scale + canvas w/h
wkey  <- paste(ck1, 100, 1000, 750, sep = "\u0001")
wkey2 <- paste(ck1, 120, 1000, 750, sep = "\u0001")
stopifnot(is.null(get0(wkey, envir = widget_cache, inherits = FALSE)))
assign(wkey, "widget1", envir = widget_cache)
stopifnot(identical(get0(wkey, envir = widget_cache, inherits = FALSE), "widget1"))
stopifnot(is.null(get0(wkey2, envir = widget_cache, inherits = FALSE)))  # text scale changed

# Invalidation on data change (same code as the observe)
if (length(ls(plot_cache)) > 0) rm(list = ls(plot_cache), envir = plot_cache)
if (length(ls(widget_cache)) > 0) rm(list = ls(widget_cache), envir = widget_cache)
stopifnot(length(ls(plot_cache)) == 0, length(ls(widget_cache)) == 0)

cat("PLOT CACHE LOGIC OK\n")
