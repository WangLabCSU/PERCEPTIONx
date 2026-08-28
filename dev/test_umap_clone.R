# test_umap_clone.R — verifies plot_clone_umap works WITHOUT ggplot2 attached,
# using the source package via devtools::load_all (proves the new
# scale_colour_manual import in NAMESPACE/plot.R fixes the worker path).
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(7)
coords <- data.frame(cell_id = paste0("C", 1:10), dim_1 = rnorm(10), dim_2 = rnorm(10))
clone_data <- data.frame(cell_id = paste0("C", 1:10), clone_id = rep(c("cl_A", "cl_B"), 5))
# --- replicate the worker's umap_clone branch (async_jobs.R run_plot) ---
common_cells <- intersect(clone_data$cell_id, coords$cell_id)
if (length(common_cells) == 0) stop("No matching cells between embedding coordinates and clone annotation.")
idx <- match(common_cells, clone_data$cell_id)
x_col <- intersect(c("dim_1", "umap_1"), names(coords))[1]
y_col <- intersect(c("dim_2", "umap_2"), names(coords))[1]
X <- coords[[x_col]][match(common_cells, coords$cell_id)]
Y <- coords[[y_col]][match(common_cells, coords$cell_id)]
n <- length(common_cells)
stopifnot(length(X) == n, length(Y) == n, length(idx) == n)
umap_data <- data.frame(X = X, Y = Y, clone_id = clone_data$clone_id[idx], stringsAsFactors = FALSE)
p <- PERCEPTIONx::plot_clone_umap(umap_data, title = "Clone Identity")
stopifnot(inherits(p, "ggplot"))
cat("plot_clone_umap OK without attached ggplot2; ggplot2 on search path:",
    any(grepl("^package:ggplot2", search())), "\n")
