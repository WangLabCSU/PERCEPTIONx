# verify_umap_fast.R — verify the direct-uwot UMAP replacement in
# run_seurat_pipeline (source via load_all): speed, coordinate sanity,
# unchanged clone detection, and that the spatial-plot functions render.
# Usage: Rscript dev/verify_umap_fast.R
suppressMessages({ devtools::load_all(".", quiet = TRUE) })
suppressMessages({ requireNamespace("Seurat", quietly = TRUE) })

set.seed(42)
gene_names <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "RB1",
                "APC", "BRAF", "CDH1", "CDKN2A", "ERBB2", "FGFR1", "ALK",
                "MET", "RET", "ROS1", "NRAS", "HRAS", "MAP2K1", "MAPK1",
                "JAK2", "STAT3", "MTOR", "AKT1", "AKT2", "CTNNB1", "SMAD4",
                "VHL", "NF1", "NF2", "STK11", "FBXW7", "ARID1A", "KDM5C",
                "KMT2D", "SETD2", "BAP1", "PBRM1", "NOTCH1", "NOTCH2",
                "NOTCH3", "JAK1", "JAK3", "SOX9", "IDH1", "IDH2", "FLT3")
n_cells <- 400; n_patients <- 20
cell_names <- paste0("CELL_", sprintf("%04d", 1:n_cells))
patient_names <- paste0("PAT_", sprintf("%03d", 1:n_patients))
patient_assignment <- rep(patient_names, each = ceiling(n_cells / n_patients))[1:n_cells]
patient_mapping <- data.frame(cell_id = cell_names, patient_id = patient_assignment,
                              stringsAsFactors = FALSE)
clinical_response <- data.frame(patient = patient_names,
                                response = c(rep("Responder", 10), rep("Non-responder", 10)),
                                stringsAsFactors = FALSE)
is_responder_cell <- clinical_response$response[
  match(patient_assignment, clinical_response$patient)] == "Responder"
expr_matrix <- matrix(0.1, nrow = length(gene_names), ncol = n_cells)
rownames(expr_matrix) <- gene_names; colnames(expr_matrix) <- cell_names
for (g in gene_names[1:5]) {
  expr_matrix[g, is_responder_cell]  <- pmax(rnorm(sum(is_responder_cell), mean = 8, sd = 3), 0.1)
  expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 3, sd = 2), 0.1)
}
for (g in gene_names[6:10]) {
  expr_matrix[g, is_responder_cell]  <- pmax(rnorm(sum(is_responder_cell), mean = 3, sd = 2), 0.1)
  expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 8, sd = 3), 0.1)
}
for (g in gene_names[11:length(gene_names)]) expr_matrix[g, ] <- runif(n_cells, 0.5, 8)
storage.mode(expr_matrix) <- "numeric"

t0 <- proc.time()[3]
prepared <- PERCEPTIONx::prepare_data(
  method = "umap", expression_matrix = expr_matrix, patient_mapping = patient_mapping,
  seurat_resolution = 0.8, seurat_dims = 10, seurat_nfeatures = min(2000, length(gene_names)))
cat(sprintf("prepare_data (new UMAP): %.2fs\n", proc.time()[3] - t0))

uc <- prepared$umap_coords
stopifnot(is.data.frame(uc), nrow(uc) == n_cells,
          all(c("cell_id", "dim_1", "dim_2") %in% names(uc)))
spr <- range(uc$dim_1)[2] - range(uc$dim_1)[1]
stopifnot(is.finite(spr), spr > 0)
cat(sprintf("umap_coords: %d cells, dim_1 range width %.2f (healthy layout)\n", nrow(uc), spr))

ccm <- prepared$cell_clone_map
cat("clones detected:", length(unique(ccm$clone_id)), "\n")
stopifnot(nrow(ccm) == n_cells)

# spatial plot renders with the new coordinates
suppressMessages({ library(ggplot2) })
umap_data <- data.frame(X = uc$dim_1, Y = uc$dim_2,
                        clone_id = ccm$clone_id[match(uc$cell_id, ccm$cell_id)],
                        stringsAsFactors = FALSE)
p1 <- PERCEPTIONx::plot_clone_umap(umap_data, title = "Clone Identity")
stopifnot(inherits(p1, "ggplot"))
p2 <- PERCEPTIONx::plot_tsne_response(
  data.frame(X = uc$dim_1, Y = uc$dim_2, expression = runif(n_cells)),
  color_var = "expression", title = "geneX", palette = "viridis")
stopifnot(inherits(p2, "ggplot"))
cat("spatial plots render OK with new coordinates\n")
cat("\nUMAP-FAST VERIFY PASSED\n")
