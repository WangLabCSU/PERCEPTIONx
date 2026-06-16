# =============================================================================
# PERCEPTION Test Pipeline
# =============================================================================
# End-to-end test script: Train -> Predict (cell + patient) -> Plot (all)
#
# Usage:
#   1. Make sure PERCEPTION package is loaded (devtools::load_all())
#   2. Modify parameters below as needed
#   3. source("test_pipeline.R")
#
# NOTE: This script downloads DepMap data (~567MB) and is NOT run by
#       testthat::test_dir(). Run it manually only when needed.
# =============================================================================

library(PERCEPTION)

# -- Parameter Settings -------------------------------------------------------

# Drug list (2 drugs for quick testing)
drug_list <- c("erlotinib", "gefitinib")

# Cancer type
cancer_type   <- "PanCan"
exclude_cancer <- "PanCan"

# Feature genes (NULL = use all DepMap genes; or provide a custom gene list)
GOI <- NULL

# Model hyperparameters
model_type     <- "glmnet"
k_features_values <- NULL  # NULL = auto-compute
ncores         <- 1

# Output directory
output_dir <- "./test_pipeline_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -- Step 0: Load DepMap Data -------------------------------------------------

message("\n", strrep("=", 60))
message("Step 0: Loading DepMap data")
message(strrep("=", 60))

# Use mirror for faster download
load_depmap(read = TRUE, mirror = TRUE)

# -- Step 1: Train Models -----------------------------------------------------

message("\n", strrep("=", 60))
message("Step 1: Training PERCEPTION models")
message(strrep("=", 60))

models <- train_models(
  drug_list        = drug_list,
  cancer_type      = cancer_type,
  exclude_cancer   = exclude_cancer,
  GOI              = GOI,
  k_features_values = k_features_values,
  model_type       = model_type,
  ncores           = ncores,
  output_dir       = output_dir
)

# -- Step 2: Evaluate Model Performance ---------------------------------------

message("\n", strrep("=", 60))
message("Step 2: Evaluating model performance")
message(strrep("=", 60))

perf <- compare_performance(models, threshold = 0.3, verbose = TRUE)

# Plot: model performance (number of drugs passing at each threshold)
pdf(file.path(output_dir, "01_model_performance.pdf"), width = 6, height = 5)
p1 <- plot_model_performance(models)
print(p1)
dev.off()
message("Saved: 01_model_performance.pdf")

# Filter significant models
sig_models <- get_significant_models(models, min_correlation = 0.3, max_pvalue = 0.05)

# -- Step 3: Cell-Level Prediction --------------------------------------------

message("\n", strrep("=", 60))
message("Step 3: Predicting drug response at cell level")
message(strrep("=", 60))

# Use excluded cell lines from DepMap (single-cell data) for prediction
cellLines_test <- get_cellLine_list(
  infunc_cancerType = cancer_type,
  infunc_drugName   = drug_list[1],
  exclude_cancer    = exclude_cancer,
  infunc_response   = get_response_matrix(drug_list[1])
)[[2]]

# Get expression matrix for test cell lines
test_cells <- DepMap$metadata_CPM_scRNA$NAME[
  DepMap$metadata_CPM_scRNA$DepMap_ID %in% cellLines_test]
expr_test <- DepMap$CPM_scRNA_CCLE_rnorm[, test_cells, drop = FALSE]

# Predict
cell_pred <- predict_drugs(models, expr_test)
message("Cell-level predictions: ", nrow(cell_pred), " cells x ", ncol(cell_pred), " drugs")
print(head(cell_pred))

# -- Step 4: Patient-Level Prediction (Simulated) -----------------------------

message("\n", strrep("=", 60))
message("Step 4: Predicting drug response at patient level (simulated)")
message(strrep("=", 60))

# Build simulated patient-clone structure
# Simulate 3 patients from test cell lines, each with 2-3 clones
set.seed(42)
patient_ids <- c("Patient_A", "Patient_B", "Patient_C")
clone_killing_list <- list()

for (pid in patient_ids) {
  n_clones <- sample(2:3, 1)
  clone_ids <- paste0(pid, "_clone", seq_len(n_clones))

  # Randomly select cells as clones
  cell_idx <- sample(seq_len(nrow(cell_pred)), n_clones, replace = TRUE)

  for (j in seq_len(n_clones)) {
    row_vals <- as.list(as.vector(cell_pred[cell_idx[j], ]))
    names(row_vals) <- colnames(cell_pred)
    clone_killing_list[[length(clone_killing_list) + 1]] <- data.frame(
      patient  = pid,
      clone_id = clone_ids[j],
      row_vals,
      check.names = FALSE
    )
  }
}

clone_killing_df <- do.call(rbind, clone_killing_list)

# Build clone_counts_df: unified columns for all clone IDs, fill 0 for missing
all_clone_ids <- unique(clone_killing_df$clone_id)
clone_counts_df <- as.data.frame(matrix(0L, nrow = length(patient_ids), ncol = length(all_clone_ids)))
colnames(clone_counts_df) <- all_clone_ids
rownames(clone_counts_df) <- patient_ids
for (i in seq_along(patient_ids)) {
  pid <- patient_ids[i]
  pid_clones <- clone_killing_df$clone_id[clone_killing_df$patient == pid]
  clone_counts_df[i, pid_clones] <- sample(50:500, length(pid_clones))
}
clone_counts_df$patients <- patient_ids

# Patient-level prediction
patient_pred <- predict_patients(clone_killing_df, clone_counts_df, mode = "weighted_average")
message("Patient-level predictions:")
print(patient_pred)

# -- Step 5: Plotting ---------------------------------------------------------

message("\n", strrep("=", 60))
message("Step 5: Generating all plots")
message(strrep("=", 60))

# -- 5a. Clone distribution --
counts_mat <- as.matrix(clone_counts_df[, all_clone_ids, drop = FALSE])
row_sums <- rowSums(counts_mat)
clone_dist_df <- data.frame(
  patients = clone_killing_df$patient,
  clones   = clone_killing_df$clone_id,
  weights  = counts_mat[
    cbind(match(clone_killing_df$patient, patient_ids),
          match(clone_killing_df$clone_id, all_clone_ids))
  ] / row_sums[match(clone_killing_df$patient, patient_ids)]
)

pdf(file.path(output_dir, "02_clone_distribution.pdf"), width = 8, height = 5)
p2 <- plot_clone_distribution(clone_dist_df)
print(p2)
dev.off()
message("Saved: 02_clone_distribution.pdf")

# -- 5b. Clone killing (lollipop) --
# Rename prediction column to comb_killing (use first drug)
drug_col <- colnames(clone_killing_df)[3]  # first drug column
clone_kill_plot_df <- clone_killing_df
clone_kill_plot_df$comb_killing <- clone_kill_plot_df[[drug_col]]
clone_kill_plot_df$weights <- clone_dist_df$weights

pdf(file.path(output_dir, "03_clone_killing.pdf"), width = 20, height = 6)
p3 <- plot_clone_killing(clone_kill_plot_df, killing_var = "comb_killing", weights_var = "weights")
print(p3)
dev.off()
message("Saved: 03_clone_killing.pdf")

# -- 5c. Response boxplot + ROC curve (using model's built-in pred vs ground truth) --
for (drug_name in names(models)) {
  model_obj <- models[[drug_name]]
  if (is.null(model_obj) || (length(model_obj) == 1 && is.na(model_obj))) next

  pred_gt <- model_obj$predVSgroundTruth$pred_gt_scRNA
  if (is.null(pred_gt) || nrow(pred_gt) == 0) next

  # Response boxplot: split by median observed value
  median_pred <- median(pred_gt$Test_pred_sc, na.rm = TRUE)
  exp_vs_pred <- data.frame(
    response = ifelse(pred_gt$Observed > median(pred_gt$Observed, na.rm = TRUE), "R", "NR"),
    predicted_killing = pred_gt$Test_pred_sc
  )
  exp_vs_pred$response <- factor(exp_vs_pred$response, levels = c("R", "NR"))

  pdf(file.path(output_dir, paste0("04_response_boxplot_", drug_name, ".pdf")),
      width = 5, height = 5)
  p4 <- plot_response_boxplot(exp_vs_pred)
  print(p4)
  dev.off()
  message("Saved: 04_response_boxplot_", drug_name, ".pdf")

  # ROC curve
  pdf(file.path(output_dir, paste0("05_roc_curve_", drug_name, ".pdf")),
      width = 5, height = 5)
  p5 <- plot_roc_curve(response = exp_vs_pred$response, predictor = exp_vs_pred$predicted_killing)
  print(p5)
  dev.off()
  message("Saved: 05_roc_curve_", drug_name, ".pdf")
}

# -- 5d. Patient response panel (composite) --
pdf(file.path(output_dir, "06_patient_response_panel.pdf"), width = 10, height = 15)
p6 <- plot_patient_response_panel(
  clone_distribution = clone_dist_df,
  clone_killing      = clone_kill_plot_df,
  exp_vs_pred        = exp_vs_pred,
  killing_col        = "comb_killing"
)
print(p6)
dev.off()
message("Saved: 06_patient_response_panel.pdf")

# -- 5e. t-SNE / UMAP visualization --
# Use test cell line expression for Seurat clustering
# Use first 100 cells to speed up
n_cells_use <- min(100, ncol(expr_test))
expr_subset <- expr_test[, 1:n_cells_use]

pdf(file.path(output_dir, "07_seurat_clustering.pdf"), width = 6, height = 5)
seurat_result <- plot_seurat_clustering(expr_subset)
print(seurat_result$umap_plot)
dev.off()
message("Saved: 07_seurat_clustering.pdf")

# t-SNE response overlay
if (!is.null(seurat_result$seurat_object)) {
  umap_coords <- as.data.frame(Seurat::Embeddings(seurat_result$seurat_object, "umap"))
  colnames(umap_coords) <- c("X", "Y")

  # Map predictions to t-SNE
  pred_for_tsne <- cell_pred[1:n_cells_use, 1]
  tsne_data <- data.frame(
    X = umap_coords$X,
    Y = umap_coords$Y,
    killing_scaled = range01(rank(-pred_for_tsne))
  )

  pdf(file.path(output_dir, "08_tsne_response.pdf"), width = 6, height = 5)
  p8 <- plot_tsne_response(tsne_data, color_var = "killing_scaled",
                           title = paste("Drug Killing -", drug_list[1]))
  print(p8)
  dev.off()
  message("Saved: 08_tsne_response.pdf")

  # t-SNE biomarker + killing side by side
  best_gene <- models[[1]]$single_best
  if (!is.null(best_gene) && best_gene %in% rownames(expr_subset)) {
    biomarker_vals <- as.numeric(expr_subset[best_gene, ])
    tsne_data$biomarker_scaled <- range01(rank(biomarker_vals))

    pdf(file.path(output_dir, "09_tsne_biomarker_killing.pdf"), width = 10, height = 5)
    p9 <- plot_tsne_biomarker_killing(tsne_data,
                                      biomarker_var = "biomarker_scaled",
                                      killing_var   = "killing_scaled",
                                      biomarker_label = paste0("Biomarker (", best_gene, ")"))
    print(p9)
    dev.off()
    message("Saved: 09_tsne_biomarker_killing.pdf")
  }
}

# -- Done ---------------------------------------------------------------------

message("\n", strrep("=", 60))
message("Pipeline complete! All outputs saved to: ", output_dir)
message(strrep("=", 60))

# List all output files
message("\nGenerated files:")
for (f in list.files(output_dir, full.names = TRUE)) {
  message("  ", f)
}
