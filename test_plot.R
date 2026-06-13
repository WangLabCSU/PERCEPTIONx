# ==============================================================================
# PERCEPTION plot.R Complete Test Script
# ==============================================================================
# This script tests all visualization functions in plot.R using real patient
# data from the PERCEPTION paper (Step2-3-5).
# ==============================================================================

# Set working directory
setwd("C:/Users/LENOVO/Desktop/PERCEPTION")

# Load required packages
required_packages <- c("ggplot2", "gridExtra", "pROC", "ggpubr", "viridis",
                       "tidyr", "dplyr", "readxl", "Seurat")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

library(ggplot2)
library(gridExtra)
library(pROC)
library(ggpubr)
library(viridis)
library(tidyr)
library(dplyr)
library(readxl)

# Load PERCEPTION functions
source('R/utils.R')
source('R/stats.R')
source('R/utils_cortest.R')
source('R/load.R')
source('R/train.R')
source('R/predict.R')
source('R/evaluate.R')
source('R/plot.R')

cat("========================================\n")
cat("PERCEPTION plot.R Test Suite\n")
cat("========================================\n\n")

# ==============================================================================
# Test 1: Multiple Myeloma Patient Response (Step3)
# ==============================================================================
cat("Test 1: Multiple Myeloma Patient Response Visualization\n")
cat("--------------------------------------------------------\n")

# Load patient data
clone_Level_z_expression <- readxl::read_xlsx(
  'Data/Supp_COhen_IdoAmit_etal/Supp_COhen_IdoAmit_etal/clone_Level_z_expression.xlsx', skip = 2)
resp <- readxl::read_xlsx(
  'Data/Supp_COhen_IdoAmit_etal/Supp_COhen_IdoAmit_etal/resp_info.xlsx')
resp <- na.omit(resp)

Clone_Counts_per_patients <- readxl::read_xlsx(
  'Data/Supp_COhen_IdoAmit_etal/Supp_COhen_IdoAmit_etal/Clone_Counts_per_patients.xlsx', sheet = 2)
Clone_Counts_per_patients <- Clone_Counts_per_patients[grep('Kydar', Clone_Counts_per_patients$...1), ]
colnames(Clone_Counts_per_patients) <- c('patients', 'c1', 'c2', 'c3')

genesUsed_toBuild <- readRDS('Data/genesUsed_toBuild.RDS')

cat("  Loaded:", nrow(resp), "patients with response data\n")
cat("  Loaded:", nrow(Clone_Counts_per_patients), "patients with clone counts\n")
cat("  Loaded:", length(genesUsed_toBuild), "feature genes\n\n")

# Preprocess expression data
rownames_for_cloneexp_matrix <- make.names(clone_Level_z_expression$Gene, unique = TRUE)
clone_Level_z_expression <- clone_Level_z_expression[, -1]
clone_Level_z_expression_df <- data.frame(clone_Level_z_expression,
                                          row.names = rownames_for_cloneexp_matrix)

# Check gene name matching
genes_available <- intersect(genesUsed_toBuild, rownames_for_cloneexp_matrix)
cat("  Genes available in both:", length(genes_available), "/", length(genesUsed_toBuild), "\n")

if (length(genes_available) < length(genesUsed_toBuild)) {
  missing_in_clone <- genesUsed_toBuild[!genesUsed_toBuild %in% rownames_for_cloneexp_matrix]
  cat("  Missing in clone expression:", length(missing_in_clone), "genes\n")
  cat("  Examples:", paste(head(missing_in_clone, 5), collapse = ", "), "\n")
}

# Use only available genes for rank normalization
clone_Level_z_expression_rnorm <- rank_normalization_mat(
  clone_Level_z_expression_df[genes_available, ])

cat("  Expression matrix dimensions:", dim(clone_Level_z_expression_rnorm), "\n\n")

# ==============================================================================
# Train models for carfilzomib and lenalidomide
# ==============================================================================
cat("  Training models for carfilzomib and lenalidomide...\n")

# Load DepMap for training
load_depmap(read = TRUE)

# Train models (using available genes only)
models_mm <- train_perception_models(
  drug_list = c("carfilzomib", "lenalidomide"),
  cancer_type = "PanCan",
  exclude_cancer = "PanCan",
  GOI = genes_available,
  ncores = 1
)

cat("  Models trained successfully\n\n")

# ==============================================================================
# Predict killing for each clone
# ==============================================================================
cat("  Predicting killing for each clone...\n")

killing_eachClone <- predict_drugs(
  model_list = models_mm,
  expr = clone_Level_z_expression_rnorm
)

killing_eachClone <- data.frame(killing_eachClone)
killing_eachClone_z <- data.frame(apply(killing_eachClone, 2, scale))
rownames(killing_eachClone_z) <- rownames(killing_eachClone)

cat("  Predictions for", nrow(killing_eachClone_z), "clones\n\n")

# Combination killing (min principle)
combination_Killing <- pmin(killing_eachClone_z$carfilzomib,
                            killing_eachClone_z$lenalidomide)
names(combination_Killing) <- rownames(killing_eachClone)

# Build clone killing dataframe (using strsplit_customv0 from utils.R)
comb_killing_df <- data.frame(
  patient = gsub('z.', '', strsplit_customv0(names(combination_Killing), '_', 1)),
  clone_id = strsplit_customv0(names(combination_Killing), '_', 2),
  comb_killing = combination_Killing
)

# ==============================================================================
# Test plot_clone_distribution
# ==============================================================================
cat("Test 1.1: plot_clone_distribution\n")

# Build clone distribution data using predict_patients helper logic
clone_cols <- setdiff(colnames(Clone_Counts_per_patients), "patients")
clone_distribution <- data.frame(t(sapply(1:nrow(Clone_Counts_per_patients), function(P) {
  patient_clones <- comb_killing_df[
    comb_killing_df$patient == Clone_Counts_per_patients$patients[P], ]$clone_id
  clone_weights <- rep(0, length(clone_cols))
  names(clone_weights) <- clone_cols
  if (length(patient_clones) > 0) {
    total_cells <- sum(Clone_Counts_per_patients[P, patient_clones])
    if (total_cells > 0) {
      existing_weights <- unlist(Clone_Counts_per_patients[P, patient_clones] / total_cells)
      clone_weights[match(names(existing_weights), names(clone_weights))] <- existing_weights
    }
  }
  clone_weights
})))
clone_distribution$patients <- Clone_Counts_per_patients$patients
clone_distribution_df <- gather(clone_distribution, clones, weights, c1:c3)
clone_distribution_df$response <- resp$response[match(clone_distribution_df$patients, resp$Patient)]

panelA <- plot_clone_distribution(clone_distribution_df, response_var = "response")
cat("  SUCCESS: Clone distribution plot created\n")
print(panelA)

# ==============================================================================
# Test plot_clone_killing
# ==============================================================================
cat("\nTest 1.2: plot_clone_killing\n")

# Add weights to clone killing dataframe
comb_killing_df$weights <- clone_distribution_df$weights[
  match(paste0(comb_killing_df$patient, '_', comb_killing_df$clone_id),
        paste0(clone_distribution_df$patients, '_', clone_distribution_df$clones))]
comb_killing_df$response <- resp$response[match(comb_killing_df$patient, resp$Patient)]
comb_killing_df$response <- factor(comb_killing_df$response, labels = c('NR', 'R'))

panelB <- plot_clone_killing(comb_killing_df,
                             killing_var = "comb_killing",
                             weights_var = "weights",
                             response_var = "response")
cat("  SUCCESS: Clone killing lollipop plot created\n")
print(panelB)

# ==============================================================================
# Test plot_response_boxplot and plot_roc_curve
# ==============================================================================
cat("\nTest 1.3: plot_response_boxplot\n")
cat("Test 1.4: plot_roc_curve\n")

# Patient-level aggregation using predict_patients
# Only pass patient, clone_id, and comb_killing columns
comb_killing_for_pred <- comb_killing_df[, c("patient", "clone_id", "comb_killing")]
patient_pred <- predict_patients(
  clone_killing_matrix = comb_killing_for_pred,
  clone_counts = Clone_Counts_per_patients,
  mode = "weighted_max"
)
most_resistant_clone_based_killing <- patient_pred[, 1]
names(most_resistant_clone_based_killing) <- Clone_Counts_per_patients$patients

# Build prediction vs observation dataframe
Exp_vs_pred_killing <- data.frame(
  resp,
  most_resistant_clone_based_killing = most_resistant_clone_based_killing[
    match(resp$Patient, names(most_resistant_clone_based_killing))]
)
Exp_vs_pred_killing$response <- factor(Exp_vs_pred_killing$response)
Exp_vs_pred_killing$response <- factor(Exp_vs_pred_killing$response,
                                       rev(levels(Exp_vs_pred_killing$response)))

panelC <- plot_response_boxplot(Exp_vs_pred_killing,
                                response_var = "response",
                                predicted_var = "most_resistant_clone_based_killing")
cat("  SUCCESS: Response boxplot created\n")
print(panelC)

panelD <- plot_roc_curve(response = Exp_vs_pred_killing$response,
                         predictor = Exp_vs_pred_killing$most_resistant_clone_based_killing)
cat("  SUCCESS: ROC curve created\n")
print(panelD)

# ==============================================================================
# Test plot_patient_response_panel (combined)
# ==============================================================================
cat("\nTest 1.5: plot_patient_response_panel (complete panel)\n")

complete_panel <- plot_patient_response_panel(
  clone_distribution = clone_distribution_df,
  clone_killing = comb_killing_df,
  exp_vs_pred = Exp_vs_pred_killing,
  response_col = "response",
  killing_col = "comb_killing",
  predicted_col = "most_resistant_clone_based_killing"
)
cat("  SUCCESS: Complete patient response panel created\n")

# Save panel
ggsave(complete_panel, filename = "test_output_mm_panel.pdf", height = 64, width = 64, limitsize =  F)
cat("  Saved to: test_output_mm_panel.pdf\n\n")

# ==============================================================================
# Test 2: Lung Cancer t-SNE Visualization (Step2)
# ==============================================================================
cat("Test 2: Lung Cancer t-SNE Visualization\n")
cat("----------------------------------------\n")

# Load lung t-SNE data
lung_tSNE <- read.csv('Data/lung_tSNE.txt', sep = '\t')
cat("  Loaded t-SNE coordinates for", nrow(lung_tSNE), "cells\n\n")

# Train model for erlotinib (EGFR inhibitor, relevant for lung cancer)
cat("  Training model for erlotinib (EGFR inhibitor)...\n")

models_lung <- train_perception_models(
  drug_list = "erlotinib",
  cancer_type = "PanCan",
  exclude_cancer = "PanCan",
  GOI = genes_available,
  ncores = 1
)

# Predict on lung scRNA data (use CPM_scRNA_CCLE_rnorm which has all cells)
lung_scRNA_rnorm <- DepMap$CPM_scRNA_CCLE_rnorm
viab_lung <- predict_drugs(
  model_list = models_lung,
  expr = lung_scRNA_rnorm
)
viab_lung <- viab_lung[, 1]  # Extract single drug vector

cat("  Predicted viability for", length(viab_lung), "cells\n")
cat("  t-SNE has", nrow(lung_tSNE), "cells\n")

# Match t-SNE cells to predicted cells
common_cells <- intersect(rownames(lung_tSNE), names(viab_lung))
cat("  Common cells:", length(common_cells), "\n\n")

# ==============================================================================
# Test plot_tsne_response
# ==============================================================================
cat("Test 2.1: plot_tsne_response\n")

# Build t-SNE data with killing overlay (only for cells with both t-SNE and prediction)
# range01 is now exported from stats.R
tsne_data <- data.frame(
  X = lung_tSNE[common_cells, "X"],
  Y = lung_tSNE[common_cells, "Y"],
  killing_scaled = range01(rank(-viab_lung[common_cells]))
)

panel_tsne <- plot_tsne_response(tsne_data,
                                 color_var = "killing_scaled",
                                 title = "Erlotinib Predicted Killing",
                                 color_label = "Predicted Viability")
cat("  SUCCESS: t-SNE response plot created\n")
print(panel_tsne)

ggsave(panel_tsne, filename = "test_output_lung_tsne.pdf", height = 32, width = 32)
cat("  Saved to: test_output_lung_tsne.pdf\n\n")

# ==============================================================================
# Test 3: Seurat Clustering (Step5)
# ==============================================================================
cat("Test 3: Seurat Clustering and UMAP\n")
cat("----------------------------------\n")

# Load lung scRNA expression
PRJNA591860 <- readRDS('Data/PRJNA591860.RDS')
cat("  Loaded expression matrix:", dim(PRJNA591860), "\n\n")

# Test run_seurat_clustering (subset for speed)
cat("Test 3.1: run_seurat_clustering\n")
cat("  (Using subset of 1000 cells for speed)\n")

# Subset for testing
subset_cells <- sample(colnames(PRJNA591860), min(1000, ncol(PRJNA591860)))
PRJNA591860_subset <- PRJNA591860[, subset_cells]

seurat_result <- run_seurat_clustering(PRJNA591860_subset,
                                       min_cells = 3,
                                       min_features = 200,
                                       nfeatures = 2000,
                                       dims = 10,
                                       resolution = 0.8)

cat("  SUCCESS: Seurat clustering completed\n")
cat("  Found", length(unique(seurat_result$cluster_ids)), "clusters\n")
print(seurat_result$umap_plot)

ggsave(seurat_result$umap_plot, filename = "test_output_seurat_umap.pdf", height = 32, width = 32)
cat("  Saved to: test_output_seurat_umap.pdf\n\n")

# ==============================================================================
# Test 4: Model Performance Plot
# ==============================================================================
cat("Test 4: Model Performance Visualization\n")
cat("---------------------------------------\n")

# Use models trained earlier
cat("Test 4.1: plot_model_performance\n")

panel_perf <- plot_model_performance(models_mm)
cat("  SUCCESS: Model performance plot created\n")
print(panel_perf)

ggsave(panel_perf, filename = "test_output_model_performance.pdf", height = 24, width = 24)
cat("  Saved to: test_output_model_performance.pdf\n\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("========================================\n")
cat("All Tests Completed Successfully!\n")
cat("========================================\n")
cat("\nGenerated output files:\n")
cat("  - test_output_mm_panel.pdf (Multiple Myeloma complete panel)\n")
cat("  - test_output_lung_tsne.pdf (Lung t-SNE visualization)\n")
cat("  - test_output_seurat_umap.pdf (Seurat UMAP)\n")
cat("  - test_output_model_performance.pdf (Model performance)\n")
cat("\nAll plot.R functions tested:\n")
cat("  ✓ plot_clone_distribution\n")
cat("  ✓ plot_clone_killing\n")
cat("  ✓ plot_response_boxplot\n")
cat("  ✓ plot_roc_curve\n")
cat("  ✓ plot_patient_response_panel\n")
cat("  ✓ plot_tsne_response\n")
cat("  ✓ run_seurat_clustering\n")
cat("  ✓ plot_model_performance\n")
cat("\n")
