# =============================================================================
# PERCEPTION Test Pipeline v2 - Simplified Workflow
# =============================================================================
# End-to-end test using the new simplified 5-step workflow:
#   1. load_depmap()       - Load DepMap reference data
#   2. train_models()      - Train drug response models
#   3. prepare_data()      - Prepare patient data (Seurat clustering + normalization)
#   4. predict_drugs()     - Predict clone-level drug response
#   5. predict_patients()  - Aggregate to patient-level response
#
# This version demonstrates the simplified user workflow where they only need to
# provide: expression matrix + patient-cell mapping (list or metadata data frame)
#
# NOTE: This script downloads DepMap data (~567MB). Run it manually only when needed.
# =============================================================================

# Load the development version of the package (use this instead of library(PERCEPTION))
devtools::load_all(".")

# -- Parameter Settings -------------------------------------------------------

# Drug list (2 drugs for quick testing)
drug_list <- c("erlotinib", "gefitinib")

# Cancer type
cancer_type   <- "PanCan"
exclude_cancer <- "PanCan"

# Feature genes (NULL = use all DepMap genes)
GOI <- NULL

# Model hyperparameters
model_type     <- "glmnet"
k_features_values <- NULL  # NULL = auto-compute
ncores         <- 1

# Output directory
output_dir <- "./test_pipeline_v2_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# STEP 1: Load DepMap Data
# =============================================================================
message("\n", strrep("=", 60))
message("STEP 1: Loading DepMap reference data")
message(strrep("=", 60))

load_depmap(read = TRUE, mirror = TRUE)

# =============================================================================
# STEP 2: Train Models
# =============================================================================
message("\n", strrep("=", 60))
message("STEP 2: Training drug response models")
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

message("Trained models for: ", paste(names(models), collapse = ", "))

# =============================================================================
# STEP 3: Prepare Patient Data (NEW SIMPLIFIED WORKFLOW)
# =============================================================================
message("\n", strrep("=", 60))
message("STEP 3: Preparing patient data")
message(strrep("=", 60))

# --- Simulate patient single-cell expression data ---
# Use DepMap single-cell data to simulate patient expression profiles
# Get test cell lines (excluded from training)
cellLines_test <- get_cellLine_list(
  infunc_cancerType = cancer_type,
  infunc_drugName   = drug_list[1],
  exclude_cancer    = exclude_cancer,
  infunc_response   = get_response_matrix(drug_list[1])
)[[2]]

# Get expression for test cells
test_cells <- DepMap$metadata_CPM_scRNA$NAME[
  DepMap$metadata_CPM_scRNA$DepMap_ID %in% cellLines_test]

# Use first 200 cells to simulate 3 patients (need enough cells for Seurat PCA)
n_cells <- 200
set.seed(42)
selected_cells <- sample(test_cells, n_cells)

# Create simulated expression matrix (genes x cells)
# Note: Using pre-normalized data for demonstration; in real use, provide raw counts
expr_matrix <- DepMap$CPM_scRNA_CCLE_rnorm[, selected_cells, drop = FALSE]

# --- Create patient-cell mapping (TWO OPTIONS DEMONSTRATED) ---
patient_ids_sim <- c("Patient_A", "Patient_B", "Patient_C")

# OPTION 1: List format (same as Running_PERCEPTION_for_new_dataset.Rmd)
sample_cell_names_list <- list(
  Patient_A = selected_cells[1:60],
  Patient_B = selected_cells[61:140],
  Patient_C = selected_cells[141:200]
)

# OPTION 2: Data frame format (more convenient for most users)
metadata_df <- data.frame(
  cell_id = selected_cells,
  patient_id = c(
    rep("Patient_A", 60),
    rep("Patient_B", 80),
    rep("Patient_C", 60)
  )
)

# --- Run prepare_data() ---
# Using Option 2 (data frame format) - simpler for users
message("Using metadata data frame format...")
prepared <- prepare_data(
  expression_matrix = expr_matrix,
  patient_mapping   = metadata_df,  # Can also use sample_cell_names_list
  genes_to_use      = NULL
)

message("\nPrepared data structure:")
message("  - clone_expression_rnorm: ", nrow(prepared$clone_expression_rnorm), " genes x ",
        ncol(prepared$clone_expression_rnorm), " clones")
message("  - clone_counts: ", nrow(prepared$clone_counts), " patients x ",
        ncol(prepared$clone_counts) - 1, " clones")
message("  - cell_clone_map: ", nrow(prepared$cell_clone_map), " cells annotated")

# Show clone distribution
message("\nClone counts per patient:")
print(prepared$clone_counts)

# =============================================================================
# STEP 4: Predict Drug Response at Clone Level
# =============================================================================
message("\n", strrep("=", 60))
message("STEP 4: Predicting clone-level drug response")
message(strrep("=", 60))

clone_pred <- predict_drugs(models, prepared$clone_expression_rnorm)

message("Clone-level predictions: ", nrow(clone_pred), " clones x ", ncol(clone_pred), " drugs")
message("\nClone predictions (first 5 clones):")
print(head(clone_pred, 5))

# =============================================================================
# STEP 5: Predict Patient-Level Response (NEW SIMPLIFIED CALL)
# =============================================================================
message("\n", strrep("=", 60))
message("STEP 5: Aggregating to patient-level response")
message(strrep("=", 60))

# NEW: Just pass clone_pred and prepared directly - no manual cbind needed!
patient_pred <- predict_patients(clone_pred, prepared, mode = "weighted_max")

message("Patient-level predictions (most resistant clone):")
print(patient_pred)

# Also try weighted_average mode
patient_pred_weighted <- predict_patients(clone_pred, prepared, mode = "weighted_average")
message("\nPatient-level predictions (weighted average):")
print(patient_pred_weighted)

# =============================================================================
# OPTIONAL: Visualization
# =============================================================================
message("\n", strrep("=", 60))
message("OPTIONAL: Generating visualizations")
message(strrep("=", 60))

# Build clone_killing_df for plotting (needed for plot functions)
clone_killing_df <- data.frame(
  patient = prepared$clone_killing_template$patient,
  clone_id = prepared$clone_killing_template$clone_id,
  clone_pred
)

# Clone distribution plot
clone_dist_df <- data.frame(
  patients = clone_killing_df$patient,
  clones   = clone_killing_df$clone_id,
  weights  = prepared$clone_counts[
    match(clone_killing_df$patient, prepared$clone_counts$patients),
    match(clone_killing_df$clone_id, colnames(prepared$clone_counts)[-1])
  ]
)
# Normalize weights per patient
for (pid in unique(clone_dist_df$patients)) {
  idx <- clone_dist_df$patients == pid
  clone_dist_df$weights[idx] <- clone_dist_df$weights[idx] / sum(clone_dist_df$weights[idx])
}

pdf(file.path(output_dir, "clone_distribution.pdf"), width = 8, height = 5)
p1 <- plot_clone_distribution(clone_dist_df)
print(p1)
dev.off()
message("Saved: clone_distribution.pdf")

# Clone killing lollipop plot (use first drug)
drug_col <- drug_list[1]
clone_kill_plot_df <- clone_killing_df
clone_kill_plot_df$comb_killing <- clone_kill_plot_df[[drug_col]]
clone_kill_plot_df$weights <- clone_dist_df$weights

pdf(file.path(output_dir, "clone_killing.pdf"), width = 12, height = 5)
p2 <- plot_clone_killing(clone_kill_plot_df, killing_var = "comb_killing", weights_var = "weights")
print(p2)
dev.off()
message("Saved: clone_killing.pdf")

# =============================================================================
# SUMMARY
# =============================================================================
message("\n", strrep("=", 60))
message("PIPELINE COMPLETE!")
message(strrep("=", 60))

message("\nThe new simplified workflow:")
message("  1. load_depmap()       -> DepMap reference data")
message("  2. train_models()      -> Drug response models")
message("  3. prepare_data()      -> Clone expression + counts (from raw scRNA + metadata)")
message("  4. predict_drugs()    -> Clone-level predictions")
message("  5. predict_patients()  -> Patient-level predictions (just pass clone_pred + prepared)")

message("\nUser inputs needed:")
message("  - expression_matrix: genes x cells (raw single-cell expression)")
message("  - patient_mapping: list OR data frame with cell_id + patient_id columns")

message("\nOutput files saved to: ", output_dir)
for (f in list.files(output_dir, full.names = TRUE)) {
  message("  ", f)
}
