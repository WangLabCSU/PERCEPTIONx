# PERCEPTION

An R package for predicting patient response and resistance to cancer treatment using single-cell transcriptomics.

## Overview

PERCEPTION (PatiEnt Response Prediction using Single-Cell Transcriptomics) is a computational framework that predicts how individual patients respond to drug treatments by leveraging both bulk and single-cell RNA sequencing data. It trains models on DepMap cell line data and applies them to patient single-cell expression profiles, enabling clone-level drug sensitivity prediction and patient-level response stratification.

This package is an R implementation of the original PERCEPTION pipeline, providing a unified and reproducible interface for model training, prediction, evaluation, and visualization.

> **Reference**: Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. *Nat Cancer* 5, 938–952 (2024). [https://doi.org/10.1038/s43018-024-00756-7](https://doi.org/10.1038/s43018-024-00756-7)
>
> **Original repository**: [https://github.com/ruppinlab/PERCEPTION](https://github.com/ruppinlab/PERCEPTION)

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("SunPast/PERCEPTION")
```

## Quick Start

### 💾 Load Data

```r
library(PERCEPTION)

# Load pre-trained models
models <- load_model("abemaciclib")

# Load DepMap reference data
DepMap <- load_depmap(read = TRUE)
```

### 🧠 Train Models

```r
# Identify available genes across expression and scRNA datasets
available_genes <- intersect(rownames(DepMap$expression_20Q4),
                             rownames(DepMap$scRNA_complete))

# Sample genes of interest
set.seed(123)
GOI_100 <- sample(available_genes, 100)

# Train a model for a single drug
models <- train_models(
  drug_list = "abemaciclib",
  cancer_type = "PanCan",
  exclude_cancer = "PanCan",
  GOI = GOI_100,
  ncores = 1
)
```

### 🎯 Predict Drug Response

```r
# Clone-level prediction (returns matrix: clones x drugs)
clone_killing <- predict_drugs(
  model_list = models,
  expr = sc_expression_rnorm
)

# Build clone_killing_matrix with patient and clone_id columns
# (clone_ids from rownames, patients extracted via strsplit_customv0)
clone_killing_df <- data.frame(
  patient = strsplit_customv0(rownames(clone_killing), "_", 1),
  clone_id = rownames(clone_killing),
  clone_killing,
  check.names = FALSE
)

# Patient-level aggregation
patient_pred <- predict_patients(
  clone_killing_matrix = clone_killing_df,
  clone_counts = clone_counts,
  mode = "weighted_max"
)
```

### 🎨 Visualize Results

```r
# t-SNE with drug response overlay
plot_tsne_response(
  tsne_data = tsne_data,
  color_var = "killing_scaled",
  title = "Drug Response"
)

# Clone distribution stacked bar chart
plot_clone_distribution(
  clone_distribution = clone_distribution,
  response_var = "response"
)

# Clone killing lollipop plot
plot_clone_killing(
  clone_killing = clone_killing,
  killing_var = "comb_killing"
)

# ROC curve
plot_roc_curve(
  response = response,
  predictor = predictor,
  smooth_curve = TRUE
)

# Response boxplot (responders vs. non-responders)
plot_response_boxplot(
  exp_vs_pred = exp_vs_pred,
  response_var = "response"
)
```

## Function Reference

### Data Loading

| Function | Description |
|----------|-------------|
| `load_depmap()` | Download and load DepMap reference datasets |
| `load_model()` | Download and load pre-trained models |
| `get_perception_mirrors()` | Get available download mirrors |
| `add_perception_mirror()` | Add a custom mirror |
| `list_perception_mirrors()` | List current mirrors |
| `reset_perception_mirrors()` | Reset to default mirrors |

### Preprocessing

| Function | Description |
|----------|-------------|
| `rank_normalization_mat()` | Rank-normalize an expression matrix |
| `range01()` | Scale a numeric vector to the 0-1 range |

### Model Training

| Function | Description |
|----------|-------------|
| `train_models()` | Full training pipeline (main entry point) |
| `get_response_matrix()` | Extract drug response data from DepMap |
| `get_cellLine_list()` | Get training/test cell line split |
| `feature_ranking_bulk()` | Rank features by correlation with drug response |
| `run_parallel_feature_ranking_bulk()` | Parallel feature ranking for multiple drugs |
| `build_on_BULK_v2()` | Build a single-drug model (glmnet or random forest) |

### Prediction

| Function | Description |
|----------|-------------|
| `predict_drugs()` | Predict drug sensitivity at clone/cell level |
| `predict_patients()` | Aggregate clone-level predictions to patient level |

### Evaluation

| Function | Description |
|----------|-------------|
| `compare_performance()` | Compare performance across model configurations |
| `get_significant_models()` | Filter models with significant stratification |
| `get_performance()` | Load pre-computed performance metrics |
| `each_patient_pseudo_bulk()` | Compute patient pseudo-bulk expression |

### Visualization

| Function | Description |
|----------|-------------|
| `plot_tsne_response()` | t-SNE/UMAP with drug response overlay |
| `plot_tsne_biomarker_killing()` | Biomarker vs. killing side-by-side on t-SNE |
| `plot_clone_distribution()` | Clone abundance stacked bar chart |
| `plot_clone_killing()` | Clone-level killing lollipop plot |
| `plot_roc_curve()` | ROC curve with AUC |
| `plot_response_boxplot()` | Responder vs. non-responder boxplot |
| `plot_model_performance()` | Model performance across thresholds |
| `plot_seurat_clustering()` | Seurat clustering and UMAP visualization |
| `plot_patient_response_panel()` | Composite patient response panel |

### Utilities

| Function | Description |
|----------|-------------|
| `err_handle()` | Error-safe evaluation (returns NA on error) |
| `stripall2match()` | Normalize strings for fuzzy matching |
| `strsplit_customv0()` | Split strings and extract elements |
| `hypergeometric_test_for_twolists()` | Hypergeometric enrichment test |
| `fdrcorr()` | FDR multiple testing correction |

## Workflow

```
DepMap Data ──► Preprocessing ──► Feature Ranking ──► Model Training
                   │                                       │
            rank_normalization_mat()              train_models()
                                                       │
Patient scRNA ──► Preprocessing ──► Clone Prediction ──► Patient Aggregation
                   │                      │                      │
            rank_normalization_mat()   predict_drugs()    predict_patients()
                                                               │
                                                    Visualization & Evaluation
                                                               │
                                              plot_roc_curve() / plot_response_boxplot()
                                              compare_performance() / get_significant_models()
```

## Data Requirements

- **DepMap reference data**: Automatically downloaded via `load_depmap()`, including bulk expression, single-cell expression, drug response (AUC), and cell line annotations.
- **Patient data**: Single-cell RNA expression matrix (genes as rows, cells as columns), rank-normalized via `rank_normalization_mat()`.
- **Clone annotations**: Mapping from cells to clones/patients, with clone proportions per patient.

## Citation

If you use this package, please cite the original PERCEPTION study:

Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. *Nat Cancer* 5, 938–952 (2024). https://doi.org/10.1038/s43018-024-00756-7

## License

MIT © PERCEPTION authors
