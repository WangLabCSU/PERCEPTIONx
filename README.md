# PERCEPTION
<!-- badges: start -->

<!-- badges: end -->

An R package for predicting patient response and resistance to cancer treatment using single-cell transcriptomics.

## 1. Overview

PERCEPTION (PERsonalized single-Cell Expression-based Planning for Treatments In ONcology) is a computational framework that predicts how individual patients respond to drug treatments by leveraging both bulk and single-cell RNA sequencing data. It trains models on DepMap cell line data and applies them to patient single-cell expression profiles, enabling clone-level drug sensitivity prediction and patient-level response stratification.

This package is an R implementation of the original PERCEPTION pipeline, providing a unified and reproducible interface for model training, prediction, evaluation, and visualization.

> **Reference**: Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. *Nat Cancer* 5, 938–952 (2024). [https://doi.org/10.1038/s43018-024-00756-7](https://doi.org/10.1038/s43018-024-00756-7)
>
> **Original repository**: [https://github.com/ruppinlab/PERCEPTION](https://github.com/ruppinlab/PERCEPTION)

## 2. Installation

Install the development version from GitHub using devtools.

```r
# install.packages("devtools")
devtools::install_github("SunPast/PERCEPTION")
```

## 3. Quick Start

### 3.1 💾 Load Data

PERCEPTION relies on DepMap reference data and optional pre-trained models. Both can be downloaded automatically with the built-in loading functions.

```r
library(PERCEPTION)

# Load pre-trained models
models <- load_model("abemaciclib")

# Load DepMap reference data
load_depmap(read = TRUE)
```

### 3.2 🧠 Train Models

Before training, identify the genes available across both bulk and single-cell expression datasets. The `train_models()` function then performs feature ranking, model building, and hyperparameter tuning in a single call.

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

### 3.3 🎯 Predict Drug Response

> **Important: Rank Normalization**
>
> PERCEPTION models are trained on **rank-normalized** expression data. If you provide your own expression data (e.g., from scRNA-seq), you **must** normalize it first using `rank_normalization_mat()`, or predictions will be unreliable.
>
> **How rank normalization works**: For each cell (column), every gene's expression value is replaced by its rank within that column, divided by the total number of genes: `x_norm = rank(x) / n`. This transforms each column into a uniform distribution over (0, 1], making the data robust to batch effects, library size differences, and outliers. Since the model coefficients capture the relationship between **relative gene expression ranks** and drug response (not absolute values), the same normalization must be applied to any new data.
>
> ```r
> # If your data is NOT already rank-normalized:
> my_expr_norm <- rank_normalization_mat(my_raw_expr)
> # Then use my_expr_norm in predict_drugs()
> ```

Prediction proceeds in two stages: first, `predict_drugs()` scores each clone's drug sensitivity from the expression matrix; then, `predict_patients()` aggregates clone-level scores into a patient-level prediction using clone proportions.

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

### 3.4 🎨 Visualize Results

PERCEPTION provides a suite of plotting functions to inspect model predictions from different perspectives: spatial (t-SNE), clonal (distribution and killing), and clinical (ROC and response stratification).

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

# ROC curve with AUC annotation
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

## 4. Function Reference

### 4.1 Data Loading

| Function | Description |
|----------|-------------|
| `load_depmap()` | Download and load DepMap reference datasets |
| `load_model()` | Download and load pre-trained models |
| `get_mirrors()` | Get available download mirrors |
| `add_mirrors()` | Add a custom mirror |
| `list_mirrors()` | List current mirrors |
| `reset_mirrors()` | Reset to default mirrors |

### 4.2 Preprocessing

| Function | Description |
|----------|-------------|
| `rank_normalization_mat()` | Rank-normalize an expression matrix |
| `range01()` | Scale a numeric vector to the 0-1 range |

### 4.3 Model Training

| Function | Description |
|----------|-------------|
| `train_models()` | Full training pipeline (main entry point) |
| `get_response_matrix()` | Extract drug response data from DepMap |
| `get_cellLine_list()` | Get training/test cell line split |
| `feature_ranking_bulk()` | Rank features by correlation with drug response |
| `run_parallel_feature_ranking_bulk()` | Parallel feature ranking for multiple drugs |
| `build_on_BULK_v2()` | Build a single-drug model (glmnet or random forest) |

### 4.4 Prediction

| Function | Description |
|----------|-------------|
| `predict_drugs()` | Predict drug sensitivity at clone/cell level |
| `predict_patients()` | Aggregate clone-level predictions to patient level |

### 4.5 Evaluation

| Function | Description |
|----------|-------------|
| `compare_performance()` | Compare performance across model configurations |
| `get_significant_models()` | Filter models with significant stratification |
| `get_performance()` | Load pre-computed performance metrics |
| `each_patient_pseudo_bulk()` | Compute patient pseudo-bulk expression |

### 4.6 Visualization

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

### 4.7 Utilities

| Function | Description |
|----------|-------------|
| `err_handle()` | Error-safe evaluation (returns NA on error) |
| `stripall2match()` | Normalize strings for fuzzy matching |
| `strsplit_customv0()` | Split strings and extract elements |
| `hypergeometric_test_for_twolists()` | Hypergeometric enrichment test |
| `fdrcorr()` | FDR multiple testing correction |

## 5. Workflow

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

## 6. Data Requirements

- **DepMap reference data**: Automatically downloaded via `load_depmap()`, including bulk expression, single-cell expression, drug response (AUC), and cell line annotations.
- **Patient data**: Single-cell RNA expression matrix (genes as rows, cells as columns), rank-normalized via `rank_normalization_mat()`.
- **Clone annotations**: Mapping from cells to clones/patients, with clone proportions per patient.

## 7. Citation

If you use this package, please cite the original PERCEPTION study:

Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. *Nat Cancer* 5, 938–952 (2024). https://doi.org/10.1038/s43018-024-00756-7

## 8. License

MIT © PERCEPTION authors
