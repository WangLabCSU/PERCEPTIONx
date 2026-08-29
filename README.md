<img src="man/figures/logo.png" class="home-logo" align="right" width="160" alt="PERCEPTIONx hex sticker">

# PERCEPTIONx

<!-- badges: start -->

[![R >= 4.1.0](https://img.shields.io/badge/R-%3E%3D%204.1.0-276DC3.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/license/mit)
[![Lifecycle: experimental](https://img.shields.io/badge/Lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![pkgdown](https://github.com/WangLabCSU/PERCEPTIONx/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/WangLabCSU/PERCEPTIONx/actions/workflows/pkgdown.yaml)
[![Last commit](https://img.shields.io/github/last-commit/WangLabCSU/PERCEPTIONx)](https://github.com/WangLabCSU/PERCEPTIONx)
[![DOI](https://img.shields.io/badge/DOI-10.1038/s43018--024--00756--7-informational)](https://doi.org/10.1038/s43018-024-00756-7)

<!-- badges: end -->

**Predicting Personalized Drug Response and Resistance from Single-Cell Tumor Transcriptomics.**

PERCEPTIONx trains drug-response models on DepMap cell-line screens and applies them to a patient's single-cell expression profile. It scores sensitivity at the clone and patient level, then turns the results into publication-ready figures — with a point-and-click Shiny web application wrapping the whole pipeline. The method is the [PERCEPTION](https://doi.org/10.1038/s43018-024-00756-7) approach (Sinha et al., *Nat Cancer* 5, 938–952, 2024).

<br clear="right">

---

## 1. Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("WangLabCSU/PERCEPTIONx")
```

---

## 2. Quick Start

The core workflow is five function calls. Load the reference data, load pre-trained models (no training needed), prepare your expression matrix, predict, and plot.

### 2.1 Load Data

```r
devtools::load_all()                      # from the package source tree
models <- load_model("abemaciclib", read = TRUE)   # pre-trained models (44 drugs)
load_depmap(read = TRUE)                  # DepMap reference (~567 MB, first run only)
```

### 2.2 Train Models (optional)

The 44 pre-trained models cover prediction out of the box. Train your own when you need other drugs or settings — `train_models()` has sensible defaults for every argument except the drug list:

```r
models <- train_models(drug_list = "erlotinib")
```

### 2.3 Predict Drug Response

`prepare_data()` clusters the cells (Seurat) and rank-normalizes the clone expression; `predict_drugs()` scores each clone; `predict_patients()` aggregates clone scores to patients (default `weighted_max`):

```r
prepared <- prepare_data(my_expression_matrix, method = "umap")

clone_pred <- predict_drugs(models, prepared$clone_expression_rnorm)
patient_pred <- predict_patients(clone_pred, prepared)
```

### 2.4 Visualize Results

Each plotting function takes the corresponding output and returns a `ggplot`:

```r
plot_tsne_response(tsne_data, color_var = "viability_scaled", title = "Drug Response")
plot_clone_distribution(clone_distribution = clone_dist, response_var = "response")
plot_clone_viability(clone_viability = clone_viability, viability_var = "comb_viability")
plot_roc_curve(response = response, predictor = predictor, smooth_curve = TRUE)
plot_response_boxplot(exp_vs_pred = exp_vs_pred, response_var = "response")
```

Every plot accepts `tooltip = TRUE` (default): with the `ggiraph` package installed, points and bars get hover tooltips (clone id, viability score, proportion, FPR/TPR). See the [R Package Tutorial](https://wanglabcsu.github.io/PERCEPTIONx/articles/pipeline.html) for the full workflow, including training, model evaluation, and the complete plotting suite.

---

## 3. Shiny Web Application

PERCEPTION-shiny wraps the whole pipeline in an interactive web dashboard. Heavy computation (training, Seurat clustering, prediction, plotting) runs in background worker processes, so the interface stays responsive while a large job runs.

```r
devtools::load_all()
run_perception_app()          # starts the app in your browser
```

The app has five tabs — Data, Train, Predict, Visualize, and Help — with a Load Demo button that generates a small synthetic dataset (49 genes x 400 cells x 20 patients) to smoke-test the whole flow. See the [Shiny App Guide](https://wanglabcsu.github.io/PERCEPTIONx/articles/shiny_app.html) for a full walkthrough.

---

## 4. Function Reference

### 4.1 Data Loading

| Function | Description |
|----------|-------------|
| `load_depmap()` | Download and load DepMap reference datasets |
| `load_model()` | Download and load pre-trained models |
| `get_mirrors()` / `add_mirrors()` / `list_mirrors()` / `reset_mirrors()` | Manage download mirrors |

### 4.2 Preprocessing

| Function | Description |
|----------|-------------|
| `prepare_data()` | Seurat clustering + rank normalization → clone-level inputs |
| `rank_normalization_mat()` | Rank-normalize an expression matrix |
| `range01()` | Scale a numeric vector to the 0-1 range |

### 4.3 Model Training

| Function | Description |
|----------|-------------|
| `train_models()` | Full training pipeline (main entry point) |
| `get_response_matrix()` | Extract drug response data from DepMap |
| `get_cellLine_list()` | Get training/test cell line split |
| `feature_ranking_bulk()` | Rank features by correlation with drug response |
| `build_on_BULK_v2()` | Build a single-drug model (glmnet or random forest) |

### 4.4 Prediction

| Function | Description |
|----------|-------------|
| `predict_drugs()` | Predict drug sensitivity at clone level |
| `predict_patients()` | Aggregate clone-level predictions to patient level |

### 4.5 Evaluation

| Function | Description |
|----------|-------------|
| `compare_performance()` | Compare performance across model configurations |
| `get_significant_models()` | Filter models with significant stratification |
| `each_patient_pseudo_bulk()` | Compute patient pseudo-bulk expression |

### 4.6 Visualization

| Function | Description |
|----------|-------------|
| `plot_tsne_response()` | t-SNE/UMAP with drug response overlay |
| `plot_tsne_biomarker_viability()` | Biomarker vs. viability side-by-side on t-SNE |
| `plot_clone_distribution()` | Clone abundance stacked bar chart |
| `plot_clone_viability()` | Clone-level viability lollipop plot |
| `plot_roc_curve()` | ROC curve with AUC |
| `plot_response_boxplot()` | Responder vs. non-responder boxplot |
| `plot_model_performance()` | Model performance across thresholds |
| `plot_model_roc()` | Validation ROC curves after training |
| `plot_seurat_clustering()` | Seurat clustering and UMAP visualization |
| `plot_patient_response_panel()` | Composite patient response panel |

---

## 5. Workflow

```
DepMap Data ──► Model Training ──► Clone Prediction ──► Patient Aggregation
                   │                    │                     │
            train_models()        predict_drugs()      predict_patients()
                                                     (weighted_max)
                                                        │
Patient scRNA ──► prepare_data()                Visualization & Evaluation
                   (clustering +                     │
                    rank normalization)        plot_*() / compare_performance()
```

---

## 6. Data Requirements

- **DepMap reference data**: downloaded automatically via `load_depmap()` (bulk expression, single-cell expression, drug response, cell line annotations).
- **Patient data**: single-cell RNA expression matrix (genes as rows, cells as columns). `prepare_data()` handles clustering and rank normalization.
- **Clinical responses (optional)**: patient-level response labels, used for validation (ROC, responder vs. non-responder boxplots).

---

## 7. Citation

If you use this package, please cite the original PERCEPTION study:

Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. *Nat Cancer* 5, 938–952 (2024). https://doi.org/10.1038/s43018-024-00756-7

---

## 8. License

MIT © PERCEPTIONx authors
