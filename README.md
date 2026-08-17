# PERCEPTIONx
<!-- badges: start -->

<!-- badges: end -->

An R package for predicting patient response and resistance to cancer treatment using single-cell transcriptomics.

## 1. Overview

PERCEPTIONx implements the PERCEPTION approach (PERsonalized single-Cell Expression-based Planning for Treatments In ONcology), a computational framework that predicts how individual patients respond to drug treatments by leveraging both bulk and single-cell RNA sequencing data. It trains models on DepMap cell line data and applies them to patient single-cell expression profiles, enabling clone-level drug sensitivity prediction and patient-level response stratification, together with an interactive visualization suite and a Shiny web application.

> **Reference**: Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. *Nat Cancer* 5, 938–952 (2024). [https://doi.org/10.1038/s43018-024-00756-7](https://doi.org/10.1038/s43018-024-00756-7)

## 2. Installation

Install the development version from GitHub using devtools.

```r
# install.packages("devtools")
devtools::install_github("WangLabCSU/PERCEPTIONx")
```

## 3. Quick Start

### 3.1 💾 Load Data

PERCEPTIONx relies on DepMap reference data and optional pre-trained models. Both can be downloaded automatically with the built-in loading functions.

```r
# From the package source tree (development mode):
devtools::load_all()

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
> PERCEPTIONx models are trained on **rank-normalized** expression data. If you provide your own expression data (e.g., from scRNA-seq), you **must** normalize it first using `rank_normalization_mat()`, or predictions will be unreliable.
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
clone_viability <- predict_drugs(
  model_list = models,
  expr = sc_expression_rnorm
)

# Build clone_viability_matrix with patient and clone_id columns
# (clone_ids from rownames, patients extracted via strsplit_customv0)
clone_viability_df <- data.frame(
  patient = strsplit_customv0(rownames(clone_viability), "_", 1),
  clone_id = rownames(clone_viability),
  clone_viability,
  check.names = FALSE
)

# Patient-level aggregation (legacy: prepared_data = clone_counts data.frame)
patient_pred <- predict_patients(
  clone_pred = clone_viability_df,
  prepared_data = clone_counts,
  mode = "weighted_max"
)

# Recommended workflow instead: pass the prepare_data() result directly
# patient_pred <- predict_patients(clone_pred, prepared)
```

### 3.4 🎨 Visualize Results

PERCEPTIONx provides a suite of plotting functions to inspect model predictions from different perspectives: spatial (t-SNE), clonal (distribution and viability), and clinical (ROC and response stratification).

```r
# t-SNE with drug response overlay
plot_tsne_response(
  tsne_data = tsne_data,
  color_var = "viability_scaled",
  title = "Drug Response"
)

# Clone distribution stacked bar chart
plot_clone_distribution(
  clone_distribution = clone_distribution,
  response_var = "response"
)

# Clone viability lollipop plot
plot_clone_viability(
  clone_viability = clone_viability,
  viability_var = "comb_viability"
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

> **Interactive tooltips (optional)**: every plotting function above accepts
> `tooltip = TRUE` (default). When the [`ggiraph`](https://cran.r-project.org/package=ggiraph)
> package is installed, points/bars get hover tooltips (clone id, viability score,
> proportion, FPR/TPR, ...). Set `tooltip = FALSE` for a plain static `ggplot`
> with the identical layout. See the package vignette §6.9 for details.

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
| `plot_tsne_biomarker_viability()` | Biomarker vs. viability side-by-side on t-SNE |
| `plot_clone_distribution()` | Clone abundance stacked bar chart |
| `plot_clone_viability()` | Clone-level viability lollipop plot |
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

### 6.1 Testing with real (large-scale) data

The built-in Shiny demo ("Load Demo" in the Data tab) only generates a small
synthetic dataset (49 genes x 400 cells x 20 patients) for smoke-testing the UI.
For a meaningful large-scale test, use:

1. **DepMap reference data** — click **Load DepMap** in the Data tab (or run
   `load_depmap()`), which downloads the full reference set (~567 MB, 15k+ genes
   x 1,000+ cell lines). This is the standard training/reference input and the
   most demanding step for memory and disk.
2. **Real patient scRNA-seq** — upload a gene x cell expression matrix. Accepted
   formats: CSV / TSV / TXT / Excel / RDS (a matrix, data.frame, or Seurat object).
   **Rank-normalize first** — either call `rank_normalization_mat()` on the matrix
   or upload raw counts and let the app normalize them during `prepare_data()`.
   The closest public example is the PERCEPTION paper's own demo: **PRJNA591860**
   (lung cancer EGFR-TKI cohort, 24 patients, Maynard et al. 2020; Zenodo
   doi:10.5281/zenodo.7860559), which ships a matching expression matrix
   (`PRJNA591860.RDS`), a cell-to-patient map, and clinical responses
   (`Sample_data_response.xlsx`). A widely used alternative is **GSE176078**
   (breast cancer scRNA-seq, 44 patients, Wu et al. 2021). Any dataset with
   10,000+ cells and multiple patients will exercise the Seurat clustering,
   prediction, and visualization steps at realistic scale.
3. **Pre-trained models** — `load_model("abemaciclib")` (or any of the 44 FDA
   drugs) avoids the cost of training and lets you go straight to prediction.

Expected runtime: with DepMap + a 44-patient scRNA cohort, Seurat clustering and
clone-level prediction take several minutes and several GB of RAM. The demo data
can be used first to verify the whole pipeline works end-to-end.

## 7. Shiny Web Application

PERCEPTION-shiny ships with an interactive web dashboard (built with Shiny) that wraps
the whole pipeline — data loading, model training, prediction, and
visualization — in a point-and-click interface.

### 7.1 Launch

```r
# From the package source tree (development mode):
devtools::load_all()
run_perception_app()          # starts the app in your browser
```

or directly from the source tree:

```r
shiny::runApp(system.file("shiny", "app", package = "PERCEPTIONx"))
```

### 7.2 Tabs

| Tab | What you can do |
|-----|-----------------|
| **Data** | Load the synthetic demo data (smoke-test), load the full DepMap reference (~567 MB), or upload your own rank-normalized single-cell matrix + clinical responses |
| **Train** | Train drug-response models (`glmnet` / `random forest`) with tunable parameters |
| **Predict** | Score clone-level and patient-level drug sensitivity for any loaded model |
| **Visualize** | Clone distribution, clone-viability lollipop, ROC curve, response boxplot, model performance, and UMAP/t-SNE overlays |
| **Help** | In-app documentation |

### 7.3 Interactive plots

All figures are rendered as **interactive SVG** (via `ggiraph`): hover any point
or bar to see a tooltip (clone id, viability score, proportion, FPR/TPR, ...).
Because the original `ggplot` object is rendered directly — never converted to
`plotly` — facet layouts and legends stay exactly as designed: no re-flow, no
overlapping labels. Static downloads are publication-quality:

| Format | Resolution |
|--------|------------|
| PNG    | 600 dpi, Cairo anti-aliased (6000 × 4200 px at default size) |
| PDF / SVG | Vector — infinitely zoomable, recommended for papers |

### 7.4 Caching & cleanup

- Downloaded DepMap data and pre-trained models are written to the R session's
  `tempdir()`, so caches are **released automatically when the app closes** —
  nothing persists on disk between sessions.
- The app cache-busts its own stylesheet on every start, so you always see the
  newest UI without manually clearing the browser cache.

## 8. Citation

If you use this package, please cite the original PERCEPTION study:

Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. *Nat Cancer* 5, 938–952 (2024). https://doi.org/10.1038/s43018-024-00756-7

## 9. License

MIT © PERCEPTIONx authors
