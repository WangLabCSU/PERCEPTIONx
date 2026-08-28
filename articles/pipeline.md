# PERCEPTIONx: Complete Pipeline Tutorial

## 1. Overview

PERCEPTIONx implements the PERCEPTION approach (PERsonalized single-Cell
Expression-based Planning for Treatments In ONcology), which predicts
how individual patients respond to drug treatments by leveraging both
bulk and single-cell RNA sequencing data. This vignette walks through
the complete pipeline:

1.  **Load Data** — Download DepMap reference data and pre-trained
    models
2.  **Train Models** — Build drug response models from cell line data
3.  **Predict** — Score drug sensitivity at clone and patient levels
4.  **Evaluate** — Assess model performance and significance
5.  **Visualize** — Generate publication-quality figures

``` r

library(PERCEPTIONx)
```

## 2. Load Data

### 2.1 DepMap Reference Data

PERCEPTIONx relies on DepMap bulk expression, single-cell expression,
and drug response (AUC) data. The
[`load_depmap()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/load_depmap.md)
function downloads and caches these datasets automatically.

``` r

# Load DepMap data (first run downloads ~500MB)
load_depmap(read = TRUE, mirror = TRUE)
```

The `mirror = TRUE` argument attempts GitHub mirrors for faster
downloads in regions with limited connectivity. You can also add custom
mirrors:

``` r

add_mirrors("https://my-mirror.example.com/PERCEPTIONx")
```

### 2.2 Pre-trained Models

For quick exploration, you can load pre-trained models without training
from scratch:

``` r

# Load a single pre-trained model
models <- load_model("abemaciclib", read = TRUE)

# Load multiple drugs
models <- load_model(c("abemaciclib", "erlotinib"), read = TRUE)
```

## 3. Train Models

### 3.1 Feature Selection

Before training, identify genes available across both bulk and
single-cell expression datasets:

``` r

available_genes <- intersect(
  rownames(DepMap$expression_20Q4),
  rownames(DepMap$scRNA_complete)
)

# Sample genes of interest
set.seed(123)
GOI_100 <- sample(available_genes, 100)
```

Setting `GOI = NULL` uses all available genes (slower but more
thorough).

### 3.2 Training

The
[`train_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/train_models.md)
function performs feature ranking, model building, and hyperparameter
tuning in a single call:

``` r

models <- train_models(
  drug_list         = c("erlotinib", "gefitinib"),
  cancer_type       = "PanCan",
  exclude_cancer    = "PanCan",
  GOI               = NULL,
  model_type        = "glmnet",
  k_features_values = NULL,  # auto-compute
  ncores            = 1,
  output_dir        = NULL   # don't save to disk
)
```

**Key parameters:**

| Parameter        | Description                                        |
|------------------|----------------------------------------------------|
| `drug_list`      | Drug name(s) to train models for                   |
| `cancer_type`    | Cancer type for training cell lines                |
| `exclude_cancer` | Cancer type to exclude (leave-one-out validation)  |
| `GOI`            | Genes of interest (NULL = use all)                 |
| `model_type`     | `"glmnet"` (elastic net) or `"rf"` (random forest) |
| `ncores`         | Number of CPU cores for parallel training          |

## 4. Predict Drug Response

> **Important: Rank Normalization**
>
> PERCEPTIONx models are trained on **rank-normalized** expression data.
> If you provide your own expression data (e.g., from scRNA-seq), you
> **must** normalize it first using
> [`rank_normalization_mat()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/rank_normalization_mat.md),
> or predictions will be unreliable.
>
> **How rank normalization works**: For each cell (column), every gene’s
> expression value is replaced by its rank within that column, divided
> by the total number of genes: `x_norm = rank(x) / n`. This transforms
> each column into a uniform distribution over (0, 1\], making the data
> robust to batch effects, library size differences, and outliers. Since
> the model coefficients capture the relationship between **relative
> gene expression ranks** and drug response (not absolute values), the
> same normalization must be applied to any new data.
>
> ``` r
>
> # If your data is NOT already rank-normalized:
> my_expr_norm <- rank_normalization_mat(my_raw_expr)
> # Then use my_expr_norm in predict_drugs()
> ```

### 4.1 Clone-Level Prediction

[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md)
scores each cell/clone’s drug sensitivity from a rank-normalized
expression matrix:

``` r

# Get test cell line expression
cellLines_test <- get_cellLine_list(
  infunc_cancerType = "PanCan",
  infunc_drugName   = "erlotinib",
  exclude_cancer    = "PanCan",
  infunc_response   = get_response_matrix("erlotinib")
)[[2]]

test_cells <- DepMap$metadata_CPM_scRNA$NAME[
  DepMap$metadata_CPM_scRNA$DepMap_ID %in% cellLines_test]
expr_test <- DepMap$CPM_scRNA_CCLE_rnorm[, test_cells, drop = FALSE]

# Predict
cell_pred <- predict_drugs(models, expr_test)
```

The result is a matrix with cells as rows and drugs as columns, where
lower values (lower viability) indicate higher drug sensitivity.

### 4.2 Patient-Level Aggregation

[`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md)
aggregates clone-level scores into patient-level predictions using clone
proportions:

``` r

# Build clone_viability_df with patient and clone_id columns
clone_viability_df <- data.frame(
  patient  = sub("_.*", "", rownames(clone_viability)),
  clone_id = rownames(clone_viability),
  clone_viability,
  check.names = FALSE
)

# Patient-level prediction
patient_pred <- predict_patients(
  clone_viability_matrix = clone_viability_df,
  clone_counts = clone_counts,
  mode = "weighted_average"
)
```

**Aggregation modes:**

| Mode               | Description                                 |
|--------------------|---------------------------------------------|
| `weighted_average` | Weighted mean of clone scores by proportion |
| `min`              | Most sensitive clone (lowest viability)     |
| `max`              | Least sensitive clone (highest viability)   |
| `weighted_max`     | Weighted maximum across clones              |

## 5. Evaluate Models

### 5.1 Performance Comparison

``` r

# Compare performance across thresholds
perf <- compare_performance(models, threshold = 0.3, verbose = TRUE)

# Plot model performance
plot_model_performance(models)
```

### 5.2 Significant Models

``` r

# Filter models with significant stratification
sig_models <- get_significant_models(
  models,
  min_correlation = 0.3,
  max_pvalue = 0.05
)
```

### 5.3 Pseudo-Bulk Expression

``` r

# Compute pseudo-bulk for each patient
pseudo_bulk <- each_patient_pseudo_bulk(
  sc_expression = sc_expression_rnorm,
  patient_clone_map = clone_mapping
)
```

## 6. Visualize Results

PERCEPTIONx provides a comprehensive suite of plotting functions:

### 6.1 t-SNE / UMAP with Drug Response

``` r

plot_tsne_response(
  tsne_data  = tsne_data,
  color_var  = "viability_scaled",
  title      = "Drug Response",
  point_size = 0.5
)
```

### 6.2 Biomarker vs. Viability Side-by-Side

``` r

plot_tsne_biomarker_viability(
  tsne_data       = tsne_data,
  biomarker_var   = "biomarker_scaled",
  viability_var     = "viability_scaled",
  biomarker_label = "Biomarker Expression",
  viability_label   = "Drug Viability"
)
```

### 6.3 Clone Distribution

``` r

plot_clone_distribution(
  clone_distribution = clone_dist_df,
  response_var       = "response"  # optional facet
)
```

### 6.4 Clone Viability Lollipop Plot

``` r

plot_clone_viability(
  clone_viability = clone_viability_df,
  viability_var   = "comb_viability",
  weights_var   = "weights"  # optional: point size by clone proportion
)
```

### 6.5 Response Boxplot

``` r

plot_response_boxplot(
  exp_vs_pred   = exp_vs_pred,
  response_var  = "response",
  predicted_var = "predicted_viability"
)
```

### 6.6 ROC Curve

``` r

plot_roc_curve(
  response     = response,
  predictor    = predictor,
  smooth_curve = TRUE
)
```

### 6.7 Seurat Clustering

``` r

result <- plot_seurat_clustering(expr_matrix)
print(result$umap_plot)
```

### 6.8 Composite Patient Response Panel

``` r

plot_patient_response_panel(
  clone_distribution = clone_dist_df,
  clone_viability      = clone_viability_df,
  exp_vs_pred        = exp_vs_pred,
  viability_col        = "comb_viability"
)
```

### 6.9 Interactive Hover Tooltips

Every plotting function accepts `tooltip = TRUE` (default). When the
[`ggiraph`](https://cran.r-project.org/package=ggiraph) package is
installed, the plot is built with interactive layers so that hovering a
point or bar shows a rich tooltip (clone id, viability score,
proportion, FPR/TPR, …):

``` r

plot_clone_viability(
  clone_viability = clone_viability_df,
  viability_var   = "comb_viability",
  weights_var   = "weights",
  tooltip       = TRUE
)
```

Set `tooltip = FALSE` to obtain a plain static `ggplot` with the
identical layout — useful when the object is passed to packages that
only accept standard `ggplot` objects. Tooltip text is auto-built from
the plot data; override it with the `tooltip_col` argument when a custom
text column already exists.

## 7. Shiny Web Application

PERCEPTION-shiny is the interactive web dashboard (built with Shiny)
that wraps the whole pipeline — data loading, model training,
prediction, and visualization — in a point-and-click interface.

Heavy computation (training, clustering, prediction, plotting) runs in
background worker processes, so the UI never blocks other users; the
main process keeps only lightweight DepMap metadata. See
`vignettes/shiny_app.Rmd` for the async architecture and the
`PERCEPTION_WORKERS` / `PERCEPTION_WORKER_IDLE_MINUTES` deployment
options.

### 7.1 Launch

``` r

library(PERCEPTIONx)
run_perception_app()          # starts the app in your browser
```

or directly from the source tree:

``` r

shiny::runApp(system.file("shiny", "app", package = "PERCEPTIONx"))
```

### 7.2 Tabs

| Tab | What you can do |
|----|----|
| **Data** | Load the synthetic demo data (smoke-test), load the full DepMap reference (~567 MB), or upload your own rank-normalized single-cell matrix + clinical responses |
| **Train** | Train drug-response models (`glmnet` / random forest) with tunable parameters |
| **Predict** | Score clone-level and patient-level drug sensitivity for any loaded model |
| **Visualize** | Clone distribution, clone-viability lollipop, ROC curve, response boxplot, model performance, and UMAP/t-SNE overlays |
| **Help** | In-app documentation |

### 7.3 Interactive plots

All figures are rendered as **interactive SVG** (via `ggiraph`): hover
any point or bar to see a tooltip (clone id, viability score,
proportion, FPR/TPR, …). Because the original `ggplot` object is
rendered directly — never converted to `plotly` — facet layouts and
legends stay exactly as designed: no re-flow, no overlapping labels.
Static downloads are publication-quality:

| Format    | Resolution                                                   |
|-----------|--------------------------------------------------------------|
| PNG       | 600 dpi, Cairo anti-aliased (6000 × 4200 px at default size) |
| PDF / SVG | Vector — infinitely zoomable, recommended for papers         |

### 7.4 Caching & cleanup

- Downloaded DepMap data and pre-trained models are written to the R
  session’s [`tempdir()`](https://rdrr.io/r/base/tempfile.html), so
  caches are **released automatically when the app closes** — nothing
  persists on disk between sessions.
- The app cache-busts its own stylesheet on every start, so you always
  see the newest UI without manually clearing the browser cache.

## 8. Complete Pipeline Script

A ready-to-run end-to-end pipeline script lives in the source repository
under `tools/` (it is not shipped with the installed package). After
cloning the repo, run:

``` r

source("tools/test_pipeline.R")
```

This script executes the full workflow (train -\> predict -\> plot) and
saves all output figures to a local directory.

## 9. Citation

If you use this package, please cite:

Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts
patient response and resistance to treatment using single-cell
transcriptomics of their tumors. *Nat Cancer* 5, 938–952 (2024).
<https://doi.org/10.1038/s43018-024-00756-7>
