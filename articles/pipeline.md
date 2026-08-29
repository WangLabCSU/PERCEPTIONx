# PERCEPTIONx: Complete Pipeline Tutorial

## 1. Overview

PERCEPTIONx predicts how individual patients respond to cancer drugs.
Models are trained on DepMap cell-line screens and applied to a
patient’s single-cell expression profile, giving clone-level viability
scores and patient-level response stratification.

The core workflow is short — five function calls:

1.  Load the DepMap reference data
2.  Load pre-trained models (no training needed)
3.  Prepare the data with
    [`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
    (clustering + rank normalization)
4.  Predict at clone and patient level
5.  Plot the results

Everything else — training new models, evaluating them, validating on
DepMap cell lines, and the full plotting suite — is optional and covered
in section 3 (Advanced Topics). If you only want prediction scores,
section 2 is all you need.

``` r

library(PERCEPTIONx)
```

------------------------------------------------------------------------

## 2. Minimal Workflow

------------------------------------------------------------------------

### 2.1 Load the DepMap reference data

[`load_depmap()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/load_depmap.md)
downloads and caches the DepMap bulk expression, single-cell expression,
drug response (AUC), and cell line annotations. The first run downloads
about 500 MB; afterwards the cache is reused.

``` r

load_depmap(read = TRUE, mirror = TRUE)
```

The `mirror = TRUE` argument uses GitHub mirrors for faster downloads in
regions with limited connectivity. You can also add custom mirrors:

``` r

add_mirrors("https://my-mirror.example.com/PERCEPTIONx")
```

------------------------------------------------------------------------

### 2.2 Load pre-trained models

PERCEPTIONx ships 44 pre-trained models for FDA-approved drugs, so you
can predict without training:

``` r

# Load a single model
models <- load_model("abemaciclib", read = TRUE)

# Or several at once
models <- load_model(c("abemaciclib", "erlotinib"), read = TRUE)
```

If you want to train models for other drugs or settings instead, see
section 3.1.

------------------------------------------------------------------------

### 2.3 Prepare the clone-level expression data

PERCEPTIONx models are trained on rank-normalized expression, and
clone-level prediction needs clones defined first.
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
does both: it clusters the cells with Seurat (UMAP or t-SNE), merges
each clone’s expression, and rank-normalizes it — one call:

``` r

prepared <- prepare_data(my_expression_matrix, method = "umap")
```

The returned list carries everything downstream needs:
`clone_expression_rnorm` (the rank-normalized clone-level matrix) and
`clone_counts` (clone abundances per patient). Rank normalization
happens inside
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
— it calls
[`rank_normalization_mat()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/rank_normalization_mat.md)
on the merged clone expression — so you do not run it yourself on this
path. See section 3.5 for the details and the `skip_clustering` option.

------------------------------------------------------------------------

### 2.4 Predict drug sensitivity

Prediction runs in two stages: clone-level scoring, then patient-level
aggregation.

#### 2.4.1 Clone-level prediction

[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md)
scores each clone against every loaded model, using the clone-level
expression from
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md).
Lower values mean lower viability, i.e. higher drug sensitivity.

``` r

clone_pred <- predict_drugs(models, prepared$clone_expression_rnorm)
```

The result is a matrix with clones as rows and drugs as columns.

#### 2.4.2 Patient-level aggregation

[`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md)
collapses clone scores to patients using the clone proportions carried
in the
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
output:

``` r

patient_pred <- predict_patients(clone_pred, prepared)
```

The default `mode = "weighted_max"` weights each clone by its proportion
and takes the most resistant clone as the patient’s score. Other
aggregation modes are documented in section 3.2.3.

------------------------------------------------------------------------

### 2.5 Visualize

Each plotting function takes the corresponding prediction output and
returns a `ggplot`. This example overlays predicted viability on a t-SNE
embedding:

``` r

plot_tsne_response(
  tsne_data = tsne_data,
  color_var = "viability_scaled",
  title     = "Drug Response"
)
```

The full set of plot functions — clone distribution, lollipop, ROC,
boxplot, UMAP overlays, and a composite panel — is listed in section
3.4.

------------------------------------------------------------------------

## 3. Advanced Topics

Everything in this section is optional. Read the parts you need; the
minimal workflow in section 2 does not depend on any of it.

------------------------------------------------------------------------

### 3.1 Train your own models

Training is only needed when you want drugs or settings not covered by
the 44 pre-trained models.
[`train_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/train_models.md)
presets every argument except the drug list, so the whole call is one
line:

``` r

models <- train_models(drug_list = "erlotinib")
```

That is all it takes — cancer type (PanCan), model type (glmnet),
cross-validation, and feature ranking are handled with sensible
defaults. The knobs you might actually touch:

| Argument | What it does | Default |
|----|----|----|
| `GOI` | Restrict features to a curated gene set | all shared genes |
| `cancer_type` / `exclude_cancer` | Cancer types to include / leave out (self-validation) | `"PanCan"` |
| `model_type` | `"glmnet"` (elastic net) or `"rf"` (random forest) | `"glmnet"` |
| `ncores` | Parallel cores | `4` |
| `output_dir` | Where fitted models are saved | `"./models"` |

One thing worth knowing: `output_dir` defaults to `"./models"`, meaning
models are written to disk. Set `output_dir = NULL` to keep them in
memory only.

The returned `models` object plugs straight into
[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md)
(section 2.4), and can be saved or re-loaded later:

``` r

saveRDS(models, "my_models.RDS")
models <- readRDS("my_models.RDS")
```

------------------------------------------------------------------------

### 3.2 Evaluate model performance

Two helpers answer how well the trained models stratify responders from
non-responders:

``` r

perf <- compare_performance(models, threshold = 0.3, verbose = TRUE)
plot_model_performance(models)

sig_models <- get_significant_models(
  models,
  min_correlation = 0.3,
  max_pvalue = 0.05
)
```

#### 3.2.1 Aggregation modes

[`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md)
offers the five clone-to-patient strategies tested in the paper (Fig.
2c):

| Mode | Description |
|----|----|
| `weighted_max` | Most resistant clone, weighted by its abundance (default; the paper’s best strategy, AUC = 0.83) |
| `weighted_average` | All clones weighted by their proportions |
| `average` | Unweighted mean across clones |
| `min` | Most sensitive clone (lowest viability) |
| `max` | Most resistant clone, unweighted |

#### 3.2.2 Pseudo-bulk expression

For bulk-level analyses, aggregate single-cell expression per patient:

``` r

pseudo_bulk <- each_patient_pseudo_bulk(
  sc_expression = sc_expression_rnorm,
  patient_clone_map = clone_mapping
)
```

------------------------------------------------------------------------

### 3.3 Validate predictions on DepMap cell lines

This is a verification exercise, not part of the prediction path: it
predicts viability for held-out DepMap cell lines and checks the scores
against known drug response. It is the pattern used to produce the
validation figures in the paper.

``` r

cellLines_test <- get_cellLine_list(
  infunc_cancerType = "PanCan",
  infunc_drugName   = "erlotinib",
  exclude_cancer    = "PanCan",
  infunc_response   = get_response_matrix("erlotinib")
)[[2]]

test_cells <- DepMap$metadata_CPM_scRNA$NAME[
  DepMap$metadata_CPM_scRNA$DepMap_ID %in% cellLines_test]
expr_test <- DepMap$CPM_scRNA_CCLE_rnorm[, test_cells, drop = FALSE]

cell_pred <- predict_drugs(models, expr_test)
```

For your own patient data, skip this block entirely and go straight to
[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md)
as in section 2.4.

------------------------------------------------------------------------

### 3.4 The plotting suite

Every plot function accepts `tooltip = TRUE` (default): when the
`ggiraph` package is installed, points and bars get hover tooltips
(clone id, viability score, proportion, FPR/TPR). Set `tooltip = FALSE`
for a plain static `ggplot` with the identical layout. Tooltip text is
auto-built from the plot data; override it with the `tooltip_col`
argument when a custom text column already exists.

``` r

# t-SNE / UMAP with drug response overlay
plot_tsne_response(tsne_data, color_var = "viability_scaled", title = "Drug Response")

# Biomarker vs. viability side-by-side
plot_tsne_biomarker_viability(
  tsne_data, biomarker_var = "biomarker_scaled", viability_var = "viability_scaled",
  biomarker_label = "Biomarker Expression", viability_label = "Drug Viability"
)

# Clone abundance stacked bar chart
plot_clone_distribution(clone_distribution = clone_dist_df, response_var = "response")

# Clone-level viability lollipop
plot_clone_viability(
  clone_viability = clone_viability_df,
  viability_var = "comb_viability",
  weights_var = "weights"          # optional: point size by clone proportion
)

# Responder vs. non-responder boxplot
plot_response_boxplot(
  exp_vs_pred = exp_vs_pred,
  response_var = "response",
  predicted_var = "predicted_viability"
)

# ROC curve with AUC
plot_roc_curve(response = response, predictor = predictor, smooth_curve = TRUE)

# Seurat clustering and UMAP
result <- plot_seurat_clustering(expr_matrix)
print(result$umap_plot)

# Composite patient response panel
plot_patient_response_panel(
  clone_distribution = clone_dist_df,
  clone_viability = clone_viability_df,
  exp_vs_pred = exp_vs_pred,
  viability_col = "comb_viability"
)
```

Interactive example:

``` r

plot_clone_viability(
  clone_viability = clone_viability_df,
  viability_var = "comb_viability",
  weights_var = "weights",
  tooltip = TRUE
)
```

------------------------------------------------------------------------

### 3.5 `prepare_data()` and Seurat clustering

The
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
function turns a raw expression matrix into the clone-level inputs used
by prediction and visualization: it clusters cells with Seurat (UMAP or
t-SNE), assigns clones, merges each clone’s expression, and
rank-normalizes it (internally calling
[`rank_normalization_mat()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/rank_normalization_mat.md)),
returning a list with `clone_expression_rnorm` and `clone_counts`.

``` r

prepared <- prepare_data(expr_matrix, method = "umap", resolution = 0.5)
clone_pred <- predict_drugs(models, prepared$clone_expression_rnorm)
patient_pred <- predict_patients(clone_pred, prepared)
```

Useful options: `resolution` (Seurat clustering resolution),
`genes_to_use` (restrict the gene space before clustering), and
`skip_clustering = TRUE` — if your data is already clone-level (one
column per clone),
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
skips clustering and only builds the normalized clone matrix and counts.

Clone identity is carried in the row names as `Patient@@Clone`;
[`parse_clone_keys()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/parse_clone_keys.md)
splits them, and
[`clone_mean_expression()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/clone_mean_expression.md)
/
[`build_clone_counts()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/build_clone_counts.md)
work on the cloned data for downstream analyses.

------------------------------------------------------------------------

### 3.6 The Shiny web application

PERCEPTION-shiny wraps the whole pipeline in a point-and-click
dashboard. Heavy computation (training, Seurat clustering, prediction,
plotting) runs in background worker processes, so the interface stays
responsive even while a large job runs. See the [Shiny App
Guide](https://wanglabcsu.github.io/PERCEPTIONx/articles/shiny_app.html)
for a full walkthrough.

``` r

library(PERCEPTIONx)
run_perception_app()          # starts the app in your browser
```

or directly from the source tree:

``` r

shiny::runApp(system.file("shiny", "app", package = "PERCEPTIONx"))
```

The app has five tabs — Data, Train, Predict, Visualize, and Help — and
ships with a Load Demo button that generates a small synthetic dataset
(49 genes x 400 cells x 20 patients) to smoke-test the whole flow.
Figures in the app are interactive SVG (ggiraph-based), with static
downloads available as 600 dpi PNG or vector PDF/SVG.

------------------------------------------------------------------------

### 3.7 End-to-end pipeline script

A ready-to-run script that executes the full workflow (train -\> predict
-\> plot) and saves all output figures lives in the source repository
under `tools/` (not shipped with the installed package). After cloning
the repo:

``` r

source("tools/test_pipeline.R")
```

------------------------------------------------------------------------

## 4. Citation

If you use this app or package, please cite both the package and the
original methodology paper:

- **Jia Ding**. PERCEPTIONx: Personalized Drug Response Prediction from
  Single-Cell Transcriptomics. R package version 0.1.0.
  <https://github.com/WangLabCSU/PERCEPTIONx>
- **Sinha, S., Vegesna, R., Mukherjee, S.** *et al.* PERCEPTION predicts
  patient response and resistance to treatment using single-cell
  transcriptomics of their tumors. *Nature Cancer* 5, 938–952 (2024).
  DOI:
  [10.1038/s43018-024-00756-7](https://doi.org/10.1038/s43018-024-00756-7)

Repository:
[github.com/WangLabCSU/PERCEPTIONx](https://github.com/WangLabCSU/PERCEPTIONx)

Feedback: <jiading682@qq.com>
