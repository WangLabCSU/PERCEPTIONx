---
name: "perceptionx-pipeline"
description: "Drives the PERCEPTIONx R package through the full drug-response prediction workflow: DepMap load, patient data preparation, model load/training, and clone/patient-level prediction. Invoke when the user asks to run, script, or automate a PERCEPTIONx analysis in R."
---

# PERCEPTIONx Pipeline

Drive the **PERCEPTIONx** R package (tumor drug-response prediction from single-cell / bulk expression, based on DepMap reference models) through a complete analysis. Use these functions and conventions instead of guessing APIs.

## Capabilities

- Load the DepMap reference dataset (567 MB) and 44 pre-trained drug-response models (glmnet, built from CCLE/DepMap)
- Prepare patient expression into rank-normalized clone-level matrices (with or without Seurat clustering)
- Train new models from user data, or reuse the pre-trained ones
- Predict drug response at clone level, then aggregate to patient level
- Evaluate model performance and export plots
- Launch the interactive Shiny GUI (`run_perception_app()`) for point-and-click use

## Installation

```r
# Main install (GitHub via ghfast proxy — faster in China)
remotes::install_git("https://ghfast.top/https://github.com/WangLabCSU/PERCEPTIONx")
# or the direct GitHub route
remotes::install_github("WangLabCSU/PERCEPTIONx")
library(PERCEPTIONx)
```

## Shiny GUI entry

```r
run_perception_app()   # interactive: opens in the browser; on Shiny Server
                       # (SHINY_PORT set) it returns the app object instead
```

The Shiny app is the same pipeline behind a UI — prefer the package functions below when scripting/automating, and `run_perception_app()` when the user wants a GUI.

## Standard workflow

```r
library(PERCEPTIONx)

# 1. (Recommended) Pin cache location so DepMap/models persist across sessions
# options(PERCEPTIONx.cache_root = "/path/to/cache")  # -> <root>/depmap, <root>/models

# 2. Load reference data (mirror = TRUE is faster in China; GitHub direct is slow)
load_depmap(read = TRUE, mirror = TRUE)            # download (if needed) + load DepMap
depmap <- get_depmap()                              # retrieve the loaded DepMap object
load_model(read = TRUE, mirror = TRUE, all = TRUE)  # load all 44 pre-trained models
# or load specific drugs: load_model("cisplatin", "paclitaxel", read = TRUE)

# 3. Prepare patient data
prepared <- prepare_data(
  method = c("umap", "tsne"),   # Seurat reduction; irrelevant when skip_clustering = TRUE
  expression_matrix = expr,     # genes as ROWS, samples (cells or clones) as COLUMNS
  patient_mapping   = mapping,  # data.frame with cell_id + patient_id columns
  skip_clustering   = FALSE     # TRUE when each column is already one clone
)
# prepared$clone_expression_rnorm -> genes x patient_clone matrix (input for predict)
# prepared$clone_counts          -> clone abundance per patient (input for aggregation)

# 4. Predict
clone_pred <- predict_drugs(model_list, prepared$clone_expression_rnorm)
patient_pred <- predict_patients(clone_pred, prepared, mode = "weighted_max", zscore = TRUE)

# 5. Evaluate & plot
perf <- compare_performance(model_list)          # CV/bulk/pseudo-bulk/scRNA metrics
p <- plot_roc_curve(response = response_vec, predictor = pred_vec)
```

## Minimal end-to-end example

Quick smoke test with the pre-trained models and simulated data (no DepMap needed):

```r
library(PERCEPTIONx)
load_model(read = TRUE, all = TRUE)              # 44 pre-trained models
set.seed(42)
expr <- matrix(rnorm(200 * 20), nrow = 200,      # 200 genes x 20 clones
               dimnames = list(paste0("G", 1:200), paste0("cl", 1:20)))
clone_pred <- predict_drugs(model_list, expr)    # clones x drugs predicted viability
head(clone_pred)
```

## Mirror management

Downloads default to GitHub (slow in China). Configure mirrors once:

```r
get_mirrors()          # current mirror list
list_mirrors()         # print them numbered
add_mirrors("https://ghfast.top/https://github.com", position = "first")
reset_mirrors()        # restore defaults
```

## Data formats (validate before running)

- **Expression matrix**: numeric, genes as rows, samples as columns. Matrix or data.frame.
- **patient_mapping**: data.frame with `cell_id` and `patient_id` columns (case-insensitive). In clone-level mode an optional `count` column (real cell number per clone) makes clone proportions accurate; without it proportions fall back to equal 1/n.
- **Clinical response**: vector/data.frame mapping each patient to response (e.g. R/NR, 0/1) — needed only for ROC/performance evaluation.

## Key parameters

| Function | Key args | Notes |
|---|---|---|
| `load_depmap` | `dest`, `read`, `mirror`, `force`, `timeout_seconds`, `retries` | `force=TRUE` deletes a stale cached RDS and re-downloads; `mirror=TRUE` tries domestic mirrors first |
| `load_model` | `...` (drug names), `all`, `read`, `mirror`, `force`, `dest` | `all=TRUE` downloads/trains all 44 models in one call and overrides `...` |
| `prepare_data` | `method`, `expression_matrix`, `patient_mapping`, `skip_clustering`, `seurat_resolution`, `seurat_dims`, `seurat_seed` | Returns list with `clone_expression_rnorm`, `clone_counts`, `cell_clone_map` |
| `train_models` | `drug_list`, `cancer_type`, `ncores`, `output_dir`, `model_type` ("glmnet"/"randomForest"), `num_folds`, `seed` | Long-running; sets `progress_cb` optionally |
| `predict_drugs` | `model_list`, `expr` | `expr` = genes x samples; returns clone-level predictions |
| `predict_patients` | `clone_pred`, `prepared_data`, `mode`, `zscore` | `mode="weighted_max"` (recommended, weighted by clone counts) or `"max"` |

## Common pitfalls

- **GitHub download slow/fails** → use `mirror = TRUE` (uses configured mirrors incl. ghfast.top); list mirrors with `get_mirrors()`, reset with `reset_mirrors()`.
- **Stale/corrupt cached RDS blocks re-download** → `force = TRUE`.
- **`get_depmap()` errors "not loaded"** → you called it before `load_depmap(read = TRUE)`.
- **ROC needs two response classes** → a single-class response cannot build a ROC curve.
- **predict_patients expects the prepared list** → pass `prepared` (from `prepare_data()`), not a raw matrix.
- **Patient aggregation**: prefer `mode = "weighted_max"` over `"max"` — it weights clones by their observed counts.

## Rules

- Never fabricate function names or arguments; call the ones listed here.
- Always validate input formats before running (rows = genes, columns = samples).
- For reproducibility set `set.seed()` before `train_models` / `prepare_data` (Seurat uses `seurat_seed`, default 42).
- Long tasks (download, Seurat clustering, multi-drug training) may take minutes — report progress to the user rather than silent waits.
