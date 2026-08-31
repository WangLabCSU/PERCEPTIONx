---
name: "perceptionx-result"
description: "Interprets PERCEPTIONx analysis outputs: clone/patient-level drug response predictions, model performance metrics, ROC curves and result plots. Invoke when the user asks to explain, visualize, or interpret PERCEPTIONx results."
---

# PERCEPTIONx Result Interpretation

Interpret and visualize outputs produced by the **PERCEPTIONx** R package. Use the functions below to inspect, compare, and plot results — do not invent outputs that the package does not produce.

## What PERCEPTIONx produces

### 1. Clone-level predictions — `predict_drugs(model_list, expr)`
- Matrix / data.frame: **rows = patient-clone (or sample), columns = drugs**
- Values = predicted viability (rank-normalized z-score scale); lower = more sensitive to the drug

### 2. Patient-level aggregation — `predict_patients(clone_pred, prepared, mode = "weighted_max", zscore = TRUE)`
- data.frame: **rows = patients, columns = drugs**
- `weighted_max` weights each clone by its observed count (recommended); `"max"` takes the raw max
- Check `prepared$clone_counts` when interpreting: patients with few clones have less certain estimates

### 3. Model performance — `compare_performance(model_list, threshold = 0.3, verbose = TRUE)`
Returns a list with:
- `perf_cv` — cross-validation performance (|correlation| = sqrt(R²)) per drug
- `perf_bulk` / `perf_pseudo_bulk` / `perf_scRNA` — held-out performance on bulk / pseudo-bulk / single-cell RNA
- `summary` — statistics for models passing the correlation threshold

Related: `get_performance(filepath)` reads metrics from a saved model RDS; `get_significant_models(model_list, min_correlation, max_pvalue)` filters models passing thresholds.

## Interpreting predictions

- **Lower predicted viability → higher drug sensitivity.** Rank drugs per patient/clone to prioritize candidates.
- **Direction check**: if higher values are expected to mean "more sensitive" for a given dataset, confirm the sign convention before reporting — the package predicts viability, so lower = better response.
- **Clone heterogeneity**: a wide spread of clone-level predictions within one patient means intra-tumoral heterogeneity — the patient-level number averages it via clone weights.
- **Model quality caveat**: predictions from models with low `compare_performance` correlation (below ~0.3) are unreliable — say so when reporting.

## Plotting functions (choose by what the user wants to see)

| Goal | Function |
|---|---|
| ROC curve + AUC badge (needs 2 response classes) | `plot_roc_curve(response, predictor, smooth_curve = TRUE)` |
| Response boxplot observed vs predicted | `plot_response_boxplot(exp_vs_pred)` |
| Clone distribution across patients (stacked) | `plot_clone_distribution(clone_distribution, response_var = ...)` |
| Clone UMAP/t-SNE colored by clone | `plot_clone_umap(tsne_data, clone_col = "clone_id")` |
| Clone viability by patient (dot plot) | `plot_clone_viability(clone_viability, viability_var = "comb_viability", drug = ...)` |
| t-SNE colored by predicted response | `plot_tsne_response(tsne_data, color_var = ...)` |
| t-SNE + biomarker + viability panel | `plot_tsne_biomarker_viability(...)` |
| Seurat clustering overview | `plot_seurat_clustering(method, expression_matrix, ...)` |
| Model performance comparison | `plot_model_performance(performance_list, ...)` / `plot_model_roc(performance_list, ...)` — both take the output of `compare_performance()` (or a named list of per-drug performance entries) |
| Per-patient response panel (3-panel: distribution + viability + obs-vs-pred) | `plot_patient_response_panel(clone_distribution, clone_viability, exp_vs_pred, response_col = "response", viability_col = "comb_viability")` |

All `plot_*` functions accept `tooltip = TRUE` when **ggiraph** is installed (interactive tooltips), plus `base_size` for font scaling.

## Saving plots

```r
export_plot_cairo(file = "roc.png", plot = p, format = "png",
                  width = 7, height = 5, res = 600)
# format: "png" / "pdf" / "tiff" / "jpeg"; res = 600 for publication-quality
```

## Workflow for a typical "interpret my results" request

1. Identify what the user has: prediction table, model list, plot object, or model RDS files.
2. Load models if needed: `load_model(read = TRUE, ...)` or `get_performance("path/to/model.RDS")`.
3. Reproduce the key numbers (best/worst drugs per patient, AUC, correlation) with the functions above.
4. Generate 1-3 targeted plots (ROC for discrimination, clone distribution for heterogeneity, response boxplot for calibration).
5. Summarize in plain language: which drugs look most promising, for which patients, and how trustworthy the models are.

## Rules

- Never claim an output the package cannot produce; verify with the listed functions first.
- Always qualify patient-level conclusions with clone-count/heterogeneity caveats.
- Report model performance alongside any prediction so the user can judge reliability.
