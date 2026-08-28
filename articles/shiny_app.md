# PERCEPTION-shiny User Guide

## PERCEPTION-shiny User Guide

**PERCEPTION-shiny** is the interactive web application (Shiny) of the
**PERCEPTIONx** R package. It wraps the complete analysis pipeline —
data loading, model training, drug-sensitivity prediction, and result
visualization — into a point-and-click interface, so you can go from a
patient single-cell expression matrix to clone-level viability scores
and patient-level response stratification without writing code.

> Methodological basis: PERCEPTION (PERsonalized single-Cell
> Expression-based Planning for Treatments In ONcology), which trains
> elastic-net models on DepMap cell-line data to predict patient
> response and resistance to drug treatment.

------------------------------------------------------------------------

### Contents

- [0. Requirements & Launch](#id_0-requirements--launch)
- [1. Interface Overview](#id_1-interface-overview)
- [2. Data Tab: Loading Data](#id_2-data-tab-loading-data)
- [3. Train Tab: Training Models
  (Optional)](#id_3-train-tab-training-models-optional)
- [4. Predict Tab: Predicting Viability
  Scores](#id_4-predict-tab-predicting-viability-scores)
- [5. Visualize Tab](#id_5-visualize-tab)
- [6. Help Tab](#id_6-help-tab)
- [7. FAQ](#id_7-faq)
- [8. Citation & Contact](#id_8-citation--contact)

------------------------------------------------------------------------

### 0. Requirements & Launch

#### 0.1 Requirements

- **R** ≥ 4.1.0
- Main dependencies: `devtools`, `shiny`, `bslib`, `Seurat`, `ggplot2`,
  `ggiraph`, `glmnet`, `caret`, `DT`, `plotly`, `waiter`, `thematic`,
  `callr`, `readxl` (the app prompts you to install any that are missing
  when a feature needs them)

#### 0.2 Launch

From the package source root, run:

``` r

devtools::load_all()          # load PERCEPTIONx from source
run_perception_app()          # launch the app (opens your browser)
```

Or run the app directory directly:

``` r

shiny::runApp("inst/shiny/app")
```

#### 0.3 Overall Flow

    DepMap reference data ──► Model training ──► Clone/patient prediction ──► Visualization & validation
              ▲                      ▲                    ▲                        ▲
          Data tab              Train tab          Predict tab               Visualize tab
     (or load pre-trained   (or skip training —
       models directly)       use the 44 pre-
                              trained models)

> For the fastest result, follow **Load Demo → Predict → Visualize**. No
> training needed.

> **Async architecture (multi-user friendly)**: all heavy computation —
> model training, Seurat clustering, prediction, and plot math — runs in
> **background worker processes**, never in the interface. The UI polls
> and shows progress, so one user’s large task never freezes other
> users.
>
> - **Training**: the standard DepMap uses a **shared master** (one
>   global background process holding a single in-memory copy of DepMap;
>   on Linux, concurrent jobs share that one copy via fork; after 12 h
>   idle it exits to release memory). Uploaded DepMaps run in **isolated
>   per-session workers** so a bad upload cannot affect others.
> - **Clustering / prediction / plots / Load Demo**: a light per-session
>   worker handles each, writing results back to files that the UI picks
>   up.
> - **Env vars** (deployment): `PERCEPTION_WORKERS` (shared-pool
>   parallelism, default 16), `PERCEPTION_WORKER_IDLE_MINUTES` (master
>   idle-exit minutes, default 720).

------------------------------------------------------------------------

### 1. Interface Overview

![PERCEPTION-shiny home page](../reference/figures/shiny-home.png)

PERCEPTION-shiny home page

The top navbar has 6 tabs: **Home, Data, Train, Predict, Visualize,
Help**. Data flows left to right: load data (Data), train/load models
(Train), predict (Predict), then visualize and validate (Visualize).

The **Home** page includes: an introduction, a four-step guide (Load
Data → Train Model → Predict → Visualize; click any step to jump to the
corresponding tab), a live data-status overview, key features, and
citation info. The **Quick Start** and **Load Demo** buttons load the
demo data in one click.

------------------------------------------------------------------------

### 2. Data Tab: Loading Data

The Data tab is the entry point. It loads four kinds of data: **demo
data / DepMap reference data / your expression matrix / clinical
responses**. Each item turns its status badge green once loaded.

> **Load Demo**
>
> Click **Load Demo** to generate a synthetic demo dataset on the fly:
> **49 genes × 400 cells × 20 patients**, automatically clustered
> (Seurat), rank-normalized, and then used to train demo models. Great
> for smoke-testing the whole flow and getting familiar with the
> interactions.

#### 2.1 Loading DepMap Reference Data (required for training)

Two ways:

1.  **Download & Load**: downloads the DepMap reference set (~567 MB,
    15k+ genes × 1,000+ cell lines) from the official mirror and loads
    it automatically. This is the standard training input and the most
    memory/disk-demanding step.
2.  **Upload a local .RDS**: if you already have the DepMap file
    (`DepMap.RDS`), browse and select it — it loads automatically.

> **About memory**: the interface process reads only DepMap **metadata**
> (gene names, drug list, component dimensions — a few hundred KB). The
> full multi-GB object is loaded by a **background worker** in a
> separate process for training only, so concurrent users do not each
> hold an 8 GB copy (see the architecture note in 0.3).

#### 2.2 Loading Models

Two ways:

1.  **Download & Load**: one-click download of the 44 FDA-approved drug
    pre-trained models (e.g. `abemaciclib`, `erlotinib`, `osimertinib`).
    Multi-selected items each show an × to remove individually
2.  **Upload a local .RDS**: select a trained model file (from
    [`train_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/train_models.md)
    or exported from the Train tab) — loads automatically

#### 2.3 Uploading an Expression Matrix (patient scRNA-seq)

- **Accepted formats**: CSV / TSV / TXT / Excel (.xlsx/.xls) / RDS
- **Accepted R objects** (RDS): a numeric matrix or data.frame (genes ×
  cells); Seurat objects are not accepted directly — export the matrix
  first
- **Orientation**: genes × cells (rows = genes, columns = cells). If the
  first column is a character column of gene names (common in Excel/CSV
  exports), it is converted to row names automatically
- **Normalization**: expression should be **rank-normalized**. If you
  upload raw counts, the app can normalize them automatically during
  clustering

#### 2.4 Uploading the Cell-to-Patient Map (Mapping)

- **Accepted formats**: CSV / TSV / TXT / Excel / RDS
- **Required columns**: `cell_id` (cell names, matching the expression
  matrix column names) and `patient_id` (the patient each cell belongs
  to). Column names are case-insensitive (`Patient`, `PATIENT` all work)
- **Named list RDS**: if the RDS is a named list of “patient → cell-name
  vector” (e.g. the paper demo’s `PRJNA591860_sample_cell_names.RDS`),
  it is converted to long format automatically; empty samples are
  dropped

#### 2.5 Uploading Clinical Responses (Response, optional but recommended)

- **Accepted formats**: CSV / TSV / TXT / Excel / RDS
- **Required columns**: `patient` (must exactly match `patient_id` in
  the Mapping) and `response` (the patient’s true clinical outcome for
  the drug)
- **Labels are normalized automatically**: `responder` / `responsive` /
  `r` / `sensitive` → **Responder**; `non-responder` / `nr` /
  `resistant` → **Non-responder**
- A check runs on upload: if the patient IDs do not overlap with the
  Mapping at all, a warning is shown
- **Why it’s needed**: only with true outcomes can you validate
  predictions (ROC curves, responder vs. non-responder boxplots). Skip
  it if you only want prediction scores

> **Example** (paper demo lung cohort PRJNA591860): label by treatment
> timepoint, keeping the groups as-is (e.g. `TN` / `RD` / `PD`) — no
> collapsing into two classes. When a binary label is needed for ROC,
> choose two groups to compare on the Visualize tab.

#### 2.6 Clustering (Seurat)

> **Don’t Forget to Run Seurat!**
>
> Clustering (Seurat) defines the cell subpopulations (“clones”) for
> downstream prediction and visualization. If you skip it, clone-level
> prediction and clone-level figures will have no input.

**Choose a clustering method**: UMAP / tSNE.

| Parameter | Description | Suggested |
|----|----|----|
| Seurat Resolution | clustering resolution; higher = finer clusters | 0.5–1.0 |
| Seurat Dims | number of PCA dimensions used | 10–30 |
| Seurat NFeatures | number of highly variable genes for clustering — fixed at 2000 in the app (not adjustable) | — |

![Clustering in the Data
tab](../reference/figures/shiny-data-clustering.png)

Clustering in the Data tab

> Notes: clustering only defines “clones”; downstream prediction is done
> at the clone level (the clone’s expression is the mean expression of
> its cells).

------------------------------------------------------------------------

### 3. Train Tab: Training Models (Optional)

**Prerequisite**: DepMap reference data must be loaded first (the top of
the page tells you what’s missing).

> If you only want to predict, **you don’t need to train** — just use
> the 44 pre-trained models from the Data tab. The Train tab is for:
> other drugs, other cancer types, custom gene sets, or inspecting model
> performance on the three validation datasets.

#### 3.1 Parameters

| Parameter | Description |
|----|----|
| **Drug Name** | free-text input; one drug per line or comma/space-separated (combination regimens supported) |
| **Cancer Type (include)** | cancer types to include, default PanCan (pan-cancer) |
| **Cancer Type (exclude)** | cancer types to exclude (to avoid self-validation) |
| **Gene Symbols** | gene list; leave empty = use all DepMap genes (recommended); paste text or upload .txt/.csv |
| **Top k Features** | keep the top-k ranked features in the model |
| **Algorithm** | elastic net `glmnet` (recommended) or random forest `rf` |
| **CPU Cores** | number of parallel cores |

#### 3.2 Interpreting the Results

- **Model Summary**: model type and hyperparameters per drug (glmnet
  alpha/lambda, or rf ntree/RMSE)
- **Performance Plot**: model performance curve (threshold
  vs. correlation)
- **Performance Metrics**: prediction–truth Pearson correlation and
  p-value on **Bulk / Pseudo-bulk / Single-cell** levels. Higher
  correlation + lower p-value = better model
- **Download Model (.RDS)**: export the trained model; re-upload it
  later on the Data tab

> Note: Performance Plot and Metrics rely on validation metrics produced
> during training and are only available for models from the Train tab /
> [`train_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/train_models.md).
> Pre-trained models
> ([`load_model()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/load_model.md))
> do not carry these fields; the app will tell you they’re not
> applicable.

------------------------------------------------------------------------

### 4. Predict Tab: Predicting Viability Scores

Select the loaded models (pre-trained from the Data tab or trained on
the Train tab) and click predict:

1.  **Clone-level**: viability score for every clone × every drug.
    Semantics: the model outputs **viability (survival)**, where
    **higher = more resistant, lower = more sensitive**. The prediction
    heatmap shows the model’s raw predictions (not normalized); the
    lollipop plot and UMAP Drug Viability use a z-score (centered at 0,
    can be negative)
2.  **Patient-level**: clone scores are aggregated to patients by clone
    proportion (default `weighted_max`), giving each patient’s
    drug-sensitivity stratification

Outputs include an interactive heatmap (clones × drugs, plotly) and
downloadable prediction tables.

![Predict tab clones×drugs
heatmap](../reference/figures/shiny-predict-heatmap.png)

Predict tab clones×drugs heatmap

------------------------------------------------------------------------

### 5. Visualize Tab

This module visualizes the prediction results, so **you must run a
prediction first**.

All figures are **interactive SVG** (built on ggiraph): hover any point
or bar to see details (clone id, viability score, proportion, FPR/TPR,
…). The in-figure toolbar has zoom disabled; export to PNG/PDF via the
download buttons at the top of the page.

#### 5.1 Clone Distribution (stacked bar chart)

![Clone distribution stacked bar
chart](../reference/figures/shiny-clone-distribution.png)

Clone distribution stacked bar chart

Shows the clone composition within each patient; one color band = one
clone (a curated palette is used for ≤ 15 clones).

> Note: clone identity depends on the data source — with global
> clustering, the same clone (color) is genuinely shared across
> patients; with clone-level input using per-patient labels
> (e.g. c1/c2/c3), **the same color across patients does not imply the
> same clone origin** — it is only a shared category label.

#### 5.2 Clone Viability (lollipop plot)

![Clone viability lollipop
plot](../reference/figures/shiny-lollipop.png)

Clone viability lollipop plot

- **Rules**: all samples, all clones, one facet per patient, one
  lollipop per clone
- **Color**: blue-white-red diverging — blue = predicted sensitive (low
  viability), red = predicted resistant (high viability). Color limits
  adapt to the data, so extreme z-scores keep a real color (never
  grey/NA).
- **Point size**: clone proportion (larger clones get bigger points)
- **Ordering**: within a patient, by proportion descending; responders
  come first when response data is present
- **Y-axis**: Predicted Viability (z-score), with a zero line as the
  reference

#### 5.3 ROC Curve

![ROC curve](../reference/figures/shiny-roc.png)

ROC curve

Uses true clinical responses (uploaded on the Data tab) against
patient-level prediction scores. AUC closer to 1 = stronger
stratification. If you predicted with multiple models, you can pick a
specific model to view.

#### 5.4 Response Boxplot (responders vs. non-responders)

![Response boxplot (R vs NR)](../reference/figures/shiny-boxplot.png)

Response boxplot (R vs NR)

Shows the distribution of prediction scores for Responder
vs. Non-responder groups, with a significance test.

#### 5.5 Model Performance

![Model performance
plot](../reference/figures/shiny-model-performance.png)

Model performance plot

This plot lives on the **Train tab** (see Performance Plot in 3.2), not
on the Visualize tab. It shows validation performance curves for trained
models (requires models from the Train tab; not applicable to
pre-trained models — the app will tell you).

#### 5.6 Spatial Plots (UMAP / t-SNE)

![Gene expression on UMAP
(SLC2A1)](../reference/figures/shiny-umap-gene.png)

Gene expression on UMAP (SLC2A1)

![Drug viability on UMAP
(erlotinib)](../reference/figures/shiny-umap-viability.png)

Drug viability on UMAP (erlotinib)

Choose a dimensionality-reduction method and a color variable:

- **Gene Expression**: per-cell expression of a single gene, winsorized
  to the 5th–95th percentiles and scaled to \[0, 1\] (range01
  normalization), continuous grey→red ramp (grey = no/low expression,
  red = high expression) — see which cell groups express the gene highly
- **Drug Viability**: each cell colored by its clone’s predicted
  viability score (“patchwork” blocks) — see which clones are predicted
  resistant (high) vs. sensitive (low)
- **Clone / Cluster**: color by clone (cluster) — see the population
  structure

> **How to read them together**
>
> If the cells with high Gene Expression are also bright (high
> viability) on the Drug Viability plot, the gene’s high expression is
> positively associated with resistance; if dark (low viability), with
> sensitivity. The two plots use different color scales — compare
> spatial patterns only, not values.
>
> For example, the two figures above give preliminary insight:
>
> 1.  **Resistance-marker clue**: the clone in the bottom-right region
>     shows high SLC2A1 expression, which overlaps exactly with the
>     region of high predicted erlotinib viability — suggesting SLC2A1
>     overexpression may mark this resistant clone.
> 2.  **Sensitive main population**: the major cell groups at the top
>     barely express SLC2A1, yet are precisely the region with the
>     lowest predicted drug viability.
>
> So in the PRJNA591860 dataset, locally abnormal SLC2A1 overexpression
> strongly points to erlotinib resistance, while most
> non-SLC2A1-expressing cells are highly drug-sensitive — preliminary
> single-cell evidence that **SLC2A1 may be a resistance target**.

#### 5.7 High-Resolution Download

| Format    | Resolution                                                    |
|-----------|---------------------------------------------------------------|
| PNG       | 600 dpi (Cairo anti-aliased), ~6000 × 4200 px at default size |
| PDF / SVG | Vector, infinitely zoomable, **recommended for publication**  |

------------------------------------------------------------------------

### 6. Help Tab

Built-in documentation covering every tab, the FAQ, and usage tips.
Check it first when you get stuck.

------------------------------------------------------------------------

### 7. FAQ

**Q1: “Maximum upload size exceeded”?** The upload limit is 1 GB. For
larger data, preprocess locally (e.g. subset genes) first; for DepMap,
prefer the “Download & Load” button (server-side download, not via the
browser).

**Q2: Model Performance says performance fields are missing?** That plot
needs validation metrics produced during training; pre-trained models
([`load_model()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/load_model.md))
don’t carry them. Train on the Train tab first, or load the demo models
on the Data tab.

**Q3: How do I construct the response labels?** Build a two-column table
`patient` / `response`. `patient` must exactly match `patient_id` in the
Mapping; `response` values are normalized automatically (responder →
Responder, etc.). Response labels keep multiple groups (e.g. TN/RD/PD) —
no timepoint collapsing; when ROC needs a binary label, pick two groups
to compare on the Visualize tab.

**Q4: Why is a viability score lower than I expected?** Scores are
relative ranks learned from DepMap (not normalized): the prediction
heatmap shows raw model values, while the lollipop plot and UMAP Drug
Viability use a z-score (centered at 0, can be negative). Higher values
mean the most resistant clone in this batch (highest viability), not an
absolute viability percentage. Clonal heterogeneity and activated
resistance pathways both raise the viability score.

**Q5: Does data persist after I close the app?** It persists. DepMap
data is cached in a persistent directory (Windows: the user data
directory; on Linux, set the `PERCEPTIONX_DEPMAP_CACHE_DIR` environment
variable), with a 12-hour unused-expiry (TTL) mechanism. Files stay on
disk after the app closes; if unused for more than 12 hours, the next
click deletes and re-downloads them.

**Q6: Why doesn’t the interface freeze during training / clustering /
prediction?** All heavy computation runs in **background worker
processes**; the interface only polls and shows progress, so you can
keep using other pages while a big task runs. If a worker ever stops
unexpectedly, the app reports “Background worker stopped” instead of
spinning forever — just submit again.

------------------------------------------------------------------------

### 8. Citation & Contact

**If you use this app/package, please cite the original methodology
paper:**

> Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts
> patient response and resistance to treatment using single-cell
> transcriptomics of their tumors. *Nature Cancer* 5, 938–952 (2024).
> DOI:
> [10.1038/s43018-024-00756-7](https://doi.org/10.1038/s43018-024-00756-7)

**Repository**:
[github.com/WangLabCSU/PERCEPTIONx](https://github.com/WangLabCSU/PERCEPTIONx)

**Feedback**: <jiading682@qq.com>

------------------------------------------------------------------------

*PERCEPTION-shiny © PERCEPTIONx authors. Screens may vary from the
version you are using.*
