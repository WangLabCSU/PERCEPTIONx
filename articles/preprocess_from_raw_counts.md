# Preprocessing raw single-cell data

## 1. When you need this vignette

[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
expects a matrix that is already at clone level: one column per clone,
or one column per cell plus a cell-to-clone table. If all you have is a
raw count matrix straight out of Cell Ranger, use
[`preprocess_sc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/preprocess_sc.md)
first. It turns

    genes x cells (raw counts, no QC, no clustering)

into

    cell_id, clone_id, patient, cluster_size_class, UMAP coordinates

which is exactly the `clone_labels` argument of
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md).

You do **not** need to know anything about single-cell analysis to use
it. Every threshold is derived from your own data.

## 2. Inputs

Two things are required:

| Argument         | What it is                                              |
|------------------|---------------------------------------------------------|
| `counts`         | Genes as rows, cells as columns. Raw counts (integers). |
| `patient_labels` | One patient ID per cell.                                |

A third is optional but important:

| Argument | When to supply it |
|----|----|
| `sample_labels` | A patient has several samples (e.g. biopsies before and after treatment). |

If the counts are still in CellRanger files,
[`read_10x()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/read_10x.md)
reads a folder, an `.h5` file, or the three files directly and returns
the matrix described above:

``` r

my_counts <- read_10x("path/to/filtered_feature_bc_matrix")
# or: read_10x("path/to/filtered_feature_bc_matrix.h5")
# or: read_10x(mtx = "matrix.mtx.gz", barcodes = "barcodes.tsv.gz",
#              features = "features.tsv.gz")
```

``` r

labels <- preprocess_sc(
  counts = my_counts,
  patient_labels = cell_to_patient  # same length as ncol(my_counts)
)
```

If your matrix is transposed,
[`preprocess_sc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/preprocess_sc.md)
stops with a message telling you to use `t(counts)` instead of silently
producing nonsense.

## 3. What it does, in order

1.  **Quality control.** For each patient separately, the gene count
    gets an upper fence at `median + nmads × MAD` and the mitochondrial
    fraction an upper fence at `min(median + nmads × MAD, max_mito)`;
    the lower gene bound is the absolute `min_genes` (default 200), not
    a MAD fence. `nmads` defaults to
    5.  Tumour tissue carries more mitochondrial reads and larger cells
        than blood, so a narrow 3-MAD fence deletes real cells: on a
        real lung squamous carcinoma raw matrix, 5 MAD recovers ~88% of
        the cells CellRanger calls versus ~78% at 3 MAD, with precision
        unchanged. A hard cap (`max_mito`, default 25%) prevents the
        fence from widening on a mostly-bad sample. Pass
        `skip_qc = TRUE` when the matrix is already curated
        (e.g. malignant cells downloaded from a paper): every cell is
        then kept, which is what the original pipeline does.
2.  **Optional doublet removal.** Only when `remove_doublets = TRUE`:
    cells that `scDblFinder` calls doublets are dropped, one patient at
    a time.
3.  **Normalisation.** `LogNormalize` with a scale factor of 10,000.
4.  **Optional exclusion of non-target cells.** Only when
    `exclude_markers` is supplied; see section 8.
5.  **Variable features and PCA.** Computed on all cells together, so
    distances mean the same thing for every patient.
6.  **Clustering.** By default all cells are clustered together, as in
    the original PERCEPTION method, so one clone can appear in several
    patients. Set `cluster_scope = "patient"` to search neighbours
    inside each patient only, which makes clone labels patient-specific.
    For a single patient the two settings are identical.
7.  **Sample integration (conditional).** With
    `cluster_scope = "patient"`, a patient that has more than one sample
    is harmonised inside that patient only. A batch effect would
    otherwise split one biological clone into two. This requires the
    optional `harmony` package; without it you get a warning and the
    step is skipped. Global clustering does not integrate.
8.  **UMAP.** For visualisation only.

## 4. Reading the QC report

``` r

labels$qc
```

| Column           | Meaning                                              |
|------------------|------------------------------------------------------|
| `n_input`        | Cells before filtering                               |
| `n_kept`         | Cells that passed                                    |
| `pct_removed`    | Percentage discarded                                 |
| `main_reason`    | Which rule removed the most cells                    |
| `mito_threshold` | The mitochondrial cap actually used for this patient |
| `n_clones`       | Clusters found                                       |

Plot it to see the fences:

``` r

plot_qc_report(
  qc_metrics = labels$qc_metrics,
  kept_cells = labels$cell_labels$cell_id,
  qc_report = labels$qc
)
```

With many patients, the default pools everyone into three panels (“all
metrics”); the per-patient fences are easier to read one metric at a
time, which gives one wrapped panel per patient:

``` r

plot_qc_report(
  qc_metrics = labels$qc_metrics,
  kept_cells = labels$cell_labels$cell_id,
  qc_report = labels$qc,
  metric = "Detected genes"
)
```

If a patient loses more than ~50% of its cells, look at the violin plot
before trusting the result. The usual causes are a genuinely bad
library, or a median that was pulled down by a large low-quality
population.

## 5. From labels to predictions

``` r

prepared <- prepare_data(
  expression_matrix = my_counts,
  clone_labels = labels$cell_labels
)

clone_pred <- predict_drugs(models, prepared$clone_expression_rnorm)
patient_pred <- predict_patients(clone_pred, prepared)
```

`prepared$label_source` records where the labels came from
(`"auto_cluster"` in this case), so a result can always be traced back
to how its clones were defined.

## 6. The three ways in

| You have | Use | Labels come from |
|----|----|----|
| A raw count matrix | [`preprocess_sc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/preprocess_sc.md) → `prepare_data(clone_labels=)` | Automatic clustering |
| A matrix plus clone labels | `prepare_data(clone_labels=)` | Your labels |
| A processed Seurat object | `prepare_data(data=)` | `Idents()` or a metadata column |

Nothing is inferred silently:
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
aborts if the assay, the layer, or the patient column cannot be
resolved, and lists the available options.

## 7. What “clone” means here

In the cell-level workflow, clone labels are **operationally defined**
by transcriptomic clustering. Cells are grouped according to
transcriptional similarity, and the resulting cluster assignments are
used as surrogate clone labels for the downstream model.

With the default global clustering a clone id can occur in several
patients: it is a transcriptional state, not a patient-specific
subpopulation. Clustering inside each patient
(`cluster_scope = "patient"`) gives patient-specific labels instead; the
two are identical for a single patient.

These labels do not necessarily represent genetically validated clonal
populations, and should not be interpreted as evidence of shared genomic
ancestry. CNV or genotype evidence is required to call true clones.

[`preprocess_sc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/preprocess_sc.md)
marks every cluster with `cluster_size_class`, which reflects **only how
many cells it contains**:

- `"large"`: at least `min_clone_cells` cells (default 10)
- `"small"`: fewer cells; the cluster is kept, not dropped

The name is deliberately not `confidence`: a cluster with 50 cells is
not more likely to be a real clone than one with 8 cells, it is just
larger.

Clusters are never removed automatically:

``` r

large_only <- labels$cell_labels[labels$cell_labels$cluster_size_class == "large", ]
```

## 8. Reducing non-target cells (optional)

This step is **not part of the original PERCEPTION method** and is off
by default. The original workflow treats every cell you pass in as a
candidate clone, so if your matrix is already tumour-cell-only, skip
this section.

Tumour tissue contains immune and stromal cells that cluster on their
own. Those clusters are not resistance clones, and keeping them dilutes
the abundance of real tumour clones. There is a second, less obvious
reason: the DepMap models are trained on cancer cell lines, so an immune
or fibroblast cluster is out-of-domain input whose predicted viability
has no training support. Even so, whether to remove them is a biological
decision, not a preprocessing rule – the tumour microenvironment is part
of the sample, and in haematological malignancies the “immune”
compartment *is* the tumour.

`exclude_markers` removes cells that strongly express a programme you
want out of the analysis:

``` r

labels <- preprocess_sc(
  counts = my_counts,
  patient_labels = cell_to_patient,
  exclude_markers = list(
    immune      = c("PTPRC", "CD3D", "CD3E", "NKG7"),
    myeloid     = c("LYZ", "CD68", "CD14"),
    fibroblast  = c("COL1A1", "DCN", "LUM"),
    endothelial = c("PECAM1", "VWF")
  ),
  exclude_min_expr = 1
)
```

Each cell is scored by the mean log-normalised expression of every set;
if the highest score reaches `exclude_min_expr`, the cell is removed.
The removed cells are returned so the decision stays auditable:

``` r

table(labels$excluded_cells$set)
```

Two warnings about interpretation:

- This is an **exclusion** rule, not a tumour-cell classifier.
  Epithelial markers (`EPCAM`, `KRT8`) do **not** imply malignancy, and
  some tumours are not epithelial at all.
- A single marker is weak evidence. Use several genes per set, and
  inspect the scores in `labels$excluded_cells` before trusting the
  filter.

Leaving `exclude_markers = NULL` (the default) keeps every cell that
passed QC.

## 9. Doublet removal (optional)

Doublets are two cells captured in one droplet. They usually pass QC – a
1:1 doublet has roughly twice the RNA of a singlet but a normal
mitochondrial fraction – so the MAD fences cannot see them. A doublet
formed from two different clones is especially problematic: it sits
between them in the embedding and can become a small cluster read as a
distinct resistant clone.

``` r

labels <- preprocess_sc(
  counts = my_counts,
  patient_labels = cell_to_patient,
  remove_doublets = TRUE
)
table(labels$doublet_cells$patient)
```

Detection uses `scDblFinder`, which is an optional dependency; without
it the step warns and is skipped, and the upper gene-count fence still
removes most high-complexity doublets. Detection is probabilistic, so
the step is off by default. On the 10x PBMC 3k reference it removes
~4-5% of cells, and roughly a quarter of the cells the gene ceiling
removes are called doublets.

## 10. Limitations

- **No ambient RNA correction and no statistical empty-droplet test.**
  Barcodes are separated by the MAD fences plus an absolute gene floor,
  not by `EmptyDrops`. On a very dilute library, low-count real cells
  can therefore be lost together with the empty droplets.
- **No cross-patient integration.** Harmony, when used, runs only
  between samples of the *same* patient (`cluster_scope = "patient"`);
  cells of two different patients are never corrected against each
  other. Global clustering puts them in one embedding but applies no
  batch correction.
- **No cell-cycle regression and no cell-type annotation.** Clusters are
  transcriptional groups, not validated clones: CNV or genotype evidence
  is required before calling clonal ancestry.
- **Doublet removal is optional and probabilistic**, and off by default.
