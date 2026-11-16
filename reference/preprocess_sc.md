# Preprocess raw single-cell data into clone labels

End-to-end preprocessing for users who only have a raw count matrix:
quality control, normalisation, variable-feature selection, PCA, and
clustering, followed by a QC report. No tissue-specific thresholds are
required; the QC fences adapt to each patient's own distribution.

## Usage

``` r
preprocess_sc(
  counts,
  patient_labels,
  sample_labels = NULL,
  nmads = 5,
  min_genes = 200,
  max_mito = 25,
  max_genes = NULL,
  skip_qc = FALSE,
  remove_doublets = FALSE,
  exclude_markers = NULL,
  exclude_min_expr = 1,
  nfeatures = 2000,
  dims = 20,
  method = c("umap", "tsne"),
  resolution = 0.8,
  min_clone_cells = 10,
  integrate = c("auto", "none", "harmony"),
  cluster_scope = c("global", "patient"),
  seed = 42,
  progress_cb = NULL
)
```

## Arguments

- counts:

  Matrix (genes as rows, cells as columns) or a `Seurat` object. Raw
  counts are expected; normalised values trigger a warning because they
  would be normalised twice.

- patient_labels:

  Character vector, one entry per cell, giving the patient a cell
  belongs to. Names are matched to `colnames(counts)`.

- sample_labels:

  Character vector, one entry per cell, giving the sample or batch a
  cell belongs to. Required only when a patient has multiple samples.
  Default = `NULL`.

- nmads:

  Numeric. Width of the MAD-based QC fence, in median absolute
  deviations. Default = 5. Tumour tissue has genuinely higher
  mitochondrial content and larger cells than blood, so a 3-MAD fence
  removes real cells: against CellRanger's own calls on a real lung
  squamous carcinoma sample, 5 MAD recovers ~88% of its cells versus
  ~78% at 3 MAD, with precision unchanged.

- min_genes:

  Integer. Minimum detected genes per cell. Default = 200.

- max_mito:

  Numeric. Hard cap on the mitochondrial fraction (%). The effective
  threshold is the smaller of this and the MAD fence. Default = 25.

- max_genes:

  Integer. Fixed upper bound on detected genes per cell. When `NULL`
  (default) the bound is estimated per patient as
  `median + nmads * MAD`, which is what removes empty droplets and
  high-complexity doublets. Set it to reproduce a fixed threshold such
  as the 2,500 genes used in most tutorials.

- skip_qc:

  Logical. When `TRUE`, the MAD fences are not applied and every cell is
  kept; the QC table then only describes the input. This reproduces the
  original PERCEPTION pipeline, which expects an already-curated matrix
  (e.g. malignant cells downloaded from a paper), where an extra
  per-patient fence only perturbs the result. Default `FALSE` – raw 10x
  data always needs the fences.

- remove_doublets:

  Logical. When `TRUE`, cells classified as doublets by scDblFinder are
  removed after QC, one patient at a time. Doublets can pass every QC
  threshold – they have normal gene counts and mitochondrial fractions –
  yet they inflate the high-complexity tail and can form a spurious
  "resistant clone" made of two fused clones. Requires the optional
  scDblFinder package; without it the step is skipped with a warning.
  Default = `FALSE`, because detection is probabilistic and should be
  opted into.

- exclude_markers:

  Optional named list of gene vectors, e.g.
  `list(immune = c("PTPRC", "CD3D"), stromal = c("COL1A1"))`. Cells
  whose strongest marker-set score reaches `exclude_min_expr` are
  removed before clustering. This is an **exclusion** rule for reducing
  non-target cells, not a tumour-cell classifier: epithelial markers do
  not imply malignancy. Default = `NULL` (keep all cells).

- exclude_min_expr:

  Numeric. Mean log-normalised expression a cell must reach in a marker
  set to be excluded. Default = 1.

- nfeatures:

  Integer. Number of variable features used for PCA. Default = 2000.

- dims:

  Integer. Number of principal components. Capped automatically for
  small datasets. Default = 20.

- method:

  Character. 2D embedding for visualisation: `"umap"` (default) or
  `"tsne"`. Clustering itself always runs on the PCA embedding; this
  only changes the coordinates returned in `cell_labels`.

- resolution:

  Numeric. Clustering resolution. Default = 0.8.

- min_clone_cells:

  Integer. Clones with fewer cells are kept but marked
  `cluster_size_class = "small"`. Default = 10.

- integrate:

  Character. Sample integration strategy within a patient: `"auto"`
  (default) applies harmony when a patient has multiple samples,
  `"none"` never integrates, `"harmony"` behaves like `"auto"`. Only
  used when `cluster_scope = "patient"`.

- cluster_scope:

  Character. `"global"` (default) clusters all cells together, matching
  the original PERCEPTION method, so one clone can span several
  patients. `"patient"` clusters inside each patient only, so clone
  labels are patient-specific. Identical for a single-patient dataset.

- seed:

  Integer. Random seed. Default = 42.

- progress_cb:

  Optional function called as `progress_cb(phase, i, n)` between stages.
  Used to drive progress bars in the Shiny app.

## Value

A list with:

- cell_labels:

  Data frame with `cell_id`, `clone_id`, `clone_key` (unique
  `"patient@clone"` key), `patient`, `cluster_size_class` (`"large"` /
  `"small"`, a cell-count category and nothing more), and
  `dim_1`/`dim_2`. This is the `clone_labels` argument of
  [`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md).

- qc:

  Data frame. Per-patient QC summary, including `n_marker_removed` when
  `exclude_markers` was used and `n_doublet_removed` when
  `remove_doublets = TRUE`.

- qc_metrics:

  Data frame. Per-cell QC metrics before filtering.

- excluded_cells:

  Data frame. Cells removed by `exclude_markers` (columns `cell_id`,
  `set`, `score`); empty when no marker filtering was requested.

- doublet_cells:

  Data frame. Cells removed by `remove_doublets` (columns `cell_id`,
  `patient`, `score`); empty when the step was not requested or no
  doublets were found.

- used_harmony:

  Logical. Whether harmony was applied.

- label_source:

  Character. Always `"auto_cluster"`.

## Details

Clustering is performed on a shared PCA embedding. By default
(`cluster_scope = "global"`) every cell is clustered together, which is
what the original PERCEPTION method does: a clone is a transcriptional
state that can appear in several patients. With
`cluster_scope = "patient"` neighbours are searched inside each patient
only, so a clone can never span two patients – useful when clone labels
are meant to be patient-specific. For a single-patient dataset the two
settings are identical.

When a patient has more than one sample (e.g. longitudinal biopsies) and
`cluster_scope = "patient"`, harmony is applied inside that patient to
prevent one biological clone from being split in two by a technical
batch effect. harmony is an optional dependency: without it, the step is
skipped with a warning.

## What "clone" means here

In the cell-level workflow, clone labels are **operationally defined**
by transcriptomic clustering. Cells are grouped by transcriptional
similarity within each patient, and the resulting cluster assignments
are used as surrogate clone labels for downstream analysis. These labels
do not necessarily represent genetically validated clonal populations
and should not be interpreted as evidence of shared genomic ancestry.
CNV- or genotype-based evidence is required to call true clones.

## Examples

``` r
if (FALSE) { # \dontrun{
  labels <- preprocess_sc(counts, patient_labels = cell_patients)
  prepared <- prepare_data(
    expression_matrix = counts,
    clone_labels = labels$cell_labels
  )
} # }
```
