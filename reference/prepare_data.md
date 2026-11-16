# Prepare patient data for PERCEPTIONx prediction

End-to-end preprocessing pipeline that takes raw single-cell expression
data and produces a rank-normalized subclone expression matrix and clone
counts table, ready for direct use with
[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md)
and
[`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md).

## Usage

``` r
prepare_data(
  method = c("umap", "tsne"),
  expression_matrix = NULL,
  patient_mapping = NULL,
  cell_col = "cell_id",
  patient_col = "patient_id",
  parse_patient = FALSE,
  patient_sep = "_",
  patient_pos = 1L,
  genes_to_use = NULL,
  seurat_resolution = 0.8,
  seurat_dims = 10L,
  seurat_nfeatures = 2000L,
  seurat_min_cells = 3L,
  seurat_min_features = 200L,
  seurat_seed = 42L,
  skip_clustering = FALSE,
  progress_cb = NULL,
  data = NULL,
  clone_labels = NULL,
  assay = NULL,
  layer = NULL,
  clone_col = NULL
)
```

## Arguments

- method:

  Character. Dimensionality reduction method passed to
  [`annotate_clones()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/annotate_clones.md).
  One of `"umap"` (default) or `"tsne"`. UMAP is faster and preserves
  global structure; t-SNE emphasizes local neighborhoods.

- expression_matrix:

  Matrix. Gene expression matrix with genes as rows and cells as
  columns.

- patient_mapping:

  List, data frame, or character vector. Patient-cell mapping in one of
  three formats:

  List format

  :   Named list where each element is a patient ID and contains a
      character vector of cell IDs. Example:
      `list(Patient_1 = c("Cell_1", "Cell_2"), Patient_2 = c("Cell_3"))`

  Data frame format

  :   Metadata with cell ID and patient ID columns. Specify column names
      via `cell_col` and `patient_col`. Example:
      `data.frame(cell_id = c("Cell_1", "Cell_2"), patient_id = c("P1", "P1"))`

  Vector format

  :   One patient ID per cell, in the column order of the expression
      matrix; a named vector is matched by cell ID instead. Use this
      when the IDs already exist as a variable, e.g. `my_patients`.

  If NULL and `clone_labels` carries a `patient` (or `patient_id`)
  column, the grouping is taken from those labels. Otherwise all cells
  are assigned to a single patient "patient1".

- cell_col:

  Character. Cell ID column name in patient_mapping data frame. Default
  = "cell_id". Only used when patient_mapping is a data frame.

- patient_col:

  Character. Patient ID column name in patient_mapping data frame.
  Default = "patient_id". Only used when patient_mapping is a data
  frame.

- parse_patient:

  Logical. If TRUE, parse patient ID from cell_col using separator.
  Default = FALSE. Auto-enabled if patient_sep or patient_pos is
  provided. Useful when cell names contain patient info (e.g.,
  "P11_M_Barcode").

- patient_sep:

  Character. Separator to split cell_col for parsing patient ID. Default
  = "\_". Providing this parameter auto-enables parse_patient.

- patient_pos:

  Integer. Position of patient ID after splitting. Default = 1 (first
  element). Providing this parameter auto-enables parse_patient.

- genes_to_use:

  Character vector. Genes to retain in the output matrix. If NULL, all
  genes in the expression matrix are used.

- seurat_resolution:

  Numeric. Clustering resolution. Default = 0.8.

- seurat_dims:

  Integer. PCA dimensions for clustering. Default = 10.

- seurat_nfeatures:

  Integer. Variable features count. Default = 2000.

- seurat_min_cells:

  Integer. Minimum cells per feature. Default = 3.

- seurat_min_features:

  Integer. Minimum features per cell. Default = 200. Auto-adjusted to
  10% of gene count if the expression matrix has fewer genes.

- seurat_seed:

  Integer. Random seed. Default = 42.

- skip_clustering:

  Logical. If TRUE, skip the Seurat clustering step and treat every
  column of `expression_matrix` as one pre-defined clone. Use this when
  you already have a clone-level expression matrix (e.g. from a
  published study). Rank normalization and clone counts are still
  applied. No UMAP/t-SNE embedding is produced in this mode.

- progress_cb:

  Optional function called as `progress_cb(phase, i, n, drug)` with the
  current pipeline stage (e.g. `"seurat-pca"`, `"mapping"`). Used to
  show live progress in the Shiny app.

- data:

  Optional `Seurat` object. When supplied, the expression matrix is read
  from `assay`/`layer`, clone labels from `clone_col` (a metadata
  column, a label vector, or `Idents()`), and patient IDs from
  `patient_col` or, when the object has no such column, from
  `patient_mapping`. Use this instead of `expression_matrix` plus
  `patient_mapping`.

- clone_labels:

  Optional data frame with `cell_id` and `clone_id` columns mapping
  every cell to a clone. Skips clustering. This is the `cell_labels`
  element returned by
  [`preprocess_sc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/preprocess_sc.md),
  and its `patient` column makes `patient_mapping` optional. Cells of
  `expression_matrix` without a label (for example cells removed by QC)
  are dropped, so the unfiltered matrix can be passed directly.

- assay:

  Character. Assay to read when `data` is a Seurat object. Default =
  `Seurat::DefaultAssay(data)`.

- layer:

  Character. Layer (Seurat v5) or slot (Seurat v4) to read. Default =
  `NULL`: the `counts` layer when the object has one, otherwise `data`.
  Set it explicitly to override, e.g. `layer = "data"` to use normalised
  values.

- clone_col:

  Character. Metadata column in `data` holding clone labels, or a vector
  of labels with one entry per cell (in matrix column order). When
  `NULL`, `Seurat::Idents(data)` is used.

## Value

A named list with:

- clone_expression_rnorm:

  Matrix. Rank-normalized clone-level expression (genes as rows,
  patient_clone as columns). Ready for
  [`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md).

- clone_counts:

  Data frame. Clone abundance per patient. Ready for
  [`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md).

- cell_clone_map:

  Data frame. Cell-to-clone mapping with columns cell_id, clone_id,
  patient, dim_1, dim_2.

- clone_viability_df_template:

  Data frame. Template with patient and clone_id columns, ready to merge
  with
  [`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md)
  output.

- umap_coords:

  Data frame. 2D embedding coordinates per cell (cell_id, dim_1, dim_2).
  Ready for
  [`plot_tsne_response()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_tsne_response.md).

- reduction_method:

  Character. The method used (`"umap"` or `"tsne"`).

## Details

The pipeline performs:

1.  Seurat clustering to define transcriptional subclones

2.  Cell-to-patient and cell-to-clone annotation

3.  Clone-level mean expression computation

4.  Rank normalization of clone expression

5.  Clone abundance table construction

## Examples

``` r
if (FALSE) { # \dontrun{
  # List format (same as Rmd)
  prepared <- prepare_data(
    expression_matrix = patient_scRNA,
    patient_mapping = cell_names_list,
    genes_to_use = GOI
  )

  # Or data frame format (from metadata)
  metadata <- data.frame(cell_id = colnames(patient_scRNA), patient_id = patient_ids)
  prepared <- prepare_data(patient_scRNA, metadata)

  # Parse patient ID from cell names (e.g., "P11_M_Barcode" -> "P11")
  metadata <- data.frame(Cell = c("P11_M_Barcode1", "P12_M_Barcode2"))
  prepared <- prepare_data(
    patient_scRNA, metadata,
    cell_col = "Cell",           # Custom column name
    parse_patient = TRUE,        # Parse from Cell column
    patient_sep = "_",           # Split by "_"
    patient_pos = 1              # Take first element
  )
  # Result: Patient IDs = "P11", "P12"

  # Use directly with prediction functions
  clone_pred <- predict_drugs(models, prepared$clone_expression_rnorm)
  patient_pred <- predict_patients(clone_pred, prepared)

  # After preprocess_sc(): the labels alone are enough, because they already
  # carry the patient of every cell. The unfiltered matrix can be passed as
  # is; cells without a label are dropped.
  labels <- preprocess_sc(counts, cell_patients)
  prepared <- prepare_data(
    expression_matrix = counts,
    clone_labels = labels$cell_labels
  )

  # Seurat object: patient IDs from metadata, or from a vector
  prepared <- prepare_data(data = so, patient_col = "patient_id")
  prepared <- prepare_data(data = so, patient_mapping = my_patients)

  # Clone labels as a vector, one per cell (same length as ncol(so))
  prepared <- prepare_data(data = so, clone_col = my_clones)
} # }
```
