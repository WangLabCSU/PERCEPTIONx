# PERCEPTIONx Patient Data Annotation Functions

Functions for annotating single-cell data with clone and patient
information, and preparing patient data for the prediction pipeline.

## Usage

``` r
run_seurat_pipeline(
  expression_matrix,
  method,
  min_cells,
  min_features,
  nfeatures,
  dims,
  resolution,
  seed
)
```
