# Compute clone-level mean expression from single-cell data

For each patient, groups cells by clone_id and computes the mean
expression of each gene within each clone. This reduces single-cell
resolution to subclone resolution, matching the original PERCEPTION
pipeline.

## Usage

``` r
clone_mean_expression(expression_matrix, cell_clone_map, patient_ids = NULL)
```

## Arguments

- expression_matrix:

  Matrix or data frame. Gene expression matrix with genes as rows and
  cells as columns.

- cell_clone_map:

  Data frame with columns `cell_id` and `clone_id`, mapping each cell to
  its clone. Typically from Seurat clustering.

- patient_ids:

  Character vector. Patient ID for each cell, same length as columns of
  expression_matrix. If NULL, all cells are assumed from one patient.

## Value

A named list of matrices, one per patient. Each matrix has genes as rows
and clone IDs as columns, with mean expression values.
