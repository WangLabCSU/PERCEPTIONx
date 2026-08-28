# Build clone abundance table from cell-clone mapping

Computes the number of cells per clone per patient, producing the
`clone_counts` data frame required by
[`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md).

## Usage

``` r
build_clone_counts(cell_clone_map, patient_ids)
```

## Arguments

- cell_clone_map:

  Data frame with columns `cell_id` and `clone_id`.

- patient_ids:

  Character vector. Patient ID for each cell, in the same order as rows
  in cell_clone_map.

## Value

A data frame with first column `patients` and remaining columns as clone
IDs with cell counts as values.
