# Compute pseudo-bulk expression for a patient

Calculates weighted average of clone-level expression to produce a
pseudo-bulk expression profile for a single patient. This is used for
evaluating model performance on pseudo-bulk data, not for prediction.

## Usage

``` r
each_patient_pseudo_bulk(
  x = 1,
  comb_viability_df,
  Clone_Counts_per_patients,
  clone_Level_z_expression_df
)
```

## Arguments

- x:

  Integer. Index of the patient. Default = 1.

- comb_viability_df:

  Data frame. Retained for backward compatibility; no longer used (clone
  ids are read from the expression column names).

- Clone_Counts_per_patients:

  Data frame. Patient-clone abundance data.

- clone_Level_z_expression_df:

  Matrix. Clone-level expression data with columns named by
  patient-clone identifiers.

## Value

Named numeric vector of pseudo-bulk expression values.
