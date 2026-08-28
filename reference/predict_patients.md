# Predict drug response at patient level

Aggregates clone-level drug response predictions to patient-level using
various weighting strategies based on clone abundance. This function
merges the former each_patient_viability (single drug) and
each_patient_viabilityv2 (multi-drug) into a unified interface.

## Usage

``` r
predict_patients(
  clone_pred,
  prepared_data,
  clone_counts = NULL,
  mode = "weighted_max",
  zscore = TRUE
)
```

## Arguments

- clone_pred:

  Matrix from
  [`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md),
  with clones as rows and drugs as columns. Row names should match clone
  column names from
  [`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md).

- prepared_data:

  List from
  [`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md),
  containing `$clone_viability_template` and `$clone_counts`.
  Alternatively, you can pass `clone_counts` directly as a data frame
  (legacy mode).

- clone_counts:

  Optional. Data frame with patients as rows and clone IDs as columns.
  Only needed if `prepared_data` is not a list from
  [`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md).

- mode:

  Character. Aggregation method:

  "weighted_max"

  :   Maximum of weighted viability across clones. Default.

  "max"

  :   Maximum viability across clones (most resistant clone)

  "weighted_average"

  :   Weighted average of clone viability by clone abundance

  "min"

  :   Minimum viability across clones (most sensitive clone)

  "average"

  :   Average of clone viability

- zscore:

  Logical. Whether to z-score scale drug columns across patients before
  aggregation. Default = TRUE. Matches the original PERCEPTION pipeline.

## Value

A data frame with patients as rows and drugs as columns, containing
aggregated viability scores.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Simple workflow: pass prepare_data() output directly
  prepared <- prepare_data(patient_scRNA)
  clone_pred <- predict_drugs(models, prepared$clone_expression)
  patient_pred <- predict_patients(clone_pred, prepared)

  # Legacy workflow: manually build clone_viability_matrix
  clone_viability_df <- data.frame(
    patient = patient_ids,
    clone_id = clone_ids,
    clone_pred
  )
  patient_pred <- predict_patients(clone_viability_df, clone_counts)
} # }
```
