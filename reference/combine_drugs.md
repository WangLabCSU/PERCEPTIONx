# Combine multi-drug predictions into a single combination score

Implements the IDA (Independent Drug Action) combination principle used
in PERCEPTION: each clone-level drug viability is first z-score scaled
across all clones, then the per-clone combined score is derived from the
drug-specific z-scores. The default rule `min` means "the most effective
drug in the combination drives the clone-level combined score" — the
exact rule hard-coded in the Shiny app's internal `combo_clone_frame()`.

## Usage

``` r
combine_drugs(
  pred,
  patient = NULL,
  clone_id = NULL,
  combine_rule = c("min", "mean", "max", "median"),
  na.rm = TRUE,
  include_single_z = FALSE
)
```

## Arguments

- pred:

  Matrix or data.frame of clone-level drug viability predictions.
  Typically the output of
  [`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md).
  Rows = clones, columns = drugs. Row names must be canonical PERCEPTION
  clone keys (`"patient@clone_id"`). Alternatively, pass a `data.frame`
  with an explicit `patient` column and a `clone_id` column.

- patient:

  Optional character vector of patient IDs, length equal to
  `nrow(pred)`. Required only when `pred` row names are NOT
  `"patient@clone_id"` keys.

- clone_id:

  Optional character vector of clone IDs, length equal to `nrow(pred)`.
  Required only when passing `patient` explicitly.

- combine_rule:

  Character scalar, per-clone aggregation after z-scoring. Must be one
  of `"min"` (IDA default), `"mean"`, `"max"`, or `"median"`.

- na.rm:

  Logical, passed through to the aggregator. Default `TRUE`.

- include_single_z:

  Logical. If `TRUE`, per-drug z-scored columns are retained in the
  output alongside `comb_viability`. Default `FALSE`.

## Value

A data.frame with columns:

- patient:

  Patient ID.

- clone_id:

  Clone ID.

- comb_viability:

  Combined viability score (already z-scored across clones, so
  [`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md)
  will not standardise it again).

- \<drug_z\>:

  (optional, when `include_single_z = TRUE`) one additional column per
  input drug, named like `<drug>.z`, containing the per-drug z-scores
  before combination.

## Details

This function is the missing glue between
[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md)
and downstream patient-level aggregation / visualisation. Its output
feeds directly into
[`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md)
(as the app's
[`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md)
already skips re-standardising the `comb_viability` column) and
[`plot_clone_viability()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_clone_viability.md).

## Pipeline


      pred     <- predict_drugs(models, rank_expr)        # clones x drugs matrix
      clone_df <- combine_drugs(pred)                    # patient / clone_id / comb_viability
      pat      <- predict_patients(clone_df, prepared_data)
      plot_clone_viability(clone_df, viability_var = "comb_viability")

The patient-level aggregated score will be identical to the one produced
by the Shiny app for the same drugs / clones.

## Caveats

The default `min` rule amplifies noise when individual drug models are
very poor — the "least bad" model can still be quite bad, but `min`
picks it as the best combination driver. For exploratory analyses you
can pass `combine_rule = "mean"` or `combine_rule = "max"` and compare.

When only one drug is present in `pred`, the combined score is simply
that drug's z-score.

## See also

[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md),
[`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md),
[`plot_clone_viability()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_clone_viability.md),
[`zscore_viability()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/zscore_viability.md)
