# Z-score scale viability values across patients

Applies base R [`scale()`](https://rdrr.io/r/base/scale.html) to each
drug column across all patients, centering to mean 0 and standard
deviation 1. This ensures comparability of predicted viability scores
between different drugs before patient-level aggregation. Matches the
original PERCEPTION pipeline.

## Usage

``` r
zscore_viability(clone_viability_df, cols = NULL)
```

## Arguments

- clone_viability_df:

  Data frame. Must have columns 'patient' and 'clone_id', plus one or
  more drug columns with predicted viability values.

- cols:

  Character vector of column names to scale. Default `NULL` scales every
  drug column. Use this to leave already-standardized columns (e.g. a
  pre-combined `comb_viability`) untouched.

## Value

A data frame with the same structure, but drug columns z-score scaled.
