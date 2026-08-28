# Compare performance of multiple trained models

Computes and summarizes performance metrics across multiple drug models.
Returns performance in cross-validation, bulk test, pseudo-bulk, and
scRNA datasets.

## Usage

``` r
compare_performance(model_list, threshold = 0.3, verbose = TRUE)
```

## Arguments

- model_list:

  A named list of trained model objects from train_models().

- threshold:

  Numeric. Minimum correlation threshold for "passing" models. Default =
  0.3.

- verbose:

  Logical. Whether to print summary statistics. Default = TRUE.

## Value

A list containing:

- perf_cv:

  Cross-validation performance (\|correlation\| = sqrt(R²); p-value not
  available) for each drug

- perf_bulk:

  Bulk test set performance for each drug

- perf_pseudo_bulk:

  Pseudo-bulk performance for each drug

- perf_scRNA:

  Single-cell RNA performance for each drug

- summary:

  Summary statistics of models passing threshold

## Examples

``` r
if (FALSE) { # \dontrun{
  models <- train_models(drug_list = c("abemaciclib", "erlotinib"))
  perf <- compare_performance(models)
} # }
```
