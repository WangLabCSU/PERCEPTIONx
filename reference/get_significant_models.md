# Get best performing models

Filters model list to return only models that meet performance criteria.
Useful for selecting significant drug models for downstream analysis.

## Usage

``` r
get_significant_models(
  model_list,
  min_correlation = 0.3,
  max_pvalue = 0.05,
  dataset = "scRNA"
)
```

## Arguments

- model_list:

  A named list of trained model objects.

- min_correlation:

  Numeric. Minimum correlation threshold. Default = 0.3.

- max_pvalue:

  Numeric. Maximum p-value threshold. Default = 0.05.

- dataset:

  Character. Which dataset to use for filtering: "scRNA" (default),
  "bulk", "pseudo_bulk", or "cv".

## Value

A filtered list containing only models meeting the criteria.

## Examples

``` r
if (FALSE) { # \dontrun{
  models <- train_models()
  significant_models <- get_significant_models(models, min_correlation = 0.3, max_pvalue = 0.05)
} # }
```
