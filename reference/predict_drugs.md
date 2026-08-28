# Predict drug response for cells or clones

Given a trained model (or list of models) and a rank-normalized
expression matrix, predicts viability scores for each cell/sample across
one or more drugs. This function merges the former viability_from_model
(single drug) and viability_in_each_dataset (multi-drug) into a unified
interface.

## Usage

``` r
predict_drugs(model_list, expr)
```

## Arguments

- model_list:

  A named list of model objects (each with a `$model` element), or a
  single model object. From
  [`train_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/train_models.md)
  or
  [`load_model()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/load_model.md).

- expr:

  Matrix or data frame. Rank-normalized expression matrix with genes as
  rows and cells/samples as columns.

## Value

A matrix with cells/samples as rows and drugs as columns, containing
predicted viability scores. Lower values indicate higher drug
sensitivity.

## Details

If a small fraction (\<= 50\\ `expr` (e.g. genes filtered out during
scRNA QC), they are imputed with the neutral rank value 0.5 and a
warning is issued. If more than half of the features are missing,
prediction stops with an error.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Single drug
  models <- load_model("erlotinib", read = TRUE)
  pred <- predict_drugs(models, expr_rnorm)

  # Multiple drugs
  models <- train_models(drug_list = c("abemaciclib", "erlotinib"),
                         cancer_type = "PanCan", exclude_cancer = "PanCan", GOI = GOI)
  pred <- predict_drugs(models, expr_rnorm)
} # }
```
