# Plot validation ROC curves for trained models

For each model and each validation dataset (bulk, pseudo-bulk,
single-cell), the observed response is stratified into the top vs bottom
50\\ (resistant vs sensitive — the PERCEPTION paper's convention,
Extended Data Fig. 4C) and a ROC curve of the predicted viability is
drawn, one curve per dataset, with the AUC annotated in the legend.
Higher AUC = the model stratifies better. This is the most informative
single-model summary after training.

## Usage

``` r
plot_model_roc(performance_list, base_size = 13)
```

## Arguments

- performance_list:

  Named list of model objects from
  [`train_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/train_models.md),
  each with a `$predVSgroundTruth` element.

- base_size:

  Numeric. Base font size. Default = 13.

## Value

A ggplot object: ROC curves per validation dataset, faceted by drug when
more than one model is provided.

## Examples

``` r
if (FALSE) { # \dontrun{
  models <- train_models(drug_list = "abemaciclib", ...)
  plot_model_roc(models)
} # }
```
