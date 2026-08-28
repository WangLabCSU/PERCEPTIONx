# Plot model performance across datasets

Visualizes the number of drugs achieving different correlation
thresholds across bulk, pseudo-bulk, and single-cell datasets.

## Usage

``` r
plot_model_performance(
  performance_list,
  threshold_range = NULL,
  base_size = 20,
  highlight_threshold = 0.3,
  tooltip = TRUE
)
```

## Arguments

- performance_list:

  Named list of model performance objects (from
  train_perception_models).

- threshold_range:

  Numeric vector. Correlation thresholds to evaluate. Default = NULL
  (auto: from 0.1 up to at least 0.6, extended to the highest observed
  correlation so the curve always reaches 0).

- base_size:

  Numeric. Base font size. Default = 20.

- highlight_threshold:

  Numeric. Threshold to highlight with vertical line. Default = 0.3.

- tooltip:

  Logical. If TRUE (default) and ggiraph is installed, points get hover
  tooltips (dataset, threshold, drug count).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
  models <- train_perception_models(c("abemaciclib", "erlotinib"), ...)
  plot_model_performance(models)
} # }
```
