# Plot predicted vs observed response boxplot

Creates a combined violin + box + jitter plot comparing predicted
viability between response groups (e.g., Responders vs Non-Responders)
with a statistical significance bracket.

## Usage

``` r
plot_response_boxplot(
  exp_vs_pred,
  response_var = "response",
  predicted_var = "predicted_viability",
  y_label = "Predicted Viability (z-score)",
  base_size = 15,
  compare_method = "wilcox.test",
  alternative = "less",
  tooltip = TRUE
)
```

## Arguments

- exp_vs_pred:

  Data frame with columns: response, predicted_viability.

- response_var:

  Character. Column name for response labels. Default = "response".

- predicted_var:

  Character. Column name for predicted values. Default =
  "predicted_viability".

- y_label:

  Character. Y-axis label. Default = "Predicted Viability (z-score)".

- base_size:

  Numeric. Base font size. Default = 15.

- compare_method:

  Character. Statistical test method. One of `"wilcox.test"` (default)
  or `"t.test"`.

- alternative:

  Character. Alternative hypothesis direction. Default = "less".

- tooltip:

  Logical. If TRUE (default) and ggiraph is installed, jittered points
  get hover tooltips.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
  exp_pred <- data.frame(
    response = factor(c("R", "NR", "R", "NR")),
    predicted_viability = c(0.8, 0.2, 0.7, 0.3)
  )
  plot_response_boxplot(exp_pred)
} # }
```
