# Plot ROC curve with AUC annotation

Generates a ROC curve from predicted vs observed response with AUC value
annotation.

## Usage

``` r
plot_roc_curve(
  response,
  predictor,
  smooth_curve = TRUE,
  base_size = 15,
  auc_digits = 3,
  title = NULL,
  tooltip = TRUE
)
```

## Arguments

- response:

  Factor or numeric. True response labels (e.g., "R"/"NR" or 0/1).

- predictor:

  Numeric. Predicted values (e.g., viability scores).

- smooth_curve:

  Logical. Whether to smooth the ROC curve. Default = TRUE.

- base_size:

  Numeric. Base font size. Default = 15.

- auc_digits:

  Integer. Number of digits for AUC display. Default = 3.

- title:

  Character. Plot title. Default = NULL.

- tooltip:

  Logical. If TRUE (default) and ggiraph is installed, curve points get
  hover tooltips (FPR / TPR).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
  response <- factor(c("R", "NR", "R", "NR", "R"))
  predictor <- c(0.8, 0.2, 0.7, 0.3, 0.9)
  plot_roc_curve(response, predictor)
} # }
```
