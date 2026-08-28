# Complete patient response visualization pipeline

Generates a comprehensive visualization panel for patient drug response
prediction, including clone distribution, clone-level viability,
response boxplot, and ROC curve. This is a convenience function that
combines multiple plot functions.

## Usage

``` r
plot_patient_response_panel(
  clone_distribution,
  clone_viability,
  exp_vs_pred,
  response_col = "response",
  viability_col = "comb_viability",
  predicted_col = "predicted_viability",
  weights_col = "weights",
  layout_matrix = NULL
)
```

## Arguments

- clone_distribution:

  Data frame. Clone weights per patient.

- clone_viability:

  Data frame. Viability scores per clone.

- exp_vs_pred:

  Data frame. Predicted vs observed response.

- response_col:

  Character. Response column name. Default = "response".

- viability_col:

  Character. Viability column name. Default = "comb_viability".

- predicted_col:

  Character. Predicted values column name. Default =
  "predicted_viability".

- weights_col:

  Character. Weights column name. Default = "weights".

- layout_matrix:

  Matrix. Layout for grid.arrange. Default = NULL (auto).

## Value

A gtable object from grid.arrange.

## Examples

``` r
if (FALSE) { # \dontrun{
  # After running prediction pipeline
  panel <- plot_patient_response_panel(
    clone_distribution = clone_dist_df,
    clone_viability = clone_kill_df,
    exp_vs_pred = response_df
  )
  ggsave(panel, filename = "patient_response.pdf", height = 15, width = 10)
} # }
```
