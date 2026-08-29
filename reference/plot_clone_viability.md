# Plot clone-level viability (lollipop plot)

Visualizes predicted drug sensitivity for each clone within patients.
Each clone is represented as a point with a stem (lollipop style).
Useful for identifying resistant clones within heterogeneous tumors.

## Usage

``` r
plot_clone_viability(
  clone_viability,
  viability_var = "comb_viability",
  weights_var = NULL,
  response_var = NULL,
  drug = NULL,
  base_size = 11,
  y_limits = c(-3, 1.2),
  viridis_scale = FALSE,
  tooltip = TRUE,
  tooltip_col = NULL
)
```

## Arguments

- clone_viability:

  Data frame with columns: patient, clone_id, and the viability column
  named by `viability_var`. Optionally also `weights_var` (clone
  proportion) and `response_var` (clinical response) columns.

- viability_var:

  Character. Column name for viability values. Default =
  "comb_viability".

- weights_var:

  Character. Optional column name for clone weights (point size).
  Default = NULL.

- response_var:

  Character. Optional column for response annotation. Default = NULL.

- drug:

  Character. Drug name, used as plot title. Default = NULL.

- base_size:

  Numeric. Base font size. Default = 11.

- y_limits:

  Numeric vector. Y-axis limits. Default = c(-3, 1.2).

- viridis_scale:

  Logical. If TRUE, uses a viridis sequential scale (dark =
  low/sensitive, yellow = high/resistant). If FALSE (default), uses the
  diverging red-blue scale (blue = low/sensitive, white = neutral, red =
  high/resistant) with data-driven limits, so every value keeps a real
  color – nothing clips to grey/NA outside a fixed window.

- tooltip:

  Logical. If TRUE (default) and ggiraph is installed, points get hover
  tooltips (clone + viability + proportion).

- tooltip_col:

  Character. Optional existing column used as the tooltip text. Default
  = NULL (auto-builds a rich tooltip).

## Value

A ggplot object.

## Details

Facet strategy is adaptive:

- \<= 12 patients: single-row grid, strips on the bottom (45 deg).

- \> 12 patients: `facet_wrap` grid, one compact panel per patient.

## Examples

``` r
if (FALSE) { # \dontrun{
  clone_kill <- data.frame(
    patient = c("P1", "P1", "P2", "P2"),
    clone_id = c("c1", "c2", "c1", "c2"),
    comb_viability = c(-0.5, 0.8, -1.2, 0.3)
  )
  plot_clone_viability(clone_kill, viability_var = "comb_viability")
} # }
```
