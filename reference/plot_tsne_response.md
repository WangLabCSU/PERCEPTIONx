# Plot UMAP with drug response overlay

Visualizes single cells in UMAP space with color overlay representing
either biomarker expression or predicted drug sensitivity.

## Usage

``` r
plot_tsne_response(
  tsne_data,
  color_var = "viability_scaled",
  title = NULL,
  color_label = "Predicted Viability",
  point_size = NULL,
  colors = NULL,
  palette = c("viridis", "diverging"),
  midpoint = 0,
  limits = NULL,
  base_size = 11,
  tooltip = TRUE,
  tooltip_col = NULL
)
```

## Arguments

- tsne_data:

  Data frame with columns: X, Y (coordinates), and optional
  biomarker/viability columns.

- color_var:

  Character. Name of the column to use for color mapping. Default =
  "viability_scaled".

- title:

  Character. Plot title. Default = NULL.

- color_label:

  Character. Legend label for color. Default = "Predicted Viability".

- point_size:

  Numeric. Point size. If NULL (default), auto-adapts to the number of
  cells to avoid overplotting.

- colors:

  Character vector. Custom gradient colors (low, mid, high) for the
  sequential palette. Default = NULL (uses built-in viridis).

- palette:

  Character. One of `"viridis"` (sequential, default) or `"diverging"`
  (blue-white-red centered at `midpoint`).

- midpoint:

  Numeric. Center value for diverging palette. Default = 0.

- limits:

  Numeric vector of length 2. Optional fixed scale limits (e.g.
  `c(0, 1)` for 0-1 expression) to pin the color at each end. Default =
  NULL (data-driven limits — recommended, so extreme values always keep
  a real color instead of clipping to grey/NA).

- base_size:

  Numeric. Base font size for theme. Default = 11.

- tooltip:

  Logical. If TRUE (default) and ggiraph is installed, points get hover
  tooltips (cell id if present, else the colored value).

- tooltip_col:

  Character. Optional existing column used as the tooltip text (e.g.
  "cell_id"). Default = NULL (auto-builds from color_var).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
  # After predicting viability for single cells
  tsne_data <- data.frame(
    X = lung_tSNE$X,
    Y = lung_tSNE$Y,
    viability_scaled = range01(rank(-viability_pred))
  )
  plot_tsne_response(tsne_data, color_var = "viability_scaled")
} # }
```
