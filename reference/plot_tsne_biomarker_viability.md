# Plot UMAP side-by-side for biomarker and viability

Creates a side-by-side comparison of biomarker expression and predicted
viability in UMAP space. Useful for visualizing correlation between
marker and response.

## Usage

``` r
plot_tsne_biomarker_viability(
  tsne_data,
  biomarker_var = "biomarker_scaled",
  viability_var = "viability_scaled",
  biomarker_label = "Biomarker Exp",
  viability_label = "Drug Viability",
  nrow = 1,
  base_size = 8
)
```

## Arguments

- tsne_data:

  Data frame with X, Y coordinates and both biomarker/viability columns.

- biomarker_var:

  Character. Column name for biomarker expression. Default =
  "biomarker_scaled".

- viability_var:

  Character. Column name for viability values. Default =
  "viability_scaled".

- biomarker_label:

  Character. Legend label for biomarker. Default = "Biomarker Exp".

- viability_label:

  Character. Legend label for viability. Default = "Drug Viability".

- nrow:

  Integer. Number of rows in arrangement. Default = 1.

- base_size:

  Numeric. Base font size. Default = 8.

## Value

A gtable object from grid.arrange.

## Examples

``` r
if (FALSE) { # \dontrun{
  tsne_data <- data.frame(
    X = lung_tSNE$X,
    Y = lung_tSNE$Y,
    biomarker_scaled = range01(rank(MDM2_expression)),
    viability_scaled = range01(rank(-viability_pred))
  )
  plot_tsne_biomarker_viability(tsne_data)
} # }
```
