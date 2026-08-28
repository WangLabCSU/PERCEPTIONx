# Plot clone distribution as stacked bar

Visualizes the proportion of each clone across patients as a stacked bar
plot. Useful for understanding tumor heterogeneity and clonal
architecture.

## Usage

``` r
plot_clone_distribution(
  clone_distribution,
  response_var = NULL,
  base_size = 15,
  tooltip = TRUE,
  tooltip_col = NULL
)
```

## Arguments

- clone_distribution:

  Data frame with columns: patients, clones, weights.

- response_var:

  Character. Optional column name for response annotation. If provided,
  facets by response. Default = NULL.

- base_size:

  Numeric. Base font size. Default = 15.

- tooltip:

  Logical. If TRUE (default) and ggiraph is installed, bar segments get
  hover tooltips (clone + proportion).

- tooltip_col:

  Character. Optional existing column used as the tooltip text. Default
  = NULL (auto-builds a rich tooltip).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
  # After computing clone weights
  clone_dist <- data.frame(
    patients = c("P1", "P1", "P1", "P2", "P2", "P2"),
    clones = c("c1", "c2", "c3", "c1", "c2", "c3"),
    weights = c(0.3, 0.5, 0.2, 0.6, 0.3, 0.1)
  )
  plot_clone_distribution(clone_dist)
} # }
```
