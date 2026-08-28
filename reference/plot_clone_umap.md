# Plot UMAP colored by clone identity

Single cells in the 2D embedding colored by clone identity — the spatial
analogue of the paper's Extended Data Fig. 8a, showing at a glance which
transcriptional subclones sit where.

## Usage

``` r
plot_clone_umap(
  tsne_data,
  clone_col = "clone_id",
  title = NULL,
  color_label = "Clone",
  point_size = NULL,
  base_size = 11,
  tooltip = TRUE,
  tooltip_col = NULL
)
```

## Arguments

- tsne_data:

  Data frame with columns: X, Y, and the clone column.

- clone_col:

  Character. Name of the column holding clone ids. Default = "clone_id".

- title:

  Character. Plot title. Default = NULL.

- color_label:

  Character. Legend label. Default = "Clone".

- point_size:

  Numeric. Point size. If NULL (default), auto-adapts to the number of
  cells to avoid overplotting.

- base_size:

  Numeric. Base font size. Default = 11.

- tooltip:

  Logical. If TRUE (default) and ggiraph is installed, points get hover
  tooltips (clone id).

- tooltip_col:

  Character. Optional existing column used as the tooltip text. Default
  = NULL (auto-builds from the clone id).

## Value

A ggplot object.
