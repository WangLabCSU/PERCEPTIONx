# Export plot to file via Cairo device

High-resolution plot export using the Cairo rendering engine for
superior anti-aliasing and cross-platform font rendering.

## Usage

``` r
export_plot_cairo(
  file,
  plot,
  format = "png",
  width = 7,
  height = 5,
  res = 600,
  draw_fun = NULL
)
```

## Arguments

- file:

  Character. Output file path.

- plot:

  A ggplot or grid object to export.

- format:

  Character. One of `"png"`, `"svg"`, or `"pdf"`. Default = `"png"`.

- width:

  Numeric. Plot width in inches. Default = 7.

- height:

  Numeric. Plot height in inches. Default = 5.

- res:

  Numeric. Output resolution (DPI) for PNG. Default = 600, minimum 96.

- draw_fun:

  Function. Optional custom draw function (e.g.
  [`gridExtra::grid.arrange`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html)
  for multi-panel layouts). If `NULL`, `print(plot)` is used.

## Value

Invisibly returns the file path. Called for its side effect of creating
the file.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
export_plot_cairo("plot.png", p)
export_plot_cairo("plot.pdf", p, format = "pdf")
} # }
```
