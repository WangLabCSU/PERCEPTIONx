# Plot QC metrics before and after filtering

Draws one violin per metric, split by whether the cell was kept, with
the per-patient thresholds drawn as horizontal lines. This is what makes
the automatic QC auditable: a user does not need to understand the MAD
rule to see whether their data were cut too aggressively.

## Usage

``` r
plot_qc_report(
  qc_metrics,
  kept_cells,
  qc_report,
  base_size = 12,
  metric = c("all", "Detected genes", "Total counts", "Mitochondrial %"),
  ncol = 6L
)
```

## Arguments

- qc_metrics:

  Data frame returned as `qc_metrics` by
  [`preprocess_sc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/preprocess_sc.md)
  (columns `cell_id`, `n_genes`, `n_counts`, `pct_mito`, `patient`).

- kept_cells:

  Character vector of cell IDs that passed QC.

- qc_report:

  Data frame returned as `qc` by
  [`preprocess_sc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/preprocess_sc.md).

- base_size:

  Numeric. Base font size. Default = 12.

- metric:

  Character. `"all"` (default) pools every patient and stacks the three
  metrics in one column – the only readable layout when a cohort has
  dozens of patients. Pass one of `"Detected genes"`, `"Total counts"`
  or `"Mitochondrial %"` to get one panel per patient (wrapped into
  rows), each with that patient's own threshold line; this is the layout
  to use when per-patient QC needs to be inspected.

- ncol:

  Integer. Number of facet columns in the per-patient layout. Default =
  6.

## Value

A `ggplot` object, faceted by metric.
