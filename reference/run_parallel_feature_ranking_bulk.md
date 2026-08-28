# Parallel feature ranking for multiple drugs

Runs feature_ranking_bulk for a list of drugs.

## Usage

``` r
run_parallel_feature_ranking_bulk(
  infunc_DrugsToUse,
  id_cancerType = "PanCan",
  infunc_exclude_cancer = "PanCan",
  infunc_GOI,
  ncores = 4
)
```

## Arguments

- infunc_DrugsToUse:

  Character vector. Drug names to rank features for.

- id_cancerType:

  Character. Cancer type. Default = "PanCan".

- infunc_exclude_cancer:

  Character. Cancer type to exclude. Default = "PanCan".

- infunc_GOI:

  Character vector. Genes of Interest.

- ncores:

  Integer. Kept for API compatibility; ignored (see NOTE above).

## Value

A list of feature ranking results, one per drug.

## Details

NOTE: feature ranking is fully vectorized (~0.15 s per drug over 15k
genes), so this deliberately runs SERIALLY. The old Windows PSOCK path
serialized the 567 MB DepMap into every worker (clusterExport) and then
still failed — workers are clean R sessions without the package's
internal helpers — before silently falling back to serial; that cost
~10x for zero benefit.
