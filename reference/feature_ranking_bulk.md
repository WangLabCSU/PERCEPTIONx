# Feature ranking for a single drug using bulk expression

Computes Pearson correlation between each gene's expression and drug
response, then ranks genes by absolute correlation. This matches the
original feature_ranking_bulk function from step0B.

## Usage

``` r
feature_ranking_bulk(
  infunc_drugName,
  infunc_cancerType = "PanCan",
  exclude_cancer = "PanCan",
  infunc_GOI,
  ncores = 4L
)
```

## Arguments

- infunc_drugName:

  Character. Name of the drug.

- infunc_cancerType:

  Character. Cancer type for training. Default = "PanCan".

- exclude_cancer:

  Character. Cancer type to exclude. Default = "PanCan".

- infunc_GOI:

  Character vector. Genes of Interest to rank.

- ncores:

  Number of cores to use for parallelization. Default = 4.

## Value

Matrix of ranked features with columns: p.value, estimate.cor (sorted by
abs(cor) descending)
