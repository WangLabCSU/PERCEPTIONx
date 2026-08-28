# Hypergeometric test for gene list overlap

Tests whether the overlap between two gene lists is statistically
significant. It answers the question: given a background set of genes
(global), a reference gene set (base_list), and a test gene set
(test_list), is the overlap between test_list and base_list more than
expected by chance?

## Usage

``` r
hypergeometric_test_for_twolists(
  test_list,
  base_list,
  global,
  lower.tail = FALSE
)
```

## Arguments

- test_list:

  Character vector. The test gene list.

- base_list:

  Character vector. The reference gene list.

- global:

  Character vector. The background gene set.

- lower.tail:

  Logical. If FALSE (default), calculates the probability of observing
  *greater than or equal to* the observed overlap (enrichment). If TRUE,
  calculates the probability of observing *less than or equal to* the
  observed overlap (depletion).

## Value

Numeric p-value.
