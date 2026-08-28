# FDR correction for multiple testing

Performs False Discovery Rate (FDR) correction using the
Benjamini-Hochberg method. This is a wrapper around
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) with method = "fdr".

## Usage

``` r
fdrcorr(test_list)
```

## Arguments

- test_list:

  A numeric vector of p-values to be corrected.

## Value

A numeric vector of FDR-adjusted p-values (q-values).
