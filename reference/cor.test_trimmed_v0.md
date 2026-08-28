# Fast correlation test (generic)

A generic function for fast correlation testing. Dispatches methods
based on the class of the first argument. This is a trimmed-down version
of the base R `cor.test` function, optimized for speed in large-scale
computations.

## Usage

``` r
cor.test_trimmed_v0(x, ...)
```

## Arguments

- x:

  First argument (typically a numeric vector).

- ...:

  Additional arguments passed to methods.

## Value

A list containing the p-value and correlation estimate. The exact
structure depends on the method used.

## See also

[`cor.test`](https://rdrr.io/r/stats/cor.test.html) for the full version
with confidence intervals.
