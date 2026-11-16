# Read a 10x Genomics count matrix

One reader for the inputs a user actually has: a CellRanger output
folder, a `.h5` file, or the three files passed explicitly. The result
is a sparse genes x cells matrix with gene symbols as row names, which
is what
[`preprocess_sc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/preprocess_sc.md)
and
[`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
expect.

## Usage

``` r
read_10x(path = NULL, mtx = NULL, barcodes = NULL, features = NULL)
```

## Arguments

- path:

  Directory holding `matrix.mtx[.gz]`, `barcodes.tsv[.gz]` and
  `genes.tsv[.gz]` / `features.tsv[.gz]`, or a `.h5` file written by
  CellRanger. Ignored when the three files are given explicitly.

- mtx:

  Path to `matrix.mtx` or `matrix.mtx.gz`.

- barcodes:

  Path to `barcodes.tsv` or `barcodes.tsv.gz`.

- features:

  Path to `genes.tsv` or `features.tsv` (either `.gz` or not).

## Value

A `dgCMatrix` with genes as rows and cells as columns.

## Examples

``` r
if (FALSE) { # \dontrun{
counts <- read_10x("path/to/filtered_feature_bc_matrix")
counts <- read_10x("path/to/filtered_feature_bc_matrix.h5")
counts <- read_10x(
  mtx = "matrix.mtx.gz",
  barcodes = "barcodes.tsv.gz",
  features = "features.tsv.gz"
)
} # }
```
