# Rank-normalize each column of a matrix

Converts each column to ranks divided by column length, producing values
in the range 0 to 1 (exclusive of 0, inclusive of 1). This is a
**crucial preprocessing step** when using your own expression data with
PERCEPTIONx.

## Usage

``` r
rank_normalization_mat(mat)
```

## Arguments

- mat:

  A numeric matrix. Rows = genes, columns = cells/samples.

## Value

A rank-normalized matrix of the same dimensions, with values between 0
and 1.

## Why rank normalization?

PERCEPTIONx models are trained on DepMap expression data that has been
rank-normalized. The model coefficients capture the relationship between
**relative gene expression ranks** and drug response, not absolute
expression values. Therefore, any new expression data fed into the model
must undergo the same rank normalization to ensure compatibility.

## How it works

For each column (cell/sample), every gene's expression value is replaced
by its rank within that column, divided by the total number of genes:
\$\$x\_{ij}^{norm} = \frac{\mathrm{rank}(x\_{ij})}{n}\$\$ where \\n\\ is
the number of rows (genes) and ties are resolved by averaging. This
transforms each column into a uniform distribution over \\(0,1\]\\,
making the data robust to batch effects, library size differences, and
outliers.

## Important

If you provide your own expression data that has **not** been
rank-normalized, predictions will be unreliable. Always run
`rank_normalization_mat()` on your raw expression matrix before passing
it to
[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md).

## Examples

``` r
raw_expr <- matrix(c(10, 200, 50, 30, 300, 80), nrow = 3, ncol = 2)
rownames(raw_expr) <- c("GENE_A", "GENE_B", "GENE_C")
norm_expr <- rank_normalization_mat(raw_expr)
```
