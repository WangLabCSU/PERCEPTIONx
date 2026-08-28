# Determine training and test cell-lines for a given drug

Identifies which cell lines should be used for training and which should
be excluded (reserved for testing) based on cancer type and single-cell
data availability.

## Usage

``` r
get_cellLine_list(
  infunc_cancerType = "PanCan",
  infunc_drugName,
  exclude_cancer = "PanCan",
  infunc_response,
  force_add_cellLines = TRUE,
  force_add_cellLines_list = NA
)
```

## Arguments

- infunc_cancerType:

  Character. Cancer type for training. Default = "PanCan".

- infunc_drugName:

  Character. Name of the drug.

- exclude_cancer:

  Character. Cancer type to exclude from training. Default = "PanCan".

- infunc_response:

  Named numeric vector. Drug response data.

- force_add_cellLines:

  Logical. Whether to force add additional cell lines. Default = TRUE.

- force_add_cellLines_list:

  Character vector. Cell line IDs to force add. Default = NA.

## Value

A list of length 2:

- \[1\] common_cellLines: Cell lines to use for training

- \[2\] cellLines2remove: Cell lines excluded (reserved for testing)
