# Get specific drug response data for cell-lines

Extracts the drug response (AUC) for a given drug from the DepMap
database. The function handles multiple screening batches by
prioritizing MTS over HTS and selecting the batch with the fewest
missing values.

## Usage

``` r
get_response_matrix(infunc_drugName)
```

## Arguments

- infunc_drugName:

  Character string. Name of the drug (e.g., "erlotinib").

## Value

A named numeric vector of AUC values for all cell lines.
