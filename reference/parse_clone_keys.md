# Parse clone keys in "Patient@Clone" format

Splits a character vector of clone identifiers (e.g. "PAT_001@0") into a
data frame with `patient` and `clone_id` columns.

## Usage

``` r
parse_clone_keys(x)
```

## Arguments

- x:

  Character vector of clone keys in `"patient@clone_id"` format.

## Value

A data frame with columns `patient` and `clone_id`.

## See also

[`build_clone_key()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/build_clone_key.md)
for the inverse operation.
