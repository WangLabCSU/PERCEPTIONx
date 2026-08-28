# Get the DepMap dataset

Retrieves the DepMap dataset from the package-level cache. This is the
recommended way to access DepMap data after loading with
`load_depmap(read = TRUE)`.

## Usage

``` r
get_depmap()
```

## Value

The DepMap list object. Errors if the data has not been loaded yet.

## Examples

``` r
if (FALSE) { # \dontrun{
load_depmap(read = TRUE)
depmap <- get_depmap()
} # }
```
