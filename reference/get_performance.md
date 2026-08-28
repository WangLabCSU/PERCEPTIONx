# Load performance metrics from a saved model file

Loads a saved model RDS file and returns only its performance portion.
The file is deserialized fully first (an RDS is read as a whole), then
the performance entries are extracted.

## Usage

``` r
get_performance(filepath, drug_names = NULL)
```

## Arguments

- filepath:

  Character. Path to a single saved model RDS file.

- drug_names:

  Character vector. Optional. Specific drug names to load. If NULL,
  loads all drugs in the file.

## Value

A list containing performance metrics for each drug model.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Load performance from a single model file
  perf <- get_performance("models/PERCEPTIONx_models_PanCan_exPanCan_20240101_120000.RDS")

  # Or use compare_performance on the result
  models <- readRDS("models/PERCEPTIONx_models_PanCan_exPanCan_20240101_120000.RDS")
  perf <- compare_performance(models)
} # }
```
