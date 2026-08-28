# Internal: predict viability for a single drug model

Extracts features from the model, matches them to expression matrix
rows, and predicts viability scores. Handles gene name format
differences (e.g., hyphens vs dots) via make.names() fallback.

## Usage

``` r
viability_from_model_internal(drug_name, model, dataset)
```

## Arguments

- drug_name:

  Character. Drug name (for error messages).

- model:

  A caret model object (from build_on_BULK_v2 output \$model).

- dataset:

  Matrix. Expression matrix with genes as rows and cells/samples as
  columns. Must be rank-normalized.

## Value

A named numeric vector of predicted viability scores.
