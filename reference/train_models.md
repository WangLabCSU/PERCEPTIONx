# Train PERCEPTIONx models for multiple drugs

This function runs the complete PERCEPTIONx training pipeline for a list
of drugs. It performs feature ranking, model training with
hyperparameter tuning across different k_features values, and selects
the best model per drug based on single-cell test performance.

## Usage

``` r
train_models(
  drug_list = NULL,
  cancer_type = "PanCan",
  exclude_cancer = "PanCan",
  GOI = NULL,
  k_features_values = NULL,
  ncores = 4,
  output_dir = "./models",
  model_type = "glmnet",
  num_folds = 3,
  num_tree = 500,
  seed = 1,
  alpha_gradient = 0.05,
  lambda_gradient = 20,
  lambda_range = c(1e-04, 1),
  cv_method = "cv",
  progress_cb = NULL
)
```

## Arguments

- drug_list:

  Character vector. Names of drugs to train models for. If NULL, uses
  the 44 FDA-approved drugs from the paper.

- cancer_type:

  Character. Cancer type for training. Default = "PanCan".

- exclude_cancer:

  Character. Cancer type to exclude from training. Default = "PanCan".

- GOI:

  Character vector. Genes of Interest to use as features.

- k_features_values:

  Numeric vector. Feature counts to test during tuning. If NULL,
  automatically calculated from expression_rnorm dimensions.

- ncores:

  Integer. Number of CPU cores for parallel feature ranking. Default =
  4.

- output_dir:

  Character. Directory to save trained model RDS file. Default =
  "./models".

- model_type:

  Character. Model type: "glmnet" or "rf". Default = "glmnet".

- num_folds:

  Integer. Number of CV folds. Default = 3.

- num_tree:

  Integer. Number of trees for random forest. Default = 500.

- seed:

  Integer. Random seed for training. Default = 1.

- alpha_gradient:

  Numeric. Step size for alpha grid in glmnet tuning. Default = 0.05.

- lambda_gradient:

  Integer. Number of lambda values in glmnet tuning grid. Default = 20.

- lambda_range:

  Numeric vector of length 2. Min and max lambda for tuning grid.
  Default = c(0.0001, 1).

- cv_method:

  Character. Cross-validation method. Default = "cv".

## Value

A named list of trained model objects, one per drug. Also saved as a
single RDS file.

## Examples

``` r
if (FALSE) { # \dontrun{
  load_depmap(read = TRUE)
  available_genes <- intersect(rownames(DepMap$expression_20Q4),
                               rownames(DepMap$scRNA_complete))
  set.seed(123)
  GOI_100 <- sample(available_genes, 100)
  train_models("abemaciclib", "PanCan", "PanCan", GOI_100, ncores = 1)
} # }
```
