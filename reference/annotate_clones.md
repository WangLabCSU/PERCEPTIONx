# Annotate cells with clone IDs via Seurat clustering

Performs Seurat clustering on a single-cell expression matrix and
returns a mapping of each cell to its cluster (clone) ID. This matches
the original PERCEPTIONx pipeline where Seurat clusters define
transcriptional subclones.

## Usage

``` r
annotate_clones(
  method = c("umap", "tsne"),
  expression_matrix,
  min_cells = 3,
  min_features = 200,
  nfeatures = 2000,
  dims = 10,
  resolution = 0.8,
  seed = 42,
  progress_cb = NULL,
  cluster_algorithm = 1,
  variable_selection = c("vst", "dispersion", "mvp")
)
```

## Arguments

- method:

  Character. Dimensionality reduction method. One of `"umap"` (default)
  or `"tsne"`. UMAP is faster and preserves global structure better;
  t-SNE emphasizes local neighborhoods.

- expression_matrix:

  Matrix. Gene expression matrix with genes as rows and cells as
  columns. Raw counts or normalized values are both accepted.

- min_cells:

  Integer. Minimum cells per feature. Default = 3.

- min_features:

  Integer. Minimum features per cell. Default = 200.

- nfeatures:

  Integer. Number of variable features. Default = 2000.

- dims:

  Integer. Number of PCA dimensions for clustering. Default = 10.

- resolution:

  Numeric. Clustering resolution. Default = 0.8.

- seed:

  Integer. Random seed for reproducibility. Default = 42.

- progress_cb:

  Optional function called with the current Seurat stage name (e.g.
  `"pca"`, `"clustering"`) between pipeline steps. Used to show live
  progress in the Shiny app.

- cluster_algorithm:

  Integer. Community-detection algorithm passed to
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html).
  Default = 1 (Louvain, no extra dependency). Use 4 (Leiden) only if the
  leidenbase package is installed.

- variable_selection:

  Character. Variable-feature selection method for
  [`Seurat::FindVariableFeatures()`](https://satijalab.org/seurat/reference/FindVariableFeatures.html).
  One of `"vst"` (default, most accurate), `"dispersion"` or `"mvp"`
  (faster on large matrices at a small cost in feature quality).

## Value

A data frame with columns: `cell_id`, `clone_id`, and `dim_1`, `dim_2`
(2D embedding coordinates for visualization).

## Examples

``` r
if (FALSE) { # \dontrun{
  cell_clone_map <- annotate_clones(patient_expression)
  cell_clone_map <- annotate_clones("tsne", patient_expression)
} # }
```
