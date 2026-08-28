# Run Seurat clustering and plot 2D embedding

Performs Seurat clustering on an expression matrix and generates a 2D
embedding visualization (UMAP or t-SNE). Useful for identifying
subclones within patient tumor samples.

## Usage

``` r
plot_seurat_clustering(
  method = c("umap", "tsne"),
  expression_matrix,
  min_cells = 3,
  min_features = 200,
  nfeatures = 2000,
  dims = 10,
  resolution = 0.8,
  seed = 1
)
```

## Arguments

- method:

  Character. Dimensionality reduction method. One of `"umap"` (default)
  or `"tsne"`.

- expression_matrix:

  Matrix. Gene expression matrix (genes as rows, cells as columns).

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

  Integer. Random seed. Default = 1.

## Value

A list containing:

- seurat_object:

  Seurat object with clustering results

- embedding_plot:

  ggplot 2D embedding visualization

- cluster_ids:

  Named vector of cluster IDs per cell

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- plot_seurat_clustering(patient_expression)
  result$embedding_plot
  result$cluster_ids
} # }
```
