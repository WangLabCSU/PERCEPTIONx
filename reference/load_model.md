# Load pre-built model of provided drugs

Downloads model files from GitHub Release if not cached locally.
Supports automatic mirror fallback for users in different regions.

## Usage

``` r
load_model(
  ...,
  dest = perception_default_model_dir(),
  read = FALSE,
  mirror = FALSE,
  mirror_url = NULL,
  timeout_seconds = 30,
  retries = 0,
  force = FALSE,
  all = FALSE
)
```

## Arguments

- ...:

  One or more drug names (e.g., "erlotinib", "gefitinib").

- dest:

  Directory to save downloaded models. Default: a cache-root derived
  path (see `options(PERCEPTIONx.cache_root)`), or "./models" when no
  cache root is set.

- read:

  Whether to read and return the downloaded model(s) as a named list.
  Default = FALSE (download only).

- mirror:

  Logical. If FALSE (default), download from GitHub directly. If TRUE,
  use mirror sites from
  [`get_mirrors()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_mirrors.md).

- mirror_url:

  Character. A specific mirror URL to use (e.g.,
  `"https://gh-proxy.com/https://github.com"`). Overrides `mirror`.

- timeout_seconds:

  Numeric, timeout for each download attempt in seconds. Default = 30.

- retries:

  Integer, number of retries for each mirror. Default = 0.

- force:

  Logical. If TRUE, delete an existing cached model file before
  downloading it again (use after an interrupted/failed download left a
  corrupt file behind). Default = FALSE (keep and reuse the cache).

- all:

  Logical. If TRUE, download all 44 pre-trained models (useful to
  pre-cache the full model set on a deployment server). Overrides any
  drug names given in `...`. Default = FALSE.

## Value

If `read = TRUE`, a named list of model objects (names = drug names). If
`read = FALSE`, invisibly returns NULL (files are saved to disk only).
The returned list can be directly passed to
[`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md),
[`compare_performance()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/compare_performance.md),
or
[`get_significant_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_significant_models.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Download and load models from GitHub directly
models <- load_model("erlotinib", "gefitinib", read = TRUE)

# Use mirror sites
models <- load_model("erlotinib", "gefitinib", read = TRUE, mirror = TRUE)

# Use a specific mirror
models <- load_model("erlotinib", read = TRUE,
                     mirror_url = "https://gh-proxy.com/https://github.com")
} # }
```
