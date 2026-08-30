# Download filtered DepMap data

Downloads the filtered DepMap RDS file for training models from GitHub
Release. The file contains bulk expression, scRNA expression, drug
response annotations, and cell line metadata required for PERCEPTIONx
model training. Supports automatic mirror fallback for users in
different regions.

## Usage

``` r
load_depmap(
  dest = perception_default_depmap_dir(),
  read = FALSE,
  mirror = FALSE,
  mirror_url = NULL,
  timeout_seconds = 300,
  retries = 1
)
```

## Arguments

- dest:

  Directory to save the downloaded file. Default: a cache-root derived
  path (see `options(PERCEPTIONx.cache_root)`), or "." when no cache
  root is set.

- read:

  Whether to read the data and assign to global environment as "DepMap".
  Default = FALSE.

- mirror:

  Logical. If FALSE (default), download from GitHub directly. If TRUE,
  use mirror sites from
  [`get_mirrors()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_mirrors.md).

- mirror_url:

  Character. A specific mirror URL to use (e.g.,
  `"https://gh-proxy.com/https://github.com"`). Overrides `mirror`.

- timeout_seconds:

  Numeric, timeout for each download attempt in seconds. Default = 300.

- retries:

  Integer, number of retries for each mirror. Default = 1.

## Value

Invisibly returns the DepMap object if read = TRUE, otherwise NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download from GitHub directly
load_depmap(read = TRUE)

# Use mirror sites
load_depmap(read = TRUE, mirror = TRUE)

# Use a specific mirror
load_depmap(read = TRUE, mirror_url = "https://gh-proxy.com/https://github.com")
} # }
```
