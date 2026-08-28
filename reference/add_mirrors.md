# Add custom download mirrors

Add custom download mirrors

## Usage

``` r
add_mirrors(urls, position = c("first", "last"))
```

## Arguments

- urls:

  Character vector of mirror URLs to add.

- position:

  Character. Where to add: "first" or "last". Default: "first".

## Value

Invisibly returns the updated mirror list.

## Examples

``` r
if (FALSE) { # \dontrun{
add_mirrors("https://my-mirror.com/https://github.com")
add_mirrors(c("https://mirror1.com/https://github.com",
              "https://mirror2.com/https://github.com"))
} # }
```
