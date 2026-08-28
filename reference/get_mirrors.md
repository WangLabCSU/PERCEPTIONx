# Get current download mirrors

Returns the current mirror list, including any user-added mirrors.
User-added mirrors are tried first by default.

## Usage

``` r
get_mirrors()
```

## Value

Character vector of mirror URLs.

## Examples

``` r
get_mirrors()
#>  [1] "https://gh-proxy.com/https://github.com"              
#>  [2] "https://ghproxy.net/https://github.com"               
#>  [3] "https://moeyy.cn/gh-proxy/https://github.com"         
#>  [4] "https://github.akams.cn/https://github.com"           
#>  [5] "http://toolwa.com/github/https://github.com"          
#>  [6] "https://v6.gh-proxy.org/https://github.com"           
#>  [7] "https://gh-proxy.org/https://github.com"              
#>  [8] "https://ghfast.top/https://github.com"                
#>  [9] "https://download.githubcdn.com?url=https://github.com"
#> [10] "https://proxy.gitwarp.top/https://github.com"         
```
