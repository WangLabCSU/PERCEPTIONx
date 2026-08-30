# Launch PERCEPTIONx Shiny Dashboard

Starts the interactive PERCEPTIONx dashboard for drug response
prediction. The dashboard provides a graphical interface for data
loading, model training, prediction, and visualization.

## Usage

``` r
run_perception_app()
```

## Value

A Shiny app object. When run interactively (no `SHINY_PORT` environment
variable set) the app is also launched in the user's browser.

## Details

Under a Shiny Server deployment the server executes this function and
expects a shiny app *object* back (a nested `runApp()` would fail), so
the function detects Shiny Server via the `SHINY_PORT` environment
variable it sets for every R worker (the same signal the shiny package
itself uses internally) and returns the app object instead of starting
its own server. This lets the same one-liner work both locally and
inside `apps/<name>/app.R`.

## Examples

``` r
if (FALSE) { # \dontrun{
library(PERCEPTIONx)
run_perception_app()
} # }
```
