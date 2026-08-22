# PERCEPTION-shiny — Shiny Web Application of the PERCEPTIONx R Package
# Theme: Deep Indigo + Warm Amber

library(shiny)
library(bslib)
library(DT)
library(plotly)
library(waiter)
library(thematic)
library(ggplot2)

# Allow large uploads: single-cell matrices / RDS reference files can far
# exceed Shiny's default 5 MB limit.
options(shiny.maxRequestSize = 1 * 1024^3)  # 1 GB (DepMap.RDS is ~567 MB)

# Auto-detect package root: works from inst/shiny/app/ -> repo root
# Falls back to installed PERCEPTIONx package if devtools not available.
pkg_root <- if (requireNamespace("devtools", quietly = TRUE)) {
  normalizePath(file.path(getwd(), "..", "..", ".."), mustWork = FALSE)
} else {
  NULL
}
if (!is.null(pkg_root) && file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
} else if (requireNamespace("PERCEPTIONx", quietly = TRUE)) {
  library(PERCEPTIONx)
} else {
  stop("Neither devtools (with repo) nor PERCEPTIONx is available. Please install PERCEPTIONx.")
}

# Source modules
source("R/shiny_helpers.R")
source("R/mod_home.R")
source("R/mod_data.R")
source("R/mod_train.R")
source("R/mod_predict.R")
source("R/mod_visualize.R")
source("R/mod_help.R")

# Theme
perception_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#165C91",
  danger = "#c13232",
  font_scale = 0.95,
  `enable-gradients` = TRUE,
  `enable-shadows` = TRUE,
  `body-bg` = "#f0f2f8",
  `body-color` = "#1e2a4a",
  `border-radius` = "0.7rem",
  `btn-border-radius` = "0.7rem"
)

ui <- page_navbar(
  theme = perception_theme,
  title = tagList(
    tags$span(class = "brand-icon",
      tags$img(src = "favicon.svg", height = "28", alt = "PERCEPTION-shiny",
               style = "vertical-align: middle; margin-right: 0.3rem;")
    ),
    tags$span(class = "brand-text", "PERCEPTION-shiny")
  ),
  id = "navbar",
  selected = "home",
  fillable = FALSE,
  header = tagList(
    use_waiter(),
    # Cache-busting: version the CSS with the app start time so browsers
    # always fetch the newest stylesheet after a restart (no stale cache).
    tags$head(tags$link(rel = "stylesheet",
                        href = paste0("styles.css?v=", format(Sys.time(), "%Y%m%d%H%M%S")))),
    tags$head(tags$link(rel = "icon", href = "favicon.svg", type = "image/svg+xml")),
    tags$head(tags$script(HTML("
      Shiny.addCustomMessageHandler('scroll-to', function(id) {
        var el = document.getElementById(id);
        if (el) el.scrollIntoView({behavior: 'smooth', block: 'start'});
      });

      // Keep the viewport still when the hidden file input receives focus on
      // 'Browse...' click (Shiny positions it at top:-1e6px, so the browser
      // would otherwise scroll the page to the top). We must NOT
      // preventDefault() here: Shiny relies on the <label> default action to
      // open the file dialog. Restoring the scroll position instead.
      $(document).on('mousedown', '.btn-file', function() {
        window.__perceptionScrollY = window.pageYOffset || document.documentElement.scrollTop;
      });
      $(document).on('click', '.btn-file', function() {
        var y = window.__perceptionScrollY || 0;
        [0, 50, 200, 500].forEach(function(ms) {
          setTimeout(function() { window.scrollTo(0, y); }, ms);
        });
      });
    ")))
  ),
  nav_spacer(),
  nav_item(tagList(
    tags$a(
      href = "https://github.com/WangLabCSU/PERCEPTIONx/blob/main/docs/PERCEPTION-shiny.md",
      target = "_blank",
      icon("book-open", class = "nav-icon"),
      title = "Tutorial"
    )
  )),
  nav_item(tagList(
    tags$a(
      href = "https://github.com/WangLabCSU/PERCEPTIONx",
      target = "_blank",
      icon("github", class = "nav-icon")
    )
  )),
  nav_panel(
    " Home",
    value = "home",
    icon = icon("house"),
    mod_home_ui("home")
  ),
  nav_panel(
    " Data",
    value = "data",
    icon = icon("database"),
    mod_data_ui("data")
  ),
  nav_panel(
    " Train",
    value = "train",
    icon = icon("brain"),
    mod_train_ui("train")
  ),
  nav_panel(
    " Predict",
    value = "predict",
    icon = icon("microscope"),
    mod_predict_ui("predict")
  ),
  nav_panel(
    " Visualize",
    value = "visualize",
    icon = icon("chart-line"),
    mod_visualize_ui("visualize")
  ),
  nav_panel(
    " Help",
    value = "help",
    icon = icon("circle-question"),
    mod_help_ui("help")
  )
)

server <- function(input, output, session) {
  # Shared reactive values
  shared <- reactiveValues(
    depmap = NULL,
    models = NULL,
    user_expr = NULL,
    user_mapping = NULL,
    user_clones = NULL,
    user_response = NULL,
    prepared_data = NULL,
    predictions = NULL,
    patient_pred = NULL,
    model_cache = list(),
    model_active = list()
  )

  # Modules — pass main session for cross-module navigation
  mod_home_server("home", shared, session)
  mod_data_server("data", shared)
  mod_train_server("train", shared, session)
  mod_predict_server("predict", shared, session)
  mod_visualize_server("visualize", shared, session)
  mod_help_server("help", session)

  # Clean up downloaded temp files when the session ends (browser tab closed)
  session$onSessionEnded(function() {
    td <- tempdir()
    if (dir.exists(td)) {
      files <- list.files(td, pattern = "\\.(RDS|rds)$", full.names = TRUE, recursive = TRUE)
      if (length(files) > 0) unlink(files, force = TRUE)
    }
  })
}

shinyApp(ui, server)
