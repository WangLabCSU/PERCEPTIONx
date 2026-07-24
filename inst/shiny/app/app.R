# PERCEPTION Shiny App — Embedded in R Package
# Theme: Deep Indigo + Warm Amber

library(shiny)
library(bslib)
library(DT)
library(plotly)
library(waiter)
library(thematic)
library(ggplot2)

# Use PERCEPTION package functions (already loaded via run_perception_app)
if (!requireNamespace("PERCEPTION", quietly = TRUE)) {
  stop("PERCEPTION package is required. Please install it first.")
}

# Source modules
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
      tags$img(src = "favicon.svg", height = "28", alt = "PERCEPTION",
               style = "vertical-align: middle; margin-right: 0.3rem;")
    ),
    tags$span(class = "brand-text", "PERCEPTION")
  ),
  id = "navbar",
  selected = "home",
  fillable = FALSE,
  header = tagList(
    use_waiter(),
    tags$head(tags$link(rel = "stylesheet", href = "styles.css")),
    tags$head(tags$link(rel = "icon", href = "favicon.svg", type = "image/svg+xml")),
    tags$head(tags$script(HTML("
      Shiny.addCustomMessageHandler('scroll-to', function(id) {
        var el = document.getElementById(id);
        if (el) el.scrollIntoView({behavior: 'smooth', block: 'start'});
      });
    ")))
  ),
  nav_spacer(),
  nav_item(tagList(
    tags$a(
      href = "https://github.com/SunPast/PERCEPTION",
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
}

shinyApp(ui, server)
