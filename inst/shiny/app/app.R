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

# Long synchronous tasks (model training, large uploads/parsing) block the
# single R thread, so the browser heartbeat cannot be answered while they run.
# The default heartbeat timeout is 60s — multi-drug training can exceed that,
# the browser drops the WebSocket, and the "hiding" messages (waiter, progress)
# are lost, leaving spinners stuck on screen. Raise it so sessions survive.
options(shiny.heartbeat.timeout = 1800)  # 30 minutes, then consider dead

# Auto-detect package root robustly, regardless of how the app is launched
# (cd inst/shiny/app && Rscript app.R, shiny::runApp("inst/shiny/app") from the
# repo root, RStudio "Run App", etc.). A candidate only counts as the LIVE repo
# if it is a SOURCE tree — DESCRIPTION plus R/*.R files. An INSTALLED package
# directory also contains DESCRIPTION (but only compiled R/*.rdb, no .R), and
# devtools::load_all() on it silently produces an EMPTY namespace (every
# PERCEPTIONx:: call then fails with "not an exported object").
is_src_root <- function(cand) {
  !is.null(cand) &&
    file.exists(file.path(cand, "DESCRIPTION")) &&
    length(list.files(file.path(cand, "R"), pattern = "\\.R$")) > 0
}
pkg_root <- NULL
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
candidates <- c(
  if (!is.null(script_dir)) file.path(script_dir, "..", "..", ".."),
  file.path(getwd(), "..", "..", ".."),  # cwd = inst/shiny/app
  file.path(getwd(), "..", ".."),        # cwd = inst/shiny
  file.path(getwd(), ".."),              # cwd = inst
  file.path(getwd())                     # cwd = repo root
)
for (cand in candidates) {
  if (is_src_root(cand)) {
    pkg_root <- normalizePath(cand)
    break
  }
}
if (!is.null(pkg_root)) {
  suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
} else if (requireNamespace("PERCEPTIONx", quietly = TRUE)) {
  warning("Repo root not detected (launched from an unusual directory); ",
          "falling back to the INSTALLED PERCEPTIONx package.")
  library(PERCEPTIONx)
} else {
  stop("Neither devtools (with repo) nor PERCEPTIONx is available. Please install PERCEPTIONx.")
}

# Background worker (callr) must load the SAME live repo code — remember where
# it is so async_jobs.R can hand the path to the worker process.
options(perception.pkg_root = pkg_root)

# Source modules
source("R/shiny_helpers.R")
source("R/async_jobs.R")
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
    tags$head(tags$link(rel = "stylesheet",
                        href = paste0("tour.css?v=", format(Sys.time(), "%Y%m%d%H%M%S")))),
    tags$head(tags$script(src = paste0("tour.js?v=", format(Sys.time(), "%Y%m%d%H%M%S")))),
    tags$head(tags$link(rel = "icon", href = "favicon.svg", type = "image/svg+xml")),
    tags$head(tags$script(HTML("
      Shiny.addCustomMessageHandler('scroll-to', function(id) {
        var el = document.getElementById(id);
        if (el) el.scrollIntoView({behavior: 'smooth', block: 'start'});
      });

      // Replace a container's innerHTML. Used by the Home page to update the
      // workflow stepper / data dashboard after the static first paint.
      Shiny.addCustomMessageHandler('set-html', function(msg) {
        var el = document.getElementById(msg.id);
        if (el) el.innerHTML = msg.html;
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
      href = "#",
      id = "tour-start",
      icon("circle-question", class = "nav-icon"),
      title = "Guided tour"
    )
  )),
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
    depmap_meta = NULL,        # lightweight DepMap metadata only (genes/drugs/dims)
    depmap_path = NULL,
    depmap_is_standard = FALSE,
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

  # NOTE: no session-end file sweep here. Shiny itself deletes a session's
  # uploaded files when it ends. A global unlink over tempdir() would delete
  # OTHER sessions' DepMap.RDS that a background worker may still need to read
  # from disk (built-in downloads all share tempdir()/DepMap.RDS), breaking
  # their training. R's tempdir() is cleaned up by the OS on process exit.

  session$onSessionEnded(function() {
    # NOTE: this callback runs OUTSIDE any reactive context. Reading a
    # reactiveValues field that was never assigned (e.g. a user who never
    # triggered any background task) errors there — wrap every read in
    # isolate() so missing fields come back as NULL instead of aborting.
    # 1. Kill THIS session's background workers. callr's supervise = TRUE only
    # kills children when the WHOLE server process exits, so without this every
    # session leaks its worker (each one holds that session's data in RAM —
    # prepared Seurat objects, prediction matrices, or even a full uploaded
    # DepMap) until the server restarts. The training MASTER is server-wide
    # and must NOT be touched; shared$train_worker is always the per-session
    # custom worker, shared$task_worker the prepare/demo/predict/plot worker.
    tw <- isolate(shared$task_worker)
    if (!is.null(tw) && inherits(tw, "r_process") && tw$is_alive()) {
      tw$kill()
    }
    cw <- isolate(shared$train_worker)
    if (!is.null(cw) && inherits(cw, "r_process") && cw$is_alive()) {
      cw$kill()
    }
    # 2. Drop this session's private job dirs (perception_tasks/<id> and, for
    # uploaded DepMap users, perception_jobs_custom/<id>). They are session-
    # scoped — no other session can be reading them — and any job that was
    # mid-flight has been aborted by the kill above.
    for (dir_key in c("task_jobs_dir", "jobs_dir")) {
      d <- isolate(shared[[dir_key]])
      if (!is.null(d) && nzchar(d) &&
          grepl("perception_tasks|perception_jobs_custom", d, fixed = FALSE)) {
        unlink(d, recursive = TRUE)
      }
    }
  })
}

shinyApp(ui, server)
