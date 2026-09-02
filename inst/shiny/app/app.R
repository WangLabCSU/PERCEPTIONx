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

# ---------------------------------------------------------------------------
# sass/bslib compile cache. On hosted deployments (Shiny Server) sass defaults
# the cache to <appDir>/app_cache/sass, but this app lives INSIDE the
# installed package, so that resolves to the read-only package install dir —
# every session then warns "cannot create dir .../app_cache (Permission
# denied)" and falls back to tempdir. Pin the cache to a writable directory
# (tempdir is always writable) so the fallback never triggers. Note: the old
# sass.cache_dir option is defunct in sass >= 0.4.1 — sass.cache now takes a
# directory path directly.
options(sass.cache = file.path(tempdir(check = TRUE), "sass-cache"))

# ---------------------------------------------------------------------------
# Static assets are INLINED into the HTML (no external file requests).
#
# The app can run two ways: locally via runApp() on the package's
# inst/shiny/app, or as a single-file launcher (apps/<name>/app.R) under
# Shiny Server. In the Shiny Server case the app directory has no www/ subdir,
# so asset URLs like "styles.css" 404 (the real files live in the installed
# package). Inlining every asset removes the www/ dependency entirely: the
# stylesheets and the tour script are embedded, and the favicon is embedded as
# a data URI (both in <link rel="icon"> and inside styles.css backgrounds).
# ---------------------------------------------------------------------------
perception_asset <- function(name) {
  paste(readLines(file.path(perception_www_dir, name), warn = FALSE), collapse = "\n")
}
perception_favicon_uri <- function() {
  svg <- perception_asset("favicon.svg")
  paste0("data:image/svg+xml;charset=utf-8,", utils::URLencode(svg, reserved = TRUE))
}

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
# Assets come from the LIVE source tree in dev, and from the installed package
# on the server. This matters: system.file() would silently serve a STALE copy
# whenever the package was installed before the latest CSS change, so "I edited
# styles.css but nothing happened" is exactly what users hit locally. In dev we
# read the working tree so CSS edits apply immediately.
perception_www_dir <- if (!is.null(pkg_root)) {
  file.path(pkg_root, "inst", "shiny", "app", "www")
} else {
  file.path(system.file("shiny", "app", package = "PERCEPTIONx"), "www")
}
if (!is.null(pkg_root)) {
  suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
} else if (requireNamespace("PERCEPTIONx", quietly = TRUE)) {
  # Deployed mode (installed package, Shiny Server): no source tree is
  # available — silently use the installed package. No warning here; this is
  # the normal, expected path for deployments.
  library(PERCEPTIONx)
} else {
  stop("Neither devtools (with repo) nor PERCEPTIONx is available. Please install PERCEPTIONx.")
}

# Background worker (callr) must load the SAME live repo code — remember where
# it is so async_jobs.R can hand the path to the worker process.
options(perception.pkg_root = pkg_root)

# Source modules into the CURRENT environment (local = TRUE). Under a local
# runApp() the app file executes in the global env, so this is a no-op change;
# but the Shiny Server path (.perception_app_object()) evaluates app.R in a
# private environment, and modules must live in that SAME environment — with
# the default source(local = FALSE) they would land in the global env and lose
# access to app.R-top-level helpers (e.g. perception_favicon_uri()).
source("R/shiny_helpers.R", local = TRUE)
source("R/async_jobs.R", local = TRUE)
source("R/mod_home.R", local = TRUE)
source("R/mod_data.R", local = TRUE)
source("R/mod_train.R", local = TRUE)
source("R/mod_predict.R", local = TRUE)
source("R/mod_visualize.R", local = TRUE)
source("R/mod_help.R", local = TRUE)

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
      tags$img(src = perception_favicon_uri(), height = "28", alt = "PERCEPTION-shiny",
               style = "vertical-align: middle; margin-right: 0.3rem;")
    ),
    tags$span(class = "brand-text", "PERCEPTION-shiny")
  ),
  id = "navbar",
  selected = "home",
  fillable = FALSE,
  header = tagList(
    use_waiter(),
    # Pre-built full-screen demo overlay, hidden by default. The Load Demo
    # buttons show it instantly from their onclick (no ~1s round trip before
    # feedback); the server hides it via 'hide-demo-overlay' when the demo
    # task finishes or fails.
    div(id = "demo-overlay", class = "demo-overlay", style = "display: none;",
      div(class = "spinner-ring"),
      h4("Preparing demo data..."),
      p(id = "demo-stage", class = "text-muted", "Running Seurat clustering...")),
    # Assets are inlined (see the perception_asset block above) — no external
    # file requests, so no www/ directory is needed on the server. The
    # favicon.svg references inside styles.css are swapped for the data URI.
    tags$head(tags$style(HTML(
      gsub("url('favicon.svg')", paste0("url(\"", perception_favicon_uri(), "\")"),
           perception_asset("styles.css"), fixed = TRUE)))),
    tags$head(tags$style(HTML(perception_asset("tour.css")))),
    tags$head(tags$script(HTML(perception_asset("tour.js")))),
    tags$head(tags$link(rel = "icon", href = perception_favicon_uri(),
                        type = "image/svg+xml")),
    tags$head(tags$script(HTML("
      Shiny.addCustomMessageHandler('scroll-to', function(id) {
        var el = document.getElementById(id);
        if (el) el.scrollIntoView({behavior: 'smooth', block: 'start'});
      });

      // Live stage text inside the demo overlay (Seurat phases).
      Shiny.addCustomMessageHandler('set-demo-stage', function(msg) {
        var el = document.getElementById('demo-stage');
        if (el && msg && msg.text) el.textContent = msg.text;
      });

      // Replace a container's innerHTML. Used by the Home page to update the
      // workflow stepper / data dashboard, and by the Data page to swap the
      // Load-Demo / Clear-Demo button after the static first paint. innerHTML
      // drops Shiny's event bindings on replaced nodes, so re-bind inputs and
      // widgets inside the container afterwards.
      //
      // unbindAll() BEFORE replacing is NOT optional: Shiny (>= 1.14) keeps a
      // GLOBAL per-id binding counter, and innerHTML alone never decrements it.
      // Without the unbind, every replace+bindAll bumps the counter and the
      // console reports the 'Duplicate input IDs were found' warning for
      // these ids (e.g. workflow stepper links, the Load-Demo button) after
      // a few state updates.
      Shiny.addCustomMessageHandler('set-html', function(msg) {
        var el = document.getElementById(msg.id);
        if (el) {
          Shiny.unbindAll(el);
          el.innerHTML = msg.html;
          Shiny.bindAll(el);
        }
      });

      // Hide the pre-built demo overlay (shown instantly by the button's
      // onclick; the server only hides it once the demo task resolves).
      Shiny.addCustomMessageHandler('hide-demo-overlay', function(msg) {
        var el = document.getElementById('demo-overlay');
        if (el) el.style.display = 'none';
      });

      // Safety net for Shiny modals: modalDialog locks body scrolling via the
      // Bootstrap 'modal-open' class. When one modal REPLACES another (e.g. a
      // demo-confirm dialog followed by a training-complete dialog) or a modal
      // is dismissed in a non-standard way, the class can be left on <body>
      // and the page can no longer scroll at all. Whenever a modal hides,
      // restore body scrolling if no modal is still open.
      $(document).on('hidden.bs.modal', function () {
        if (document.querySelectorAll('.modal.show, .modal.in').length === 0) {
          document.body.classList.remove('modal-open');
          document.body.style.removeProperty('overflow');
          document.body.style.removeProperty('padding-right');
        }
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
      icon("circle-info", class = "nav-icon"),
      title = "Guided tour"
    )
  )),
  nav_item(tagList(
    tags$a(
      href = "https://wanglabcsu.github.io/PERCEPTIONx/articles/shiny_app.html",
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

  # Pre-warm the per-session task worker right away: the process spawn +
  # package load (~5 s) happens in the background while the user is still on
  # the landing page, so the FIRST plot/demo/predict does not pay it. If it
  # fails (no real risk), tasks re-spawn lazily as before.
  tryCatch(ensure_session_worker(shared), error = function(e) NULL)
  # Then kick off a warm-up task in that worker (fire-and-forget) so the
  # heavy Seurat / caret / glmnet namespaces are already loaded when the user
  # first runs Seurat clustering / the demo — that ~3-4 s load no longer sits
  # inside the first prepare/demo task.
  tryCatch(submit_session_task(shared, "warmup", list()), error = function(e) NULL)

  # Modules — pass main session for cross-module navigation
  mod_home_server("home", shared, session)
  mod_data_server("data", shared)
  mod_train_server("train", shared, session)
  mod_predict_server("predict", shared, session)
  mod_visualize_server("visualize", shared, session)
  mod_help_server("help", session, shared)

  # NOTE: no session-end file sweep here. Shiny itself deletes a session's
  # uploaded files when it ends. A global sweep over the cache/temp areas
  # would risk deleting files that other sessions' background workers still
  # need (the shared DepMap cache, an in-progress download, or a queued job),
  # breaking their training. Per-session job dirs are removed below; R's
  # tempdir() is cleaned up by the OS on process exit, and the DepMap cache
  # on disk is NEVER auto-expired (a pre-downloaded copy must be reused by
  # every user, not deleted after idle time).

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
    # NOTE: deliberately NO closeAllConnections() here. It was added to silence
    # "closing unused connection N" warnings in the server log, but it corrupts
    # the process's R connection table — any later connection use (a refreshed
    # session in the same process, or the shutdown sequence itself) then fails
    # with "invalid connection", e.g. callr::r_bg failing to start workers
    # ("Failed to start the background worker after 3 attempts"). Log warnings
    # are cosmetic; a broken session is not.
  })
}

shinyApp(ui, server)
