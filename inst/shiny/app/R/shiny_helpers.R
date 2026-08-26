# Shared helpers for the Shiny app — sourced BEFORE the modules in app.R.
# Anything used by more than one module lives here so the entry points cannot
# drift apart (e.g. the Home "Load Demo" button and the Data-tab demo button).

# ---------------------------------------------------------------------------
# DepMap metadata extraction
#
# The main Shiny process NEVER holds the full multi-GB DepMap object anymore
# (training runs in the background master). Everything the UI needs — gene
# list, drug list, lineages, component dimensions, a few column names — is a
# few hundred KB, extracted here. Built-in downloads cache the metadata to a
# sidecar file so later sessions read kilobytes instead of gigabytes.
# ---------------------------------------------------------------------------
extract_depmap_meta <- function(depmap_path, cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file) && file.exists(depmap_path)) {
    # Only trust a cache that is newer than the source DepMap.RDS itself:
    # if the RDS was replaced (re-download / different version) the cached
    # genes/drugs would be silently stale. A corrupt cache falls through to
    # a full re-extraction instead of erroring. (Both files must exist —
    # file.mtime() returns NA for a missing file and the comparison would
    # then produce if(NA) instead of a clean fall-through.)
    if (file.mtime(cache_file) >= file.mtime(depmap_path)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (!is.null(cached)) return(cached)
    }
  }
  DepMap <- readRDS(depmap_path)
  if (!is.list(DepMap)) {
    stop("The file is not a list. Expected a DepMap.RDS from the PERCEPTIONx release (or an object with the same structure).")
  }
  required_fields <- c("secondary_prism", "secondary_screen_drugAnnotation",
                       "expression_rnorm", "scRNA_complete", "scRNA_subset_rnorm",
                       "CPM_scRNA_CCLE_rnorm", "annotation_20Q4",
                       "metadata_CPM_scRNA", "expression_20Q4")
  missing_fields <- setdiff(required_fields, names(DepMap))
  if (length(missing_fields) > 0) {
    stop("The RDS is missing required components: ",
         paste(missing_fields, collapse = ", "),
         ". Re-save it with these components, or use the built-in 'Download & Load' (known-good).")
  }
  nms <- names(DepMap)
  components <- lapply(nms, function(nm) {
    obj <- DepMap[[nm]]
    if (is.matrix(obj) || is.data.frame(obj)) {
      list(nrow = nrow(obj), ncol = ncol(obj),
           cols_preview = head(colnames(obj), 12))
    } else {
      list(nrow = NA_integer_, ncol = NA_integer_, cols_preview = character(0))
    }
  })
  names(components) <- nms
  lineages <- sort(unique(as.character(DepMap$annotation_20Q4$lineage)))
  lineages <- lineages[!is.na(lineages) & nzchar(lineages)]
  meta <- list(
    loaded = TRUE,
    genes  = rownames(DepMap$expression_rnorm),
    drugs  = unique(as.character(DepMap$secondary_screen_drugAnnotation$CommonName)),
    lineages = lineages,
    components = components
  )
  if (!is.null(cache_file)) saveRDS(meta, cache_file)
  rm(DepMap); gc()
  meta
}

# ---------------------------------------------------------------------------
# Workflow stepper HTML (Home page)
#
# Pure renderer: builds the pipeline stepper + hint from a plain list of
# shared state (reactiveValuesToList(shared) or the initial empty list).
# It is called BOTH at UI build time (static first paint — no server round
# trip, so the section shows instantly) and on the server whenever shared
# state changes (HTML is pushed to the browser via the 'set-html' handler).
#
# @param shared A list mirroring the app's shared reactiveValues.
# @param ns     Namespace function for the actionLink ids (identity when
#               no namespacing is needed).
# ---------------------------------------------------------------------------
workflow_stepper_html <- function(shared, ns = identity) {
  data_ok    <- !is.null(shared$prepared_data) ||
                (!is.null(shared$user_expr) && !is.null(shared$user_clones))
  train_ok   <- !is.null(shared$models) && length(shared$models) > 0
  predict_ok <- !is.null(shared$predictions)
  done <- c(data_ok, train_ok, predict_ok, FALSE)

  # "Active" = the first incomplete step (the user's next action).
  active_idx <- which(!done)[1]
  if (is.na(active_idx)) active_idx <- 4L

  # Red "blocked" state — used SPARINGLY (red reads as error). Only for
  # Train when the user has data but no DepMap reference and no models:
  # training is genuinely impossible until one of those is provided.
  depmap_ok <- !is.null(shared$depmap_meta)
  blocked2 <- data_ok && !train_ok && !depmap_ok

  step_defs <- list(
    list(id = "wf_1", title = "Load Data"),
    list(id = "wf_2", title = "Train"),
    list(id = "wf_3", title = "Predict"),
    list(id = "wf_4", title = "Visualize")
  )
  status_done <- c("Data ready", "Model ready", "Predictions ready", "Complete")
  status_next <- c("Start here", "Next step", "Next step", "Next step")
  status_todo <- c("Pending",  "Pending",  "Pending",  "Pending")

  step_html <- lapply(seq_along(step_defs), function(i) {
    st <- if (done[i]) "done"
          else if (i == 2L && blocked2) "blocked"
          else if (i == active_idx) "active"
          else "todo"
    label <- if (done[i]) status_done[i]
             else if (i == 2L && blocked2) "Needs DepMap"
             else if (i == active_idx) status_next[i]
             else status_todo[i]
    node <- if (done[i]) icon("check", style = "font-size: 1.05rem;")
            else if (i == 2L && blocked2) icon("triangle-exclamation", style = "font-size: 1rem;")
            else tags$span(i)
    tags$div(class = paste("wf-step", st),
      actionLink(ns(step_defs[[i]]$id), NULL, class = "wf-step-link",
        div(class = "wf-node", node),
        div(class = "wf-step-title", step_defs[[i]]$title),
        div(class = "wf-step-status", label)
      )
    )
  })

  # Connector lines between steps — filled once the NEXT step is done.
  lines <- lapply(seq_len(length(step_defs) - 1), function(i)
    div(class = paste("wf-line", if (done[i + 1]) "filled" else "")))

  children <- list()
  for (i in seq_along(step_defs)) {
    children[[length(children) + 1]] <- step_html[[i]]
    if (i < length(step_defs)) children[[length(children) + 1]] <- lines[[i]]
  }

  hint <- if (blocked2) {
    "Data is ready, but training needs the DepMap reference — load it in the Data tab (or load pre-trained models) first."
  } else switch(active_idx,
    "1" = "Start: load data — or click Load Demo to try the whole pipeline instantly.",
    "2" = "Data ready -> next: train (or load) a drug response model.",
    "3" = "Model ready -> next: run clone- and patient-level prediction.",
    "4" = "Predictions ready -> go to Visualize to generate plots.",
    "All steps complete -> everything is ready. Explore the tabs!"
  )

  tagList(
    div(class = "wf-steps", children),
    div(class = if (blocked2) "wf-hint wf-hint-warn" else "wf-hint",
        icon(if (blocked2) "triangle-exclamation" else "circle-info"), hint)
  )
}

# ---------------------------------------------------------------------------
# Data status dashboard HTML (Home page)
#
# Same pattern as workflow_stepper_html(): pure renderer used both for the
# static first paint and for subsequent server-side updates.
# ---------------------------------------------------------------------------
data_dashboard_html <- function(shared) {
  # DepMap stats come from the lightweight metadata (the main process never
  # holds the full multi-GB object; the training master does).
  dm <- shared$depmap_meta
  depmap_cell_lines <- if (!is.null(dm) && !is.null(dm$components$expression_rnorm))
    dm$components$expression_rnorm$ncol else NULL
  depmap_drugs <- if (!is.null(dm) && !is.null(dm$components$secondary_prism))
    dm$components$secondary_prism$nrow else NULL
  depmap_sc_models <- if (!is.null(dm) && !is.null(dm$components$scRNA_complete))
    dm$components$scRNA_complete$ncol else NULL

  items <- list(
    list(name = "DepMap Reference", icon_name = "database",
         loaded = !is.null(dm),
         detail = if (!is.null(dm))
           paste0(if (!is.null(depmap_cell_lines)) paste0(depmap_cell_lines, " cell lines, "),
                  if (!is.null(depmap_drugs)) paste0(depmap_drugs, " drugs"),
                  if (!is.null(depmap_sc_models)) paste0(", ", depmap_sc_models, " scRNA-seq profiles"))
         else "Not loaded"),
    list(name = "Expression Matrix", icon_name = "table",
         loaded = !is.null(shared$user_expr),
         detail = if (!is.null(shared$user_expr))
           paste0(nrow(shared$user_expr), " genes x ", ncol(shared$user_expr), " clones")
         else "Not loaded"),
    list(name = "Clone Map (Seurat)", icon_name = "layer-group",
         loaded = !is.null(shared$prepared_data),
         detail = if (!is.null(shared$prepared_data))
           paste0(ncol(shared$prepared_data$clone_expression_rnorm), " clones, ",
                  nrow(shared$prepared_data$clone_counts), " patients")
         else "Not loaded"),
    list(name = "Clinical Response", icon_name = "heartbeat",
         loaded = !is.null(shared$user_response),
         detail = if (!is.null(shared$user_response))
           paste0(nrow(shared$user_response), " patients")
         else "Not loaded"),
    list(name = "Drug Models", icon_name = "cube",
         loaded = !is.null(shared$models),
         detail = if (!is.null(shared$models))
           paste0(length(shared$models), " model(s): ",
                  paste(names(shared$models), collapse = ", "))
         else "Not loaded"),
    list(name = "Predictions", icon_name = "microscope",
         loaded = !is.null(shared$predictions),
         detail = if (!is.null(shared$predictions))
           paste0(nrow(shared$predictions), " clones x ", ncol(shared$predictions), " drugs")
         else "Not loaded")
  )

  tagList(
    div(class = "dashboard-grid",
      lapply(items, function(item) {
        div(class = paste("dashboard-card", if (item$loaded) "card-loaded" else "card-empty"),
          div(class = "dashboard-card-icon", icon(item$icon_name)),
          div(class = "dashboard-card-body",
            div(class = "dashboard-card-name", item$name),
            div(class = "dashboard-card-detail", item$detail)
          ),
          div(class = "dashboard-card-status",
            span(class = paste("status-dot", if (item$loaded) "green" else "gray"))
          )
        )
      })
    )
  )
}

# ---------------------------------------------------------------------------
# Demo data pipeline (single source of truth)
#
# Builds the demo expression matrix, runs Seurat clustering via prepare_data(),
# trains the two demo drug models, and populates `shared`. Both the Home page
# "Load Demo" button and the Data-tab "Load Demo Data" button call this.
#
# @param shared     The app's shared reactiveValues.
# @param on_success Optional function called after everything is loaded
#                   (e.g. mark the demo active, navigate to a tab).
#
# @return Invisibly, a list with the loaded counts (genes/cells/patients/clones).
# ---------------------------------------------------------------------------
run_demo_pipeline <- function(shared, session, on_success = NULL) {
  # The whole demo (structured data + Seurat clustering + 2 demo models) runs
  # in the per-session background worker (task "demo", see async_jobs.R), so
  # the main Shiny thread never blocks other users while it prepares.
  w <- Waiter$new(
    html = tagList(
      div(class = "spinner-ring"),
      h4("Preparing demo data..."),
      p(class = "text-muted", "Running Seurat clustering")
    ),
    color = "rgba(255,255,255,0.85)"
  )
  w$show()
  jobid <- tryCatch(
    submit_session_task(shared, "demo", list()),
    error = function(e) {
      w$hide()
      showNotification(paste("Demo failed to start:", e$message), type = "error", duration = 10)
      NULL
    }
  )
  if (is.null(jobid)) return(invisible(NULL))
  poll_task(shared, session, jobid,
    on_done = function(res) {
      w$hide()
      shared$user_response <- res$user_response
      shared$user_mapping  <- res$user_mapping
      shared$user_expr     <- res$user_expr       # stays CELL-level for spatial plots
      shared$prepared_data <- res$prepared_data
      shared$user_clones   <- res$user_clones
      shared$models        <- res$models
      shared$model_cache   <- res$models
      shared$model_active  <- setNames(rep(TRUE, length(res$models)), names(res$models))
      n_clones <- ncol(res$prepared_data$clone_expression_rnorm)
      if (!is.null(on_success)) on_success()
      showNotification(paste0("Demo data loaded: ", nrow(res$user_expr), " genes x ",
                              ncol(res$user_expr), " cells, ", nrow(res$user_response),
                              " patients, ", n_clones, " clones detected. Models ready!"),
                       type = "message", duration = 6)
    },
    on_error = function(msg) {
      w$hide()
      showNotification(paste("Demo data preparation error:", msg), type = "error", duration = 10)
    })
  invisible(NULL)
}
