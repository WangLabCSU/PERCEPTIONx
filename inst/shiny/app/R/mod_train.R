# Training Module
mod_train_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        div(class = "section-header",
          icon("brain"),
          h4("Model Training")
        ),
        div(class = "info-box",
          icon("info-circle"), " Train drug response models on DepMap bulk expression data. ",
          "Models are built using elastic net (glmnet) or random forest with cross-validated hyperparameter tuning. ",
          "The trained model is then refined on single-cell expression profiles."
        )
      )
    ),

    # Prerequisites check
    fluidRow(style = "margin-top: 1rem;",
      column(12,
        uiOutput(ns("prereq_check"))
      )
    ),

    fluidRow(style = "margin-top: 0.5rem;",
      # Parameters Panel — align-self:flex-start keeps it at natural height
      # when the results column grows after training (see styles.css).
      column(5, class = "train-config-col",
        div(class = "card train-params-card animate-fade-in-up",
          div(class = "card-header",
            icon("sliders-h"), " Parameters"
          ),
          div(class = "card-body",

            # Section: Drug & Cancer Type
            div(class = "param-section",
              h6(class = "param-section-title", icon("tag"), " Drug & Cancer Type"),
              textAreaInput(ns("drug"), "Drug Names (one per line or comma-separated)",
                            value = "abemaciclib", rows = 3, height = "90px",
                            width = "100%",
                            placeholder = "e.g.\nabemaciclib\nerlotinib\nAny drug with response data in your DepMap upload can be trained."),
              uiOutput(ns("drug_hint")),
              selectizeInput(ns("cancer_type"), "Cancer Type (include)",
                             choices = NULL, selected = "PanCan", width = "100%",
                             options = list(maxItems = 1, placeholder = "Select cancer type")),
              selectizeInput(ns("exclude_cancer"), "Cancer Type (exclude)",
                             choices = NULL, selected = "PanCan", width = "100%",
                             options = list(maxItems = 1, placeholder = "Select cancer type to exclude"))
            ),

            # Section: Genes of Interest
            div(class = "param-section",
              h6(class = "param-section-title", icon("dna"), " Genes of Interest"),
              div(class = "info-box", style = "margin-bottom: 0.6rem; font-size: 0.78rem; padding: 0.6rem 0.8rem;",
                icon("info-circle"),
                " Leave empty to use ", strong("all genes"), " from DepMap data (recommended for most users). ",
                "Or enter specific gene symbols to restrict feature selection to a subset."
              ),
              textAreaInput(ns("goi"), "Gene Symbols (optional — empty = all genes)",
                            placeholder = "Enter gene symbols, one per line or comma-separated.\nLeave empty to use all genes.",
                            rows = 8,
                            height = "120px",
                            width = "100%"),
              fileInput(ns("goi_file"), "Or Upload Gene List (.txt, .csv)",
                        accept = c(".txt", ".csv")),
              numericInput(ns("k_features"), "Top k Features",
                           value = 100, min = 10, max = 5000, step = 10,
                           width = "100%")
            ),

            # Section: Model Configuration
            div(class = "param-section",
              h6(class = "param-section-title", icon("cog"), " Model Configuration"),
              selectInput(ns("model_type"), "Algorithm",
                          choices = c("Elastic Net (glmnet)" = "glmnet", "Random Forest" = "rf"),
                          selected = "glmnet", width = "100%"),
              selectInput(ns("ncores"), "CPU Cores",
                          choices = 1:max(1L, min(4L, parallel::detectCores(), na.rm = TRUE)),
                          selected = 1, width = "100%"),
              tags$small(class = "text-muted",
                "Maximum 4 cores — hard cap for shared-server safety.")
            ),

            div(class = "train-action-row",
              actionButton(ns("train"), "Start Training",
                           class = "btn-primary", icon = icon("play")),
              actionButton(ns("reset"), "Reset",
                           class = "btn-reset btn-sm", icon = icon("undo"))
            )
          )
        )
      ),

      # Results Panel
      column(7, class = "train-results-col",
        div(class = "card",
          div(class = "card-header",
            icon("chart-bar"), " Training Results"
          ),
          div(class = "card-body",
            # Progress
            uiOutput(ns("progress")),

            # Model download — at the top so it never sits next to / competes
            # with the plot image downloads inside the tabs below.
            uiOutput(ns("download_btn")),

            # Model summary
            uiOutput(ns("model_summary")),

            # Performance metrics — p-value & correlation
            uiOutput(ns("perf_metrics")),

            # Performance plots — validation ROC (default) + threshold curve.
            # Rendered dynamically only after training completes: keeps the
            # heavy plotly.js bundle OFF the initial page (so entering the
            # module is instant) and shows no empty frames before a run.
            uiOutput(ns("perf_plots"))
          )
        )
      )
    )
  )
}

mod_train_server <- function(id, shared, main_session) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    trained <- reactiveVal(NULL)

    # Prerequisites check
    output$prereq_check <- renderUI({
      missing <- c()
      if (is.null(shared$depmap_meta)) missing <- c(missing, "DepMap Data")
      if (length(missing) > 0) {
        div(class = "info-box", style = "border-left-color: var(--accent);",
          icon("exclamation-triangle"), " Missing prerequisites: ",
          strong(paste(missing, collapse = ", ")),
          ". Please go to the Data tab first to load the required data.",
          actionButton(ns("go_data"), "Go to Data", class = "btn-accent btn-sm", style = "margin-left: 0.5rem;")
        )
      } else {
        div(class = "info-box", style = "border-left-color: var(--success);",
          icon("check-circle"), " All prerequisites loaded. Ready to train!"
        )
      }
    })

    observeEvent(input$go_data, {
      bslib::nav_select("navbar", selected = "data", session = main_session)
    })

    # Populate drug & cancer type choices when DepMap metadata is loaded
    observe({
      meta <- shared$depmap_meta
      if (!is.null(meta)) {
        # Drug hint: how many drugs have response data in the loaded DepMap.
        # Any of them can be trained — the 44 FDA-approved list is only a
        # recommended subset, not a limitation.
        available <- meta$drugs
        output$drug_hint <- renderUI({
          if (length(available) == 0) return(NULL)
          tags$small(class = "text-muted",
            paste0(length(available), " drugs with response data in the loaded DepMap — ",
                   "any of them can be trained, beyond the 44 recommended ones. ",
                   "Matching is case-insensitive.")
          )
        })

        # Cancer type choices: PanCan + unique lineages from DepMap annotation
        cancer_choices <- c("PanCan", meta$lineages)
        updateSelectizeInput(session, "cancer_type", choices = cancer_choices, selected = "PanCan", server = TRUE)
        updateSelectizeInput(session, "exclude_cancer", choices = cancer_choices, selected = "PanCan", server = TRUE)
      }
    })

    # Parse GOI input
    goi_parsed <- reactive({
      genes <- NULL
      if (!is.null(input$goi_file)) {
        file <- input$goi_file
        if (grepl("\\.csv$", file$name)) {
          df <- read.csv(file$datapath, stringsAsFactors = FALSE)
          genes <- as.character(df[[1]])
        } else {
          genes <- readLines(file$datapath)
        }
      } else if (nchar(trimws(input$goi)) > 0) {
        # NB: use explicit separator chars — "\s" inside a character class is
        # not treated as whitespace by R's default regex engine.
        genes <- unlist(strsplit(input$goi, "[,\r\n\t ]+"))
        genes <- trimws(genes)
        genes <- genes[nchar(genes) > 0]
      }
      genes
    })

    # Parse drug text input (one per line or comma/space separated)
    drug_parsed <- reactive({
      if (is.null(input$drug)) return(character(0))
      d <- unlist(strsplit(input$drug, "[,\r\n\t ]+"))
      d <- trimws(d)
      d <- d[nchar(d) > 0]
      unique(d)
    })

    # Waiter overlay + in-place progress updater. Stored in a reactiveVal so
    # BOTH the job-submission handler and the async polling observer can show,
    # update, and hide the same overlay.
    train_waiter <- reactiveVal(NULL)
    set_train_overlay <- function(stage, detail, value) {
      pct <- round(100 * max(0, min(1, value)))
      session$sendCustomMessage("set-html", list(id = ns("trn_stage"), html = stage))
      session$sendCustomMessage("set-html", list(id = ns("trn_detail"), html = detail))
      session$sendCustomMessage("set-html",
        list(id = ns("trn_bar"),
             html = paste0("<div class='train-progress-bar' style='width: ", pct, "%;'></div>")))
    }

    # Train
    observeEvent(input$train, {
      # Validate
      if (is.null(shared$depmap_meta)) {
        showNotification("Please load DepMap data first (Data tab)", type = "error")
        return()
      }
      # GOI is optional — if empty, use all genes from DepMap
      goi <- goi_parsed()
      if (is.null(goi) || length(goi) == 0) {
        goi <- shared$depmap_meta$genes
        showNotification(paste0("Using all ", length(goi), " genes from DepMap data"), type = "message")
      }
      drugs <- drug_parsed()
      if (length(drugs) == 0) {
        showNotification("Please enter at least one drug name", type = "error")
        return()
      }
      # Filter to drugs that actually have response data in the loaded DepMap.
      # Matching is case- and punctuation-insensitive (same rule as
      # train_models()'s stripall2match): "carfilzomib" matches "Carfilzomib".
      if (!is.null(shared$depmap_meta$drugs)) {
        available <- shared$depmap_meta$drugs
        norm_key <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(x)))
        matched <- match(norm_key(drugs), norm_key(available))
        missing <- drugs[is.na(matched)]
        if (length(missing) > 0) {
          showNotification(paste0("Not found in DepMap response data (skipped): ",
                                  paste(missing, collapse = ", "),
                                  ". Drugs outside DepMap (e.g. daratumumab) cannot be trained here — ",
                                  "load a pre-trained model instead (Predict tab)."),
                           type = "warning", duration = 8)
          drugs <- available[matched[!is.na(matched)]]
        }
        if (length(drugs) == 0) {
          showNotification("No valid drug names left — check spellings", type = "error")
          return()
        }
      }

      # ONE overlay layer: spinner + stage text + progress bar live together
      # inside the same waiter overlay. Stored in train_waiter() so the async
      # polling observer can keep updating it and finally hide it.
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Training model..."),
          p(id = ns("trn_stage"), class = "text-muted", "Submitting background job..."),
          p(id = ns("trn_detail"), class = "text-muted",
            style = "font-size: 0.82rem; opacity: 0.75;",
            "Elastic net tuning + single-cell refinement"),
          # NB: the id sits on the TRACK. set-html replaces the bar inside it
          # with a fresh .train-progress-bar whose width % is therefore relative
          # to the full track. (Putting the id on the bar itself caused a bar
          # nested inside a 5%-wide bar — the visible fill never moved.)
          div(id = ns("trn_bar"), class = "train-progress-track",
            div(class = "train-progress-bar", style = "width: 5%;")
          )
        ),
        color = "rgba(255,255,255,0.85)"
      )
      train_waiter(w)
      w$show()

      # Training runs in a SEPARATE R process (callr worker, see async_jobs.R):
      # the main Shiny thread is never blocked, so the UI / other users stay
      # responsive. We submit a job, then the polling observer below tracks it.
      tryCatch({
        params <- list(
          drug_list         = drugs,
          cancer_type       = input$cancer_type,
          exclude_cancer    = input$exclude_cancer,
          GOI               = goi,
          k_features_values = input$k_features,
          model_type        = input$model_type,
          ncores            = input$ncores,
          output_dir        = tempdir()
        )
        jobid <- submit_train_job(shared, params)
        shared$active_job <- jobid
        shared$active_job_drugs <- drugs
        showNotification(sprintf("Training job %s running in the background — you can keep using the app.",
                                 jobid),
                         type = "message", duration = 6)
      }, error = function(e) {
        w$hide()
        train_waiter(NULL)
        showNotification(paste("Failed to start training job:", e$message), type = "error", duration = 10)
      })
    })

    # Poll the background training job every second; drive the overlay and
    # finalize when the worker reports done / error. Also detects a DEAD
    # worker (master OOM-killed / custom worker crashed): without this the
    # job would spin "queued" forever because status.txt never appears.
    observe({
      jobid <- shared$active_job
      if (is.null(jobid)) return()
      worker_alive <- if (isTRUE(shared$depmap_is_standard)) {
        m <- getOption("perception.master")
        is.null(m) || m$is_alive()   # NULL = not yet spawned = still recoverable
      } else {
        w <- shared$train_worker
        is.null(w) || w$is_alive()
      }
      if (!worker_alive) {
        shared$active_job <- NULL
        shared$active_job_drugs <- NULL
        if (!is.null(train_waiter())) { train_waiter()$hide(); train_waiter(NULL) }
        showNotification("The background training worker stopped unexpectedly (it may have run out of memory or been killed). Please try again.",
                         type = "error", duration = 12)
        return()
      }
      st <- tryCatch(read_job_state(shared, jobid),
                     error = function(e) list(status = "error", message = conditionMessage(e)))
      w <- train_waiter()

      if (st$status == "done") {
        job_dir <- file.path(shared$jobs_dir, jobid)
        result <- tryCatch(readRDS(file.path(job_dir, "result.rds")), error = function(e) NULL)
        requested_drugs <- shared$active_job_drugs %||% character(0)
        shared$active_job <- NULL
        shared$active_job_drugs <- NULL
        # Result is in memory now — drop the on-disk job dir (models can be
        # large) so the shared jobs pool does not accumulate GBs over time.
        unlink(job_dir, recursive = TRUE)
        if (!is.null(w)) { w$hide(); train_waiter(NULL) }
        if (is.null(result)) {
          showNotification("Background job finished but produced no result (see console).",
                           type = "error", duration = 10)
          return()
        }
        trained(result)
        shared$models <- result
        if (is.null(shared$model_cache)) shared$model_cache <- list()
        for (nm in names(result)) {
          shared$model_cache[[nm]] <- result[[nm]]
          shared$model_active[[nm]] <- TRUE
        }
        failed_drugs <- setdiff(requested_drugs, names(result))
        if (length(failed_drugs) > 0) {
          showNotification(paste0("No model produced for: ",
                                  paste(failed_drugs, collapse = ", "),
                                  " (feature ranking or model building failed — see console warnings)."),
                           type = "warning", duration = 10)
        }
        showNotification("Training completed successfully!", type = "message")
        return()
      }
      if (st$status == "error") {
        shared$active_job <- NULL
        shared$active_job_drugs <- NULL
        unlink(file.path(shared$jobs_dir, jobid), recursive = TRUE)
        if (!is.null(w)) { w$hide(); train_waiter(NULL) }
        showNotification(paste("Training job failed:", st$message), type = "error", duration = 12)
        return()
      }
      # running / queued → keep the overlay in sync
      if (!is.null(w)) {
        if (st$status == "queued") {
          set_train_overlay("Waiting for background worker",
                            "Loading DepMap in the worker process (first run only)...", 0.02)
        } else if (identical(st$phase, "rank")) {
          set_train_overlay("Feature ranking", "Elastic net tuning + single-cell refinement", 0.05)
        } else if (identical(st$phase, "done")) {
          set_train_overlay("All drugs trained", "Saving model files...", 1)
        } else {
          n <- if (is.null(st$n) || !is.finite(st$n) || st$n <= 0) 1 else st$n
          i <- if (is.null(st$i) || is.na(st$i)) 0 else st$i
          set_train_overlay(sprintf("Training %d/%d: %s", i, n, st$drug),
                            "Elastic net tuning + single-cell refinement",
                            max(0.05, (i - 1) / n))
        }
      }
      invalidateLater(1000, session)
    })

    # Reset
    observeEvent(input$reset, {
      trained(NULL)
      updateTextAreaInput(session, "drug", value = "abemaciclib")
      updateSelectizeInput(session, "cancer_type", selected = "PanCan")
      updateSelectizeInput(session, "exclude_cancer", selected = "PanCan")
      updateTextAreaInput(session, "goi", value = "")
    })

    # Progress UI
    output$progress <- renderUI({
      if (is.null(trained())) {
        div(class = "text-muted", style = "text-align: center; padding: 3rem;",
          icon("brain", style = "font-size: 3rem; opacity: 0.15; display: block; margin-bottom: 0.8rem;"),
          "Enter drug names, configure parameters, and click ",
          strong("Start Training"),
          br(),
          span(style = "font-size: 0.82rem; opacity: 0.7;", "Check the prerequisites above for available data")
        )
      } else {
        div(class = "status-badge loaded",
          span(class = "status-dot green"),
          "Training complete"
        )
      }
    })

    # Model Summary
    output$model_summary <- renderUI({
      req(trained())
      models <- trained()
      drug_names <- names(models)

      summary_rows <- lapply(drug_names, function(d) {
        m <- models[[d]]
        # m is a build_on_BULK_v2 output list with $model (caret train object)
        if (!is.null(m$model) && inherits(m$model, "train")) {
          model_method <- m$model$method
          if (model_method == "glmnet") {
            bt <- m$model$bestTune
            info <- paste0("glmnet, alpha=", round(bt$alpha, 3),
                          ", lambda=", signif(bt$lambda, 4))
          } else if (model_method == "rf") {
            fm <- m$model$finalModel
            best_rmse <- min(m$model$results$RMSE, na.rm = TRUE)
            info <- paste0("rf, ntree=", fm$ntree,
                          ", RMSE=", signif(best_rmse, 4))
          } else {
            info <- paste("method:", model_method)
          }
        } else if (!is.null(m$model)) {
          info <- paste("model class:", paste(class(m$model), collapse = ", "))
        } else {
          info <- "No model object"
        }
        tags$tr(
          tags$td(strong(d)),
          tags$td(info)
        )
      })

      tagList(
        h6(style = "font-weight: 700; color: var(--primary-dark);", icon("info-circle"), " Model Summary"),
        tags$table(class = "table table-sm",
          tags$thead(tags$tr(tags$th("Drug"), tags$th("Parameters"))),
          tags$tbody(summary_rows)
        )
      )
    })

    # Validation ROC (default view after training). Rendered via plotly for
    # interactivity. The AUC card layers are plotly-safe now: all text uses
    # fixed colors and the color swatches are point markers (no fill-mapped
    # rects / colored text that the converter renders unreliably).
    # Keep the ROC/performance ggplot objects so the tab download buttons can
    # export crisp static images (export_plot_cairo) without needing orca.
    roc_gg <- reactiveVal(NULL)
    perf_gg <- reactiveVal(NULL)

    output$perf_roc_plot <- renderPlotly({
      req(trained())
      p <- PERCEPTIONx::plot_model_roc(trained(), base_size = 13)
      roc_gg(p)
      # plot_model_roc drops a drug only when its held-out response data
      # cannot be split into two classes at all — surface the reason instead
      # of leaving a silent gap.
      dropped <- attr(p, "dropped_drugs")
      if (length(dropped) > 0) {
        why <- attr(p, "dropped_reasons")
        detail <- if (length(why)) {
          paste(names(why), paste0("(", why, ")"), collapse = "; ")
        } else {
          paste(dropped, collapse = ", ")
        }
        showNotification(paste0("ROC not drawn for: ", detail,
                                " — held-out cell lines have too few or non-varied drug responses in your DepMap data."),
                         type = "warning", duration = 12)
      }
      ggplotly(p) %>%
        layout(font = list(family = "Inter, sans-serif", size = 12))
    })

    # Performance Plot
    output$perf_plot <- renderPlotly({
      req(trained())
      models <- trained()
      # tooltip = FALSE: return a plain ggplot (ggiraph geoms are not
      # convertible by ggplotly and would drop the point layers with warnings)
      p <- plot_model_performance(models, base_size = 13, tooltip = FALSE)
      perf_gg(p)
      ggplotly(p, tooltip = c("x", "y", "colour")) %>%
        layout(font = list(family = "Inter, sans-serif", size = 12))
    })

    output$download_roc_png <- downloadHandler(
      filename = function() paste0("validation_roc_", format(Sys.Date(), "%Y%m%d"), ".png"),
      content = function(file) {
        p <- roc_gg()
        if (is.null(p)) {
          showNotification("Generate the ROC plot first.", type = "warning")
          return()
        }
        PERCEPTIONx::export_plot_cairo(file, p, format = "png", width = 12, height = 10, res = 600)
      }
    )
    output$download_roc_pdf <- downloadHandler(
      filename = function() paste0("validation_roc_", format(Sys.Date(), "%Y%m%d"), ".pdf"),
      content = function(file) {
        p <- roc_gg()
        if (is.null(p)) {
          showNotification("Generate the ROC plot first.", type = "warning")
          return()
        }
        PERCEPTIONx::export_plot_cairo(file, p, format = "pdf", width = 12, height = 10)
      }
    )
    output$download_perf_png <- downloadHandler(
      filename = function() paste0("performance_curve_", format(Sys.Date(), "%Y%m%d"), ".png"),
      content = function(file) {
        p <- perf_gg()
        if (is.null(p)) {
          showNotification("Generate the performance plot first.", type = "warning")
          return()
        }
        PERCEPTIONx::export_plot_cairo(file, p, format = "png", width = 10, height = 7, res = 600)
      }
    )

    # Performance plots wrapper — only rendered after training completes so the
    # results card starts with just the hint and no empty plot frames. The
    # plotly.js bundle is therefore NOT part of the initial page (the module
    # opens instantly); it loads on demand when the first plot renders.
    output$perf_plots <- renderUI({
      req(trained())
      # With several drugs the ROC panels are stacked vertically (one per
      # drug; 2 columns once there are 5+), so grow the plot height per row
      # instead of squashing them into the single-drug 340px frame.
      n_models <- length(trained())
      roc_ncol <- if (n_models <= 4) 1 else 2
      roc_h <- if (n_models > 1) min(1500, 150 + 240 * ceiling(n_models / roc_ncol)) else 340
      div(class = "viz-plot-wrapper",
        tabsetPanel(
          tabPanel("Validation ROC",
            plotlyOutput(ns("perf_roc_plot"), height = paste0(roc_h, "px")),
            tags$small(class = "text-muted", style = "display: block; margin-top: 0.3rem;",
              "ROC of the predicted viability in stratifying the top vs bottom 50% observed response by rank (paper's convention), one curve per validation dataset (bulk / pseudo-bulk / single-cell) with AUC annotated in a box. 0.5 = random."),
            div(style = "margin-top: 0.4rem;",
              downloadButton(ns("download_roc_png"), "PNG", class = "btn-outline-primary btn-sm"),
              downloadButton(ns("download_roc_pdf"), "PDF", class = "btn-outline-primary btn-sm")
            )
          ),
          tabPanel("Performance Curve",
            plotlyOutput(ns("perf_plot"), height = "320px"),
            tags$small(class = "text-muted", style = "display: block; margin-top: 0.3rem;",
              "For each threshold, the proportion of trained drugs whose predicted-observed Pearson correlation exceeds it (per dataset). Most informative when several drugs are trained together; with a single drug the curve is a simple step."),
            div(style = "margin-top: 0.4rem;",
              downloadButton(ns("download_perf_png"), "PNG", class = "btn-outline-primary btn-sm")
            )
          )
        )
      )
    })

    # Performance Metrics — display p-value & correlation
    output$perf_metrics <- renderUI({
      req(trained())
      models <- trained()
      drug_names <- names(models)

      rows <- lapply(drug_names, function(d) {
        m <- models[[d]]
        # Helper to safely extract cor & p-value from named vector or data.frame
        extract_metrics <- function(perf) {
          if (is.null(perf)) return(list(cor = NA, p = NA))
          if (is.data.frame(perf)) {
            list(cor = perf$estimate.cor[1], p = perf$p.value[1])
          } else if (is.numeric(perf)) {
            # Named vector from unlist(cor.test(...))
            cor_val <- if ("estimate.cor" %in% names(perf)) perf["estimate.cor"] else
                       if ("cor" %in% names(perf)) perf["cor"] else NA
            p_val <- if ("p.value" %in% names(perf)) perf["p.value"] else NA
            list(cor = unname(cor_val), p = unname(p_val))
          } else {
            list(cor = NA, p = NA)
          }
        }

        bulk <- extract_metrics(m$performance_in_bulk)
        pseudo <- extract_metrics(m$performance_in_pseudo_bulk)
        sc <- extract_metrics(m$performance_in_scRNA)

        fmt <- function(v, digits = 3) {
          if (is.na(v)) return("—")
          format(round(v, digits), nsmall = digits)
        }

        # p-value: always show the number (scientific notation when very
        # small), never an abbreviated "< 0.001".
        fmt_p <- function(v) {
          if (is.na(v)) return("—")
          if (v < 0.001) format(v, scientific = TRUE, digits = 2)
          else format(round(v, 3), nsmall = 3)
        }

        tags$tr(
          tags$td(strong(d)),
          tags$td(fmt(bulk$cor), tags$span(class = "text-muted", style = "font-size:0.72rem; display:block;", "p = ", fmt_p(bulk$p))),
          tags$td(fmt(pseudo$cor), tags$span(class = "text-muted", style = "font-size:0.72rem; display:block;", "p = ", fmt_p(pseudo$p))),
          tags$td(fmt(sc$cor), tags$span(class = "text-muted", style = "font-size:0.72rem; display:block;", "p = ", fmt_p(sc$p)))
        )
      })

      tagList(
        h6(style = "font-weight: 700; color: var(--primary-dark); margin-top: 1rem;",
            icon("ruler"), " Performance Metrics (Pearson correlation)"),
        div(class = "table-responsive",
          tags$table(class = "table table-sm",
            tags$thead(tags$tr(
              tags$th("Drug"),
              tags$th("Bulk"),
              tags$th("Pseudo-bulk"),
              tags$th("Single-cell")
            )),
            tags$tbody(rows)
          )
        ),
        tags$small(class = "text-muted",
          "Correlation between predicted and observed drug response. Higher correlation + lower p-value = better model."
        )
      )
    })

    # Download
    output$download_btn <- renderUI({
      req(trained())
      div(style = "margin-top: 1rem;",
        downloadButton(ns("download_model"), "Download Model (.RDS)",
                       class = "btn-outline-primary btn-sm")
      )
    })

    output$download_model <- downloadHandler(
      filename = function() {
        drug_str <- paste(drug_parsed(), collapse = "_")
        if (nchar(drug_str) > 80) drug_str <- substr(drug_str, 1, 80)
        paste0(drug_str, "_model_", format(Sys.Date(), "%Y%m%d"), ".RDS")
      },
      content = function(file) {
        saveRDS(trained(), file)
      }
    )
  })
}
