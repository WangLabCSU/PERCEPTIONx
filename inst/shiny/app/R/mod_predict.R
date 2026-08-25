# Prediction Module

# Empty-state card annotations. Rendered STATICALLY in the UI: the results
# cards live in a hidden nav tab, so server-rendered uiOutput content arrives
# only after the tab is shown / an event fires (same trap as the Home page
# stepper). The server just clears — or restores — them via the 'set-html'
# handler once predictions exist.
clone_desc_html <- function() {
  div(class = "card-desc",
    icon("fire", style = "font-size: 2rem; opacity: 0.3; display: block; margin-bottom: 0.35rem;"),
    "Predicted viability of every clone per drug (higher = more resistant).",
    br(),
    "Select model & expression sources, then click ",
    strong("Run Prediction"), "."
  )
}

patient_desc_html <- function() {
  div(class = "card-desc",
    icon("user", style = "font-size: 2rem; opacity: 0.3; display: block; margin-bottom: 0.35rem;"),
    "Clone viabilities aggregated per patient (weighted_max).",
    br(),
    "Appears here after running prediction."
  )
}

mod_predict_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        div(class = "section-header",
          icon("microscope"),
          h4("Drug Response Prediction")
        ),
        div(class = "info-box",
          icon("info-circle"), " Predict drug sensitivity at clone level using trained models, ",
          "then aggregate to patient-level response using weighted averaging strategies. ",
          strong("Viability semantics:"), " higher viability = more resistant (less sensitive); ",
          "lower viability = more sensitive (stronger drug killing)."
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
      # Configuration
      column(4, class = "predict-config-col",
        div(class = "card",
          div(class = "card-header",
            icon("cog"), " Configuration"
          ),
          div(class = "card-body",
            # Model source
            div(class = "param-section",
              h6(class = "param-section-title", icon("cube"), " Model Source"),
              radioButtons(ns("model_source"), NULL,
                           choices = c("From Training/Demo" = "trained", "Upload .RDS" = "upload"),
                           selected = "trained", inline = TRUE),
              conditionalPanel(
                condition = paste0("input['", ns("model_source"), "'] == 'upload'"),
                fileInput(ns("model_upload"), "Upload Model (.RDS)",
                          accept = c(".rds", ".RDS"))
              )
            ),

            # Expression source
            div(class = "param-section",
              h6(class = "param-section-title", icon("table"), " Expression Data"),
              radioButtons(ns("expr_source"), NULL,
                           choices = c("From Data Tab" = "loaded", "Upload New" = "upload"),
                           selected = "loaded", inline = TRUE),
              conditionalPanel(
                condition = paste0("input['", ns("expr_source"), "'] == 'upload'"),
                fileInput(ns("expr_upload"), "Upload Expression",
                          accept = c(".csv", ".rds", ".RDS"))
              )
            ),

            # Aggregation mode
            div(class = "param-section",
              h6(class = "param-section-title", icon("layer-group"), " Patient Aggregation"),
              selectInput(ns("agg_mode"), "Aggregation Mode",
                          choices = c(
                            "Weighted Max (recommended)" = "weighted_max",
                            "Weighted Average" = "weighted_average",
                            "Min (Most Sensitive)" = "min",
                            "Max (Most Resistant)" = "max"
                          ),
                          selected = "weighted_max"),
              div(class = "info-box", style = "margin-top: 0.5rem; font-size: 0.78rem;",
                strong("weighted_max"), ": top N most resistant clones (recommended)",
                br(), strong("weighted_average"), ": weighted average across all clones",
                br(), strong("min/max"), ": most sensitive/resistant clone"
              )
            ),

            hr(),

            actionButton(ns("predict"), "Run Prediction",
                         class = "btn-primary", icon = icon("play"))
          )
        )
      ),

      # Results
      column(8, class = "predict-results-col",
        # Clone-level heatmap
        div(class = "card",
          div(class = "card-header",
            icon("fire"), " Clone-Level Predictions"
          ),
          div(class = "card-body",
            div(id = ns("clone_desc"), clone_desc_html()),
            uiOutput(ns("clone_heatmap_area"))
          )
        ),

        # Patient-level table
        div(class = "card", style = "margin-top: 1rem;",
          div(class = "card-header",
            icon("user"), " Patient-Level Predictions"
          ),
          div(class = "card-body",
            div(id = ns("patient_desc"), patient_desc_html()),
            uiOutput(ns("patient_table_area")),
            uiOutput(ns("patient_download"))
          )
        )
      )
    )
  )
}

mod_predict_server <- function(id, shared, main_session) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Prerequisites check
    output$prereq_check <- renderUI({
      missing <- c()
      if (is.null(shared$models)) missing <- c(missing, "Trained Model")
      if (is.null(shared$prepared_data) && is.null(shared$user_expr)) {
        missing <- c(missing, "Expression Data")
      }
      if (is.null(shared$prepared_data) && is.null(shared$user_clones)) {
        missing <- c(missing, "Clone Annotation (run Seurat clustering in Data tab)")
      }
      if (length(missing) > 0) {
        div(class = "info-box", style = "border-left-color: var(--accent);",
          icon("exclamation-triangle"), " Missing prerequisites: ",
          strong(paste(missing, collapse = ", ")),
          ". Please load data first (Data tab) or use 'Load Demo Data'.",
          actionButton(ns("go_data"), "Go to Data", class = "btn-accent btn-sm", style = "margin-left: 0.5rem;")
        )
      } else {
        div(class = "info-box", style = "border-left-color: var(--success);",
          icon("check-circle"), " All prerequisites loaded. Ready to predict!"
        )
      }
    })

    observeEvent(input$go_data, {
      bslib::nav_select("navbar", selected = "data", session = main_session)
    })

    # Get model
    current_model <- reactive({
      if (input$model_source == "trained") {
        shared$models
      } else {
        if (!is.null(input$model_upload)) {
          readRDS(input$model_upload$datapath)
        } else {
          NULL
        }
      }
    })

    # Get expression — require prepare_data() output (clone-level, rank-normalized).
    # Never fall back to the raw cell-level user_expr: it is not rank-normalized
    # and would produce unreliable predictions.
    current_expr <- reactive({
      if (input$expr_source == "loaded") {
        if (!is.null(shared$prepared_data)) {
          shared$prepared_data$clone_expression_rnorm
        } else {
          NULL
        }
      } else {
        if (!is.null(input$expr_upload)) {
          file <- input$expr_upload
          if (grepl("\\.rds$|\\.RDS$|\\.Rds$", file$name)) {
            readRDS(file$datapath)
          } else {
            as.matrix(read.csv(file$datapath, row.names = 1, check.names = FALSE))
          }
        } else {
          NULL
        }
      }
    })

    # Run prediction — submitted to the background per-session worker. The
    # worker runs predict_drugs() + predict_patients() and ships the results
    # back via result.rds; the UI (clone_pred/patient_pred reactiveVals below)
    # updates when the job completes.
    clone_pred <- reactiveVal(NULL)
    patient_pred <- reactiveVal(NULL)

    observeEvent(input$predict, {
      model <- current_model()
      expr <- current_expr()

      if (is.null(model)) {
        showNotification("No model available. Train or upload a model first.", type = "error")
        return()
      }
      if (is.null(expr)) {
        showNotification("No expression data available. Load or upload data first.", type = "error")
        return()
      }

      # Pre-check: gene overlap between model features and expression data
      # (check ALL models, not just the first one)
      extract_features <- function(m) {
        if (inherits(m, "train")) {
          feats <- setdiff(colnames(m$trainingData), ".outcome")
          if (length(feats) > 0) return(feats)
          feats <- colnames(m$trainingData)
          if (length(feats) > 1) return(feats[-1])
        }
        if (!is.null(m$preProcess) && !is.null(m$preProcess$mean)) {
          return(names(m$preProcess$mean))
        }
        if (!is.null(m$finalModel) && !is.null(m$finalModel$beta)) {
          return(rownames(m$finalModel$beta))
        }
        character(0)
      }
      # Collect features from all models
      all_features <- unique(unlist(lapply(names(model), function(drug) {
        extract_features(model[[drug]]$model)
      })))
      if (length(all_features) > 0) {
        overlap <- intersect(all_features, rownames(expr))
        # Try make.names() fallback (hyphens vs dots)
        if (length(overlap) == 0) {
          expr_rownames_made <- make.names(rownames(expr), unique = TRUE)
          features_made <- make.names(all_features, unique = TRUE)
          overlap <- intersect(features_made, expr_rownames_made)
        }
        if (length(overlap) == 0) {
          showNotification(
            paste0("Gene name mismatch: none of the model's ", length(all_features),
                   " features found in your expression data. ",
                   "Check gene name formats (hyphens vs dots) between training and prediction data."),
            type = "error", duration = 12)
          return()
        }
        if (length(overlap) < length(all_features) * 0.5) {
          showNotification(
            paste0("Note: ", length(overlap), "/", length(all_features),
                   " model features found in data. Predictions may be less reliable."),
            type = "warning", duration = 8)
        }
      }

      w <- Waiter$new(
        html = tagList(div(class = "spinner-ring"), h4("Predicting..."),
                       p(class = "text-muted", "Running predict_drugs() in background worker")),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()

      # Build the (small) patient-level input. Standard path: pass only the
      # template + counts (not the whole prepared_data). Legacy path: pass
      # user_clones and let the worker aggregate.
      patient_input <- NULL
      legacy_clones <- NULL
      if (!is.null(shared$prepared_data)) {
        patient_input <- list(
          clone_viability_template = shared$prepared_data$clone_viability_template,
          clone_counts = shared$prepared_data$clone_counts
        )
      } else {
        legacy_clones <- shared$user_clones
      }

      jobid <- tryCatch(
        submit_session_task(shared, "predict", list(
          model_list    = model,
          expr          = expr,
          patient_input = patient_input,
          legacy_clones = legacy_clones,
          mode          = input$agg_mode
        )),
        error = function(e) { w$hide(); NULL }
      )
      if (is.null(jobid)) return()

      poll_task(shared, session, jobid,
        on_done = function(res) {
          w$hide()
          clone_pred(res$clone_pred)
          patient_pred(res$patient_pred)
          shared$predictions <- res$clone_pred
          shared$patient_pred <- res$patient_pred
          showNotification("Prediction complete (background worker)", type = "message")
        },
        on_error = function(msg) {
          w$hide()
          showNotification(paste("Prediction error:", msg), type = "error", duration = 10)
        })
    })

    # Patient-level aggregation now happens inside the background worker
    # (task "predict", see async_jobs.R); patient_pred() is a reactiveVal set
    # when the job completes.

    # Clone heatmap
    output$clone_heatmap <- renderPlotly({
      req(clone_pred())
      mat <- clone_pred()

      plotly::plot_ly(
        x = colnames(mat),
        y = rownames(mat),
        z = mat,
        type = "heatmap",
        colorscale = list(
          c(0, 1),
          c("rgb(58, 79, 138)", "rgb(232, 145, 58)")
        ),
        colorbar = list(title = "Viability"),
        hovertemplate = "Clone: %{y}<br>Drug: %{x}<br>Viability: %{z:.3f}<extra></extra>"
      ) %>%
        layout(
          xaxis = list(title = "Drug", tickangle = 45),
          yaxis = list(title = "Clone"),
          font = list(family = "Inter, sans-serif", size = 11),
          margin = list(l = 120, b = 80)
        )
    })

    # Heatmap area — only rendered after a prediction run, so no empty
    # heatmap frame shows before the user clicks Run Prediction. Height grows
    # with the number of clones so many rows are not squashed.
    output$clone_heatmap_area <- renderUI({
      req(clone_pred())
      mat <- clone_pred()
      n_rows <- nrow(mat)
      h <- if (is.null(n_rows) || n_rows <= 0) 400
           else min(700, max(240, 140 + 16 * n_rows))
      plotlyOutput(ns("clone_heatmap"), height = paste0(h, "px"))
    })

    # Card annotations (icon + description + hint) are statically pre-rendered
    # in the UI (visible on first paint). Clear them once predictions exist,
    # and restore them if a later run produces no result again.
    observe({
      has_res <- !is.null(clone_pred())
      session$sendCustomMessage("set-html",
        list(id = ns("clone_desc"),
             html = if (has_res) "" else as.character(clone_desc_html())))
      session$sendCustomMessage("set-html",
        list(id = ns("patient_desc"),
             html = if (has_res) "" else as.character(patient_desc_html())))
    })

    # Patient table
    output$patient_table <- renderDT({
      req(patient_pred())
      df <- as.data.frame(patient_pred())
      df$patient <- rownames(df)
      df <- df[, c("patient", setdiff(names(df), "patient"))]

      datatable(df,
                options = list(
                  pageLength = 15,
                  dom = "ftip",
                  columnDefs = list(list(className = "dt-center", targets = "_all"))
                ),
                rownames = FALSE,
                class = "display") %>%
        formatRound(names(df)[sapply(df, is.numeric)], 4)
    })

    # Patient table area — only rendered once patient predictions exist.
    output$patient_table_area <- renderUI({
      req(patient_pred())
      DTOutput(ns("patient_table"))
    })

    # Patient download
    output$patient_download <- renderUI({
      req(patient_pred())
      div(style = "margin-top: 0.75rem;",
        downloadButton(ns("download_pred"), "Download Predictions (.csv)",
                       class = "btn-outline-primary btn-sm")
      )
    })

    output$download_pred <- downloadHandler(
      filename = function() {
        paste0("predictions_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        df <- as.data.frame(patient_pred())
        df$patient <- rownames(df)
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}