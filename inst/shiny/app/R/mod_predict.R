# Prediction Module

# Empty-state card annotations. Rendered STATICALLY in the UI so the results
# cards show immediately on first paint (no empty flash waiting for a server
# round trip). The server just clears — or restores — them via the 'set-html'
# handler once predictions exist.
clone_desc_html <- function() {
  div(
    class = "card-desc",
    icon(
      "fire",
      style = "font-size: 2rem; opacity: 0.3; display: block; margin-bottom: 0.35rem;"
    ),
    "Predicted viability of every clone per drug (higher = more resistant).",
    br(),
    "Select model & expression sources, then click ",
    strong("Run Prediction"),
    "."
  )
}

patient_desc_html <- function() {
  div(
    class = "card-desc",
    icon(
      "user",
      style = "font-size: 2rem; opacity: 0.3; display: block; margin-bottom: 0.35rem;"
    ),
    "Clone viabilities aggregated per patient (weighted_max).",
    br(),
    "Appears here after running prediction."
  )
}

# Publication-style heatmap of clone x drug viability. Kept as a top-level
# function (not inline in a reactive) so it can be rendered and exported
# without a live Shiny session.
clone_heatmap_plot <- function(mat, base_size = 11) {
  mat <- as.matrix(mat)
  clones <- rownames(mat)
  drugs <- colnames(mat)
  long <- data.frame(
    clone = factor(rep(clones, times = length(drugs)), levels = rev(clones)),
    drug = factor(rep(drugs, each = length(clones)), levels = drugs),
    viability = as.vector(mat),
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$drug, y = .data$clone, fill = .data$viability)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      name = "Viability",
      limits = c(0, 1),
      oob = scales::squish
    ) +
    ggplot2::labs(
      title = "Predicted viability by clone and drug",
      subtitle = sprintf(
        "%d clone%s × %d drug%s · higher = more resistant",
        length(clones),
        if (length(clones) == 1L) "" else "s",
        length(drugs),
        if (length(drugs) == 1L) "" else "s"
      ),
      caption = "Values are model predictions, not measured IC50s.",
      x = "Drug",
      y = "Clone"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.4, colour = "grey35"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = ggplot2::rel(1.05)
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        colour = "grey30",
        size = ggplot2::rel(0.82)
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        colour = "grey40",
        size = ggplot2::rel(0.7)
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

mod_predict_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        12,
        div(
          class = "section-header",
          icon("microscope"),
          h4("Drug Response Prediction")
        ),
        div(
          class = "info-box",
          icon("info-circle"),
          " Predict drug sensitivity at clone level using trained models, ",
          "then aggregate to patient-level response using weighted averaging strategies. ",
          strong("Viability semantics:"),
          " higher viability = more resistant (less sensitive); ",
          "lower viability = more sensitive (stronger drug killing)."
        )
      )
    ),

    # Prerequisites check
    fluidRow(
      style = "margin-top: 1rem;",
      column(12, uiOutput(ns("prereq_check")))
    ),

    fluidRow(
      style = "margin-top: 0.5rem;",
      # Configuration
      column(
        4,
        class = "predict-config-col",
        div(
          class = "card",
          div(class = "card-header", icon("cog"), " Configuration"),
          div(
            class = "card-body",
            # Sources — single source of truth: models & expression come from
            # the Data tab (management-integrated), never uploaded here.
            div(
              class = "param-section",
              h6(class = "param-section-title", icon("link"), " Sources"),
              p(
                class = "param-note",
                "Prediction runs on the ",
                strong("active models"),
                " and the expression data prepared in the Data tab. Load, upload ",
                "or train models there — every uploaded .RDS appears in the ",
                "management list automatically."
              )
            ),

            # Aggregation mode
            div(
              class = "param-section",
              h6(
                class = "param-section-title",
                icon("layer-group"),
                " Patient Aggregation"
              ),
              selectInput(
                ns("agg_mode"),
                "Aggregation Mode",
                choices = c(
                  "Weighted Max (recommended)" = "weighted_max",
                  "Weighted Average" = "weighted_average",
                  "Min (Most Sensitive)" = "min",
                  "Max (Most Resistant)" = "max"
                ),
                selected = "weighted_max"
              ),
              p(
                class = "param-note",
                strong("weighted_max"),
                ": highest clone viability, weighted by clone proportion ",
                "(recommended). ",
                strong("weighted_average"),
                ": weighted mean across clones. ",
                strong("min / max"),
                ": the most sensitive / most resistant clone."
              )
            ),

            hr(),

            actionButton(
              ns("predict"),
              "Run Prediction",
              class = "btn-primary",
              icon = icon("play")
            )
          )
        )
      ),

      # Results
      column(
        8,
        class = "predict-results-col",
        # Clone-level heatmap
        div(
          class = "card",
          div(class = "card-header", icon("fire"), " Clone-Level Predictions"),
          div(
            class = "card-body",
            div(id = ns("clone_desc"), clone_desc_html()),
            uiOutput(ns("clone_heatmap_area"))
          )
        ),

        # Patient-level table
        div(
          class = "card",
          style = "margin-top: 1rem;",
          div(
            class = "card-header",
            icon("user"),
            " Patient-Level Predictions"
          ),
          div(
            class = "card-body",
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

    # Guards against double-submitting the predict task (rapid double clicks).
    pred_busy <- reactiveVal(FALSE)

    # Prerequisites check
    output$prereq_check <- renderUI({
      missing <- c()
      if (is.null(shared$models)) {
        missing <- c(missing, "Trained Model")
      }
      if (is.null(shared$prepared_data) && is.null(shared$user_expr)) {
        missing <- c(missing, "Expression Data")
      }
      if (is.null(shared$prepared_data) && is.null(shared$user_clones)) {
        missing <- c(
          missing,
          "Clone Annotation (run Seurat clustering in Data tab)"
        )
      }
      if (length(missing) > 0) {
        div(
          class = "info-box",
          style = "border-left-color: var(--accent);",
          icon("exclamation-triangle"),
          " Missing prerequisites: ",
          strong(paste(missing, collapse = ", ")),
          ". Please load data first (Data tab) or use 'Load Demo Data'.",
          actionButton(
            ns("go_data"),
            "Go to Data",
            class = "btn-accent btn-sm",
            style = "margin-left: 0.5rem;"
          )
        )
      } else {
        div(
          class = "info-box",
          style = "border-left-color: var(--success);",
          icon("check-circle"),
          " All prerequisites loaded. Ready to predict!"
        )
      }
    })

    observeEvent(input$go_data, {
      bslib::nav_select("navbar", selected = "data", session = main_session)
    })

    # Model source: always the models the user loaded/activated on the Data tab
    # (Pre-trained download, .RDS upload, or models trained on the Train tab).
    current_model <- reactive(shared$models)

    # Expression: always the clone-level, rank-normalized data prepared on the
    # Data tab. The raw cell-level user_expr is NOT rank-normalized and would
    # produce unreliable predictions.
    current_expr <- reactive({
      if (is.null(shared$prepared_data)) {
        NULL
      } else {
        shared$prepared_data$clone_expression_rnorm
      }
    })

    # Run prediction — submitted to the background per-session worker. The
    # worker runs predict_drugs() + predict_patients() and ships the results
    # back via result.rds; the UI (clone_pred/patient_pred reactiveVals below)
    # updates when the job completes.
    clone_pred <- reactiveVal(NULL)
    patient_pred <- reactiveVal(NULL)

    # The module-local results are mirrors of shared$predictions / $patient_pred
    # (set on every completed job). When the app clears them — Clear Demo,
    # uploading new data, deactivating models — the local copies must follow,
    # or the Predict tab keeps showing stale results that contradict the
    # global state (e.g. "Predictions Not loaded" in the sidebar).
    observe({
      if (is.null(shared$predictions) && !is.null(clone_pred())) {
        clone_pred(NULL)
      }
      if (is.null(shared$patient_pred) && !is.null(patient_pred())) {
        patient_pred(NULL)
      }
    })

    observeEvent(input$predict, {
      model <- current_model()
      expr <- current_expr()

      # Completely empty state: point the user somewhere actionable instead of
      # a generic "no model" / "no data" pair.
      if (is.null(model) && is.null(expr)) {
        showNotification(
          "No model or expression data loaded yet — use 'Load Demo Data' or upload your data in the Data tab first.",
          type = "warning",
          duration = 8
        )
        return()
      }
      if (is.null(model)) {
        showNotification(
          "No model available. Train or upload a model first.",
          type = "error",
          duration = 8
        )
        return()
      }
      if (is.null(expr)) {
        showNotification(
          "No expression data available. Load or upload data first.",
          type = "error",
          duration = 8
        )
        return()
      }
      if (pred_busy()) {
        showNotification(
          "A prediction is already running",
          type = "warning",
          duration = 5
        )
        return()
      }

      # Pre-check: gene overlap between model features and expression data
      # (check ALL models, not just the first one)
      # NOTE: this runs BEFORE pred_busy(TRUE) — a blocked prediction must not
      # leave the busy flag set, or the next click would wrongly report
      # "A prediction is already running".
      extract_features <- function(m) {
        if (inherits(m, "train")) {
          feats <- setdiff(colnames(m$trainingData), ".outcome")
          if (length(feats) > 0) {
            return(feats)
          }
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
          showModal(modalDialog(
            title = "Data incompatibility — prediction blocked",
            tags$div(
              div(
                style = "text-align: center; margin-bottom: 0.6rem;",
                icon(
                  "exclamation-triangle",
                  style = "color: #e08a3c; font-size: 1.8rem;"
                )
              ),
              p(paste0(
                "None of the model's ",
                length(all_features),
                " features were found in the expression data (",
                nrow(expr),
                " genes available)."
              )),
              p(
                "The model was trained on the DepMap gene space; the current ",
                "expression data shares no genes with it. Check gene name formats ",
                "(hyphens vs dots, e.g. MIF-1 vs MIF.1) between training and ",
                "prediction data."
              )
            ),
            footer = modalButton("OK"),
            easyClose = TRUE
          ))
          return()
        }
        if (length(overlap) < length(all_features) * 0.5) {
          showNotification(
            paste0(
              "Note: ",
              length(overlap),
              "/",
              length(all_features),
              " model features found in data. Predictions may be less reliable."
            ),
            type = "warning",
            duration = 8
          )
        }
      }

      # All validation passed — only now claim the busy flag so a blocked run
      # never leaves it set.
      pred_busy(TRUE)

      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Predicting..."),
          p(class = "text-muted", "Running prediction...")
        ),
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
        submit_session_task(
          shared,
          "predict",
          list(
            model_list = model,
            expr = expr,
            patient_input = patient_input,
            legacy_clones = legacy_clones,
            mode = input$agg_mode
          )
        ),
        error = function(e) {
          pred_busy(FALSE)
          w$hide()
          NULL
        }
      )
      if (is.null(jobid)) {
        return()
      }

      poll_task(
        shared,
        session,
        jobid,
        on_done = function(res) {
          pred_busy(FALSE)
          w$hide()
          clone_pred(res$clone_pred)
          patient_pred(res$patient_pred)
          shared$predictions <- res$clone_pred
          shared$patient_pred <- res$patient_pred
          showNotification("Prediction complete", type = "message")
        },
        on_error = function(msg) {
          pred_busy(FALSE)
          w$hide()
          showNotification(
            paste("Prediction error:", msg),
            type = "error",
            duration = 10
          )
        }
      )
    })

    # Patient-level aggregation now happens inside the background worker
    # (task "predict", see async_jobs.R); patient_pred() is a reactiveVal set
    # when the job completes.

    # Clone heatmap — a static ggplot rather than plotly. This is the figure
    # that ends up in a paper, so it gets the same treatment as the QC report
    # (title, subtitle, caption, no web chrome) and a 300 dpi PNG export.
    clone_heatmap_height <- function() {
      mat <- clone_pred()
      if (is.null(mat)) {
        return(420L)
      }
      as.integer(min(780L, max(320L, 190L + 20L * nrow(mat))))
    }

    clone_heatmap_obj <- shiny::reactive({
      mat <- clone_pred()
      shiny::req(mat)
      clone_heatmap_plot(mat)
    })

    output$clone_heatmap <- renderPlot(
      clone_heatmap_obj(),
      height = function() clone_heatmap_height(),
      res = 110
    )

    # Heatmap area — only rendered after a prediction run, so no empty
    # heatmap frame shows before the user clicks Run Prediction. Height grows
    # with the number of clones so many rows are not squashed.
    # NOTE: return an EXPLICIT empty div when there are no results. A bare
    # req() failure inside renderUI leaves the previous output in place
    # (Shiny keeps stale content on silent errors), which is exactly how the
    # old heatmap survived "Clear Demo Data".
    output$clone_heatmap_area <- renderUI({
      mat <- clone_pred()
      if (is.null(mat)) {
        return(div())
      }
      # Demo models are SIMULATED (trained on random noise, not real drug
      # response data) but carry real drug names — the results must carry a
      # visible banner whenever any of them appear in the prediction.
      demo <- intersect(colnames(mat), demo_models(shared$models))
      banner <- if (length(demo) > 0) {
        div(
          class = "info-box",
          style = "border-left-color: var(--warning); margin-bottom: 0.6rem;",
          icon("exclamation-triangle"),
          strong("Demo simulation: "),
          "Demo models are SIMULATED placeholders, not trained on real drug response data. ",
          "These predictions demonstrate the workflow only and must NOT be interpreted as real drug response."
        )
      } else {
        NULL
      }
      tagList(
        banner,
        div(
          class = "figure-toolbar",
          downloadButton(
            ns("clone_heatmap_download"),
            "Download PNG",
            class = "btn-outline-primary qc-download-btn",
            icon = icon("download")
          )
        ),
        plotOutput(
          ns("clone_heatmap"),
          height = paste0(clone_heatmap_height(), "px")
        )
      )
    })

    output$clone_heatmap_download <- downloadHandler(
      filename = function() {
        paste0("clone_viability_", format(Sys.Date(), "%Y%m%d"), ".png")
      },
      content = function(file) {
        p <- clone_heatmap_obj()
        ggplot2::ggsave(
          file,
          plot = p,
          width = 8.5,
          height = max(4, clone_heatmap_height() / 100),
          dpi = 300,
          bg = "white",
          limitsize = FALSE
        )
      }
    )

    # Card annotations (icon + description + hint) are statically pre-rendered
    # in the UI (visible on first paint). Clear them once predictions exist,
    # and restore them if a later run produces no result again.
    observe({
      has_res <- !is.null(clone_pred())
      session$sendCustomMessage(
        "set-html",
        list(
          id = ns("clone_desc"),
          html = if (has_res) "" else as.character(clone_desc_html())
        )
      )
      session$sendCustomMessage(
        "set-html",
        list(
          id = ns("patient_desc"),
          html = if (has_res) "" else as.character(patient_desc_html())
        )
      )
    })

    # Patient table
    output$patient_table <- renderDT({
      req(patient_pred())
      df <- as.data.frame(patient_pred())
      df$patient <- rownames(df)
      df <- df[, c("patient", setdiff(names(df), "patient"))]

      datatable(
        df,
        options = list(
          pageLength = 15,
          dom = "ftip",
          columnDefs = list(list(className = "dt-center", targets = "_all"))
        ),
        rownames = FALSE,
        class = "display"
      ) %>%
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
      div(
        style = "margin-top: 0.75rem;",
        downloadButton(
          ns("download_pred"),
          "Download Predictions (.csv)",
          class = "btn-outline-primary btn-sm"
        )
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
