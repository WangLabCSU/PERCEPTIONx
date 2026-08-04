# Prediction Module
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
          "then aggregate to patient-level response using weighted averaging strategies."
        )
      )
    ),

    # Prerequisites check
    fluidRow(style = "margin-top: 1rem;",
      column(12,
        uiOutput(ns("prereq_check"))
      )
    ),

    fluidRow(style = "margin-top: 0.5rem; align-items: flex-start;",
      # Configuration
      column(4,
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
                            "Min (Most Resistant)" = "min",
                            "Max (Most Sensitive)" = "max"
                          ),
                          selected = "weighted_max"),
              div(class = "info-box", style = "margin-top: 0.5rem; font-size: 0.78rem;",
                strong("weighted_max"), ": top N most resistant clones (recommended)",
                br(), strong("weighted_average"), ": weighted average across all clones",
                br(), strong("min/max"), ": most resistant/sensitive clone"
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
            uiOutput(ns("clone_pred_status")),
            plotlyOutput(ns("clone_heatmap"), height = "400px")
          )
        ),

        # Patient-level table
        div(class = "card", style = "margin-top: 1rem;",
          div(class = "card-header",
            icon("user"), " Patient-Level Predictions"
          ),
          div(class = "card-body",
            DTOutput(ns("patient_table")),
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

    # Get expression — prefer prepared_data$clone_expression_rnorm
    current_expr <- reactive({
      if (input$expr_source == "loaded") {
        if (!is.null(shared$prepared_data)) {
          shared$prepared_data$clone_expression_rnorm
        } else {
          shared$user_expr
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

    # Run prediction — Stage 1: clone-level via predict_drugs()
    clone_pred <- eventReactive(input$predict, {
      model <- current_model()
      expr <- current_expr()

      if (is.null(model)) {
        showNotification("No model available. Train or upload a model first.", type = "error")
        return(NULL)
      }
      if (is.null(expr)) {
        showNotification("No expression data available. Load or upload data first.", type = "error")
        return(NULL)
      }

      w <- Waiter$new(
        html = tagList(div(class = "spinner-ring"), h4("Predicting..."),
                       p(class = "text-muted", "Running predict_drugs()")),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()

      tryCatch({
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
            w$hide()
            showNotification(
              paste0("Gene name mismatch: none of the model's ", length(all_features),
                     " features found in your expression data. ",
                     "Check gene name formats (hyphens vs dots) between training and prediction data."),
              type = "error", duration = 12)
            return(NULL)
          }
          if (length(overlap) < length(all_features) * 0.5) {
            showNotification(
              paste0("Note: ", length(overlap), "/", length(all_features),
                     " model features found in data. Predictions may be less reliable."),
              type = "warning", duration = 8)
          }
        }

        result <- PERCEPTION::predict_drugs(
          model_list = model,
          expr = expr
        )
        w$hide()
        showNotification("Clone-level prediction complete (predict_drugs)", type = "message")
        result
      }, error = function(e) {
        w$hide()
        showNotification(paste("Prediction error:", e$message), type = "error", duration = 10)
        NULL
      })
    })

    # Stage 2: patient-level via predict_patients()
    patient_pred <- reactive({
      req(clone_pred())
      tryCatch({
        clone_pred_mat <- clone_pred()  # matrix: clones x drugs

        if (!is.null(shared$prepared_data)) {
          # Standard path: use prepare_data() output directly
          result <- PERCEPTION::predict_patients(
            clone_pred_mat,
            shared$prepared_data,
            mode = input$agg_mode
          )
        } else {
          # Fallback: manually build from user_clones (legacy)
          clone_data <- shared$user_clones
          clone_rows <- rownames(clone_pred_mat)
          clone_to_patient <- split(clone_data$patient, clone_data$clone_id)
          clone_to_patient <- lapply(clone_to_patient, unique)

          clone_killing_list <- lapply(clone_rows, function(cl) {
            pat <- clone_to_patient[[cl]]
            if (is.null(pat) || length(pat) == 0) return(NULL)
            if (length(pat) > 1) pat <- pat[1]
            df <- data.frame(
              patient = pat,
              clone_id = cl,
              stringsAsFactors = FALSE
            )
            for (drug in colnames(clone_pred_mat)) {
              df[[drug]] <- clone_pred_mat[cl, drug]
            }
            df
          })
          clone_killing_df <- do.call(rbind, clone_killing_list)

          if (is.null(clone_killing_df) || nrow(clone_killing_df) == 0) {
            showNotification("No matching clones between prediction and annotation", type = "error")
            return(NULL)
          }

          clone_counts <- as.data.frame.matrix(
            table(clone_data$patient, clone_data$clone_id)
          )
          clone_counts$patients <- rownames(clone_counts)

          result <- PERCEPTION::predict_patients(
            clone_killing_df,
            clone_counts,
            mode = input$agg_mode
          )
        }

        shared$patient_pred <- result
        shared$predictions <- clone_pred_mat
        result
      }, error = function(e) {
        showNotification(paste("Patient prediction error:", e$message), type = "error")
        NULL
      })
    })

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

    output$clone_pred_status <- renderUI({
      if (is.null(clone_pred())) {
        div(class = "text-muted", style = "text-align: center; padding: 2rem;",
          icon("arrow-left"), " Configure parameters and click Run Prediction"
        )
      } else {
        tagList(
          span(class = "status-badge loaded",
            span(class = "status-dot green"),
            paste(nrow(clone_pred()), "clones x", ncol(clone_pred()), "drugs")
          )
        )
      }
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