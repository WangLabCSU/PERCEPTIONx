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
      # Parameters Panel
      column(5,
        div(class = "card animate-fade-in-up",
          div(class = "card-header",
            icon("sliders-h"), " Parameters"
          ),
          div(class = "card-body",

            # Section: Drug & Cancer Type
            div(class = "param-section",
              h6(class = "param-section-title", icon("tag"), " Drug & Cancer Type"),
              selectizeInput(ns("drug"), "Drug Name (select one or more for combination)",
                             choices = NULL, selected = "abemaciclib", multiple = TRUE,
                             options = list(maxItems = NULL, placeholder = "Select drug(s)",
                                            plugins = list("remove_button"))),
              selectizeInput(ns("cancer_type"), "Cancer Type (include)",
                             choices = NULL, selected = "PanCan",
                             options = list(maxItems = 1, placeholder = "Select cancer type")),
              selectizeInput(ns("exclude_cancer"), "Cancer Type (exclude)",
                             choices = NULL, selected = "PanCan",
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
                            height = "120px"),
              fileInput(ns("goi_file"), "Or Upload Gene List (.txt, .csv)",
                        accept = c(".txt", ".csv")),
              numericInput(ns("k_features"), "Top k Features",
                           value = 100, min = 10, max = 5000, step = 10)
            ),

            # Section: Model Configuration
            div(class = "param-section",
              h6(class = "param-section-title", icon("cog"), " Model Configuration"),
              selectInput(ns("model_type"), "Algorithm",
                          choices = c("Elastic Net (glmnet)" = "glmnet", "Random Forest" = "rf"),
                          selected = "glmnet"),
              numericInput(ns("ncores"), "CPU Cores",
                           value = 1, min = 1, max = parallel::detectCores(), step = 1)
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
      column(7,
        div(class = "card animate-fade-in-up delay-1",
          div(class = "card-header",
            icon("chart-bar"), " Training Results"
          ),
          div(class = "card-body",
            # Progress
            uiOutput(ns("progress")),

            # Model summary
            uiOutput(ns("model_summary")),

            # Performance metrics — p-value & correlation
            uiOutput(ns("perf_metrics")),

            # Performance plot
            div(class = "viz-plot-wrapper",
              plotlyOutput(ns("perf_plot"), height = "320px")
            ),

            # Download
            uiOutput(ns("download_btn"))
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
      if (is.null(shared$depmap)) missing <- c(missing, "DepMap Data")
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

    # Populate drug & cancer type choices when DepMap is loaded
    observe({
      depmap <- shared$depmap
      if (!is.null(depmap)) {
        # Drug choices: 44 FDA-approved drugs (recommended subset from train_models default)
        fda_drugs <- c(
          "abemaciclib", "afatinib", "axitinib", "azacitidine", "cladribine",
          "clofarabine", "cobimetinib", "dabrafenib", "dasatinib", "daunorubicin",
          "decitabine", "docetaxel", "doxorubicin", "epirubicin", "erlotinib",
          "etoposide", "gefitinib", "gemcitabine", "homoharringtonine", "ibrutinib",
          "icotinib", "ixabepilone", "lapatinib", "lenvatinib", "midostaurin",
          "niraparib", "osimertinib", "paclitaxel", "palbociclib", "ponatinib",
          "romidepsin", "sunitinib", "temsirolimus", "teniposide", "thioguanine",
          "topotecan", "trametinib", "vandetanib", "vemurafenib", "vinblastine",
          "vincristine", "vindesine", "vinflunine", "vinorelbine"
        )
        # Filter to drugs that actually exist in DepMap
        if (!is.null(depmap$secondary_screen_drugAnnotation)) {
          available <- unique(depmap$secondary_screen_drugAnnotation$CommonName)
          fda_drugs <- intersect(fda_drugs, available)
        }
        updateSelectizeInput(session, "drug", choices = fda_drugs, selected = "abemaciclib", server = TRUE)

        # Cancer type choices: PanCan + unique lineages from DepMap annotation
        cancer_choices <- "PanCan"
        if (!is.null(depmap$annotation_20Q4) && !is.null(depmap$annotation_20Q4$lineage)) {
          lineages <- as.character(depmap$annotation_20Q4$lineage)
          lineages <- sort(unique(lineages))
          lineages <- lineages[!is.na(lineages) & nchar(lineages) > 0]
          cancer_choices <- c("PanCan", lineages)
        }
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
        genes <- unlist(strsplit(input$goi, "[,\\s\\n]+"))
        genes <- trimws(genes)
        genes <- genes[nchar(genes) > 0]
      }
      genes
    })

    # Train
    observeEvent(input$train, {
      # Validate
      if (is.null(shared$depmap)) {
        showNotification("Please load DepMap data first (Data tab)", type = "error")
        return()
      }
      # GOI is optional — if empty, use all genes from DepMap
      goi <- goi_parsed()
      if (is.null(goi) || length(goi) == 0) {
        goi <- rownames(shared$depmap$expression_rnorm)
        showNotification(paste0("Using all ", length(goi), " genes from DepMap data"), type = "message")
      }
      if (is.null(input$drug) || length(input$drug) == 0) {
        showNotification("Please select at least one drug", type = "error")
        return()
      }

      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Training model..."),
          p(class = "text-muted", "This may take a few seconds")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()

      tryCatch({
        result <- PERCEPTIONx::train_models(
          drug_list = input$drug,
          cancer_type = input$cancer_type,
          exclude_cancer = input$exclude_cancer,
          GOI = goi,
          k_features_values = input$k_features,
          model_type = input$model_type,
          ncores = input$ncores,
          output_dir = tempdir()
        )
        trained(result)
        shared$models <- result
        # Also cache permanently and mark active
        if (is.null(shared$model_cache)) shared$model_cache <- list()
        for (nm in names(result)) {
          shared$model_cache[[nm]] <- result[[nm]]
          shared$model_active[[nm]] <- TRUE
        }
        w$hide()
        showNotification("Training completed successfully!", type = "message")
      }, error = function(e) {
        w$hide()
        showNotification(paste("Training error:", e$message), type = "error")
      })
    })

    # Reset
    observeEvent(input$reset, {
      trained(NULL)
      updateSelectizeInput(session, "drug", selected = "abemaciclib")
      updateSelectizeInput(session, "cancer_type", selected = "PanCan")
      updateSelectizeInput(session, "exclude_cancer", selected = "PanCan")
      updateTextAreaInput(session, "goi", value = "")
    })

    # Progress UI
    output$progress <- renderUI({
      if (is.null(trained())) {
        div(class = "text-muted", style = "text-align: center; padding: 2rem;",
          icon("hourglass-start"), " No training run yet"
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

    # Performance Plot
    output$perf_plot <- renderPlotly({
      req(trained())
      models <- trained()
      p <- plot_model_performance(models, base_size = 13)
      ggplotly(p, tooltip = c("x", "y", "colour")) %>%
        layout(font = list(family = "Inter, sans-serif", size = 12))
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
          if (v < 0.001) return("< 0.001")
          format(round(v, digits), nsmall = digits)
        }

        tags$tr(
          tags$td(strong(d)),
          tags$td(fmt(bulk$cor), tags$span(class = "text-muted", style = "font-size:0.72rem; display:block;", "p = ", fmt(bulk$p))),
          tags$td(fmt(pseudo$cor), tags$span(class = "text-muted", style = "font-size:0.72rem; display:block;", "p = ", fmt(pseudo$p))),
          tags$td(fmt(sc$cor), tags$span(class = "text-muted", style = "font-size:0.72rem; display:block;", "p = ", fmt(sc$p)))
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
        drug_str <- paste(input$drug, collapse = "_")
        paste0(drug_str, "_model_", format(Sys.Date(), "%Y%m%d"), ".RDS")
      },
      content = function(file) {
        saveRDS(trained(), file)
      }
    )
  })
}
