# Home Page Module — Redesigned
mod_home_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Hero Section — Lightweight ERV-Atlas-Style Banner
    div(class = "hero",
      div(class = "hero-inner",
        div(class = "hero-text",
          div(class = "hero-title-row",
            h1(class = "hero-title", "PERCEPTION-shiny")
          ),
          p(class = "hero-tagline",
            "Precision oncology from single-cell transcriptomics \u2014 ",
            "trained on DepMap, applied to patient tumors."
          ),
          p(class = "hero-desc",
            "PERCEPTION (PERsonalized single-Cell Expression-based Planning for Treatments In ONcology) ",
            "predicts patient response and resistance to cancer treatment using single-cell expression profiles. ",
            "It trains elastic-net models on DepMap cell-line data, then applies them to patient single-cell data ",
            "for clone-level drug sensitivity and patient-level response stratification."
          ),
          div(class = "hero-actions",
            actionButton(ns("go_data"), "Quick Start", class = "btn-hero-primary", icon = icon("rocket")),
            actionButton(ns("go_demo"), "Load Demo", class = "btn-hero-secondary", icon = icon("flask")),
            tags$a(href = "https://github.com/WangLabCSU/PERCEPTIONx/blob/main/docs/PERCEPTION-shiny.md",
                   target = "_blank", class = "btn-hero-secondary",
                   icon("book-open"), " Tutorial")
          )
        ),
        div(class = "hero-brand",
          div(class = "hero-logo",
            tags$img(src = "favicon.svg", alt = "PERCEPTION-shiny logo")
          )
        )
      )
    ),

    # Stat Cards — aligned with hero-inner (max-width 1200px, centered)
    fluidRow(style = "padding: 1.5rem 0 0.5rem;",
      div(class = "stat-grid",
        div(class = "stat-card stat-blue animate-fade-in-up delay-1",
          div(class = "stat-icon", icon("dna")),
          div(class = "stat-value", "15K+"),
          div(class = "stat-label", "DepMap Genes")
        ),
        div(class = "stat-card stat-teal animate-fade-in-up delay-2",
          div(class = "stat-icon", icon("table-cells")),
          div(class = "stat-value", "Single-Cell"),
          div(class = "stat-label", "Resolution")
        ),
        div(class = "stat-card stat-amber animate-fade-in-up delay-3",
          div(class = "stat-icon", icon("pills")),
          div(class = "stat-value", "44"),
          div(class = "stat-label", "Pretrained Drug Models")
        )
      )
    ),

    hr(class = "section-divider"),

    # Guided Workflow — live pipeline stepper (state derives from shared data)
    fluidRow(style = "padding: 0 2rem;",
      div(class = "section-header animate-fade-in",
        icon("route"),
        h4("Pipeline — Where You Are"),
        p(class = "text-muted animate-fade-in", style = "font-size: 0.88rem; margin: 0.2rem 0 0 auto;",
          "Status updates live as you load data, train and predict. Click any step to jump to it.")
      ),
      div(class = "card animate-fade-in-up",
        uiOutput(ns("workflow_stepper"))
      )
    ),

    hr(class = "section-divider"),

    # Data Status Dashboard — real-time overview
    div(class = "data-status-section",
      div(class = "section-header animate-fade-in",
        icon("dashboard"),
        h4("Data Status"),
        p(class = "text-muted animate-fade-in", style = "font-size: 0.88rem; margin: 0.2rem 0 0 auto;",
          "Real-time overview of loaded data.")
      ),
      uiOutput(ns("data_dashboard"))
    ),

    hr(class = "section-divider"),

    # Key Features
    fluidRow(style = "padding: 0 2rem 1rem;",
      div(class = "section-header animate-fade-in",
        icon("star"),
        h4("Key Features")
      ),
      div(class = "feature-grid",
        div(class = "feature-card card-red animate-fade-in-up delay-1",
          div(class = "feature-icon icon-red", icon("cloud-arrow-down")),
          h5("DepMap Integration"),
          p("Seamlessly load and cache DepMap cell line expression and drug response data for model training.")
        ),
        div(class = "feature-card card-purple animate-fade-in-up delay-2",
          div(class = "feature-icon icon-purple", icon("code-branch")),
          h5("Elastic Net & Random Forest"),
          p("Train drug response models using glmnet (elastic net) or random forest with cross-validated hyperparameter tuning.")
        ),
        div(class = "feature-card card-blue animate-fade-in-up delay-3",
          div(class = "feature-icon icon-blue", icon("dna")),
          h5("Clone-Level Prediction"),
          p("Predict drug sensitivity at single-cell clone resolution, then aggregate to patient-level response.")
        ),
        div(class = "feature-card card-amber animate-fade-in-up delay-4",
          div(class = "feature-icon icon-amber", icon("palette")),
          h5("Rich Visualization"),
          p("Generate UMAP plots, clone distribution stacks, ROC curves, and comprehensive patient response panels.")
        )
      )
    ),

    # Citation & Repository
    fluidRow(style = "padding: 0 2rem 3rem;",
      div(class = "section-header animate-fade-in",
        icon("book-open"),
        h4("Citation & Repository")
      ),
      div(class = "citation-box animate-fade-in",
        div(class = "citation-row",
          div(class = "citation-col",
            div(class = "citation-item",
              h6(icon("file-lines"), " Original Paper"),
              p(
                strong("Sinha, S., Vegesna, R., Mukherjee, S."),
                " et al. PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. ",
                em("Nat Cancer"), " 5, 938-952 (2024). ",
                br(),
                a("DOI: 10.1038/s43018-024-00756-7",
                  href = "https://doi.org/10.1038/s43018-024-00756-7",
                  target = "_blank")
              )
            ),
            div(class = "citation-item",
              h6(icon("code-branch"), " Original Codebase"),
              p(
                icon("github"),
                " ",
                a("github.com/dm-lab-04/PERCEPTION",
                  href = "https://github.com/dm-lab-04/PERCEPTION",
                  target = "_blank"),
                tags$small(" (original Python implementation by Sinha et al.)",
                           class = "text-muted")
              )
            )
          ),
          div(class = "citation-col",
            div(class = "citation-item",
              h6(icon("box-open"), " R Package Repository"),
              p(
                icon("github"),
                " ",
                a("github.com/WangLabCSU/PERCEPTIONx",
                  href = "https://github.com/WangLabCSU/PERCEPTIONx",
                  target = "_blank"),
                br(),
                tags$small("Install with: ", code("devtools::install_github('WangLabCSU/PERCEPTIONx')"),
                           class = "text-muted")
              )
            )
          )
        )
      )
    )
  )
}

mod_home_server <- function(id, shared, main_session) {
  moduleServer(id, function(input, output, session) {
    # Navigation buttons — use main session for cross-module nav
    observeEvent(input$go_data, bslib::nav_select("navbar", selected = "data", session = main_session))
    observeEvent(input$wf_1, bslib::nav_select("navbar", selected = "data", session = main_session))
    observeEvent(input$wf_2, bslib::nav_select("navbar", selected = "train", session = main_session))
    observeEvent(input$wf_3, bslib::nav_select("navbar", selected = "predict", session = main_session))
    observeEvent(input$wf_4, bslib::nav_select("navbar", selected = "visualize", session = main_session))

    # Live workflow stepper — step states are derived directly from the shared
    # data so the pipeline tracker updates the moment anything is loaded.
    output$workflow_stepper <- renderUI({
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
      depmap_ok <- !is.null(shared$depmap)
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
    })

    # Load Demo Data — full pipeline, via the shared helper (same code path
    # as the Data-tab "Load Demo Data" button; see shiny_helpers.R).
    observeEvent(input$go_demo, {
      run_demo_pipeline(shared, on_success = function()
        bslib::nav_select("navbar", selected = "data", session = main_session))
    })

    # Data Status Dashboard — real-time overview
    output$data_dashboard <- renderUI({
      # DepMap stats (only computed when data is loaded)
      depmap_cell_lines <- if (!is.null(shared$depmap) && !is.null(shared$depmap$expression_rnorm))
        ncol(shared$depmap$expression_rnorm) else NULL
      depmap_drugs <- if (!is.null(shared$depmap) && !is.null(shared$depmap$secondary_prism))
        nrow(shared$depmap$secondary_prism) else NULL
      depmap_sc_models <- if (!is.null(shared$depmap) && !is.null(shared$depmap$scRNA_complete))
        ncol(shared$depmap$scRNA_complete) else NULL

      items <- list(
        list(name = "DepMap Reference", icon_name = "database",
             loaded = !is.null(shared$depmap),
             detail = if (!is.null(shared$depmap))
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
    })
  })
}
