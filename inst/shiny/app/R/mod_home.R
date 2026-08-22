# Home Page Module — Redesigned
mod_home_ui <- function(id) {
  ns <- NS(id)
  # Empty state for the FIRST paint: the stepper + dashboard are rendered
  # statically here (instant, no server round trip) and only get updated
  # later via the 'set-html' message when shared state changes.
  empty_state <- list()
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
        div(id = ns("workflow_stepper"),
            workflow_stepper_html(empty_state, ns))
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
      div(id = ns("data_dashboard"),
          data_dashboard_html(empty_state))
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
    ns <- session$ns
    # Navigation buttons — use main session for cross-module nav
    observeEvent(input$go_data, bslib::nav_select("navbar", selected = "data", session = main_session))
    observeEvent(input$wf_1, bslib::nav_select("navbar", selected = "data", session = main_session))
    observeEvent(input$wf_2, bslib::nav_select("navbar", selected = "train", session = main_session))
    observeEvent(input$wf_3, bslib::nav_select("navbar", selected = "predict", session = main_session))
    observeEvent(input$wf_4, bslib::nav_select("navbar", selected = "visualize", session = main_session))

    # Live workflow stepper — step states are derived directly from the shared
    # data so the pipeline tracker updates the moment anything is loaded.
    # Rendered once statically in the UI; this observer pushes updates when
    # shared state changes (no initial server round trip on first paint).
    observe({
      state <- reactiveValuesToList(shared)
      session$sendCustomMessage("set-html", list(
        id   = ns("workflow_stepper"),
        html = as.character(workflow_stepper_html(state, ns))
      ))
    })

    # Load Demo Data — full pipeline, via the shared helper (same code path
    # as the Data-tab "Load Demo Data" button; see shiny_helpers.R).
    observeEvent(input$go_demo, {
      run_demo_pipeline(shared, on_success = function()
        bslib::nav_select("navbar", selected = "data", session = main_session))
    })

    # Data Status Dashboard — real-time overview
    observe({
      state <- reactiveValuesToList(shared)
      session$sendCustomMessage("set-html", list(
        id   = ns("data_dashboard"),
        html = as.character(data_dashboard_html(state))
      ))
    })
  })
}
