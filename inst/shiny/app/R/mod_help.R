# Help Module — Reference Documentation
mod_help_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        div(class = "section-header",
          icon("circle-question"),
          h4("Reference Documentation")
        ),
        div(class = "info-box",
          icon("book"), " Detailed reference for data formats, parameters, and interpretation. ",
          "For a quick getting-started guide, visit the ", strong("Home"), " page."
        )
      )
    ),

    fluidRow(style = "margin-top: 1rem;",
      # Left sidebar: section navigation
      column(3,
        div(class = "card help-nav-card",
          div(class = "card-header",
            icon("list-ul"), " Sections"
          ),
          div(class = "card-body", style = "padding: 0.6rem 0.5rem !important;",
            div(class = "help-nav-list",
              actionButton(ns("nav_data_req"), "Data Requirements", class = "help-nav-btn"),
              actionButton(ns("nav_training"), "Training Parameters", class = "help-nav-btn"),
              actionButton(ns("nav_prediction"), "Aggregation Modes", class = "help-nav-btn"),
              actionButton(ns("nav_viz_guide"), "Visualization Guide", class = "help-nav-btn"),
              actionButton(ns("nav_faq"), "FAQ", class = "help-nav-btn"),
              actionButton(ns("nav_citation"), "Citation", class = "help-nav-btn")
            )
          )
        )
      ),

      # Main content
      column(9,
        # Data Requirements
        div(id = ns("section_data_req"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("database")),
            div(
              h5("Data Requirements"),
              p(class = "help-subtitle", "What data you need and how to prepare it")
            )
          ),
          div(class = "help-data-grid",
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--primary);", icon("layer-group")),
              strong("DepMap Data"),
              p("Automatically downloaded via ", code("load_depmap()"), ". Includes bulk expression, single-cell expression, and drug response (AUC) data from DepMap."),
              tags$span(class = "status-badge loaded", "Auto-loaded")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--accent);", icon("table")),
              strong("Expression Matrix"),
              p("A gene x cell single-cell expression matrix. Raw counts or normalized values. Rank normalization is applied automatically during Seurat clustering."),
              tags$span(class = "status-badge unloaded", "User upload")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--success);", icon("users")),
              strong("Patient-Cell Mapping"),
              p("CSV with columns: ", code("cell_id"), " and ", code("patient_id"), ". Maps each cell to its patient. Clones are auto-detected via Seurat clustering."),
              tags$span(class = "status-badge unloaded", "User upload")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--primary-light);", icon("stethoscope")),
              strong("Clinical Response"),
              p("(Optional) CSV with columns: ", code("patient"), ", ", code("response"), " (Responder/Non-responder). Required for ROC curves and boxplots."),
              tags$span(class = "status-badge unloaded", "Optional")
            )
          )
        ),

        # Training Parameters
        div(id = ns("section_training"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("sliders-h")),
            div(
              h5("Training Parameters"),
              p(class = "help-subtitle", "Configure your model training")
            )
          ),
          div(class = "help-param-grid",
            div(class = "help-param-item",
              div(class = "help-param-name", icon("tag"), " Drug Name"),
              p("Must match a drug in DepMap secondary_prism data (e.g. abemaciclib, erlotinib). Use the dropdown on the Data tab for the full list of 44 supported drugs.")
            ),
            div(class = "help-param-item",
              div(class = "help-param-name", icon("filter"), " Cancer Type"),
              p("Filter cell lines by cancer type. Use ", code("PanCan"), " for all cancer types, or specify e.g. ", code("Breast"), ", ", code("Lung"), ".")
            ),
            div(class = "help-param-item",
              div(class = "help-param-name", icon("ban"), " Exclude Cancer"),
              p("Exclude specific cancer types from training. Set to ", code("PanCan"), " to exclude none.")
            ),
            div(class = "help-param-item",
              div(class = "help-param-name", icon("dna"), " Genes of Interest"),
              p("Gene symbols for feature selection. Leave ", strong("empty to use all genes"), " from DepMap (~15K). Should overlap with your scRNA data for prediction.")
            ),
            div(class = "help-param-item",
              div(class = "help-param-name", icon("sort-numeric-down"), " k_features"),
              p("Number of top-ranked features (by Pearson correlation) to use. Default: 100. Higher = more features but slower.")
            ),
            div(class = "help-param-item",
              div(class = "help-param-name", icon("cog"), " Model Type"),
              p(code("glmnet"), " (elastic net, recommended) with 5-fold CV, or ", code("rf"), " (random forest).")
            )
          )
        ),

        # Aggregation Modes
        div(id = ns("section_prediction"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("layer-group")),
            div(
              h5("Patient Aggregation Modes"),
              p(class = "help-subtitle", "How clone-level predictions become patient-level")
            )
          ),
          div(class = "help-data-grid",
            div(class = "help-data-card", style = "border-left: 3px solid var(--success);",
              div(class = "help-data-icon", style = "color: var(--success);", icon("star")),
              strong("weighted_max"),
              tags$span(class = "status-badge loaded", "Recommended"),
              p("Weighted average of the N most resistant clones, where N is determined by clone frequency.")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--primary);", icon("balance-scale")),
              strong("weighted_average"),
              p("Weighted average across all clones, weighted by clone proportion.")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--error);", icon("arrow-down")),
              strong("min"),
              p("Takes the most resistant clone's prediction (pessimistic).")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--success-light);", icon("arrow-up")),
              strong("max"),
              p("Takes the most sensitive clone's prediction (optimistic).")
            )
          )
        ),

        # Visualization Guide
        div(id = ns("section_viz_guide"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("chart-line")),
            div(
              h5("Visualization Guide"),
              p(class = "help-subtitle", "Available plot types and their data requirements")
            )
          ),
          div(class = "help-viz-grid",
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("chart-bar")),
              div(
                strong("Clone Distribution"),
                p("Stacked bar chart showing clone proportions per patient."),
                tags$span(class = "status-badge unloaded", "Needs: Clones")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("chart-line")),
              div(
                strong("Clone Killing Lollipop"),
                p("Lollipop chart of predicted viability per clone."),
                tags$span(class = "status-badge unloaded", "Needs: Predictions + Clones")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("wave-square")),
              div(
                strong("ROC Curve"),
                p("Receiver Operating Characteristic curve for response prediction."),
                tags$span(class = "status-badge unloaded", "Needs: Patient Pred + Response")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("box-open")),
              div(
                strong("Response Boxplot"),
                p("Boxplot comparing predicted viability between responders and non-responders."),
                tags$span(class = "status-badge unloaded", "Needs: Patient Pred + Response")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("gauge-high")),
              div(
                strong("Model Performance"),
                p("Cross-validation curve from glmnet model."),
                tags$span(class = "status-badge unloaded", "Needs: Trained Model")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("people-group")),
              div(
                strong("Patient Response Panel"),
                p("Comprehensive panel combining boxplot, ROC, and statistics."),
                tags$span(class = "status-badge unloaded", "Needs: Patient Pred + Response")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("braille")),
              div(
                strong("UMAP Gene Expression"),
                p("UMAP plot colored by selected gene expression."),
                tags$span(class = "status-badge unloaded", "Needs: Predictions + Expression + UMAP")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("arrows-left-right")),
              div(
                strong("UMAP Drug Killing"),
                p("UMAP plot colored by predicted drug killing."),
                tags$span(class = "status-badge unloaded", "Needs: Predictions + UMAP")
              )
            )
          )
        ),

        # FAQ
        div(id = ns("section_faq"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("comment-dots")),
            div(
              h5("Frequently Asked Questions"),
              p(class = "help-subtitle", "Common questions and answers")
            )
          ),
          div(class = "help-faq-list",
            div(class = "help-faq-item",
              div(class = "help-faq-q", icon("question-circle"), " How long does training take?"),
              div(class = "help-faq-a",
                "Training a single drug model typically takes 30 seconds to 5 minutes depending on the number of features and cores used. Using all genes (~15K) with k_features=100 is the default and works well."
              )
            ),
            div(class = "help-faq-item",
              div(class = "help-faq-q", icon("question-circle"), " Can I use my own bulk expression data?"),
              div(class = "help-faq-a",
                "Currently, PERCEPTION trains on DepMap data. You can use the trained model to predict on your own single-cell data after rank normalization."
              )
            ),
            div(class = "help-faq-item",
              div(class = "help-faq-q", icon("question-circle"), " What if my gene names don't match?"),
              div(class = "help-faq-a",
                "Gene symbols must match between training and prediction data. The app automatically intersects genes between DepMap and your expression matrix. Check for hyphen vs dot formatting (e.g. ", code("BRCA1"), " vs ", code("BRCA.1"), ")."
              )
            ),
            div(class = "help-faq-item",
              div(class = "help-faq-q", icon("question-circle"), " How do I preprocess my data?"),
              div(class = "help-faq-a",
                "Use ", code("rank_normalization_mat()"), " to rank-normalize your expression matrix. This is the same preprocessing applied to DepMap data."
              )
            ),
            div(class = "help-faq-item",
              div(class = "help-faq-q", icon("question-circle"), " Should I specify Genes of Interest?"),
              div(class = "help-faq-a",
                "Leave it ", strong("empty"), " to use all ~15K genes from DepMap (recommended for most users). Specify a subset only if you want to restrict feature selection to known biomarkers."
              )
            ),
            div(class = "help-faq-item",
              div(class = "help-faq-q", icon("question-circle"), " Download fails for DepMap or models?"),
              div(class = "help-faq-a",
                "Enable the ", strong("Mirror"), " checkbox (default ON) to use mirror download. If still failing, download the files manually from GitHub Releases and use the file upload option."
              )
            )
          )
        ),

        # Citation
        div(id = ns("section_citation"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("book-open")),
            div(
              h5("Citation"),
              p(class = "help-subtitle", "If you use PERCEPTION, please cite")
            )
          ),
          div(class = "citation-box",
            "Sinha, S., Vegesna, R., Mukherjee, S. et al. PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. ",
            em("Nat Cancer"), " 5, 938-952 (2024). ",
            a("DOI: 10.1038/s43018-024-00756-7",
              href = "https://doi.org/10.1038/s43018-024-00756-7",
              target = "_blank"),
            br(), br(),
            icon("github"),
            a(" github.com/SunPast/PERCEPTION",
              href = "https://github.com/SunPast/PERCEPTION",
              target = "_blank")
          )
        )
      )
    )
  )
}

mod_help_server <- function(id, main_session) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Section navigation scroll — use namespaced IDs
    observeEvent(input$nav_data_req, {
      session$sendCustomMessage("scroll-to", ns("section_data_req"))
    })
    observeEvent(input$nav_training, {
      session$sendCustomMessage("scroll-to", ns("section_training"))
    })
    observeEvent(input$nav_prediction, {
      session$sendCustomMessage("scroll-to", ns("section_prediction"))
    })
    observeEvent(input$nav_viz_guide, {
      session$sendCustomMessage("scroll-to", ns("section_viz_guide"))
    })
    observeEvent(input$nav_faq, {
      session$sendCustomMessage("scroll-to", ns("section_faq"))
    })
    observeEvent(input$nav_citation, {
      session$sendCustomMessage("scroll-to", ns("section_citation"))
    })
  })
}
