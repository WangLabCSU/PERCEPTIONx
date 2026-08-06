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

    # Guided Workflow — How to Use
    fluidRow(style = "padding: 0 2rem;",
      div(class = "section-header animate-fade-in",
        icon("route"),
        h4("How to Use — Follow These Steps")
      ),
      p(class = "text-muted animate-fade-in", style = "font-size: 0.88rem; margin-bottom: 1.2rem;",
        "Click any step below to jump directly to that section. Data flows left to right."),
      div(class = "guided-steps",
        div(class = "step-card animate-fade-in-up delay-1",
          div(class = "step-head",
            div(class = "step-number", "1"),
            div(class = "step-icon", icon("upload"))
          ),
          div(class = "step-title", "Load Data"),
          div(class = "step-desc", "Load DepMap reference data or upload your own expression matrix and clone annotations."),
          actionButton(ns("go_step1"), NULL, class = "step-go", icon = icon("arrow-right"), label = "Go to Data")
        ),
        div(class = "step-card animate-fade-in-up delay-2",
          div(class = "step-head",
            div(class = "step-number", "2"),
            div(class = "step-icon", icon("gear"))
          ),
          div(class = "step-title", "Train Model"),
          div(class = "step-desc", "Train a drug response model on DepMap bulk data using elastic net or random forest."),
          actionButton(ns("go_step2"), NULL, class = "step-go", icon = icon("arrow-right"), label = "Go to Train")
        ),
        div(class = "step-card animate-fade-in-up delay-3",
          div(class = "step-head",
            div(class = "step-number", "3"),
            div(class = "step-icon", icon("flask"))
          ),
          div(class = "step-title", "Predict"),
          div(class = "step-desc", "Apply the model to single-cell data for clone-level and patient-level drug sensitivity."),
          actionButton(ns("go_step3"), NULL, class = "step-go", icon = icon("arrow-right"), label = "Go to Predict")
        ),
        div(class = "step-card animate-fade-in-up delay-4",
          div(class = "step-head",
            div(class = "step-number", "4"),
            div(class = "step-icon", icon("chart-line"))
          ),
          div(class = "step-title", "Visualize"),
          div(class = "step-desc", "Generate UMAP, ROC, clone distribution, and patient response panel plots."),
          actionButton(ns("go_step4"), NULL, class = "step-go", icon = icon("arrow-right"), label = "Go to Visualize")
        )
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
    observeEvent(input$go_step1, bslib::nav_select("navbar", selected = "data", session = main_session))
    observeEvent(input$go_step2, bslib::nav_select("navbar", selected = "train", session = main_session))
    observeEvent(input$go_step3, bslib::nav_select("navbar", selected = "predict", session = main_session))
    observeEvent(input$go_step4, bslib::nav_select("navbar", selected = "visualize", session = main_session))

    # Load Demo Data — Full pipeline coverage
    observeEvent(input$go_demo, {
      set.seed(42)
      gene_names <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "RB1",
                       "APC", "BRAF", "CDH1", "CDKN2A", "ERBB2", "FGFR1", "ALK",
                       "MET", "RET", "ROS1", "NRAS", "HRAS", "MAP2K1", "MAPK1",
                       "JAK2", "STAT3", "MTOR", "AKT1", "AKT2", "CTNNB1", "SMAD4",
                       "VHL", "NF1", "NF2", "STK11", "FBXW7", "ARID1A", "KDM5C",
                       "KMT2D", "SETD2", "BAP1", "PBRM1", "NOTCH1", "NOTCH2",
                       "NOTCH3", "JAK1", "JAK3", "SOX9", "IDH1", "IDH2", "FLT3")
      n_cells <- 400
      n_patients <- 20
      cell_names <- paste0("CELL_", sprintf("%04d", 1:n_cells))
      patient_names <- paste0("PAT_", sprintf("%03d", 1:n_patients))

      # 1. Patient-Cell mapping + Clinical response (must precede expression generation)
      patient_assignment <- rep(patient_names, each = ceiling(n_cells / n_patients))[1:n_cells]
      patient_mapping <- data.frame(
        cell_id = cell_names,
        patient_id = patient_assignment,
        stringsAsFactors = FALSE
      )
      shared$user_mapping <- patient_mapping

      response_labels <- c(rep("Responder", 10), rep("Non-responder", 10))
      clinical_response <- data.frame(
        patient = patient_names,
        response = response_labels,
        stringsAsFactors = FALSE
      )
      shared$user_response <- clinical_response

      # 2. STRUCTURED expression matrix so biomarker plots show meaningful correlation.
      # Two drug-biomarker groups with OPPOSITE patterns:
      #   - abemaciclib biomarkers (genes 1-5): HIGH in responders, LOW in non-responders
      #   - erlotinib biomarkers (genes 6-10): LOW in responders, HIGH in non-responders
      # Other genes are background noise.
      is_responder_cell <- clinical_response$response[match(patient_assignment, clinical_response$patient)] == "Responder"
      abemaciclib_markers <- gene_names[1:5]   # TP53, BRCA1, EGFR, MYC, KRAS
      erlotinib_markers   <- gene_names[6:10]  # PIK3CA, PTEN, RB1, APC, BRAF

      expr_matrix <- matrix(0.1, nrow = length(gene_names), ncol = n_cells)
      rownames(expr_matrix) <- gene_names
      colnames(expr_matrix) <- cell_names

      for (g in abemaciclib_markers) {
        expr_matrix[g, is_responder_cell]  <- pmax(rnorm(sum(is_responder_cell),  mean = 12, sd = 2), 0.1)
        expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 2,  sd = 1), 0.1)
      }
      for (g in erlotinib_markers) {
        expr_matrix[g, is_responder_cell]  <- pmax(rnorm(sum(is_responder_cell),  mean = 2,  sd = 1), 0.1)
        expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 12, sd = 2), 0.1)
      }
      for (g in setdiff(gene_names, c(abemaciclib_markers, erlotinib_markers))) {
        expr_matrix[g, ] <- runif(n_cells, 0.5, 8)
      }
      storage.mode(expr_matrix) <- "numeric"

      # 3. Run prepare_data (Seurat clustering + rank normalization)
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Preparing demo data..."),
          p(class = "text-muted", "Running Seurat clustering")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      tryCatch({
        prepared <- PERCEPTIONx::prepare_data(
          method = "umap",
          expression_matrix = expr_matrix,
          patient_mapping = patient_mapping,
          seurat_resolution = 0.8,
          seurat_dims = 10,
          seurat_nfeatures = min(2000, length(gene_names))
        )
        shared$prepared_data <- prepared
        shared$user_expr <- expr_matrix  # Keep cell-level expression for spatial gene plots
        shared$user_clones <- prepared$cell_clone_map
        w$hide()

        # 4. Train models on STRUCTURED training data that mirrors the demo expression
        if (requireNamespace("caret", quietly = TRUE) && requireNamespace("glmnet", quietly = TRUE)) {
          make_structured_training <- function(marker_genes, direction, n_train = 160) {
            x_train <- matrix(0, nrow = length(gene_names), ncol = n_train)
            rownames(x_train) <- gene_names
            half <- n_train %/% 2
            responder_like <- seq_len(half)
            nonresponder_like <- seq.int(half + 1, n_train)

            for (g in marker_genes) {
              x_train[g, responder_like]    <- pmax(rnorm(half, mean = 12, sd = 2), 0.1)
              x_train[g, nonresponder_like] <- pmax(rnorm(half, mean = 2,  sd = 1), 0.1)
            }
            for (g in setdiff(gene_names, marker_genes)) {
              x_train[g, ] <- runif(n_train, 0.5, 8)
            }
            y_train <- direction * colMeans(x_train[marker_genes, , drop = FALSE])
            y_train <- y_train + rnorm(n_train, sd = 0.3)
            as.data.frame(cbind(y = y_train, t(x_train)))
          }

          make_drug_model <- function(drug_name, seed, marker_genes, direction) {
            set.seed(seed)
            train_df <- make_structured_training(marker_genes, direction, n_train = 160)
            caret_model <- caret::train(
              y ~ ., data = train_df,
              method = "glmnet",
              trControl = caret::trainControl(method = "cv", number = 3),
              tuneLength = 3
            )
            obj <- list(
              model = caret_model,
              performance_in_scRNA = data.frame(estimate.cor = c(0.45, 0.38), p.value = c(0.001, 0.01)),
              performance_in_bulk = data.frame(estimate.cor = c(0.55, 0.48), p.value = c(0.0001, 0.001)),
              performance_in_pseudo_bulk = data.frame(estimate.cor = c(0.50, 0.42), p.value = c(0.0005, 0.005)),
              predVSgroundTruth = list(pred_gt_scRNA = data.frame(Observed = rnorm(20), Test_pred_sc = rnorm(20))),
              single_best = marker_genes[1]
            )
            attr(obj, "drug_name") <- drug_name
            obj
          }

          shared$models <- list(
            abemaciclib = make_drug_model("abemaciclib", 101, abemaciclib_markers, direction = +1),
            erlotinib   = make_drug_model("erlotinib",   202, erlotinib_markers,   direction = -1)
          )
          shared$model_cache <- shared$models
          shared$model_active <- list(abemaciclib = TRUE, erlotinib = TRUE)
        }

        n_clones <- ncol(prepared$clone_expression_rnorm)
        showNotification(paste0("Demo data loaded! ", length(gene_names), " genes × ", n_cells,
                                " cells, ", n_patients, " patients, ", n_clones,
                                " clones. Explore all tabs!"),
                         type = "message", duration = 6)
        bslib::nav_select("navbar", selected = "data", session = main_session)
      }, error = function(e) {
        w$hide()
        showNotification(paste("Demo data preparation error:", e$message), type = "error", duration = 10)
      })
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
