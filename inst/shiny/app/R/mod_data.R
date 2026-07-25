# Data Loading Module
mod_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        div(class = "section-header",
          icon("database"),
          h4("Data Management")
        )
      )
    ),

    # Demo Data Banner
    fluidRow(
      column(12,
        div(class = "demo-banner animate-fade-in-up",
          div(class = "demo-icon", icon("flask")),
          div(class = "demo-text",
            strong("Try it out! "),
            "Click the button below to load demo data (50 genes × 400 cells, 20 patients) and explore all features without uploading anything."
          ),
          actionButton(ns("load_demo"), "Load Demo Data", class = "btn-demo btn-sm", icon = icon("play"))
        )
      )
    ),

    fluidRow(class = "data-management-row",
      # DepMap Data
      column(6,
        div(class = "card animate-fade-in-up delay-1",
          div(class = "card-header",
            icon("database"), " DepMap Data"
          ),
          div(class = "card-body",
            p(class = "text-muted", style = "font-size: 0.86rem;",
              "Load DepMap reference datasets including bulk expression, single-cell expression, drug response (AUC), and cell line annotations. This is a filtered version derived from the original DepMap release used in the PERCEPTION article, with unused tables and objects removed for efficiency."),
            tags$small(class = "text-muted", style = "display: block; margin-top: 0.3rem;",
              "To download manually, visit ",
              tags$a(href = "https://github.com/SunPast/PERCEPTION/releases/tag/depmap",
                     target = "_blank", "GitHub Release", style = "color: var(--primary); text-decoration: underline;"),
              "."
            ),
            hr(),
            div(class = "inline-form-row", style = "display: flex; gap: 0.5rem; align-items: center; flex-wrap: wrap;",
              actionButton(ns("load_depmap"), "Download & Load",
                           class = "btn-primary btn-sm", icon = icon("download")),
              checkboxInput(ns("depmap_mirror"), "Use mirror", value = TRUE)
            ),
            div(style = "margin-top: 0.6rem; border-top: 1px dashed var(--border); padding-top: 0.6rem;",
              tags$small(class = "text-muted", "Or load a pre-downloaded DepMap.RDS:"),
              div(style = "display: flex; gap: 0.5rem; align-items: center; margin-top: 0.3rem;",
                fileInput(ns("depmap_file"), NULL, accept = c(".RDS", ".rds"), width = "100%"),
                actionButton(ns("load_depmap_local"), "Load",
                             class = "btn-sm", style = "background:#869791; color:#fff; border-color:#869791; white-space:nowrap;")
              )
            ),
            div(style = "margin-top: 0.75rem;",
              uiOutput(ns("depmap_status"))
            )
          )
        )
      ),

      # Pre-trained Models
      column(6,
        div(class = "card animate-fade-in-up delay-2",
          div(class = "card-header",
            icon("cube"), " Pre-trained Models"
          ),
          div(class = "card-body",
            p(class = "text-muted", style = "font-size: 0.86rem;",
              "Load pre-trained drug response models from the PERCEPTION GitHub Release repository. 44 models are available, each trained on DepMap bulk expression with Elastic Net regression and 5-fold cross-validation. Models are cached locally after first download."),
            tags$small(class = "text-muted", style = "display: block; margin-top: 0.3rem;",
              "To download manually, visit ",
              tags$a(href = "https://github.com/SunPast/PERCEPTION/releases/tag/models-v1",
                     target = "_blank", "GitHub Release", style = "color: var(--primary); text-decoration: underline;"),
              "."
            ),
            hr(),
            div(class = "inline-form-row", style = "display: flex; gap: 0.5rem; align-items: center; flex-wrap: wrap;",
              tags$label(class = "control-label", style = "margin:0; font-size:0.88rem; font-weight:600; white-space:nowrap; line-height:38px;", "Drug"),
              selectizeInput(ns("model_name"), label = NULL, width = "300px",
                             choices = c("abemaciclib", "afatinib", "axitinib", "azacitidine", "cladribine",
                                         "clofarabine", "cobimetinib", "dabrafenib", "dasatinib", "daunorubicin",
                                         "decitabine", "docetaxel", "doxorubicin", "epirubicin", "erlotinib",
                                         "etoposide", "gefitinib", "gemcitabine", "homoharringtonine", "ibrutinib",
                                         "icotinib", "ixabepilone", "lapatinib", "lenvatinib", "midostaurin",
                                         "niraparib", "osimertinib", "paclitaxel", "palbociclib", "ponatinib",
                                         "romidepsin", "sunitinib", "temsirolimus", "teniposide", "thioguanine",
                                         "topotecan", "trametinib", "vandetanib", "vemurafenib", "vinblastine",
                                         "vincristine", "vindesine", "vinflunine", "vinorelbine"),
                             selected = "abemaciclib",
                             multiple = TRUE,
                             options = list(placeholder = "Select one or more drugs...", maxOptions = 50)),
              actionButton(ns("load_model"), "Download & Load",
                           class = "btn-primary", icon = icon("download")),
              checkboxInput(ns("model_mirror"), "Mirror", value = TRUE)
            ),
            div(style = "margin-top: 0.6rem; border-top: 1px dashed var(--border); padding-top: 0.6rem;",
              tags$small(class = "text-muted", "Or load a pre-downloaded model .RDS:"),
              div(style = "display: flex; gap: 0.5rem; align-items: center; margin-top: 0.3rem;",
                fileInput(ns("model_file"), NULL, accept = c(".RDS", ".rds"), width = "100%"),
                actionButton(ns("load_model_local"), "Load",
                             class = "btn-sm", style = "background:#869791; color:#fff; border-color:#869791; white-space:nowrap;")
              )
            ),
            div(style = "margin-top: 0.75rem;",
              uiOutput(ns("model_status"))
            )
          )
        )
      )
    ),

    # --- Loaded Models Management (placed after Data Management, before Upload Your Data) ---
    fluidRow(style = "margin-top: 1.5rem;",
      column(12,
        div(class = "card animate-fade-in-up",
          div(class = "card-header",
            icon("boxes-stacked"), " Loaded Models Management",
            tags$span(style = "margin-left: auto; font-size: 0.78rem; color: var(--text-muted);",
              span(class = "status-dot green", style = "display: inline-block; vertical-align: middle; margin-right: 0.2rem;"), "Active  ",
              span(class = "status-dot gray", style = "display: inline-block; vertical-align: middle; margin-right: 0.2rem;"), "Inactive"
            )
          ),
          div(class = "card-body", style = "padding: 0.8rem 1.2rem !important;",
            uiOutput(ns("models_management"))
          )
        )
      )
    ),

    # Upload Your Data Section
    fluidRow(style = "margin-top: 1.5rem;",
      column(12,
        div(class = "section-header",
          icon("upload"),
          h4("Upload Your Data")
        )
      )
    ),

    fluidRow(
      # Expression Matrix
      column(4,
        div(class = "card animate-fade-in-up delay-1",
          div(class = "card-header",
            icon("table"), " Expression Matrix"
          ),
          div(class = "card-body",
            p(class = "text-muted", style = "font-size: 0.82rem;",
              "Single-cell expression matrix (genes × cells). Raw counts or normalized values. CSV or RDS format."),
            div(class = "data-format-hint",
              icon("table"), " Format: rows = genes, columns = cells",
              tags$pre(class = "data-format-example", "        CELL_001  CELL_002  CELL_003
TP53      2.1       0.0       5.3
BRCA1     0.0       1.8       3.2
EGFR      4.7       2.1       0.0")
            ),
            fileInput(ns("expr_file"), "Upload Expression",
                      accept = c(".csv", ".rds", ".RDS"), width = "100%"),
            uiOutput(ns("expr_status"))
          )
        )
      ),

      # Patient-Cell Mapping
      column(4,
        div(class = "card animate-fade-in-up delay-2",
          div(class = "card-header",
            icon("users"), " Patient-Cell Mapping"
          ),
          div(class = "card-body",
            p(class = "text-muted", style = "font-size: 0.82rem;",
              "CSV with columns: ", tags$code("cell_id"), " and ", tags$code("patient_id"),
              ". Maps each cell to its patient. Clones will be auto-detected via Seurat clustering."),
            div(class = "data-format-hint",
              icon("users"), " Format: cell_id + patient_id",
              tags$pre(class = "data-format-example", "cell_id    patient_id
CELL_001   PAT_001
CELL_002   PAT_001
CELL_003   PAT_002")
            ),
            fileInput(ns("mapping_file"), "Upload Mapping",
                      accept = c(".csv"), width = "100%"),
            uiOutput(ns("mapping_status"))
          )
        )
      ),

      # Clinical Response
      column(4,
        div(class = "card animate-fade-in-up delay-3",
          div(class = "card-header",
            icon("heartbeat"), " Clinical Response"
          ),
          div(class = "card-body",
            p(class = "text-muted", style = "font-size: 0.82rem;",
              "Patient response data for evaluation. CSV with columns: patient, response (Responder/Non-responder)."),
            div(class = "data-format-hint",
              icon("heartbeat"), " Format: patient + response",
              tags$pre(class = "data-format-example", "patient    response
PAT_001    Responder
PAT_002    Non-responder
PAT_003    Responder")
            ),
            fileInput(ns("response_file"), "Upload Response",
                      accept = c(".csv"), width = "100%"),
            uiOutput(ns("response_status"))
          )
        )
      )
    ),

    # Seurat Clustering & Preprocessing
    fluidRow(style = "margin-top: 1rem;",
      column(12,
        div(class = "card animate-fade-in-up",
          div(class = "card-header",
            icon("shapes"), " Clone Detection & Preprocessing"
          ),
          div(class = "card-body seurat-body",
            # Left: description
            div(class = "seurat-left",
              p(class = "text-muted", style = "font-size: 0.85rem; line-height: 1.5; margin-bottom: 0.6rem;",
                "Run Seurat clustering to automatically detect transcriptional subclones, ",
                "compute clone-level mean expression, rank-normalize the data, and build ",
                "the clone abundance table required for prediction."
              ),
              p(class = "text-muted", style = "font-size: 0.8rem; line-height: 1.5; margin-bottom: 0;",
                strong("Method"), " — UMAP (default) preserves global structure and is faster with large datasets. ",
                "t-SNE emphasizes fine local neighborhoods and may reveal finer substructure at the cost of speed.",
                br(),
                strong("Resolution"), " controls clustering granularity — higher values produce more clones (finer subclones); ",
                "lower values produce fewer, broader clones. Default 0.8 suits most datasets.",
                br(),
                strong("PCA Dims"), " sets the number of principal components used for clustering — ",
                "higher values capture more biological signal but may include noise. Default 10 is standard for scRNA-seq."
              ),
              div(class = "info-box", style = "margin-top: 0.6rem; margin-bottom: 0; font-size: 0.8rem; padding: 0.5rem 0.7rem;",
                icon("info-circle"),
                "Requires Expression Matrix and Patient-Cell Mapping loaded first. ",
                "Seurat will cluster cells into subclones — you do not need to provide clone annotations manually."
              )
            ),
            # Right: controls
            div(class = "seurat-right",
              selectInput(ns("seurat_method"), "Reduction Method",
                          choices = c("UMAP" = "umap", "t-SNE" = "tsne"),
                          selected = "umap", width = "100%"),
              div(class = "seurat-params",
                numericInput(ns("seurat_resolution"), "Resolution",
                             value = 0.8, min = 0.1, max = 2, step = 0.1, width = "100%"),
                numericInput(ns("seurat_dims"), "PCA Dims",
                             value = 10, min = 2, max = 30, step = 1, width = "100%")
              ),
              actionButton(ns("run_seurat"), "Run Seurat Clustering",
                           class = "btn-primary seurat-run-btn",
                           icon = icon("wand-magic-sparkles"),
                           disabled = "disabled",
                           title = "Load expression matrix and patient-cell mapping first"),
              tags$script(HTML(paste0("
Shiny.addCustomMessageHandler('seurat-btn-state-", ns("run_seurat"), "', function(msg) {
  var btn = document.getElementById('", ns("run_seurat"), "');
  if (!btn) return;
  if (msg.enabled) {
    btn.removeAttribute('disabled');
    btn.removeAttribute('title');
  } else {
    btn.setAttribute('disabled', 'disabled');
    btn.setAttribute('title', 'Load expression matrix and patient-cell mapping first');
  }
});
")))
            )
          ),
          div(class = "seurat-status-bar",
            uiOutput(ns("seurat_status"))
          )
        )
      )
    ),

    # Data Overview & Preview (merged)
    fluidRow(style = "margin-top: 1.5rem;",
      column(12,
        div(class = "card animate-fade-in",
          div(class = "card-header",
            icon("eye"), " Data Overview & Preview"
          ),
          div(class = "card-body",
            # Status badges at top
            uiOutput(ns("status_overview")),
            hr(),
            tabsetPanel(
              tabPanel("Expression", DTOutput(ns("expr_preview"))),
              tabPanel("Clone Map", DTOutput(ns("clone_preview"))),
              tabPanel("Response", DTOutput(ns("response_preview"))),
              tabPanel("DepMap", DTOutput(ns("depmap_preview")))
            )
          )
        )
      )
    )
  )
}

mod_data_server <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- Enable/disable Run Seurat button based on required data (#10) ---
    observe({
      ready <- !is.null(shared$user_expr) && !is.null(shared$user_mapping)
      session$sendCustomMessage(
        paste0("seurat-btn-state-", ns("run_seurat")),
        list(enabled = ready)
      )
    })

    # --- Load Demo Data ---
    observeEvent(input$load_demo, {
      set.seed(42)
      gene_names <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "RB1",
                       "APC", "BRAF", "CDH1", "CDKN2A", "ERBB2", "FGFR1", "ALK",
                       "MET", "RET", "ROS1", "NRAS", "HRAS", "MAP2K1", "MAPK1",
                       "JAK2", "STAT3", "MTOR", "AKT1", "AKT2", "CTNNB1", "SMAD4",
                       "VHL", "NF1", "NF2", "STK11", "FBXW7", "ARID1A", "KDM5C",
                       "KMT2D", "SETD2", "BAP1", "PBRM1", "NOTCH1", "NOTCH2",
                       "NOTCH3", "JAK1", "JAK3", "SOX9", "IDH1", "IDH2", "FLT3")
      n_cells <- 400
      n_patients <- 25
      cell_names <- paste0("CELL_", sprintf("%04d", 1:n_cells))
      patient_names <- paste0("PAT_", sprintf("%03d", 1:n_patients))

      # Clinical response — defined FIRST so expression can be structured around it
      # 13 Responders + 12 Non-responders for meaningful box plots
      response_labels <- c(rep("Responder", 13), rep("Non-responder", 12))
      clinical_response <- data.frame(
        patient = patient_names,
        response = response_labels,
        stringsAsFactors = FALSE
      )
      shared$user_response <- clinical_response

      # Patient-Cell mapping
      patient_assignment <- rep(patient_names, each = ceiling(n_cells / n_patients))[1:n_cells]
      patient_mapping <- data.frame(
        cell_id = cell_names,
        patient_id = patient_assignment,
        stringsAsFactors = FALSE
      )
      shared$user_mapping <- patient_mapping

      # Build STRUCTURED expression so biomarker plots show meaningful correlation.
      # Two drug-biomarker groups with OPPOSITE patterns:
      #   - abemaciclib biomarkers (genes 1-5): HIGH in responders, LOW in non-responders
      #   - erlotinib biomarkers (genes 6-10): LOW in responders, HIGH in non-responders
      # Other genes are background noise.
      # This guarantees: when models trained on similarly structured data predict
      # on this expression, predictions correlate with the biomarker genes.
      is_responder_cell <- clinical_response$response[match(patient_assignment, clinical_response$patient)] == "Responder"
      abemaciclib_markers <- gene_names[1:5]   # TP53, BRCA1, EGFR, MYC, KRAS
      erlotinib_markers   <- gene_names[6:10]  # PIK3CA, PTEN, RB1, APC, BRAF
      noise_genes         <- gene_names[11:length(gene_names)]

      expr_matrix <- matrix(0.1, nrow = length(gene_names), ncol = n_cells)
      rownames(expr_matrix) <- gene_names
      colnames(expr_matrix) <- cell_names

      # abemaciclib markers: moderately HIGH in responders, LOW in non-responders
      for (g in abemaciclib_markers) {
        expr_matrix[g, is_responder_cell] <- pmax(rnorm(sum(is_responder_cell), mean = 8, sd = 3), 0.1)
        expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 3, sd = 2), 0.1)
      }
      # erlotinib markers: OPPOSITE pattern (LOW in responders, HIGH in non-responders)
      for (g in erlotinib_markers) {
        expr_matrix[g, is_responder_cell] <- pmax(rnorm(sum(is_responder_cell), mean = 3, sd = 2), 0.1)
        expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 8, sd = 3), 0.1)
      }
      # Noise genes: random uniform (no response signal)
      for (g in noise_genes) {
        expr_matrix[g, ] <- runif(n_cells, 0.5, 8)
      }
      storage.mode(expr_matrix) <- "numeric"
      shared$user_expr <- expr_matrix

      # Run prepare_data to get clone annotation, rank-normalized expression, clone counts
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
        prepared <- PERCEPTION::prepare_data(
          method = "umap",
          expression_matrix = expr_matrix,
          patient_mapping = patient_mapping,
          seurat_resolution = 0.8,
          seurat_dims = 10,
          seurat_nfeatures = min(2000, length(gene_names))
        )
        shared$prepared_data <- prepared
        shared$user_clones <- prepared$cell_clone_map
        # Keep shared$user_expr as cell-level expression (NOT clone-level) for biomarker plots
        w$hide()

        # Train models on STRUCTURED training data that mirrors the demo expression
        # pattern, so model features (informative genes) drive predictions.
        if (requireNamespace("caret", quietly = TRUE) && requireNamespace("glmnet", quietly = TRUE)) {
          # Generate structured training data with the same response-correlated pattern
          make_structured_training <- function(marker_genes, direction, n_train = 160) {
            # n_train samples: half "responder-like", half "non-responder-like"
            x_train <- matrix(0, nrow = length(gene_names), ncol = n_train)
            rownames(x_train) <- gene_names
            half <- n_train %/% 2
            responder_like <- seq_len(half)
            nonresponder_like <- seq.int(half + 1, n_train)

            # Marker genes follow the SAME pattern as demo expression
            for (g in marker_genes) {
              x_train[g, responder_like] <- pmax(rnorm(half, mean = 8, sd = 3), 0.1)
              x_train[g, nonresponder_like] <- pmax(rnorm(half, mean = 3, sd = 2), 0.1)
            }
            # Noise genes
            for (g in setdiff(gene_names, marker_genes)) {
              x_train[g, ] <- runif(n_train, 0.5, 8)
            }

            # y = signed mean of marker gene expression:
            #   direction = +1: high marker expr → high killing (responder sensitive)
            #   direction = -1: high marker expr → low killing (responder resistant)
            y_train <- direction * colMeans(x_train[marker_genes, , drop = FALSE])
            y_train <- y_train + rnorm(n_train, sd = 0.3)  # small noise

            train_df <- as.data.frame(cbind(y = y_train, t(x_train)))
            train_df
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

          # abemaciclib: high marker expression → high killing (responders sensitive)
          # erlotinib: high marker expression → low killing (responders resistant — opposite pattern)
          shared$models <- list(
            abemaciclib = make_drug_model("abemaciclib", 101, abemaciclib_markers, direction = +1),
            erlotinib   = make_drug_model("erlotinib",   202, erlotinib_markers,   direction = -1)
          )
          shared$model_cache <- shared$models
          shared$model_active <- list(abemaciclib = TRUE, erlotinib = TRUE)
        }

        n_clones <- ncol(prepared$clone_expression_rnorm)
        showNotification(paste0("Demo data loaded: ", length(gene_names), " genes × ", n_cells,
                                " cells, ", n_patients, " patients, ", n_clones,
                                " clones detected. Models and predictions ready!"),
                         type = "message", duration = 6)
      }, error = function(e) {
        w$hide()
        showNotification(paste("Demo data preparation error:", e$message), type = "error", duration = 10)
      })
    })

    # --- Load DepMap (download) ---
    observeEvent(input$load_depmap, {
      use_mirror <- isTRUE(input$depmap_mirror)
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Loading DepMap data..."),
          p(class = "text-muted", "This may take a few minutes (567MB download)")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      tryCatch({
        shared$depmap <- PERCEPTION::load_depmap(read = TRUE, mirror = use_mirror,
                                                   timeout_seconds = 600, retries = 2)
        w$hide()
        showNotification("DepMap data loaded successfully", type = "message")
      }, error = function(e) {
        w$hide()
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
        showNotification("Tip: Try enabling 'Use mirror' or download DepMap.RDS manually and use the 'Load' button below.", type = "warning", duration = 15)
      })
    })

    # --- Load DepMap (local file) ---
    observeEvent(input$load_depmap_local, {
      req(input$depmap_file)
      file <- input$depmap_file
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Reading DepMap.RDS..."),
          p(class = "text-muted", "This may take a few seconds")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      tryCatch({
        DepMap <- readRDS(file$datapath)
        depmap_env <- getFromNamespace(".depmap_env", "PERCEPTION")
        assign("DepMap", DepMap, envir = depmap_env)
        do.call("assign", list("DepMap", DepMap, envir = .GlobalEnv))
        shared$depmap <- DepMap
        w$hide()
        showNotification(paste("DepMap loaded from file:", length(DepMap), "datasets"), type = "message")
      }, error = function(e) {
        w$hide()
        showNotification(paste("Error reading file:", e$message), type = "error")
      })
    })

    output$depmap_status <- renderUI({
      if (is.null(shared$depmap)) {
        tagList(span(class = "status-badge unloaded", span(class = "status-dot gray"), "Not loaded"))
      } else {
        n_datasets <- length(shared$depmap)
        tagList(
          span(class = "status-badge loaded", span(class = "status-dot green"), "Loaded"),
          br(), br(),
          tags$small(class = "text-muted", paste(n_datasets, "datasets available"))
        )
      }
    })

    # --- Load Model (download) ---
    observeEvent(input$load_model, {
      req(input$model_name)
      drug_list <- trimws(input$model_name)
      drug_list <- drug_list[nzchar(drug_list)]
      if (length(drug_list) == 0) return()
      use_mirror <- isTRUE(input$model_mirror)
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4(sprintf("Loading %d model(s)...", length(drug_list))),
          p(class = "text-muted", "This may take a few seconds per drug")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      tryCatch({
        if (is.null(shared$models)) shared$models <- list()
        loaded_count <- 0
        failed_drugs <- c()
        for (drug in drug_list) {
          tryCatch({
            new_model <- PERCEPTION::load_model(drug, read = TRUE, mirror = use_mirror,
                                                 timeout_seconds = 120, retries = 2)
            shared$models[[drug]] <- new_model[[drug]]
            shared$model_cache[[drug]] <- new_model[[drug]]
            shared$model_active[[drug]] <- TRUE
            loaded_count <- loaded_count + 1
          }, error = function(e) {
            failed_drugs <<- c(failed_drugs, paste0(drug, " (", conditionMessage(e), ")"))
          })
        }
        w$hide()
        if (loaded_count > 0) {
          showNotification(sprintf("Loaded %d model(s): %s", loaded_count,
                                   paste(drug_list[1:min(loaded_count, length(drug_list))], collapse = ", ")),
                           type = "message", duration = 5)
        }
        if (length(failed_drugs) > 0) {
          showNotification(paste("Failed to load:", paste(failed_drugs, collapse = "; ")),
                           type = "error", duration = 10)
        }
      }, error = function(e) {
        w$hide()
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })

    # --- Load Model (local file) ---
    observeEvent(input$load_model_local, {
      req(input$model_file)
      file <- input$model_file
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Reading model..."),
          p(class = "text-muted", "This may take a few seconds")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      tryCatch({
        model_obj <- readRDS(file$datapath)
        drug_name <- tools::file_path_sans_ext(basename(file$name))
        if (!is.null(attr(model_obj, "drug_name"))) {
          drug_name <- attr(model_obj, "drug_name")
        }
        if (is.null(shared$models)) shared$models <- list()
        shared$models[[drug_name]] <- model_obj
        shared$model_cache[[drug_name]] <- model_obj
        shared$model_active[[drug_name]] <- TRUE
        w$hide()
        showNotification(paste("Model loaded from file:", drug_name), type = "message")
      }, error = function(e) {
        w$hide()
        showNotification(paste("Error reading file:", e$message), type = "error")
      })
    })

    output$model_status <- renderUI({
      if (is.null(shared$models)) {
        span(class = "status-badge unloaded", span(class = "status-dot gray"), "No model loaded")
      } else {
        n_drugs <- length(shared$models)
        tagList(
          span(class = "status-badge loaded", span(class = "status-dot green"), "Loaded"),
          br(), br(),
          tags$small(class = "text-muted", paste(n_drugs, "drug model(s)"))
        )
      }
    })

    # --- Models Management: status badges + toggle buttons ---
    output$models_management <- renderUI({
      cache <- shared$model_cache
      if (is.null(cache) || length(cache) == 0) {
        div(class = "text-muted", style = "text-align: center; padding: 1.5rem;",
          icon("inbox", style = "font-size: 1.8rem; opacity: 0.3; display: block; margin-bottom: 0.5rem;"),
          "No models in cache. Use the Pre-trained Models card above to load a model,",
          "or train one in the Train tab."
        )
      } else {
        drug_names <- names(cache)
        cards <- lapply(drug_names, function(d) {
          btn_id <- paste0("toggle_", gsub("[^A-Za-z0-9]", "_", d))
          is_active <- isTRUE(shared$model_active[[d]])
          state_class <- if (is_active) "active" else "inactive"
          dot_class <- if (is_active) "green" else "gray"
          label <- if (is_active) "Active" else "Inactive"
          div(class = "model-mgmt-card",
            div(class = "model-mgmt-status",
              strong(d)
            ),
            actionButton(ns(btn_id), label,
                         class = paste("btn-model-toggle btn-sm", state_class),
                         icon = tags$span(class = "status-dot", class = dot_class,
                                          style = "width:8px;height:8px;border-radius:50%;display:inline-block;"))
          )
        })
        div(class = "model-mgmt-grid", cards)
      }
    })

    # Toggle observers for all 44 drugs
    drug_list <- c("abemaciclib", "afatinib", "axitinib", "azacitidine", "cladribine",
                   "clofarabine", "cobimetinib", "dabrafenib", "dasatinib", "daunorubicin",
                   "decitabine", "docetaxel", "doxorubicin", "epirubicin", "erlotinib",
                   "etoposide", "gefitinib", "gemcitabine", "homoharringtonine", "ibrutinib",
                   "icotinib", "ixabepilone", "lapatinib", "lenvatinib", "midostaurin",
                   "niraparib", "osimertinib", "paclitaxel", "palbociclib", "ponatinib",
                   "romidepsin", "sunitinib", "temsirolimus", "teniposide", "thioguanine",
                   "topotecan", "trametinib", "vandetanib", "vemurafenib", "vinblastine",
                   "vincristine", "vindesine", "vinflunine", "vinorelbine")
    for (d in drug_list) {
      local({
        drug_local <- d
        btn_id <- paste0("toggle_", gsub("[^A-Za-z0-9]", "_", drug_local))
        observeEvent(input[[btn_id]], {
          cat(file = stderr(), paste0("[DEBUG] toggle fired for: ", drug_local, "\n"))
          is_active <- isTRUE(shared$model_active[[drug_local]])
          if (is_active) {
            # Deactivate: remove from active models but keep in cache
            shared$models[[drug_local]] <- NULL
            shared$model_active[[drug_local]] <- FALSE
            showNotification(paste("Model deactivated:", drug_local, "(kept in cache)"), type = "message")
          } else {
            # Activate: restore from cache
            if (!is.null(shared$model_cache[[drug_local]])) {
              if (is.null(shared$models)) shared$models <- list()
              shared$models[[drug_local]] <- shared$model_cache[[drug_local]]
              shared$model_active[[drug_local]] <- TRUE
              showNotification(paste("Model activated:", drug_local), type = "message")
            }
          }
        }, ignoreInit = TRUE)
      })
    }

    # --- Upload Expression ---
    observeEvent(input$expr_file, {
      file <- input$expr_file
      tryCatch({
        if (grepl("\\.rds$|\\.RDS$|\\.Rds$", file$name)) {
          mat <- readRDS(file$datapath)
        } else {
          mat <- as.matrix(read.csv(file$datapath, row.names = 1, check.names = FALSE))
        }
        # Ensure numeric (Seurat LogMap requires numeric, not integer)
        if (is.integer(mat)) mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat),
                                           dimnames = dimnames(mat))
        shared$user_expr <- mat
        shared$prepared_data <- NULL  # Reset prepared data when expression changes
        showNotification(paste("Expression matrix loaded:", nrow(mat), "genes x", ncol(mat), "cells"), type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })

    output$expr_status <- renderUI({
      if (is.null(shared$user_expr)) return(NULL)
      tagList(
        span(class = "status-badge loaded", span(class = "status-dot green"), "Loaded"),
        tags$small(class = "text-muted",
          paste0(nrow(shared$user_expr), " genes x ", ncol(shared$user_expr), " cells")
        )
      )
    })

    # --- Upload Patient-Cell Mapping ---
    observeEvent(input$mapping_file, {
      file <- input$mapping_file
      tryCatch({
        shared$user_mapping <- read.csv(file$datapath, stringsAsFactors = FALSE)
        shared$prepared_data <- NULL  # Reset prepared data when mapping changes
        showNotification(paste("Patient-cell mapping loaded:", nrow(shared$user_mapping), "cells"), type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })

    output$mapping_status <- renderUI({
      if (is.null(shared$user_mapping)) return(NULL)
      n_patients <- length(unique(shared$user_mapping$patient_id))
      tagList(
        span(class = "status-badge loaded", span(class = "status-dot green"), "Loaded"),
        tags$small(class = "text-muted",
          paste0(nrow(shared$user_mapping), " cells, ", n_patients, " patients")
        )
      )
    })

    # --- Run Seurat Clustering (prepare_data) ---
    observeEvent(input$run_seurat, {
      req(shared$user_expr, shared$user_mapping)
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Running Seurat clustering..."),
          p(class = "text-muted", "Detecting subclones and preparing data")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      tryCatch({
        prepared <- PERCEPTION::prepare_data(
          method = input$seurat_method,
          expression_matrix = shared$user_expr,
          patient_mapping = shared$user_mapping,
          seurat_resolution = input$seurat_resolution,
          seurat_dims = input$seurat_dims
        )
        shared$prepared_data <- prepared
        shared$user_clones <- prepared$cell_clone_map
        # Keep shared$user_expr as cell-level expression for biomarker plots
        w$hide()
        showNotification(paste0("Clustering complete: ", ncol(prepared$clone_expression_rnorm),
                                " clones detected across ", nrow(prepared$clone_counts), " patients"),
                         type = "message", duration = 5)
      }, error = function(e) {
        w$hide()
        showNotification(paste("Clustering error:", e$message), type = "error", duration = 10)
      })
    })

    output$seurat_status <- renderUI({
      if (is.null(shared$prepared_data)) {
        tags$span(class = "status-badge unloaded",
          span(class = "status-dot gray"),
          "Not run yet — upload expression + mapping, then click Run Seurat Clustering"
        )
      } else {
        pd <- shared$prepared_data
        tagList(
          span(class = "status-badge loaded", span(class = "status-dot green"), "Complete"),
          tags$small(class = "text-muted",
            paste0(ncol(pd$clone_expression_rnorm), " clones, ",
                   nrow(pd$clone_counts), " patients, ",
                   nrow(pd$clone_expression_rnorm), " genes (rank-normalized)")
          )
        )
      }
    })

    # --- Upload Clinical Response ---
    observeEvent(input$response_file, {
      file <- input$response_file
      tryCatch({
        shared$user_response <- read.csv(file$datapath, stringsAsFactors = FALSE)
        showNotification(paste("Clinical response loaded:", nrow(shared$user_response), "patients"), type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })

    output$response_status <- renderUI({
      if (is.null(shared$user_response)) return(NULL)
      tagList(
        span(class = "status-badge loaded", span(class = "status-dot green"), "Loaded"),
        tags$small(class = "text-muted", paste(nrow(shared$user_response), "patients"))
      )
    })

    # --- Status Overview ---
    output$status_overview <- renderUI({
      items <- list(
        list(name = "Expression Matrix", loaded = !is.null(shared$user_expr)),
        list(name = "Patient-Cell Mapping", loaded = !is.null(shared$user_mapping)),
        list(name = "Clone Map (Seurat)", loaded = !is.null(shared$prepared_data)),
        list(name = "Clinical Response", loaded = !is.null(shared$user_response)),
        list(name = "Trained Model", loaded = !is.null(shared$models)),
        list(name = "DepMap Data", loaded = !is.null(shared$depmap))
      )

      tagList(
        div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 0.8rem;",
          lapply(items, function(item) {
            div(class = if (item$loaded) "status-badge loaded" else "status-badge unloaded",
              style = "justify-content: center; padding: 0.5rem 0.8rem; font-size: 0.82rem;",
              span(class = if (item$loaded) "status-dot green" else "status-dot gray"),
              strong(item$name)
            )
          })
        )
      )
    })

    # --- Data Previews ---
    output$depmap_preview <- renderDT({
      req(shared$depmap)
      summary_df <- data.frame(
        Dataset = names(shared$depmap),
        Rows = sapply(shared$depmap, function(x) nrow(x)),
        Cols = sapply(shared$depmap, function(x) ncol(x))
      )
      datatable(summary_df, options = list(pageLength = 10, dom = "tp"),
                rownames = FALSE, class = "display")
    })

    output$expr_preview <- renderDT({
      req(shared$user_expr)
      datatable(shared$user_expr[, 1:min(10, ncol(shared$user_expr)), drop = FALSE],
                options = list(pageLength = 10, dom = "tp", scrollX = TRUE),
                class = "display")
    })

    output$clone_preview <- renderDT({
      req(shared$user_clones)
      datatable(shared$user_clones, options = list(pageLength = 10, dom = "tp"),
                rownames = FALSE, class = "display")
    })

    output$response_preview <- renderDT({
      req(shared$user_response)
      datatable(shared$user_response, options = list(pageLength = 10, dom = "tp"),
                rownames = FALSE, class = "display")
    })

  })
}
