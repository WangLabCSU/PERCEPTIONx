# Shared helpers for the Shiny app — sourced BEFORE the modules in app.R.
# Anything used by more than one module lives here so the entry points cannot
# drift apart (e.g. the Home "Load Demo" button and the Data-tab demo button).

# ---------------------------------------------------------------------------
# Demo data pipeline (single source of truth)
#
# Builds the demo expression matrix, runs Seurat clustering via prepare_data(),
# trains the two demo drug models, and populates `shared`. Both the Home page
# "Load Demo" button and the Data-tab "Load Demo Data" button call this.
#
# @param shared     The app's shared reactiveValues.
# @param on_success Optional function called after everything is loaded
#                   (e.g. mark the demo active, navigate to a tab).
#
# @return Invisibly, a list with the loaded counts (genes/cells/patients/clones).
# ---------------------------------------------------------------------------
run_demo_pipeline <- function(shared, on_success = NULL) {
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

  # Clinical response — defined FIRST so expression can be structured around it.
  response_labels <- c(rep("Responder", 10), rep("Non-responder", 10))
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

  # STRUCTURED expression so biomarker plots show meaningful correlation.
  # Two drug-biomarker groups with OPPOSITE patterns:
  #   - abemaciclib biomarkers (genes 1-5): HIGH in responders, LOW in non-responders
  #   - erlotinib biomarkers (genes 6-10): LOW in responders, HIGH in non-responders
  is_responder_cell <- clinical_response$response[match(patient_assignment, clinical_response$patient)] == "Responder"
  abemaciclib_markers <- gene_names[1:5]
  erlotinib_markers   <- gene_names[6:10]
  noise_genes         <- gene_names[11:length(gene_names)]

  expr_matrix <- matrix(0.1, nrow = length(gene_names), ncol = n_cells)
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- cell_names

  for (g in abemaciclib_markers) {
    expr_matrix[g, is_responder_cell] <- pmax(rnorm(sum(is_responder_cell), mean = 8, sd = 3), 0.1)
    expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 3, sd = 2), 0.1)
  }
  for (g in erlotinib_markers) {
    expr_matrix[g, is_responder_cell] <- pmax(rnorm(sum(is_responder_cell), mean = 3, sd = 2), 0.1)
    expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 8, sd = 3), 0.1)
  }
  for (g in noise_genes) {
    expr_matrix[g, ] <- runif(n_cells, 0.5, 8)
  }
  storage.mode(expr_matrix) <- "numeric"
  shared$user_expr <- expr_matrix

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
    shared$user_clones <- prepared$cell_clone_map
    # shared$user_expr stays CELL-level (not clone-level) for spatial gene plots.
    w$hide()

    # Train demo models on structured data mirroring the demo expression.
    if (requireNamespace("caret", quietly = TRUE) && requireNamespace("glmnet", quietly = TRUE)) {
      make_structured_training <- function(marker_genes, direction, n_train = 160) {
        x_train <- matrix(0, nrow = length(gene_names), ncol = n_train)
        rownames(x_train) <- gene_names
        half <- n_train %/% 2
        responder_like <- seq_len(half)
        nonresponder_like <- seq.int(half + 1, n_train)

        for (g in marker_genes) {
          x_train[g, responder_like] <- pmax(rnorm(half, mean = 8, sd = 3), 0.1)
          x_train[g, nonresponder_like] <- pmax(rnorm(half, mean = 3, sd = 2), 0.1)
        }
        for (g in setdiff(gene_names, marker_genes)) {
          x_train[g, ] <- runif(n_train, 0.5, 8)
        }
        # y = signed mean marker expression, interpreted as VIABILITY
        #   (high = resistant, low = sensitive).
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
        abemaciclib = make_drug_model("abemaciclib", 101, abemaciclib_markers, direction = -1),
        erlotinib   = make_drug_model("erlotinib",   202, erlotinib_markers,   direction = +1)
      )
      shared$model_cache <- shared$models
      shared$model_active <- list(abemaciclib = TRUE, erlotinib = TRUE)
    }

    n_clones <- ncol(prepared$clone_expression_rnorm)
    if (!is.null(on_success)) on_success()
    showNotification(paste0("Demo data loaded: ", length(gene_names), " genes × ", n_cells,
                            " cells, ", n_patients, " patients, ", n_clones,
                            " clones detected. Models ready!"),
                     type = "message", duration = 6)
    invisible(list(genes = length(gene_names), cells = n_cells,
                   patients = n_patients, clones = n_clones))
  }, error = function(e) {
    w$hide()
    showNotification(paste("Demo data preparation error:", e$message), type = "error", duration = 10)
    invisible(NULL)
  })
}
