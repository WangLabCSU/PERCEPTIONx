# Shared helpers for the Shiny app — sourced BEFORE the modules in app.R.
# Anything used by more than one module lives here so the entry points cannot
# drift apart (e.g. the Home "Load Demo" button and the Data-tab demo button).

# ---------------------------------------------------------------------------
# Workflow stepper HTML (Home page)
#
# Pure renderer: builds the pipeline stepper + hint from a plain list of
# shared state (reactiveValuesToList(shared) or the initial empty list).
# It is called BOTH at UI build time (static first paint — no server round
# trip, so the section shows instantly) and on the server whenever shared
# state changes (HTML is pushed to the browser via the 'set-html' handler).
#
# @param shared A list mirroring the app's shared reactiveValues.
# @param ns     Namespace function for the actionLink ids (identity when
#               no namespacing is needed).
# ---------------------------------------------------------------------------
workflow_stepper_html <- function(shared, ns = identity) {
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
}

# ---------------------------------------------------------------------------
# Data status dashboard HTML (Home page)
#
# Same pattern as workflow_stepper_html(): pure renderer used both for the
# static first paint and for subsequent server-side updates.
# ---------------------------------------------------------------------------
data_dashboard_html <- function(shared) {
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
}

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
