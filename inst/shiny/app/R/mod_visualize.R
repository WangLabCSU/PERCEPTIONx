# Visualization Module — Enhanced
mod_visualize_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        div(class = "section-header",
          icon("chart-line"),
          h4("Visualization")
        )
      )
    ),

    # Prerequisites Check
    fluidRow(
      column(12,
        uiOutput(ns("prereq_check"))
      )
    ),

    fluidRow(
      # Sidebar - Plot Selection + Data Status
      column(3, class = "viz-sidebar-col",
        div(class = "card animate-fade-in-up viz-plot-type-card",
          div(class = "card-header",
            icon("palette"), " Plot Type"
          ),
          div(class = "card-body viz-plot-type",
            radioButtons(ns("plot_type"), "",
                         choices = c(
                           "Clone Distribution" = "clone_dist",
                           "Clone Viability Lollipop" = "clone_kill",
                           "ROC Curve" = "roc",
                           "Response Boxplot" = "boxplot"
                         ),
                         selected = "clone_dist"),

            conditionalPanel(
              condition = paste0("input['", ns("plot_type"), "'] == 'clone_kill' || input['", ns("plot_type"), "'] == 'boxplot' || input['", ns("plot_type"), "'] == 'roc'"),
              selectizeInput(ns("drug_name_common"), "Drug Name",
                             choices = NULL, selected = NULL, width = "100%",
                             options = list(maxItems = 1, placeholder = "Select a drug"))
            ),

            # ROC comparison pair — shown ONLY when the ROC plot type is
            # selected (and the response has >2 groups; see server side).
            conditionalPanel(
              condition = paste0("input['", ns("plot_type"), "'] == 'roc'"),
              uiOutput(ns("roc_pair_picker"))
            ),

            actionButton(ns("generate"), "Generate Plot", width = "100%",
                         class = "btn-primary btn-sm", icon = icon("wand-magic-sparkles")),

            hr(),

            div(class = "viz-advanced-section",
              h6(class = "viz-advanced-title", icon("map"), " Spatial Visualizations"),
              p(class = "viz-advanced-desc", textOutput(ns("spatial_desc"), inline = TRUE)),

              radioButtons(ns("plot_type_advanced"), NULL,
                           choices = c(
                             "Clone Identity" = "umap_clone",
                             "Drug Viability" = "umap_viability",
                             "Gene Expression" = "umap_gene"
                           ),
                           selected = character(0)),

              conditionalPanel(
                condition = paste0("input['", ns("plot_type_advanced"), "'] == 'umap_gene'"),
                selectizeInput(ns("umap_gene"), "Gene",
                               choices = NULL, selected = NULL, width = "100%",
                               options = list(maxItems = 1, placeholder = "Select a gene"))
              ),
              conditionalPanel(
                condition = paste0("input['", ns("plot_type_advanced"), "'] == 'umap_viability'"),
                selectizeInput(ns("umap_drug"), "Drug",
                               choices = NULL, selected = NULL, width = "100%",
                               options = list(maxItems = 1, placeholder = "Select a drug"))
              ),
              conditionalPanel(
                condition = paste0("input['", ns("plot_type_advanced"), "'] == 'umap_gene' || input['", ns("plot_type_advanced"), "'] == 'umap_viability' || input['", ns("plot_type_advanced"), "'] == 'umap_clone'"),
                actionButton(ns("generate_advanced"), "Generate Spatial Plot", width = "100%",
                             class = "btn-primary btn-sm", icon = icon("wand-magic-sparkles"))
              )
            ),
            
            br(),

            # Plot size controls — at the bottom so the spatial plots above
            # are not crowded out.
            div(class = "viz-size-controls",
              tags$label(class = "viz-size-label", icon("arrows-alt"), " Plot Size"),
              div(class = "viz-size-row",
                div(class = "viz-size-col",
                  tags$span(class = "viz-size-mini-label", "Width"),
                  numericInput(ns("plot_width"), NULL, value = 100, min = 50, max = 200, step = 10)
                ),
                div(class = "viz-size-col",
                  tags$span(class = "viz-size-mini-label", "Height"),
                  numericInput(ns("plot_height"), NULL, value = 100, min = 50, max = 200, step = 10)
                )
              ),
              div(class = "viz-size-row",
                sliderInput(ns("plot_text_size"), "Text size", min = 60, max = 160, value = 100, step = 10, post = "%", ticks = TRUE)
              )
            )
          )
        ),

        # Data Availability Summary — fixed 300px
        div(class = "card animate-fade-in-up delay-1 viz-data-status-card",
          div(class = "card-header",
            icon("clipboard-check"), " Data Status"
          ),
          div(class = "card-body", style = "font-size: 0.82rem; padding: 0.7rem 1rem !important; overflow-y: auto;",
            uiOutput(ns("data_status"))
          )
        )
      ),

      # Main Plot Area — aligned to top of container
      column(9, class = "viz-main-col",
        div(class = "card animate-fade-in-up delay-1 viz-output-card",
          div(class = "card-header",
            icon("image"), " Output",
            tags$span(class = "viz-output-header-info",
              uiOutput(ns("plot_info_inline")),
              tags$span(style = "margin-left: 0.5rem;",
                uiOutput(ns("plot_download_inline"))
              )
            )
          ),
          div(class = "card-body",
            uiOutput(ns("plot_status")),
            ggiraph::girafeOutput(ns("main_plot"))
          )
        ),
        # Plot explanation card
        uiOutput(ns("plot_explanation"))
      )
    )
  )
}

mod_visualize_server <- function(id, shared, main_session) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    current_plot <- reactiveVal(NULL)
    # Track which plot type was last generated (defined here, before use).
    current_plot_type <- reactiveVal(NULL)

    # Map response labels to canonical short classes: R/NR spellings collapse
    # to R/NR; ANY other label (e.g. longitudinal time points TN/RD/PD) is kept
    # as-is (uppercased) so multi-group response data is not silently dropped.
    label_resp <- function(x) {
      y <- tolower(trimws(as.character(x)))
      y[y %in% c("responder", "response", "responsive", "r", "sensitive", "sensitivity")] <- "R"
      y[y %in% c("non-responder", "nonresponder", "non responder", "non-responsive",
                 "nonresponsive", "nr", "resistant", "resistance", "progressor",
                 "progression", "non", "non_responder")] <- "NR"
      keep <- !y %in% c("R", "NR")
      y[keep] <- toupper(y[keep])
      y
    }

    # Distinct response groups in the loaded clinical response (order of
    # appearance): "R"/"NR" for two-class data, or TN/RD/PD for longitudinal
    # time-point data, etc.
    response_groups <- reactive({
      cr <- shared$user_response
      if (is.null(cr) || is.null(cr$response)) return(character(0))
      grps <- unique(label_resp(cr$response))
      grps[!is.na(grps) & nzchar(trimws(grps)) & grps != "NA"]
    })

    # ROC comparison pair picker: only shown when the response has >2 groups,
    # since a ROC needs exactly two classes. Defaults to PD vs RD (progressed
    # vs residual/responding) matching the paper's Fig. 4b semantics; TN
    # (treatment-naive baseline) is excluded from that comparison.
    output$roc_pair_picker <- renderUI({
      grps <- response_groups()
      if (length(grps) <= 2) return(NULL)
      def <- if (all(c("PD", "RD") %in% grps)) c("PD", "RD") else grps[1:2]
      tagList(
        tags$small(class = "text-muted", style = "display: block; margin: 0.4rem 0 0.2rem;",
          icon("info-circle"), " Response has ", length(grps), " groups — pick the two for the ROC:"),
        div(class = "viz-size-row",
          div(class = "viz-size-col",
            # selectize = FALSE: a native <select> keeps a constant width,
            # unlike the selectize widget whose pill makes the box change
            # width before/after a selection.
            selectInput(ns("roc_group_a"), "ROC Group A",
                        choices = grps, selected = def[1], width = "100%",
                        selectize = FALSE)),
          div(class = "viz-size-col",
            selectInput(ns("roc_group_b"), "ROC Group B",
                        choices = grps, selected = def[2], width = "100%",
                        selectize = FALSE))
        )
      )
    })

    # ---------------------------------------------------------------------
    # Combination (multi-drug) helpers — replicate the paper's pipeline
    # (Fig. 2b–e; Methods "Testing prediction strategies for multiple
    # myeloma", code Step3_Figure3.Rmd):
    #   clone level: z-score EACH drug's viability across ALL clones (global),
    #                then per-clone combination = pmin over drugs (IDA
    #                principle: the single most effective drug dominates).
    #   patient level (strategy 5, the one the paper chose): most-resistant
    #                clone weighted by its abundance = max(comb * weight).
    # ---------------------------------------------------------------------
    combo_clone_frame <- function() {
      pred_mat <- shared$predictions
      if (is.null(pred_mat) || ncol(pred_mat) < 1) {
        stop("No clone-level predictions to combine.")
      }
      pred_mat <- as.matrix(pred_mat)
      # Global z-score per drug across all clones (sd == 0 columns stay 0).
      z_mat <- pred_mat
      for (j in seq_len(ncol(pred_mat))) {
        col <- pred_mat[, j]
        s <- stats::sd(col, na.rm = TRUE)
        z_mat[, j] <- if (is.na(s) || s == 0) rep(0, nrow(pred_mat)) else
          (col - mean(col, na.rm = TRUE)) / s
      }
      comb <- apply(z_mat, 1, min, na.rm = TRUE)

      # Attach patient / clone_id (row order follows the prediction matrix,
      # which matches the clone_viability_template built by prepare_data()).
      if (!is.null(shared$prepared_data$clone_viability_template)) {
        tmpl <- shared$prepared_data$clone_viability_template
        clone_viability_df <- data.frame(patient = tmpl$patient,
                                         clone_id = tmpl$clone_id,
                                         comb_viability = unname(comb),
                                         stringsAsFactors = FALSE)
      } else {
        parsed <- PERCEPTIONx::parse_clone_keys(names(comb))
        clone_viability_df <- data.frame(patient = parsed$patient,
                                         clone_id = parsed$clone_id,
                                         comb_viability = unname(comb),
                                         stringsAsFactors = FALSE)
      }

      # Clone abundance = real cell-count proportion per patient.
      clone_data <- shared$user_clones
      if (!is.null(clone_data) && nrow(clone_data) > 0) {
        clone_viability_df$weights <- vapply(seq_len(nrow(clone_viability_df)), function(i) {
          pat <- clone_viability_df$patient[i]
          cl  <- clone_viability_df$clone_id[i]
          n_p <- sum(clone_data$patient == pat, na.rm = TRUE)
          if (n_p == 0) NA_real_ else
            sum(clone_data$patient == pat & clone_data$clone_id == cl, na.rm = TRUE) / n_p
        }, numeric(1))
      } else {
        clone_viability_df$weights <- NA_real_
      }

      if (!is.null(shared$user_response)) {
        clone_viability_df$response <- label_resp(shared$user_response$response[
          match(clone_viability_df$patient, shared$user_response$patient)
        ])
      }
      clone_viability_df
    }

    combo_patient_frame <- function(clone_viability_df) {
      patients <- unique(clone_viability_df$patient)
      scores <- vapply(patients, function(p) {
        sub <- clone_viability_df[clone_viability_df$patient == p, ]
        sub <- sub[!is.na(sub$comb_viability) & !is.na(sub$weights), ]
        if (nrow(sub) == 0) return(NA_real_)
        max(sub$comb_viability * sub$weights)  # most-resistant clone, weighted
      }, numeric(1))
      data.frame(patient = patients, combination = scores, stringsAsFactors = FALSE)
    }


    # Reactive plot size (px) — driven by user's Width/Height size controls.
    # Used by the static renderPlot (crisp original ggplot, no plotly overlap).
    plot_size <- reactive({
      w_pct <- if (is.null(input$plot_width)) 100 else as.numeric(input$plot_width)
      h_pct <- if (is.null(input$plot_height)) 100 else as.numeric(input$plot_height)
      list(w = round(1000 * w_pct / 100), h = round(750 * h_pct / 100))
    })

    # Reactive text size scale (50-160%)
    text_scale <- reactive({
      if (is.null(input$plot_text_size)) 100 else as.numeric(input$plot_text_size) / 100
    })

    # Reactive reduction method label (UMAP or t-SNE)
    reduction_method <- reactive({
      m <- shared$prepared_data$reduction_method
      if (is.null(m)) "umap" else m
    })
    reduction_label <- reactive({
      m <- reduction_method()
      if (m == "none") "no embedding" else if (m == "tsne") "t-SNE" else "UMAP"
    })

    # Normalize embedding coordinate columns — old prepare_data() returns
    # umap_1/umap_2, new returns dim_1/dim_2. Accept either.
    get_embedding_xy <- function(coords, cell_ids) {
      if (!is.data.frame(coords)) stop("umap_coords is not a data frame")
      if (nrow(coords) == 0L) stop("umap_coords has 0 rows — re-run Seurat clustering")
      if (!"cell_id" %in% names(coords)) stop("umap_coords has no 'cell_id' column")
      x_col <- intersect(c("dim_1", "umap_1"), names(coords))[1]
      y_col <- intersect(c("dim_2", "umap_2"), names(coords))[1]
      if (is.na(x_col) || is.na(y_col))
        stop("Embedding coordinate columns not found in umap_coords (names: ",
             paste(names(coords), collapse = ", "), ")")
      idx <- match(cell_ids, coords$cell_id)
      if (all(is.na(idx)))
        stop("0 matching cell IDs between query and umap_coords")
      list(X = coords[[x_col]][idx], Y = coords[[y_col]][idx])
    }

    # Safe scale: 0—1 using percentiles. Falls back to min-max if range01 fails.
    safe_range01 <- function(x) {
      if (length(x) == 0L) stop("empty input to safe_range01")
      r <- tryCatch(PERCEPTIONx::range01(x), error = function(e) NULL)
      if (is.null(r) || length(r) != length(x)) {
        rng <- range(x, na.rm = TRUE)
        if (rng[2] == rng[1]) return(rep(0.5, length(x)))
        (x - rng[1]) / (rng[2] - rng[1])
      } else {
        r
      }
    }

    output$spatial_desc <- renderText({
      paste0("Uses ", reduction_label(), " coordinates from Seurat clustering. No extra files needed.")
    })

    # Prerequisites Check
    output$prereq_check <- renderUI({
      missing <- c()
      if (is.null(shared$predictions) && is.null(shared$models) &&
          is.null(shared$user_clones) && is.null(shared$user_expr)) {
        missing <- c(missing, "No data loaded")
      }

      if (length(missing) > 0) {
        div(class = "info-box animate-fade-in", style = "border-left-color: var(--accent); margin-bottom: 0.4rem;",
          icon("exclamation-triangle"), " ",
          strong("No visualization data available."),
          " Please load data and run predictions first. ",
          actionButton(ns("go_data"), "Go to Data", class = "btn-accent btn-sm", style = "margin-left: 0.5rem;")
        )
      } else {
        div(class = "info-box animate-fade-in", style = "border-left-color: var(--success); margin-bottom: 0.4rem;",
          icon("check-circle"), " Data available for visualization. Select a plot type and click Generate."
        )
      }
    })

    # Data Status Summary
    output$data_status <- renderUI({
      status_items <- list(
        list(label = "Expression", cond = !is.null(shared$user_expr)),
        list(label = "Clones", cond = !is.null(shared$user_clones)),
        list(label = "Model", cond = !is.null(shared$models)),
        list(label = "Predictions", cond = !is.null(shared$predictions)),
        list(label = "Clinical Response", cond = !is.null(shared$user_response))
      )

      tagList(
        lapply(status_items, function(item) {
          div(class = "status-row",
            span(class = paste("status-dot", if (item$cond) "green" else "gray")),
            sprintf(" %s: ", item$label),
            if (item$cond) tags$span(class = "status-badge loaded", "Loaded")
            else tags$span(class = "status-badge unloaded", "Missing")
          )
        })
      )
    })

    # Navigation
    observeEvent(input$go_data, bslib::nav_select("navbar", selected = "data", session = main_session))

    # Plot type requirements info
    plot_requirements <- list(
      clone_dist = c("user_clones"),
      clone_kill = c("predictions", "user_clones"),
      roc = c("patient_pred", "user_response"),
      boxplot = c("patient_pred", "user_response"),
      umap_gene = c("predictions", "user_expr"),
      umap_viability = c("predictions"),
      umap_clone = c("user_clones")
    )

    plot_labels <- list(
      clone_dist = "Clone Distribution",
      clone_kill = "Clone Viability Lollipop",
      roc = "ROC Curve",
      boxplot = "Response Boxplot",
      umap_gene = "Gene Expression",
      umap_viability = "Drug Viability",
      umap_clone = "Clone Identity"
    )
    # For display, prefix with method label
    spatial_plot_label <- function(pt) {
      paste(reduction_label(), plot_labels[[pt]])
    }

    # Populate drug choices from trained models (or predictions as fallback).
    # "Combination" is prepended and used as the DEFAULT selection: it applies
    # the paper's IDA principle (per-clone combination viability = pmin of the
    # z-scored single-drug viabilities) and the most-resistant-clone weighted
    # aggregation (weighted_max) for patient-level scores.
    observe({
      drug_choices <- NULL
      if (!is.null(shared$models)) {
        drug_choices <- names(shared$models)
      } else if (!is.null(shared$predictions)) {
        drug_choices <- colnames(shared$predictions)
      }
      if (!is.null(drug_choices) && length(drug_choices) > 0) {
        combo_choices <- c("Combination" = "Combination",
                           setNames(drug_choices, drug_choices))
        updateSelectizeInput(session, "drug_name_common",
                             choices = combo_choices, selected = "Combination",
                             server = TRUE)
        # Spatial UMAP plots remain per-drug only (no combination there).
        updateSelectizeInput(session, "umap_drug", choices = drug_choices, server = TRUE)
      }
    })

    # Populate gene choices when expression data is available
    observe({
      expr_mat <- shared$user_expr
      if (!is.null(expr_mat)) {
        updateSelectizeInput(session, "umap_gene", choices = rownames(expr_mat), server = TRUE)
      }
    })

    # Generate Plot
    observeEvent(input$generate, {
      pt <- input$plot_type
      reqs <- plot_requirements[[pt]]

      # Combination mode derives patient-level scores from the clone-level
      # predictions + clone annotation, so boxplot/ROC no longer need the
      # per-drug patient_pred matrix.
      is_combo <- identical(input$drug_name_common, "Combination")
      if (is_combo && pt %in% c("boxplot", "roc")) {
        reqs <- c("predictions", "user_clones", "user_response")
      }

      # Check requirements
      missing <- c()
      if ("user_clones" %in% reqs && is.null(shared$user_clones)) missing <- c(missing, "Clone Annotation")
      if ("predictions" %in% reqs && is.null(shared$predictions)) missing <- c(missing, "Predictions")
      if ("patient_pred" %in% reqs && is.null(shared$patient_pred)) missing <- c(missing, "Patient Predictions")
      if ("user_response" %in% reqs && is.null(shared$user_response)) missing <- c(missing, "Clinical Response")
      if ("models" %in% reqs && is.null(shared$models)) missing <- c(missing, "Trained Model")

      if (length(missing) > 0) {
        showNotification(paste("Missing data for this plot:", paste(missing, collapse = ", ")), type = "warning", duration = 5)
        return()
      }

      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Generating plot..."),
          p(class = "text-muted", "This may take a few seconds")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()

      tryCatch({
        p <- switch(pt,

          "clone_dist" = {
            # Build clone_distribution data frame: patients, clones, weights
            clone_data <- shared$user_clones
            clone_dist_list <- lapply(unique(clone_data$patient), function(p) {
              p_cells <- clone_data[clone_data$patient == p, ]
              p_clones <- unique(p_cells$clone_id)
              n_p_cells <- nrow(p_cells)
              # Display label: when clone ids are patient-qualified (e.g.
              # "Kydar01_c1"), strip the known patient prefix so clones appear
              # as shared categories (c1/c2/c3) across patients, matching the
              # paper's clone distribution figures. Any other id scheme (e.g.
              # Seurat cluster ids "0","1") is kept unchanged.
              display_clones <- vapply(p_clones, function(cl) {
                if (!startsWith(cl, p)) return(cl)
                rest <- sub("^[\\._-]", "", substring(cl, nchar(p) + 1))
                if (nchar(rest) == 0) cl else rest
              }, character(1))
              data.frame(
                patients = p,
                clones = display_clones,
                weights = sapply(p_clones, function(cl) sum(p_cells$clone_id == cl) / n_p_cells),
                stringsAsFactors = FALSE
              )
            })
            clone_distribution <- do.call(rbind, clone_dist_list)
            # Add response column for faceting — keeps every group (R/NR, or
            # TN/RD/PD longitudinal time points); patients without a valid
            # response are dropped.
            if (!is.null(shared$user_response)) {
              clone_distribution$response <- shared$user_response$response[
                match(clone_distribution$patients, shared$user_response$patient)
              ]
              clone_distribution$response <- label_resp(clone_distribution$response)
              clone_distribution <- clone_distribution[!is.na(clone_distribution$response), ]
            }
            plot_clone_distribution <- PERCEPTIONx::plot_clone_distribution
            p <- plot_clone_distribution(clone_distribution, response_var = "response")
          },

          "clone_kill" = {
            # Build clone_viability data frame: patient, clone_id, comb_viability, weights
            clone_data <- shared$user_clones
            pred_mat <- shared$predictions
            drug <- if (nchar(input$drug_name_common) > 0) input$drug_name_common else colnames(pred_mat)[1]

            if (is_combo) {
              # Combination clone lollipop (paper Fig. 2b): z-score each drug
              # across all clones, then per-clone pmin (IDA principle). The
              # comb_viability column already holds the z-scored combination,
              # so it is NOT rescaled again below.
              combo_df <- combo_clone_frame()
              combo_df <- combo_df[!is.na(combo_df$response), ]
              if (nrow(combo_df) == 0) {
                stop("No clones with a valid R/NR response for the combination plot.")
              }
              PERCEPTIONx::plot_clone_viability(combo_df, viability_var = "comb_viability",
                                                weights_var = "weights", response_var = "response",
                                                drug = "Combination")
            } else {
            # Use clone_viability_template from prepared_data if available (most reliable)
            if (!is.null(shared$prepared_data$clone_viability_template)) {
              tmpl <- shared$prepared_data$clone_viability_template
              clone_viability_df <- data.frame(
                patient = tmpl$patient,
                clone_id = tmpl$clone_id,
                comb_viability = pred_mat[rownames(pred_mat), drug],
                stringsAsFactors = FALSE
              )
              # Clone abundance weights (cell-count proportion) so single-drug
              # lollipops show the same point-size legend as Combination mode.
              if (!is.null(clone_data) && nrow(clone_data) > 0) {
                clone_viability_df$weights <- vapply(seq_len(nrow(clone_viability_df)), function(i) {
                  pat <- clone_viability_df$patient[i]
                  cl  <- clone_viability_df$clone_id[i]
                  n_p <- sum(clone_data$patient == pat, na.rm = TRUE)
                  if (n_p == 0) NA_real_ else
                    sum(clone_data$patient == pat & clone_data$clone_id == cl, na.rm = TRUE) / n_p
                }, numeric(1))
              }
            } else {
              # Fallback: parse rownames
              parsed <- PERCEPTIONx::parse_clone_keys(rownames(pred_mat))
              pred_clone_ids <- parsed$clone_id
              pred_patients  <- parsed$patient

              clone_kill_list <- lapply(unique(clone_data$patient), function(p) {
                p_clones <- unique(clone_data$clone_id[clone_data$patient == p])
                p_clones <- intersect(p_clones, pred_clone_ids)
                if (length(p_clones) == 0) return(NULL)
                n_p_cells <- sum(clone_data$patient == p)
                pred_rows <- which(pred_clone_ids %in% p_clones & pred_patients == p)
                if (length(pred_rows) == 0) {
                  pred_rows <- match(p_clones, pred_clone_ids)
                }
                data.frame(
                  patient = p,
                  clone_id = pred_clone_ids[pred_rows],
                  comb_viability = pred_mat[pred_rows, drug],
                  weights = sapply(pred_clone_ids[pred_rows], function(cl) sum(clone_data$clone_id == cl) / n_p_cells),
                  stringsAsFactors = FALSE
                )
              })
              clone_viability_df <- do.call(rbind, clone_kill_list)
            }
            if (is.null(clone_viability_df) || nrow(clone_viability_df) == 0) {
              stop("No matching clones between prediction matrix and clone annotation.")
            }
            if (!is.null(shared$user_response)) {
              clone_viability_df$response <- label_resp(shared$user_response$response[
                match(clone_viability_df$patient, shared$user_response$patient)
              ])
              clone_viability_df <- clone_viability_df[!is.na(clone_viability_df$response), ]
            }
            # Match the paper's lollipop: z-score predicted viability across all
            # clones (per drug), so most (sensitive) clones fall below the zero
            # line and resistant outliers point upward — same as the paper.
            clone_viability_df$comb_viability <- as.numeric(scale(clone_viability_df$comb_viability))
            p_result <- PERCEPTIONx::plot_clone_viability(clone_viability_df, viability_var = "comb_viability",
                                          weights_var = "weights", response_var = "response",
                                          drug = drug)
            p_result
            }
          },

          "roc" = {
            cr <- shared$user_response
            if (is_combo) {
              # Combination ROC (paper Fig. 2e): patient-level combination score
              # = most-resistant clone weighted by abundance (weighted_max).
              combo_df <- combo_clone_frame()
              pat_df <- combo_patient_frame(combo_df)
              rv <- label_resp(cr$response[match(pat_df$patient, cr$patient)])
              predictor_vec <- pat_df$combination
            } else {
              pp <- shared$patient_pred
              drug <- if (nchar(input$drug_name_common) > 0) input$drug_name_common else colnames(pp)[1]
              rv <- label_resp(cr$response[match(rownames(pp), cr$patient)])
              predictor_vec <- pp[[drug]]
            }
            # Drop patients without a valid response — NA must never become a class
            # (that would bias the ROC / AUC).
            keep <- !is.na(rv) & !is.na(predictor_vec)
            rv <- rv[keep]
            predictor_vec <- predictor_vec[keep]
            # ROC is binary: use the two groups directly when there are exactly
            # two; otherwise use the pair selected in roc_pair_picker.
            grps <- response_groups()
            if (length(grps) > 2) {
              a <- input$roc_group_a
              b <- input$roc_group_b
              if (is.null(a) || !(a %in% grps)) a <- if ("PD" %in% grps) "PD" else grps[1]
              if (is.null(b) || !(b %in% grps)) b <- if ("RD" %in% grps) "RD" else grps[2]
              sel <- unique(c(a, b))
            } else {
              sel <- grps
            }
            response_vec <- factor(rv, levels = sel)
            keep2 <- !is.na(response_vec)
            response_vec <- response_vec[keep2]
            predictor_vec <- predictor_vec[keep2]
            # Auto-disable smoothing when sample size is small (< 10 patients)
            smooth <- length(response_vec) >= 10
            title <- if (is_combo) "Combination" else drug
            if (length(grps) > 2) title <- paste0(title, " — ", paste(sel, collapse = " vs "))
            PERCEPTIONx::plot_roc_curve(response = response_vec, predictor = predictor_vec,
                           smooth_curve = smooth, title = title)
          },

          "boxplot" = {
            cr <- shared$user_response
            if (is_combo) {
              # Combination response boxplot (paper Fig. 2d): patient-level
              # combination score = most-resistant clone weighted by abundance.
              combo_df <- combo_clone_frame()
              pat_df <- combo_patient_frame(combo_df)
              rv <- label_resp(cr$response[match(pat_df$patient, cr$patient)])
              predictor_vec <- pat_df$combination
            } else {
              pp <- shared$patient_pred
              drug <- if (nchar(input$drug_name_common) > 0) input$drug_name_common else colnames(pp)[1]
              rv <- label_resp(cr$response[match(rownames(pp), cr$patient)])
              predictor_vec <- pp[[drug]]
            }
            # Drop patients without a valid response; keep ALL response groups
            # (2 for R/NR data, 3 for longitudinal time points like TN/RD/PD).
            keep <- !is.na(rv) & !is.na(predictor_vec)
            rv <- rv[keep]
            predictor_vec <- predictor_vec[keep]
            exp_vs_pred <- data.frame(
              response = factor(rv, levels = unique(rv)),
              predicted_viability = predictor_vec,
              stringsAsFactors = FALSE
            )
            # plot_response_boxplot has no 'title' parameter — use ggplot2::labs() after
            p <- PERCEPTIONx::plot_response_boxplot(exp_vs_pred)
            p <- p + ggplot2::ggtitle(if (is_combo) "Combination" else drug)
            p
          }
        )

        if (!is.null(p)) {
          current_plot(p)
          current_plot_type(pt)
          w$hide()
          showNotification(paste0(plot_labels[[pt]], " generated successfully"), type = "message")
        } else {
          w$hide()
          showNotification("Selected plot type is not available with current data", type = "warning")
        }
      }, error = function(e) {
        w$hide()
        showNotification(paste("Plot error:", e$message), type = "error", duration = 8)
      })
    })

    # Generate Advanced Plot (UMAP spatial visualizations)
    observeEvent(input$generate_advanced, {
      pt <- input$plot_type_advanced
      if (is.null(pt) || nchar(pt) == 0) {
        showNotification("Select a spatial plot type first.", type = "warning")
        return()
      }

      # Check for 2D embedding coordinates from prepare_data()
      umap_coords <- shared$prepared_data$umap_coords
      if (is.null(umap_coords)) {
        showNotification("No 2D embedding available (clustering was skipped or not run). Re-run with Seurat clustering to enable spatial plots.", type = "warning", duration = 8)
        return()
      }

      clone_data <- shared$user_clones
      if (is.null(clone_data)) {
        showNotification("No clone annotation found. Run Seurat clustering first.", type = "warning")
        return()
      }

      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Generating plot..."),
          p(class = "text-muted", "This may take a few seconds")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()

      tryCatch({
        # Build common embedding data frame
        common_cells <- intersect(clone_data$cell_id, umap_coords$cell_id)
        if (length(common_cells) == 0) {
          w$hide()
          showNotification(paste0("No matching cells between ", reduction_label(), " coordinates and clone annotation."), type = "error")
          return()
        }

        p <- switch(pt,

          "umap_gene" = {
            gene <- input$umap_gene
            if (is.null(gene) || nchar(gene) == 0) {
              showNotification("Select a gene first.", type = "warning")
              NULL
            } else {
              expr_mat <- shared$user_expr
              if (is.null(expr_mat) || !(gene %in% rownames(expr_mat))) {
                showNotification(paste0("Gene '", gene, "' not found in expression matrix."), type = "error")
                NULL
              } else {
                cell_expr <- as.numeric(expr_mat[gene, common_cells])
                xy <- get_embedding_xy(umap_coords, common_cells)
                n <- length(common_cells)
                stopifnot(length(xy$X) == n, length(xy$Y) == n, length(cell_expr) == n)
                umap_data <- data.frame(
                  X = xy$X,
                  Y = xy$Y,
                  expression = scale(cell_expr)[, 1],
                  row.names = common_cells
                )
                PERCEPTIONx::plot_tsne_response(umap_data, color_var = "expression",
                                                title = gene, color_label = "Expression (z-score)",
                                                palette = "diverging", midpoint = 0,
                                                base_size = 11)
              }
            }
          },

          "umap_viability" = {
            pred_mat <- shared$predictions
            if (is.null(pred_mat)) {
              showNotification("No clone-level predictions found. Run predictions first.", type = "warning")
              NULL
            } else {
              drug <- if (!is.null(input$umap_drug) && nchar(input$umap_drug) > 0) input$umap_drug else colnames(pred_mat)[1]

              # Build Patient@@CloneID key for each cell to look up its
              # clone-level viability value (correct per-patient-per-clone)
              cell_keys <- PERCEPTIONx::build_clone_key(clone_data$patient, clone_data$clone_id)
              pred_keys <- rownames(pred_mat)
              cell_viability <- setNames(pred_mat[match(cell_keys, pred_keys), drug],
                                       clone_data$cell_id)

              kill_common <- intersect(names(cell_viability), umap_coords$cell_id)
              kill_common <- kill_common[!is.na(cell_viability[kill_common])]
              if (length(kill_common) == 0) {
                showNotification(paste0("No matching cells between ", reduction_label(), " coordinates and prediction data."), type = "error")
                NULL
              } else {
                raw_vals <- cell_viability[kill_common]
                scaled_vals <- safe_range01(raw_vals)
                xy <- get_embedding_xy(umap_coords, kill_common)
                # Defensive: ensure all columns have the same length
                n <- length(kill_common)
                stopifnot(length(xy$X) == n, length(xy$Y) == n, length(scaled_vals) == n)
                umap_data <- data.frame(
                  X = xy$X,
                  Y = xy$Y,
                  viability_scaled = scaled_vals,
                  row.names = kill_common
                )
                PERCEPTIONx::plot_tsne_response(umap_data, color_var = "viability_scaled",
                                                title = drug, color_label = "Predicted Viability",
                                                palette = "viridis", base_size = 11)
              }
            }
          },

          "umap_clone" = {
            # Clone identity on the embedding (paper Extended Data Fig. 8a):
            # each cell colored by its clone, no prediction needed.
            idx <- match(common_cells, clone_data$cell_id)
            xy <- get_embedding_xy(umap_coords, common_cells)
            n <- length(common_cells)
            stopifnot(length(xy$X) == n, length(xy$Y) == n, length(idx) == n)
            umap_data <- data.frame(
              X = xy$X,
              Y = xy$Y,
              clone_id = clone_data$clone_id[idx],
              stringsAsFactors = FALSE
            )
            PERCEPTIONx::plot_clone_umap(umap_data, title = "Clone Identity")
          },

          NULL
        )

        if (!is.null(p)) {
          current_plot(p)
          current_plot_type(pt)
          w$hide()
          showNotification(paste0(spatial_plot_label(pt), " generated successfully"), type = "message")
        } else {
          w$hide()
        }
      }, error = function(e) {
        w$hide()
        showNotification(paste("Plot error:", e$message), type = "error", duration = 8)
      })
    })

    # --- Interactive rendering via ggiraph (hover tooltips) ---
    # Renders the ORIGINAL ggplot as an SVG widget. Because the plot is never
    # converted to plotly, facet/strip layouts stay exactly as designed (no
    # re-flow, no overlap). Hovering a point shows a rich tooltip. Canvas is
    # controlled with CSS (.girafe.html-widget) so it grows with the actual
    # plot height, is centred, and can never overflow the card.
    output$main_plot <- ggiraph::renderGirafe({
      req(current_plot())
      p_obj <- current_plot()
      scale <- text_scale()
      sz <- plot_size()

      # Lollipop: stretch the canvas height so every facet gets more vertical
      # room (many panels sharing one row otherwise squeeze text/ticks).
      if (identical(current_plot_type(), "clone_kill")) {
        sz$h <- round(sz$h * 1.3)
      }

      # Objects ggiraph cannot render (grid.arrange gtable / plotly widget):
      # fall back to an explanatory static-looking panel.
      if (inherits(p_obj, "gtable") || inherits(p_obj, "grob") || inherits(p_obj, "gTree") ||
          inherits(p_obj, "plotly") || inherits(p_obj, "htmlwidget")) {
        p_fallback <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = "This plot type cannot be shown interactively.",
                            size = 5, color = "#5a6a8a") +
          ggplot2::theme_void()
        p_obj <- p_fallback
      }

      # Apply text size scale relative to the plot's OWN base font size (not a
      # hardcoded value) so 100% -> base and >100% grows monotonically.
      base <- p_obj$theme$text$size
      if (is.null(base) || length(base) == 0 || !is.finite(base)) base <- 11

      p_scaled <- if (scale != 1) {
        p_obj +
          theme(text = element_text(size = base * scale),
                plot.title = element_text(size = (base + 1) * scale, hjust = 0, vjust = 0,
                                          face = "plain", margin = margin(b = 6)),
                plot.subtitle = element_text(size = base * scale),
                axis.title = element_text(size = base * scale),
                axis.text = element_text(size = (base - 1) * scale),
                legend.text = element_text(size = (base - 2) * scale),
                legend.title = element_text(size = (base - 1) * scale),
                strip.text = element_text(size = (base - 1) * scale))
      } else {
        p_obj
      }

      # width_svg/height_svg set the SVG viewBox aspect ratio (canvas shape).
      ggiraph::girafe(
        ggobj = p_scaled,
        width_svg = max(sz$w / 100, 4),
        height_svg = max(sz$h / 100, 3),
        options = list(
          ggiraph::opts_sizing(rescale = TRUE, width = 1),
          ggiraph::opts_tooltip(
            use_cursor = TRUE, opacity = 0.92, offx = 14, offy = 14,
            css = "background-color:#1e2a4a;color:#ffffff;border-radius:5px;
                   padding:7px 11px;font-size:12px;line-height:1.45;
                   box-shadow:0 3px 10px rgba(0,0,0,0.30);"
          ),
          ggiraph::opts_hover(
            css = "cursor:pointer;fill-opacity:1;stroke:#1e2a4a;stroke-width:1.5px;"
          ),
          ggiraph::opts_zoom(min = 1, max = 1)
        )
      )
    })

    output$plot_status <- renderUI({
      if (is.null(current_plot())) {
        div(class = "text-muted", style = "text-align: center; padding: 3rem;",
          icon("chart-area", style = "font-size: 3rem; opacity: 0.15; display: block; margin-bottom: 0.8rem;"),
          "Select a plot type, configure parameters, and click ",
          strong("Generate Plot"),
          br(),
          span(style = "font-size: 0.82rem; opacity: 0.7;", "Check Data Status panel for available data")
        )
      } else {
        NULL
      }
    })

    # Inline plot info in header (replaces separate Plot Details card)
    output$plot_info_inline <- renderUI({
      req(current_plot())
      pt <- current_plot_type()
      if (is.null(pt)) return(NULL)
      # Spatial plots get method prefix
      label <- if (pt %in% c("umap_gene", "umap_viability", "umap_clone")) spatial_plot_label(pt) else plot_labels[[pt]]
      tags$span(class = "viz-info-inline",
        strong(label),
        tags$span(class = "text-muted", style = "margin-left: 0.5rem; font-size: 0.78rem;",
          format(Sys.time(), "%Y-%m-%d %H:%M")
        )
      )
    })

    # Plot explanation card
    plot_explanations <- list(
      clone_dist = list(
        title = "Clone Distribution",
        desc = "Shows clone proportions within each patient. With global clustering, the same clone (color) genuinely spans patients. With clone-level input using per-patient labels (e.g. c1/c2/c3), the same color across patients does not imply the same clone origin — it is a shared category label.",
        requires = "Clone annotation (Seurat clustering or clone-level input)"
      ),
      clone_kill = list(
        title = "Clone Viability Lollipop",
        desc = "Displays predicted drug viability scores for each clone within each patient. Taller = higher viability = more resistant; shorter = lower viability = more sensitive. Use this to identify which subclones are most affected by treatment.",
        requires = "Clone-level Predictions + Clone Annotation"
      ),
      roc = list(
        title = "ROC Curve",
        desc = "Receiver Operating Characteristic curve evaluating how well the predicted viability score discriminates Responders (R) from Non-Responders (NR). AUC closer to 1.0 indicates better prediction accuracy.",
        requires = "Patient-level Predictions + Clinical Response"
      ),
      boxplot = list(
        title = "Response Boxplot",
        desc = "Patient-level predicted viability scores per clinical response group. Two groups (R/NR): p-value is a one-sided Wilcoxon test asking whether responders are predicted more sensitive (lower viability). Three or more groups (e.g. TN/RD/PD): p-value is a global Kruskal-Wallis test asking whether any group differs. Non-parametric tests keep the p-value valid even with the small patient sample sizes inherent to clinical response data.",
        requires = "Patient-level Predictions + Clinical Response"
      ),
      umap_gene = list(
        title = "Gene Expression",
        desc = "Shows the expression level of a selected gene across all single cells in the 2D embedding space. Color gradient indicates expression intensity — brighter colors = higher expression. Use this to examine how a gene's expression pattern relates to the transcriptional subclone landscape.",
        requires = "Clone-level Predictions + Expression Matrix + 2D embedding coordinates"
      ),
      umap_viability = list(
        title = "Drug Viability",
        desc = "Shows the predicted drug viability score across all single cells in the 2D embedding space. Color gradient indicates viability — brighter = higher viability (more resistant), darker = lower viability (more sensitive). Use this to identify which regions of the embedding (i.e., which subclones) are most affected by a given drug.",
        requires = "Clone-level Predictions + 2D embedding coordinates (auto-generated by Seurat)"
      ),
      umap_clone = list(
        title = "Clone Identity",
        desc = "Shows where each clone sits in the 2D embedding space, with every cell colored by its clone (paper Extended Data Fig. 8a style). Use this to see the spatial layout of each transcriptional subclone and how it relates to the rest of the tumor.",
        requires = "Clone annotation + 2D embedding coordinates (auto-generated by Seurat)"
      )
    )

    output$plot_explanation <- renderUI({
      # Follow the currently SELECTED plot type (spatial radio takes priority
      # when active), so the annotation updates the moment the radio changes —
      # not only after Generate is clicked. A selected-but-not-yet-generated
      # type previews its own explanation.
      pt <- if (!is.null(input$plot_type_advanced) && nzchar(input$plot_type_advanced)) {
        input$plot_type_advanced
      } else {
        input$plot_type
      }
      if (is.null(pt) || nchar(pt) == 0) return(NULL)

      info <- plot_explanations[[pt]]
      if (is.null(info)) return(NULL)

      # Prefix title for spatial plots
      display_title <- if (pt %in% c("umap_gene", "umap_viability", "umap_clone"))
        paste(reduction_label(), info$title) else info$title

      # Clone-level input without a count column falls back to equal weights;
      # say so explicitly so users do not misread the proportions.
      equal_weight_note <- NULL
      if (pt == "clone_dist" && !is.null(shared$prepared_data) &&
          identical(shared$prepared_data$reduction_method, "none") &&
          !isTRUE(shared$mapping_has_count)) {
        equal_weight_note <- div(class = "info-box",
          style = "border-left-color: var(--warning, #d97706); margin: 0.6rem 0 0.8rem 0; font-size: 0.82rem; line-height: 1.5;",
          icon("triangle-exclamation"),
          " Clone-level input without a count column: proportions are equal (1/n per clone). Add a count column with real cell numbers in the mapping file to show true proportions."
        )
      }

      div(class = "card viz-explanation-card",
        div(class = "card-header",
          icon("circle-info"), " About This Plot"
        ),
        div(class = "card-body",
          h6(strong(display_title)),
          p(class = "text-muted", style = "font-size: 0.85rem; line-height: 1.5;", info$desc),
          equal_weight_note,
          tags$span(class = "viz-explanation-req",
            icon("clipboard-check", style = "font-size: 0.75rem;"),
            tags$span(style = "font-size: 0.78rem; font-weight: 600;", "Requires: "),
            tags$span(style = "font-size: 0.78rem;", info$requires)
          )
        )
      )
    })

    # Inline download buttons in header
    output$plot_download_inline <- renderUI({
      req(current_plot())
      tagList(
        downloadButton(ns("download_png"), "PNG", class = "btn-outline-primary btn-sm"),
        downloadButton(ns("download_pdf"), "PDF", class = "btn-outline-primary btn-sm")
      )
    })

    # A short, safe file stem for the CURRENTLY generated plot (not the
    # selected radio, which can differ for spatial plots).
    download_stem <- function() {
      pt <- current_plot_type()
      if (is.null(pt) || nchar(pt) == 0) pt <- input$plot_type
      paste0(pt, "_", format(Sys.Date(), "%Y%m%d"))
    }

    output$download_png <- downloadHandler(
      filename = function() paste0(download_stem(), ".png"),
      content = function(file) {
        p_obj <- current_plot()
        if (inherits(p_obj, "plotly") || inherits(p_obj, "htmlwidget")) {
          # Plotly object — use orca if available, else export via webshot
          tryCatch(
            plotly::orca(p_obj, file, scale = 2),
            error = function(e) {
              # Fallback: htmlwidget screenshot
              htmlwidgets::saveWidget(p_obj, tempfile(fileext = ".html"))
              showNotification("Plotly plot saved as interactive HTML (orca not available for PNG export)", type = "message")
            }
          )
        } else {
          PERCEPTIONx::export_plot_cairo(file, p_obj, format = "png", width = 10, height = 7, res = 600)
        }
      }
    )

    output$download_pdf <- downloadHandler(
      filename = function() paste0(download_stem(), ".pdf"),
      content = function(file) {
        p_obj <- current_plot()
        if (inherits(p_obj, "plotly") || inherits(p_obj, "htmlwidget")) {
          tryCatch(
            plotly::orca(p_obj, file),
            error = function(e) {
              showNotification("PDF export for plotly plots requires orca. Saved as HTML instead.", type = "warning")
            }
          )
        } else {
          PERCEPTIONx::export_plot_cairo(file, p_obj, format = "pdf", width = 10, height = 7)
        }
      }
    )

  })
}
