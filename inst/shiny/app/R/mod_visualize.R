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
                           "Clone Killing Lollipop" = "clone_kill",
                           "ROC Curve" = "roc",
                           "Response Boxplot" = "boxplot",
                           "Model Performance" = "model_perf"
                         ),
                         selected = "clone_dist"),

            conditionalPanel(
              condition = paste0("input['", ns("plot_type"), "'] == 'clone_kill' || input['", ns("plot_type"), "'] == 'boxplot' || input['", ns("plot_type"), "'] == 'roc'"),
              selectizeInput(ns("drug_name_common"), "Drug Name",
                             choices = NULL, selected = NULL,
                             options = list(maxItems = 1, placeholder = "Select a drug"))
            ),

            actionButton(ns("generate"), "Generate Plot",
                         class = "btn-primary btn-sm", icon = icon("wand-magic-sparkles")),

            # Plot size controls
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
            ),

            hr(),

            div(class = "viz-advanced-section",
              h6(class = "viz-advanced-title", icon("map"), " Spatial Visualizations"),
              p(class = "viz-advanced-desc", "Uses UMAP coordinates from Seurat clustering. No extra files needed."),

              radioButtons(ns("plot_type_advanced"), NULL,
                           choices = c(
                             "UMAP Gene Expression" = "umap_gene",
                             "UMAP Drug Killing" = "umap_killing"
                           ),
                           selected = character(0)),

              conditionalPanel(
                condition = paste0("input['", ns("plot_type_advanced"), "'] == 'umap_gene'"),
                selectizeInput(ns("umap_gene"), "Gene",
                               choices = NULL, selected = NULL,
                               options = list(maxItems = 1, placeholder = "Select a gene"))
              ),
              conditionalPanel(
                condition = paste0("input['", ns("plot_type_advanced"), "'] == 'umap_killing'"),
                selectizeInput(ns("umap_drug"), "Drug",
                               choices = NULL, selected = NULL,
                               options = list(maxItems = 1, placeholder = "Select a drug"))
              ),
              conditionalPanel(
                condition = paste0("input['", ns("plot_type_advanced"), "'] == 'umap_gene' || input['", ns("plot_type_advanced"), "'] == 'umap_killing'"),
                actionButton(ns("generate_advanced"), "Generate UMAP Plot",
                             class = "btn-primary btn-sm", icon = icon("wand-magic-sparkles"))
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
            plotlyOutput(ns("main_plot"), height = "750px"),
            tags$script(HTML(sprintf("
              (function() {
                var plotId = '%s';
                // Height handler
                var hName = '%s';
                if (!window[hName]) {
                  window[hName] = true;
                  Shiny.addCustomMessageHandler(hName, function(h) {
                    var el = document.getElementById(plotId);
                    if (el) {
                      el.style.height = h + 'px';
                      var widgets = el.querySelectorAll('.html-widget, .plotly.html-widget');
                      widgets.forEach(function(w) {
                        w.style.height = h + 'px';
                        w.style.maxHeight = h + 'px';
                      });
                      window.dispatchEvent(new Event('resize'));
                    }
                  });
                }
                // Width handler
                var wName = '%s';
                if (!window[wName]) {
                  window[wName] = true;
                  Shiny.addCustomMessageHandler(wName, function(w) {
                    var el = document.getElementById(plotId);
                    if (el) {
                      el.style.width = w + '%';
                      var widgets = el.querySelectorAll('.html-widget, .plotly.html-widget');
                      widgets.forEach(function(wd) {
                        wd.style.width = w + '%';
                        wd.style.maxWidth = w + '%';
                      });
                      window.dispatchEvent(new Event('resize'));
                    }
                  });
                }
              })();
            ", ns("main_plot"), ns("viz_update_plot_height"), ns("viz_update_plot_width"))))
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

    # Reactive plot height (in px) — driven by user's Width/Height size controls
    # Apply via session$sendCustomMessage + JS handler (no widget re-render)
    observeEvent(input$plot_height, {
      h_pct <- if (is.null(input$plot_height)) 100 else as.numeric(input$plot_height)
      h <- round(750 * h_pct / 100)
      session$sendCustomMessage(ns("viz_update_plot_height"), as.character(h))
    }, ignoreInit = TRUE)

    # Reactive plot width (in %) — driven by user's Width size control
    observeEvent(input$plot_width, {
      w_pct <- if (is.null(input$plot_width)) 100 else as.numeric(input$plot_width)
      session$sendCustomMessage(ns("viz_update_plot_width"), as.character(w_pct))
    }, ignoreInit = TRUE)

    # Reactive text size scale (50-160%)
    text_scale <- reactive({
      if (is.null(input$plot_text_size)) 100 else as.numeric(input$plot_text_size) / 100
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
      model_perf = c("models"),
      umap_gene = c("predictions", "user_expr"),
      umap_killing = c("predictions")
    )

    plot_labels <- list(
      clone_dist = "Clone Distribution",
      clone_kill = "Clone Killing Lollipop",
      roc = "ROC Curve",
      boxplot = "Response Boxplot",
      model_perf = "Model Performance",
      umap_gene = "UMAP Gene Expression",
      umap_killing = "UMAP Drug Killing"
    )

    # Populate drug choices from trained models (or predictions as fallback)
    observe({
      drug_choices <- NULL
      if (!is.null(shared$models)) {
        drug_choices <- names(shared$models)
      } else if (!is.null(shared$predictions)) {
        drug_choices <- colnames(shared$predictions)
      }
      if (!is.null(drug_choices) && length(drug_choices) > 0) {
        updateSelectizeInput(session, "drug_name_common", choices = drug_choices, server = TRUE)
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
        p <- suppressWarnings(switch(pt,

          "clone_dist" = {
            # Build clone_distribution data frame: patients, clones, weights
            clone_data <- shared$user_clones
            clone_dist_list <- lapply(unique(clone_data$patient), function(p) {
              p_cells <- clone_data[clone_data$patient == p, ]
              p_clones <- unique(p_cells$clone_id)
              n_p_cells <- nrow(p_cells)
              data.frame(
                patients = p,
                clones = p_clones,
                weights = sapply(p_clones, function(cl) sum(p_cells$clone_id == cl) / n_p_cells),
                stringsAsFactors = FALSE
              )
            })
            clone_distribution <- do.call(rbind, clone_dist_list)
            # Add response column for faceting
            if (!is.null(shared$user_response)) {
              clone_distribution$response <- shared$user_response$response[
                match(clone_distribution$patients, shared$user_response$patient)
              ]
            }
            plot_clone_distribution <- PERCEPTION::plot_clone_distribution
            p <- plot_clone_distribution(clone_distribution, response_var = "response")
          },

          "clone_kill" = {
            # Build clone_killing data frame: patient, clone_id, comb_killing, weights
            clone_data <- shared$user_clones
            pred_mat <- shared$predictions
            drug <- if (nchar(input$drug_name_common) > 0) input$drug_name_common else colnames(pred_mat)[1]

            # Use clone_killing_template from prepared_data if available (most reliable)
            if (!is.null(shared$prepared_data$clone_killing_template)) {
              tmpl <- shared$prepared_data$clone_killing_template
              clone_killing_df <- data.frame(
                patient = tmpl$patient,
                clone_id = tmpl$clone_id,
                comb_killing = pred_mat[rownames(pred_mat), drug],
                stringsAsFactors = FALSE
              )
            } else {
              # Fallback: parse rownames
              pred_clone_ids <- sapply(strsplit(rownames(pred_mat), "@@"), `[`, 2)
              pred_patients <- sapply(strsplit(rownames(pred_mat), "@@"), `[`, 1)

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
                  comb_killing = pred_mat[pred_rows, drug],
                  weights = sapply(pred_clone_ids[pred_rows], function(cl) sum(clone_data$clone_id == cl) / n_p_cells),
                  stringsAsFactors = FALSE
                )
              })
              clone_killing_df <- do.call(rbind, clone_kill_list)
            }
            if (is.null(clone_killing_df) || nrow(clone_killing_df) == 0) {
              stop("No matching clones between prediction matrix and clone annotation.")
            }
            if (!is.null(shared$user_response)) {
              clone_killing_df$response <- shared$user_response$response[
                match(clone_killing_df$patient, shared$user_response$patient)
              ]
            }
            p_result <- PERCEPTION::plot_clone_killing(clone_killing_df, killing_var = "comb_killing",
                                          weights_var = "weights", response_var = "response",
                                          drug = drug)
            p_result
          },

          "roc" = {
            # Build exp_vs_pred data frame for ROC
            pp <- shared$patient_pred
            cr <- shared$user_response
            drug <- if (nchar(input$drug_name_common) > 0) input$drug_name_common else colnames(pp)[1]
            response_vec <- cr$response[match(rownames(pp), cr$patient)]
            # Convert to R/NR factor
            response_vec <- ifelse(tolower(response_vec) %in% c("responder", "r"), "R", "NR")
            response_vec <- factor(response_vec, levels = c("R", "NR"))
            predictor_vec <- pp[[drug]]
            # Auto-disable smoothing when sample size is small (< 10 patients)
            n_pts <- sum(!is.na(response_vec) & !is.na(predictor_vec))
            smooth <- n_pts >= 10
            PERCEPTION::plot_roc_curve(response = response_vec, predictor = predictor_vec,
                           smooth_curve = smooth, title = drug)
          },

          "boxplot" = {
            # Build exp_vs_pred data frame
            pp <- shared$patient_pred
            cr <- shared$user_response
            drug <- if (nchar(input$drug_name_common) > 0) input$drug_name_common else colnames(pp)[1]
            response_vec <- cr$response[match(rownames(pp), cr$patient)]
            response_vec <- ifelse(tolower(response_vec) %in% c("responder", "r"), "R", "NR")
            exp_vs_pred <- data.frame(
              response = factor(response_vec, levels = c("R", "NR")),
              predicted_killing = pp[[drug]],
              stringsAsFactors = FALSE
            )
            # plot_response_boxplot has no 'title' parameter — use ggplot2::labs() after
            p <- PERCEPTION::plot_response_boxplot(exp_vs_pred)
            p <- p + ggplot2::ggtitle(drug)
            p
          },

          "model_perf" = {
            PERCEPTION::plot_model_performance(shared$models)
          },

          NULL
        ))

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

      # Check for UMAP coordinates from prepare_data()
      umap_coords <- shared$prepared_data$umap_coords
      if (is.null(umap_coords)) {
        showNotification("UMAP coordinates not found. Run Seurat clustering in Data module first.", type = "warning", duration = 8)
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
        # Build common UMAP data frame
        common_cells <- intersect(clone_data$cell_id, umap_coords$cell_id)
        if (length(common_cells) == 0) {
          w$hide()
          showNotification("No matching cells between UMAP coordinates and clone annotation.", type = "error")
          return()
        }

        p <- suppressWarnings(switch(pt,

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
                umap_data <- data.frame(
                  X = umap_coords$umap_1[match(common_cells, umap_coords$cell_id)],
                  Y = umap_coords$umap_2[match(common_cells, umap_coords$cell_id)],
                  expression = scale(cell_expr)[, 1],
                  row.names = common_cells
                )
                PERCEPTION::plot_tsne_response(umap_data, color_var = "expression",
                                                title = gene, color_label = "Expression",
                                                point_size = 0.8, base_size = 11)
              }
            }
          },

          "umap_killing" = {
            pred_mat <- shared$predictions
            if (is.null(pred_mat)) {
              showNotification("No clone-level predictions found. Run predictions first.", type = "warning")
              NULL
            } else {
              drug <- if (!is.null(input$umap_drug) && nchar(input$umap_drug) > 0) input$umap_drug else colnames(pred_mat)[1]

              pred_clone_ids <- sapply(strsplit(rownames(pred_mat), "@@"), `[`, 2)
              cell_pred <- setNames(pred_mat[, drug], pred_clone_ids)
              cell_killing <- cell_pred[clone_data$clone_id]
              names(cell_killing) <- clone_data$cell_id

              kill_common <- intersect(names(cell_killing), umap_coords$cell_id)
              if (length(kill_common) == 0) {
                showNotification("No matching cells between UMAP coordinates and prediction data.", type = "error")
                NULL
              } else {
                umap_data <- data.frame(
                  X = umap_coords$umap_1[match(kill_common, umap_coords$cell_id)],
                  Y = umap_coords$umap_2[match(kill_common, umap_coords$cell_id)],
                  killing_scaled = scale(cell_killing[kill_common])[, 1],
                  row.names = kill_common
                )
                PERCEPTION::plot_tsne_response(umap_data, color_var = "killing_scaled",
                                                title = drug, color_label = "Predicted Killing",
                                                point_size = 0.8, base_size = 11)
              }
            }
          },

          NULL
        ))

        if (!is.null(p)) {
          current_plot(p)
          current_plot_type(pt)
          w$hide()
          showNotification(paste0(plot_labels[[pt]], " generated successfully"), type = "message")
        } else {
          w$hide()
        }
      }, error = function(e) {
        w$hide()
        showNotification(paste("Plot error:", e$message), type = "error", duration = 8)
      })
    })

    output$main_plot <- renderPlotly({
      req(current_plot())
      p_obj <- current_plot()
      scale <- text_scale()
      # If already a plotly object (e.g. combined UMAP plots), return as-is
      if (inherits(p_obj, "plotly") || inherits(p_obj, "htmlwidget")) {
        p_obj
      } else {
        # Apply text size scale to ggplot's base font size
        # Use theme() addition (not theme_bw()) to preserve existing title/legend customization
        p_scaled <- if (scale != 1) {
          p_obj +
            theme(text = element_text(size = 11 * scale),
                  plot.title = element_text(size = 12 * scale, hjust = 0, vjust = 0, face = "plain",
                                            margin = margin(b = 6)),
                  plot.subtitle = element_text(size = 11 * scale),
                  axis.title = element_text(size = 11 * scale),
                  axis.text = element_text(size = 10 * scale),
                  legend.text = element_text(size = 9 * scale),
                  legend.title = element_text(size = 10 * scale),
                  strip.text = element_text(size = 10 * scale))
        } else {
          p_obj
        }
        ggplotly(p_scaled, tooltip = c("x", "y", "label")) %>%
          layout(
            font = list(family = "Inter, sans-serif", size = 12 * scale, color = "#1e2a4a"),
            paper_bgcolor = "transparent",
            plot_bgcolor = "transparent",
            xaxis = list(gridcolor = "#eef0f6", zerolinecolor = "#dfe3ee"),
            yaxis = list(gridcolor = "#eef0f6", zerolinecolor = "#dfe3ee"),
            legend = list(
              bgcolor = "rgba(255,255,255,0.8)",
              bordercolor = "#dfe3ee",
              borderwidth = 1,
              font = list(size = 11, color = "#1e2a4a")
            )
          )
      }
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
      tags$span(class = "viz-info-inline",
        strong(plot_labels[[pt]]),
        tags$span(class = "text-muted", style = "margin-left: 0.5rem; font-size: 0.78rem;",
          format(Sys.time(), "%Y-%m-%d %H:%M")
        )
      )
    })

    # Plot explanation card
    plot_explanations <- list(
      clone_dist = list(
        title = "Clone Distribution",
        desc = "Shows the proportion of each transcriptional subclone within every patient's cell population. Each bar represents a clone, colored by patient. This helps identify which subclones dominate each patient's tumor.",
        requires = "Clone Annotation (from Seurat clustering)"
      ),
      clone_kill = list(
        title = "Clone Killing Lollipop",
        desc = "Displays predicted drug killing scores for each clone within each patient. Lollipop height indicates sensitivity — taller = more sensitive to the drug. Use this to identify which subclones are most affected by treatment.",
        requires = "Clone-level Predictions + Clone Annotation"
      ),
      roc = list(
        title = "ROC Curve",
        desc = "Receiver Operating Characteristic curve evaluating how well the predicted killing score discriminates Responders (R) from Non-Responders (NR). AUC closer to 1.0 indicates better prediction accuracy.",
        requires = "Patient-level Predictions + Clinical Response"
      ),
      boxplot = list(
        title = "Response Boxplot",
        desc = "Compares predicted killing scores between Responders and Non-Responders. The p-value (Wilcoxon test) indicates whether the difference is statistically significant. A lower p-value suggests the model captures clinically meaningful differences.",
        requires = "Patient-level Predictions + Clinical Response"
      ),
      model_perf = list(
        title = "Model Performance",
        desc = "Summarizes cross-validation performance of the drug response models across bulk, pseudo-bulk, and single-cell datasets. Shows the number of drugs achieving various correlation thresholds between predicted and observed drug response.",
        requires = "Trained Models"
      ),
      umap_gene = list(
        title = "UMAP Gene Expression",
        desc = "Shows the expression level of a selected gene across all single cells in UMAP space. Color gradient indicates expression intensity — brighter colors = higher expression. Use this to examine how a gene's expression pattern relates to the transcriptional subclone landscape.",
        requires = "Clone-level Predictions + Expression Matrix + UMAP coordinates"
      ),
      umap_killing = list(
        title = "UMAP Drug Killing",
        desc = "Shows the predicted drug killing score across all single cells in UMAP space. Color gradient indicates sensitivity — red = sensitive, blue = resistant. Use this to identify which regions of the UMAP (i.e., which subclones) are most affected by a given drug.",
        requires = "Clone-level Predictions + UMAP coordinates (auto-generated by Seurat)"
      )
    )

    # Track which plot type was last generated
    current_plot_type <- reactiveVal(NULL)

    output$plot_explanation <- renderUI({
      pt <- current_plot_type()
      if (is.null(pt) || nchar(pt) == 0 || is.null(current_plot())) return(NULL)

      info <- plot_explanations[[pt]]
      if (is.null(info)) return(NULL)

      div(class = "card viz-explanation-card",
        div(class = "card-header",
          icon("circle-info"), " About This Plot"
        ),
        div(class = "card-body",
          h6(strong(info$title)),
          p(class = "text-muted", style = "font-size: 0.85rem; line-height: 1.5;", info$desc),
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

    output$download_png <- downloadHandler(
      filename = function() paste0(input$plot_type, "_", format(Sys.Date(), "%Y%m%d"), ".png"),
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
          ggplot2::ggsave(file, p_obj, device = "png", width = 10, height = 7, dpi = 300)
        }
      }
    )

    output$download_pdf <- downloadHandler(
      filename = function() paste0(input$plot_type, "_", format(Sys.Date(), "%Y%m%d"), ".pdf"),
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
          ggplot2::ggsave(file, p_obj, device = "pdf", width = 10, height = 7)
        }
      }
    )

  })
}
