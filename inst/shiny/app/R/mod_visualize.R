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
            uiOutput(ns("main_plot_host"))
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
    # Guards against double-submitting a plot task (rapid double clicks).
    viz_busy <- reactiveVal(FALSE)

    # --- In-memory plot cache ---
    # plot_cache   : request signature (plot type + parameters) -> finished ggplot
    # widget_cache : plot signature + text scale + canvas size -> rendered
    #                girafe widget, so returning to a plot reuses the SVG
    #                instead of paying the fixed ~1.4 s conversion again.
    # Both are dropped whenever the underlying source data is replaced/cleared.
    # current_plot_key records which signature the currently shown plot came
    # from, so the widget cache can be keyed stably per plot.
    plot_cache <- new.env(parent = emptyenv())
    widget_cache <- new.env(parent = emptyenv())
    current_plot_key <- reactiveVal(NULL)

    # When ALL source data disappears (e.g. "Clear Demo Data" on the Data
    # tab), drop the cached plot so the Output area collapses back to the
    # placeholder shown on first launch instead of a stale blank block. It
    # expands again once new data is loaded and a new plot is generated.
    observe({
      no_data <- is.null(shared$predictions) && is.null(shared$models) &&
                 is.null(shared$user_clones) && is.null(shared$user_expr)
      if (no_data && !is.null(current_plot())) {
        current_plot(NULL)
        current_plot_type(NULL)
        current_plot_key(NULL)
      }
    })

    # Any source-object change (demo load/clear, uploads, new predictions)
    # makes every cached plot stale — drop both caches.
    observe({
      invisible(list(shared$predictions, shared$user_clones, shared$user_expr,
                     shared$user_response, shared$patient_pred,
                     shared$models, shared$prepared_data))
      if (length(ls(plot_cache)) > 0) rm(list = ls(plot_cache), envir = plot_cache)
      if (length(ls(widget_cache)) > 0) rm(list = ls(widget_cache), envir = widget_cache)
    })

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
          icon("info-circle"), " Response has ", length(grps), " groups. Pick two for ROC:"),
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

    # Generate Plot — computation runs in the background per-session worker
    # (task "plot", see async_jobs.R); the returned ggplot is drawn by
    # renderGirafe below, so the main thread never blocks on plot math.
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

      # Cache hit: the identical request was generated before — reuse it
      # instantly (no background task, no waiter overlay).
      ck <- paste(pt, is_combo, input$drug_name_common,
                  input$roc_group_a, input$roc_group_b, sep = "\u0001")
      cached <- get0(ck, envir = plot_cache, inherits = FALSE)
      if (!is.null(cached)) {
        current_plot(cached$plot)
        current_plot_type(pt)
        current_plot_key(ck)
        showNotification(paste0(plot_labels[[pt]], " loaded from cache"), type = "message", duration = 3)
        return()
      }

      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Generating plot..."),
          p(class = "text-muted", "Drawing plot...")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      if (viz_busy()) {
        showNotification("A plot is already being generated", type = "warning", duration = 5)
        return()
      }
      viz_busy(TRUE)
      w$show()

      pd <- shared$prepared_data
      jobid <- tryCatch(
        submit_session_task(shared, "plot", list(
          plot_type = pt,
          data = list(
            predictions   = shared$predictions,
            user_clones   = shared$user_clones,
            user_response = shared$user_response,
            patient_pred  = shared$patient_pred,
            # Chart plots never need the expression matrix or embedding
            # coordinates — carrying them made params.rds serialize/deserialize
            # slow (lollipop / distribution / ROC / box).
            prepared = list(
              umap_coords = if (!is.null(pd) && pt %in% c("umap_gene", "umap_viability", "umap_clone")) pd$umap_coords else NULL,
              clone_viability_template = if (!is.null(pd)) pd$clone_viability_template else NULL
            )
          ),
          params = list(
            drug = if (!is.null(input$drug_name_common) && nchar(input$drug_name_common) > 0) input$drug_name_common else "",
            combo = is_combo,
            roc_group_a = input$roc_group_a,
            roc_group_b = input$roc_group_b
          )
        )),
        error = function(e) { viz_busy(FALSE); w$hide(); NULL }
      )
      if (is.null(jobid)) return()

      poll_task(shared, session, jobid,
        on_done = function(res) {
          viz_busy(FALSE)
          w$hide()
          if (!is.null(res$plot)) {
            assign(ck, list(plot = res$plot, message = res$message), envir = plot_cache)
            current_plot(res$plot)
            current_plot_type(pt)
            current_plot_key(ck)
            showNotification(paste0(plot_labels[[pt]], " generated successfully"), type = "message")
          } else {
            msg <- res$message
            if (is.null(msg) || !nzchar(msg)) msg <- "Selected plot type is not available with current data"
            showNotification(msg, type = "warning")
          }
        },
        on_error = function(msg) {
          viz_busy(FALSE)
          w$hide()
          showNotification(paste("Plot error:", msg), type = "error", duration = 8)
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

      # Cache hit: reuse the identical spatial plot instantly.
      ck <- paste(pt, input$umap_gene, input$umap_drug, sep = "\u0001")
      cached <- get0(ck, envir = plot_cache, inherits = FALSE)
      if (!is.null(cached)) {
        current_plot(cached$plot)
        current_plot_type(pt)
        current_plot_key(ck)
        showNotification(paste0(spatial_plot_label(pt), " loaded from cache"), type = "message", duration = 3)
        return()
      }

      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Generating plot..."),
          p(class = "text-muted", "Drawing plot...")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      if (viz_busy()) {
        showNotification("A plot is already being generated", type = "warning", duration = 5)
        return()
      }
      viz_busy(TRUE)
      w$show()

      # Gene Expression: ship ONLY the selected gene's per-cell expression
      # vector (a small named numeric) instead of the whole expression matrix
      # — a full matrix can be hundreds of MB, and serializing it into
      # params.rds took tens of seconds. NULL when the gene is invalid; the
      # worker then reports "not found".
      gene_vec <- NULL
      if (pt == "umap_gene" && !is.null(input$umap_gene) && nchar(input$umap_gene) > 0 &&
          !is.null(shared$user_expr) && input$umap_gene %in% rownames(shared$user_expr)) {
        gene_vec <- setNames(as.numeric(shared$user_expr[input$umap_gene, ]),
                             colnames(shared$user_expr))
      }

      pd <- shared$prepared_data
      jobid <- tryCatch(
        submit_session_task(shared, "plot", list(
          plot_type = pt,
          data = list(
            predictions   = shared$predictions,
            user_clones   = shared$user_clones,
            user_response = shared$user_response,
            patient_pred  = shared$patient_pred,
            # No plot ships the raw expression matrix anymore — umap_gene
            # gets only the requested gene's vector, umap_viability/clone
            # need none of it. umap_coords (n_cells x 2) stays small.
            umap_gene_expr = gene_vec,
            prepared = list(
              umap_coords = if (!is.null(pd) && pt %in% c("umap_gene", "umap_viability", "umap_clone")) pd$umap_coords else NULL,
              clone_viability_template = if (!is.null(pd)) pd$clone_viability_template else NULL
            )
          ),
          params = list(
            drug = "",
            combo = FALSE,
            roc_group_a = NULL,
            roc_group_b = NULL,
            umap_gene = input$umap_gene,
            umap_drug = input$umap_drug
          )
        )),
        error = function(e) { viz_busy(FALSE); w$hide(); NULL }
      )
      if (is.null(jobid)) return()

      poll_task(shared, session, jobid,
        on_done = function(res) {
          viz_busy(FALSE)
          w$hide()
          if (!is.null(res$plot)) {
            assign(ck, list(plot = res$plot, message = res$message), envir = plot_cache)
            current_plot(res$plot)
            current_plot_type(pt)
            current_plot_key(ck)
            showNotification(paste0(spatial_plot_label(pt), " generated successfully"), type = "message")
          } else {
            msg <- res$message
            if (is.null(msg) || !nzchar(msg)) msg <- "Selected plot type is not available with current data"
            showNotification(msg, type = "warning")
          }
        },
        on_error = function(msg) {
          viz_busy(FALSE)
          w$hide()
          showNotification(paste("Plot error:", msg), type = "error", duration = 8)
        })
    })

    # --- Interactive rendering via ggiraph (hover tooltips) ---
    # Renders the ORIGINAL ggplot as an SVG widget. Because the plot is never
    # converted to plotly, facet/strip layouts stay exactly as designed (no
    # re-flow, no overlap). Hovering a point shows a rich tooltip. Canvas is
    # controlled with CSS (.girafe.html-widget) so it grows with the actual
    # plot height, is centred, and can never overflow the card.
    #
    # Text-scale the plot (relative to its OWN base font size, so 100% -> base
    # and >100% grows monotonically), and replace gtable/grob/plotly objects
    # — which ggiraph cannot render — with an explanatory static panel.
    scaled_plot <- reactive({
      req(current_plot())
      p_obj <- current_plot()
      if (inherits(p_obj, "gtable") || inherits(p_obj, "grob") || inherits(p_obj, "gTree") ||
          inherits(p_obj, "plotly") || inherits(p_obj, "htmlwidget")) {
        p_obj <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = "This plot type cannot be shown interactively.",
                            size = 5, color = "#5a6a8a") +
          ggplot2::theme_void()
      }
      scale <- text_scale()
      base <- p_obj$theme$text$size
      if (is.null(base) || length(base) == 0 || !is.finite(base)) base <- 11
      if (scale == 1) return(p_obj)
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
    })

    # Conditional host: mount the girafe output ONLY while a plot exists.
    # Without this, Shiny hides an empty widget with `visibility: hidden`,
    # which still occupies its full layout space — leaving a large blank
    # block after "Clear Demo Data". Removing the element entirely makes the
    # output area collapse to just the placeholder, like on first launch.
    output$main_plot_host <- renderUI({
      if (is.null(current_plot())) NULL else ggiraph::girafeOutput(ns("main_plot"))
    })

    output$main_plot <- ggiraph::renderGirafe({
      req(current_plot())
      sz <- plot_size()
      # Lollipop: stretch the canvas height so every facet gets more vertical
      # room (many panels sharing one row otherwise squeeze text/ticks).
      if (identical(current_plot_type(), "clone_kill")) {
        sz$h <- round(sz$h * 1.3)
      }
      # Widget cache: the girafe SVG conversion is the fixed ~1.4 s cost.
      # Reuse the finished widget whenever the same plot is shown with the
      # same text scale and canvas size (e.g. switching back to a cached
      # plot) instead of converting again.
      wkey <- paste(current_plot_key(), text_scale(), sz$w, sz$h, sep = "\u0001")
      cached_w <- get0(wkey, envir = widget_cache, inherits = FALSE)
      if (!is.null(cached_w)) return(cached_w)
      # width_svg/height_svg set the SVG viewBox aspect ratio (canvas shape).
      w <- ggiraph::girafe(
        ggobj = scaled_plot(),
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
      assign(wkey, w, envir = widget_cache)
      w
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
        desc = "Patient-level predicted viability (z-score) per clinical response group. Two groups (R/NR): p-value is a one-sided Wilcoxon test asking whether responders are predicted more sensitive (lower viability). Three or more groups (e.g. TN/RD/PD): pairwise two-sided Wilcoxon tests with BH (FDR) multiple-testing correction — the adj. p on each bracket tells you which specific pairs differ. Non-parametric tests keep the p-value valid even with the small patient sample sizes inherent to clinical response data.",
        requires = "Patient-level Predictions + Clinical Response"
      ),
      umap_gene = list(
        title = "Gene Expression",
        desc = "Shows the expression level of a selected gene (0-1 normalized, 5th-95th percentile) across all single cells in the 2D embedding space. Grey = no/low expression, red = high expression (grey-to-red sequential ramp). Use this to examine how a gene's expression pattern relates to the transcriptional subclone landscape.",
        requires = "Clone-level Predictions + Expression Matrix + 2D embedding coordinates"
      ),
      umap_viability = list(
        title = "Drug Viability",
        desc = "Shows the predicted drug viability (z-score: 0 = cohort mean across clones, higher = more resistant) across all single cells in the 2D embedding space. Red = high viability (resistant), blue = low viability (sensitive), white = neutral — the same diverging red-blue ramp as the lollipop, with data-driven limits so extreme values keep a real color. Use this to identify which regions of the embedding (i.e., which subclones) are most affected by a given drug.",
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
        # current_plot() always holds a ggplot here (ggiraph conversion happens
        # only at display time), so export straight via the Cairo backend —
        # no orca / headless browser needed.
        PERCEPTIONx::export_plot_cairo(file, current_plot(), format = "png",
                                       width = 10, height = 7, res = 600)
      }
    )

    output$download_pdf <- downloadHandler(
      filename = function() paste0(download_stem(), ".pdf"),
      content = function(file) {
        PERCEPTIONx::export_plot_cairo(file, current_plot(), format = "pdf",
                                       width = 10, height = 7)
      }
    )

  })
}
