# Visualization Module — "Gallery Cockpit"
#
# Redesigned around a single user question: "what can I draw, and how?"
#  * A data-readiness bar at the top says at a glance what is loaded.
#  * All 7 plot types are visible at once as a clickable gallery. Ready cards
#    are active; cards whose data is missing are greyed out and say exactly
#    what they need (and let you jump to the tab that provides it).
#  * Clicking a card draws the plot IMMEDIATELY — there is no Generate button,
#    and switching plot types never requires a confirm step. Parameter changes
#    (drug / gene / ROC pair) redraw the current plot automatically.
#  * The plot explanation follows the currently selected card (no more stale
#    "About this plot" that sticks to the spatial radio), and the generated-at
#    timestamp is shown in UTC.
#
# The heavy lifting (background worker submit, plot_cache/widget_cache, girafe
# rendering, PNG/PDF download) is unchanged from the previous layout.

mod_visualize_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(column(12,
      div(class = "section-header", icon("chart-line"), h4("Visualization"))
    )),

    # Data readiness bar — one full-width row.
    fluidRow(column(12, uiOutput(ns("readiness_bar")))),

    fluidRow(
      # LEFT column: plot gallery + parameters.
      column(4, class = "viz-sidebar-col",
        div(class = "card animate-fade-in-up viz-gallery-card",
          div(class = "card-header",
            icon("images"), " Plot Gallery"
          ),
          div(class = "card-body",
            uiOutput(ns("gallery_grid")),
            tags$small(class = "text-muted", style = "display: block; margin-top: 0.5rem;",
              icon("lightbulb"), " Grey cards are missing data — hover to see what they need.")
          )
        ),
        # Parameters appear here once a card is clicked.
        uiOutput(ns("param_panel"))
      ),

      # RIGHT column: output + about.
      column(8, class = "viz-main-col",
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
        uiOutput(ns("plot_explanation"))
      )
    )
  )
}

mod_visualize_server <- function(id, shared, main_session) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    current_plot <- reactiveVal(NULL)
    # The currently SELECTED plot type (single source of truth — replaces the
    # old two radio groups). Set when a gallery card is clicked.
    viz_selected <- reactiveVal(NULL)
    current_plot_type <- reactiveVal(NULL)
    viz_busy <- reactiveVal(FALSE)

    # --- In-memory plot cache ---
    # plot_cache   : request signature (plot type + parameters) -> finished ggplot
    # widget_cache : plot signature + text scale + canvas size -> rendered
    #                girafe widget, so returning to a plot reuses the SVG
    #                instead of paying the fixed ~1.4 s conversion again.
    plot_cache <- new.env(parent = emptyenv())
    widget_cache <- new.env(parent = emptyenv())
    current_plot_key <- reactiveVal(NULL)

    # When ALL source data disappears (e.g. "Clear Demo Data"), drop the
    # cached plot and the selection so the Output area collapses back to the
    # placeholder instead of a stale blank block.
    observe({
      no_data <- is.null(shared$predictions) && is.null(shared$models) &&
                 is.null(shared$user_clones) && is.null(shared$user_expr)
      if (no_data && !is.null(current_plot())) {
        current_plot(NULL)
        current_plot_type(NULL)
        current_plot_key(NULL)
        viz_selected(NULL)
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

    response_groups <- reactive({
      cr <- shared$user_response
      if (is.null(cr) || is.null(cr$response)) return(character(0))
      grps <- unique(label_resp(cr$response))
      grps[!is.na(grps) & nzchar(trimws(grps)) & grps != "NA"]
    })

    # Reactive plot size / text scale / reduction label (unchanged).
    plot_size <- reactive({
      w_pct <- if (is.null(input$plot_width)) 100 else as.numeric(input$plot_width)
      h_pct <- if (is.null(input$plot_height)) 100 else as.numeric(input$plot_height)
      list(w = round(1000 * w_pct / 100), h = round(750 * h_pct / 100))
    })
    text_scale <- reactive({
      if (is.null(input$plot_text_size)) 100 else as.numeric(input$plot_text_size) / 100
    })
    reduction_method <- reactive({
      m <- shared$prepared_data$reduction_method
      if (is.null(m)) "umap" else m
    })
    reduction_label <- reactive({
      m <- reduction_method()
      if (m == "none") "no embedding" else if (m == "tsne") "t-SNE" else "UMAP"
    })

    # ---- Gallery metadata (single source for cards, requirements, labels) ----
    gallery_items <- list(
      list(id = "clone_dist",  icon = "chart-pie",   title = "Clone Distribution",
           tagline = "Clone proportions within each patient",
           req = c("user_clones"), spatial = FALSE),
      list(id = "clone_kill",  icon = "chart-column", title = "Clone Viability Lollipop",
           tagline = "Predicted viability per clone",
           req = c("predictions", "user_clones"), spatial = FALSE),
      list(id = "roc",         icon = "chart-line",  title = "ROC Curve",
           tagline = "How well scores separate responders from non-responders",
           req = c("patient_pred", "user_response"), spatial = FALSE),
      list(id = "boxplot",     icon = "chart-gantt", title = "Response Boxplot",
           tagline = "Predicted viability by clinical response group",
           req = c("patient_pred", "user_response"), spatial = FALSE),
      list(id = "umap_clone",  icon = "map",         title = "Clone Identity",
           tagline = "Where each clone sits in the embedding",
           req = c("user_clones"), spatial = TRUE),
      list(id = "umap_viability", icon = "flask",    title = "Drug Viability",
           tagline = "Predicted drug viability across cells",
           req = c("predictions"), spatial = TRUE),
      list(id = "umap_gene",   icon = "dna",         title = "Gene Expression",
           tagline = "Expression of one gene across cells",
           req = c("predictions", "user_expr"), spatial = TRUE)
    )

    plot_requirements <- setNames(
      lapply(gallery_items, function(it) it$req),
      vapply(gallery_items, function(it) it$id, character(1))
    )
    plot_labels <- setNames(
      lapply(gallery_items, function(it) it$title),
      vapply(gallery_items, function(it) it$id, character(1))
    )

    spatial_plot_label <- function(pt) {
      paste(reduction_label(), plot_labels[[pt]])
    }
    plot_label <- function(pt) {
      if (isTRUE(pt %in% c("umap_gene", "umap_viability", "umap_clone"))) spatial_plot_label(pt) else plot_labels[[pt]]
    }

    # Requirements check helpers
    req_missing <- function(reqs) {
      missing <- c()
      if ("user_clones" %in% reqs && is.null(shared$user_clones)) missing <- c(missing, "Clone Annotation")
      if ("predictions" %in% reqs && is.null(shared$predictions)) missing <- c(missing, "Predictions")
      if ("patient_pred" %in% reqs && is.null(shared$patient_pred)) missing <- c(missing, "Patient Predictions")
      if ("user_response" %in% reqs && is.null(shared$user_response)) missing <- c(missing, "Clinical Response")
      if ("models" %in% reqs && is.null(shared$models)) missing <- c(missing, "Trained Model")
      if ("user_expr" %in% reqs && is.null(shared$user_expr)) missing <- c(missing, "Expression")
      missing
    }
    # boxplot/ROC are satisfiable through the Combination route (predictions +
    # clones + response) even without per-drug patient_pred.
    card_ready <- function(pt) {
      if (pt %in% c("boxplot", "roc")) {
        return(!is.null(shared$user_response) &&
                 ((!is.null(shared$patient_pred)) ||
                  (!is.null(shared$predictions) && !is.null(shared$user_clones))))
      }
      length(req_missing(plot_requirements[[pt]])) == 0
    }
    effective_reqs <- function(pt) {
      reqs <- plot_requirements[[pt]]
      if (pt %in% c("boxplot", "roc") &&
          identical(input$drug_name_common, "Combination")) {
        reqs <- c("predictions", "user_clones", "user_response")
      }
      reqs
    }

    # ---- Data readiness bar ----
    output$readiness_bar <- renderUI({
      items <- list(
        list(key = "Expression",  ok = !is.null(shared$user_expr),     tab = "data"),
        list(key = "Clones",      ok = !is.null(shared$user_clones),   tab = "data"),
        list(key = "Response",    ok = !is.null(shared$user_response), tab = "data"),
        list(key = "Models",      ok = !is.null(shared$models),        tab = "train"),
        list(key = "Predictions", ok = !is.null(shared$predictions),   tab = "predict")
      )
      div(class = "viz-readiness-bar",
        span(class = "viz-readiness-title", icon("clipboard-check"), " Data readiness:"),
        lapply(items, function(it) {
          div(class = paste("readiness-pill", if (it$ok) "ok" else "missing"),
            span(class = paste("readiness-dot", if (it$ok) "green" else "gray")),
            it$key,
            if (it$ok) tags$span(class = "readiness-badge", "ready")
            else tags$button(
              class = "readiness-go",
              onclick = sprintf("Shiny.setInputValue('%s', '%s');",
                                ns("go_tab"), it$tab),
              "load \u2192")
          )
        })
      )
    })
    observeEvent(input$go_tab, {
      if (nzchar(input$go_tab)) bslib::nav_select("navbar", selected = input$go_tab, session = main_session)
    })

    # ---- Plot gallery grid ----
    output$gallery_grid <- renderUI({
      cards <- lapply(gallery_items, function(it) {
        ready <- card_ready(it$id)
        missing <- req_missing(plot_requirements[[it$id]])
        div(class = paste("viz-card", if (ready) "ready" else "missing"),
          title = if (ready) paste(it$title, "- click to draw") else
            paste0("Missing: ", paste(missing, collapse = ", ")),
          # Hovering a card previews its full description in the About panel
          # (without drawing); clicking draws it.
          onmouseover = sprintf("Shiny.setInputValue('%s', '%s');",
                                ns("preview_plot"), it$id),
          onmouseout  = sprintf("Shiny.setInputValue('%s', '');",
                                ns("preview_plot")),
          onclick = sprintf("Shiny.setInputValue('%s', '%s');",
                            ns("pick_plot"), it$id),
          div(class = "viz-card-icon", icon(it$icon)),
          div(class = "viz-card-title", it$title),
          div(class = "viz-card-tagline", it$tagline),
          div(class = "viz-card-status",
            if (ready) tags$span(class = "viz-card-ready", icon("check-circle"), " Ready")
            else tags$span(class = "viz-card-missing",
              icon("circle-exclamation"),
              paste("Needs:", paste(missing, collapse = ", ")))
          )
        )
      })
      div(class = "viz-gallery-grid", cards)
    })

    # Hover preview: hovering a gallery card shows that plot's full
    # description in the About panel WITHOUT drawing it — users learn what a
    # plot is before committing to it. Moving off clears back to the selected
    # plot (or to nothing).
    hover_plot <- reactiveVal(NULL)
    observeEvent(input$preview_plot, {
      hp <- input$preview_plot
      if (is.null(hp) || !nzchar(hp)) hover_plot(NULL) else hover_plot(hp)
    })

    # Clicking a card selects it and draws immediately. When the required
    # data is missing, do NOT select it — otherwise the parameter panel would
    # appear for a plot that cannot be drawn.
    observeEvent(input$pick_plot, {
      pt <- input$pick_plot
      if (is.null(pt) || !nzchar(pt)) return()
      missing <- req_missing(effective_reqs(pt))
      if (length(missing) > 0) {
        showNotification(paste0(plot_labels[[pt]], " needs: ",
                                paste(missing, collapse = ", "),
                                ". Load them first."), type = "warning", duration = 6)
        return()
      }
      viz_selected(pt)
      generate_plot(pt)
    })

    # ROC pair picker (only meaningful when response has >2 groups).
    output$roc_pair_picker <- renderUI({
      grps <- response_groups()
      if (length(grps) <= 2) return(NULL)
      def <- if (all(c("PD", "RD") %in% grps)) c("PD", "RD") else grps[1:2]
      tagList(
        tags$small(class = "text-muted", style = "display: block; margin: 0.4rem 0 0.2rem;",
          icon("info-circle"), " Response has ", length(grps), " groups. Pick two for ROC:"),
        div(class = "viz-size-row",
          div(class = "viz-size-col",
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

    # ---- Parameter panel for the selected plot ----
    # NOTE: choices are built INLINE here instead of via updateSelectizeInput.
    # The panel is rendered on demand (only after a card is clicked), so a
    # server-side update fired while the input does not exist yet would be
    # silently dropped — the dropdowns would come up empty. Building the
    # choices directly in renderUI keeps them correct on first paint and
    # re-freshes them whenever shared data changes (renderUI re-evaluates).
    output$param_panel <- renderUI({
      pt <- viz_selected()
      if (is.null(pt)) return(NULL)
      has_drug  <- pt %in% c("clone_kill", "boxplot", "roc", "umap_viability")
      has_gene  <- pt == "umap_gene"
      has_roc   <- pt == "roc"

      drug_choices <- NULL
      if (!is.null(shared$models)) drug_choices <- names(shared$models)
      else if (!is.null(shared$predictions)) drug_choices <- colnames(shared$predictions)
      combo_choices <- if (length(drug_choices) > 0) {
        c("Combination" = "Combination", setNames(drug_choices, drug_choices))
      } else NULL
      gene_choices <- if (!is.null(shared$user_expr)) rownames(shared$user_expr) else NULL

      tagList(
        div(class = "card animate-fade-in-up viz-param-card",
          div(class = "card-header", icon("sliders"), " Plot Parameters",
            tags$span(class = "text-muted", style = "font-size: 0.78rem; margin-left: 0.6rem;",
              "Changes redraw the plot automatically.")
          ),
          div(class = "card-body",
            fluidRow(
              if (has_drug && pt != "umap_viability") column(6,
                selectizeInput(ns("drug_name_common"), "Drug",
                               choices = combo_choices, selected = "Combination",
                               width = "100%",
                               options = list(maxItems = 1, placeholder = "Select a drug"))
              ),
              if (pt == "umap_viability") column(6,
                selectizeInput(ns("umap_drug"), "Drug",
                               choices = if (length(drug_choices) > 0) setNames(drug_choices, drug_choices) else NULL,
                               selected = if (length(drug_choices) > 0) drug_choices[1] else NULL,
                               width = "100%",
                               options = list(maxItems = 1, placeholder = "Select a drug"))
              ),
              if (has_gene) column(6,
                selectizeInput(ns("umap_gene"), "Gene",
                               choices = gene_choices,
                               selected = if (!is.null(gene_choices) && length(gene_choices) > 0) gene_choices[1] else NULL,
                               width = "100%",
                               options = list(maxItems = 1, placeholder = "Select a gene"))
              ),
              if (has_roc) column(12, uiOutput(ns("roc_pair_picker"))),
              column(12,
                div(class = "viz-size-controls",
                  tags$label(class = "viz-size-label", icon("arrows-alt"), " Plot Size"),
                  div(class = "viz-size-row",
                    div(class = "viz-size-col",
                      tags$span(class = "viz-size-mini-label", "Width"),
                      numericInput(ns("plot_width"), NULL, value = 100, min = 50, max = 200, step = 10)),
                    div(class = "viz-size-col",
                      tags$span(class = "viz-size-mini-label", "Height"),
                      numericInput(ns("plot_height"), NULL, value = 100, min = 50, max = 200, step = 10))
                  ),
                  div(class = "viz-size-row",
                    sliderInput(ns("plot_text_size"), "Text size",
                                min = 60, max = 160, value = 100, step = 10,
                                post = "%", ticks = TRUE)
                  )
                )
              )
            )
          )
        )
      )
    })

    # Parameter changes redraw the current plot automatically (debounced so
    # fast typing/sliding does not fire a plot per keystroke).
    auto_redraw <- reactive({
      list(pt = viz_selected(),
           drug = input$drug_name_common,
           gene = input$umap_gene,
           udrug = input$umap_drug,
           ra = input$roc_group_a,
           rb = input$roc_group_b)
    }) |> debounce(300)
    observeEvent(auto_redraw(), {
      pt <- viz_selected()
      if (is.null(pt)) return()
      missing <- req_missing(effective_reqs(pt))
      if (length(missing) > 0) return()
      generate_plot(pt)
    })

    # ---- Core: generate a plot (normal or spatial) in the background ----
    generate_plot <- function(pt) {
      is_spatial <- pt %in% c("umap_gene", "umap_viability", "umap_clone")

      # Drug: the parameter panel may not have initialized on the FIRST card
      # click, so input$drug_name_common can still be NULL then — fall back to
      # the panel's default "Combination" so the plot matches what the
      # dropdown shows (otherwise the first drug leaks in as a silent default).
      drug_sel <- if (!is.null(input$drug_name_common) && nzchar(input$drug_name_common)) {
        input$drug_name_common
      } else "Combination"
      is_combo <- identical(drug_sel, "Combination")

      # umap_gene: the parameter panel may not have rendered on the FIRST card
      # click, so input$umap_gene can still be NULL/empty then — fall back to
      # the first gene instead of drawing an empty "Select a gene first" plot.
      gene_sel <- if (pt == "umap_gene") {
        if (!is.null(input$umap_gene) && nzchar(input$umap_gene)) {
          input$umap_gene
        } else if (!is.null(shared$user_expr) && nrow(shared$user_expr) > 0) {
          rownames(shared$user_expr)[1]
        } else ""
      } else ""

      # Cache hit: identical request was generated before — reuse instantly.
      ck <- if (is_spatial) paste(pt, gene_sel, input$umap_drug, sep = "\u0001")
            else paste(pt, is_combo, drug_sel,
                       input$roc_group_a, input$roc_group_b, sep = "\u0001")
      cached <- get0(ck, envir = plot_cache, inherits = FALSE)
      if (!is.null(cached)) {
        current_plot(cached$plot)
        current_plot_type(pt)
        current_plot_key(ck)
        showNotification(paste0(plot_label(pt), " loaded from cache"), type = "message", duration = 3)
        return()
      }

      # Busy guard: another plot is mid-flight; skip silently. A subsequent
      # parameter change / card click will regenerate anyway, so there is no
      # need to nag the user while fast-switching plots.
      if (viz_busy()) return()
      viz_busy(TRUE)
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Generating plot..."),
          p(class = "text-muted", "Drawing plot...")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()

      pd <- shared$prepared_data
      if (is_spatial) {
        umap_coords <- pd$umap_coords
        if (is.null(umap_coords)) {
          viz_busy(FALSE); w$hide()
          showNotification("No 2D embedding available (clustering was skipped or not run). Re-run with Seurat clustering to enable spatial plots.", type = "warning", duration = 8)
          return()
        }
        if (is.null(shared$user_clones)) {
          viz_busy(FALSE); w$hide()
          showNotification("No clone annotation found. Run Seurat clustering first.", type = "warning")
          return()
        }
        # umap_gene: ship only the selected gene's per-cell vector (small),
        # never the whole expression matrix.
        gene_vec <- NULL
        if (pt == "umap_gene" && nzchar(gene_sel) &&
            !is.null(shared$user_expr) && gene_sel %in% rownames(shared$user_expr)) {
          gene_vec <- setNames(as.numeric(shared$user_expr[gene_sel, ]),
                               colnames(shared$user_expr))
        }
        jobid <- tryCatch(
          submit_session_task(shared, "plot", list(
            plot_type = pt,
            data = list(
              # Spatial plots need only the slices listed here — carrying the
              # full predictions/response matrices made params.rds serialize
              # slowly (Gene Expression felt slow for no reason).
              predictions   = if (pt == "umap_viability") shared$predictions else NULL,
              user_clones   = shared$user_clones,
              user_response = NULL,
              patient_pred  = NULL,
              umap_gene_expr = gene_vec,
              prepared = list(
                umap_coords = umap_coords,
                clone_viability_template = if (!is.null(pd)) pd$clone_viability_template else NULL
              )
            ),
            params = list(
              drug = "",
              combo = FALSE,
              roc_group_a = NULL,
              roc_group_b = NULL,
              umap_gene = gene_sel,
              umap_drug = input$umap_drug
            )
          )),
          error = function(e) { viz_busy(FALSE); w$hide(); NULL }
        )
      } else {
        jobid <- tryCatch(
          submit_session_task(shared, "plot", list(
            plot_type = pt,
            data = list(
              predictions   = shared$predictions,
              user_clones   = shared$user_clones,
              user_response = shared$user_response,
              patient_pred  = shared$patient_pred,
              prepared = list(
                umap_coords = NULL,
                clone_viability_template = if (!is.null(pd)) pd$clone_viability_template else NULL
              )
            ),
            params = list(
              drug = drug_sel,
              combo = is_combo,
              roc_group_a = input$roc_group_a,
              roc_group_b = input$roc_group_b
            )
          )),
          error = function(e) { viz_busy(FALSE); w$hide(); NULL }
        )
      }
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
            showNotification(paste0(plot_label(pt), " generated successfully"), type = "message")
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
    }

    # --- Interactive rendering via ggiraph (unchanged) ---
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

    output$main_plot_host <- renderUI({
      if (is.null(current_plot())) NULL else ggiraph::girafeOutput(ns("main_plot"))
    })

    output$main_plot <- ggiraph::renderGirafe({
      req(current_plot())
      sz <- plot_size()
      if (identical(current_plot_type(), "clone_kill")) {
        sz$h <- round(sz$h * 1.3)
      }
      wkey <- paste(current_plot_key(), text_scale(), sz$w, sz$h, sep = "\u0001")
      cached_w <- get0(wkey, envir = widget_cache, inherits = FALSE)
      if (!is.null(cached_w)) return(cached_w)
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
          "Click a card in the Plot Gallery to draw it.",
          br(),
          span(style = "font-size: 0.82rem; opacity: 0.7;",
            "Check the data readiness bar above for what is loaded")
        )
      } else {
        NULL
      }
    })

    # Inline plot info in header: title + generated-at time (UTC).
    output$plot_info_inline <- renderUI({
      req(current_plot())
      pt <- current_plot_type()
      if (is.null(pt)) return(NULL)
      utc_now <- format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M UTC")
      tags$span(class = "viz-info-inline",
        strong(plot_label(pt)),
        tags$span(class = "text-muted", style = "margin-left: 0.5rem; font-size: 0.78rem;",
          "Generated at ", utc_now)
      )
    })

    # Compact "About this plot" — hovering a gallery card previews its
    # description here; otherwise it follows the selected plot.
    output$plot_explanation <- renderUI({
      pt <- if (!is.null(hover_plot())) hover_plot() else viz_selected()
      if (is.null(pt) || !nzchar(pt)) return(NULL)
      is_preview <- !is.null(hover_plot()) && !identical(hover_plot(), viz_selected())
      item <- gallery_items[[match(pt, vapply(gallery_items, function(x) x$id, ""))]]
      if (is.null(item)) return(NULL)
      info <- list(
        clone_dist = "Shows clone proportions within each patient. With global clustering, the same clone (color) genuinely spans patients. With clone-level input using per-patient labels (e.g. c1/c2/c3), the same color across patients does not imply the same clone origin — it is a shared category label.",
        clone_kill = "Displays predicted drug viability scores for each clone within each patient. Taller = higher viability = more resistant; shorter = lower viability = more sensitive. Use this to identify which subclones are most affected by treatment.",
        roc = "Receiver Operating Characteristic curve evaluating how well the predicted viability score discriminates Responders (R) from Non-Responders (NR). AUC closer to 1.0 indicates better prediction accuracy.",
        boxplot = "Patient-level predicted viability (z-score) per clinical response group. Two groups (R/NR): p-value is a one-sided Wilcoxon test asking whether responders are predicted more sensitive (lower viability). Three or four groups (e.g. TN/RD/PD): pairwise two-sided Wilcoxon tests with BH (FDR) correction — the adj. p on each bracket tells you which pairs differ. Five or more groups: a single Kruskal-Wallis omnibus test. Non-parametric tests keep the p-value valid even with small patient sample sizes.",
        umap_gene = "Expression level of a selected gene (0-1 normalized, 5th-95th percentile) across all single cells in the 2D embedding. Grey = no/low expression, red = high expression. Use this to examine how a gene's expression relates to the transcriptional subclone landscape.",
        umap_viability = "Predicted drug viability (z-score: 0 = cohort mean across clones, higher = more resistant) across all single cells in the 2D embedding. Red = resistant, blue = sensitive, white = neutral. Use this to identify which regions of the embedding (i.e., which subclones) are most affected by a given drug.",
        umap_clone = "Shows where each clone sits in the 2D embedding space, with every cell colored by its clone. Use this to see the spatial layout of each transcriptional subclone and how it relates to the rest of the tumor."
      )
      requires_map <- list(
        clone_dist = "Clone annotation (Seurat clustering or clone-level input)",
        clone_kill = "Clone-level Predictions + Clone Annotation",
        roc = "Patient-level Predictions + Clinical Response",
        boxplot = "Patient-level Predictions + Clinical Response",
        umap_gene = "Clone-level Predictions + Expression Matrix + 2D embedding coordinates",
        umap_viability = "Clone-level Predictions + 2D embedding coordinates (auto-generated by Seurat)",
        umap_clone = "Clone annotation + 2D embedding coordinates (auto-generated by Seurat)"
      )
      equal_weight_note <- NULL
      if (pt == "clone_dist" && !is.null(shared$prepared_data) &&
          identical(shared$prepared_data$reduction_method, "none") &&
          !isTRUE(shared$mapping_has_count)) {
        equal_weight_note <- div(class = "info-box",
          style = "border-left-color: var(--warning, #d97706); margin: 0.6rem 0 0.8rem 0; font-size: 0.82rem; line-height: 1.5;",
          icon("triangle-exclamation"),
          " Clone-level input without a count column: proportions are equal (1/n per clone). Add a count column with real cell numbers in the mapping file to show true proportions.")
      }
      div(class = "card viz-explanation-card",
        div(class = "card-header",
          icon("circle-info"), " About This Plot",
          if (is_preview) tags$span(class = "viz-preview-tag",
            icon("eye"), " preview — click the card to draw")
        ),
        div(class = "card-body",
          h6(strong(if (item$spatial) spatial_plot_label(pt) else item$title)),
          p(class = "text-muted", style = "font-size: 0.85rem; line-height: 1.5; margin-bottom: 0.4rem;",
            info[[pt]]),
          tags$span(class = "viz-explanation-req",
            icon("clipboard-check", style = "font-size: 0.75rem;"),
            tags$span(style = "font-size: 0.78rem; font-weight: 600;", "Requires: "),
            tags$span(style = "font-size: 0.78rem;", requires_map[[pt]])
          ),
          equal_weight_note
        )
      )
    })

    # Inline download buttons in header (unchanged)
    output$plot_download_inline <- renderUI({
      req(current_plot())
      tagList(
        downloadButton(ns("download_png"), "PNG", class = "btn-outline-primary btn-sm"),
        downloadButton(ns("download_pdf"), "PDF", class = "btn-outline-primary btn-sm")
      )
    })

    output$download_png <- downloadHandler(
      filename = function() paste0("perceptionx_plot_", gsub("[-: ]", "", format(Sys.time(), tz = "UTC")), ".png"),
      content = function(file) {
        p <- scaled_plot()
        sz <- plot_size()
        if (inherits(p, "gtable") || inherits(p, "grob") || inherits(p, "gTree") ||
            inherits(p, "plotly") || inherits(p, "htmlwidget")) {
          writeLines("This plot type cannot be exported as PNG.", file)
          return()
        }
        tryCatch({
          ragg::agg_png(file, width = sz$w, height = sz$h, res = 100)
          print(p)
          dev.off()
        }, error = function(e) {
          if (grDevices::dev.cur() > 1) grDevices::dev.off()
          png(file, width = sz$w, height = sz$h, res = 100)
          print(p)
          dev.off()
        })
      }
    )

    output$download_pdf <- downloadHandler(
      filename = function() paste0("perceptionx_plot_", gsub("[-: ]", "", format(Sys.time(), tz = "UTC")), ".pdf"),
      content = function(file) {
        p <- scaled_plot()
        if (inherits(p, "gtable") || inherits(p, "grob") || inherits(p, "gTree") ||
            inherits(p, "plotly") || inherits(p, "htmlwidget")) {
          writeLines("This plot type cannot be exported as PDF.", file)
          return()
        }
        grDevices::pdf(file, width = 12, height = 9)
        print(p)
        grDevices::dev.off()
      }
    )

    # Pre-render every viz output so switching to this tab NEVER shows a
    # blank page: by default Shiny suspends hidden outputs and renders them
    # on demand (the tab would look empty for a request round-trip or two).
    # These outputs are cheap when empty (gallery grid, status placeholder,
    # NULL param panel), so computing them eagerly costs nothing and makes
    # the tab appear fully rendered the moment it is activated.
    for (out_name in c("readiness_bar", "gallery_grid", "param_panel",
                       "plot_status", "main_plot_host", "plot_info_inline",
                       "plot_download_inline", "plot_explanation")) {
      outputOptions(output, out_name, suspendWhenHidden = FALSE)
    }
  })
}
