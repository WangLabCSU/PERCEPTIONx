# Help Module — Reference Documentation
mod_help_ui <- function(id) {
  ns <- NS(id)
  # Keep the in-app citation in sync with the installed package version
  # (falls back to the current release when running from source).
  pkg_ver <- tryCatch(as.character(utils::packageVersion("PERCEPTIONx")),
                      error = function(e) "0.3.0")
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
              actionButton(ns("nav_overview"), "Walkthroughs", class = "help-nav-btn"),
              actionButton(ns("nav_drugs"), "Supported Drugs", class = "help-nav-btn"),
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
        # Walkthroughs — pick Demo or Advanced; each opens a centered,
        # slide-style guided tour (game-tutorial style prev/next paging).
        div(id = ns("section_overview"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("graduation-cap")),
            div(
              h5("Guided Walkthroughs"),
              p(class = "help-subtitle", "New to the app? Pick a tour below. It opens as an interactive slideshow in the center of the screen.")
            )
          ),

          div(class = "tw-entry-grid",
            div(class = "tw-entry tw-entry-demo",
              div(class = "tw-entry-top",
                div(class = "tw-entry-ic", icon("wand-magic-sparkles")),
                div(class = "tw-entry-title",
                  strong("Try the demo"),
                  span("~1 minute \u00b7 no files needed")
                )
              ),
              p("Load the built-in demo data, run the predictions, then draw and export the plots - nothing to upload or train."),
              actionButton(ns("tour_demo"), "Start the Demo Tour", class = "btn-tw btn-tw-demo", width = "100%")
            ),
            div(class = "tw-entry tw-entry-adv",
              div(class = "tw-entry-top",
                div(class = "tw-entry-ic", icon("flask-vial")),
                div(class = "tw-entry-title",
                  strong("Analyze your own data"),
                  span("~2 minutes \u00b7 uploads required")
                )
              ),
              p("Walk through uploading expression / mapping / response files, (optional) training, prediction and how to read the plots."),
              actionButton(ns("tour_advanced"), "Start the Advanced Tour", class = "btn-tw btn-tw-adv", width = "100%")
            )
          ),

          p(class = "text-muted", style = "font-size: 0.8rem; margin-top: 0.6rem; margin-bottom: 0;",
            icon("book"), " Detailed reference for every parameter and format sits in the sections below.")
        ),

        # Supported Drugs - searchable catalogue with the pre-trained 44 marked
        div(id = ns("section_drugs"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("pills")),
            div(
              h5("Supported Drugs"),
              p(class = "help-subtitle", "Find a compound, see if a model is ready, or train it yourself")
            )
          ),
          p("This is the drug screen shipped with our DepMap data (PRISM): ",
            strong("1448 compounds"), ", most of them research tool molecules rather than clinical drugs. ",
            "The ", strong("44 pre-trained"), " ones are clinically relevant and come with ready-made models. ",
            "Any other compound can be trained on the Train tab once DepMap is loaded."),
          div(class = "help-drug-legend",
            span(class = "help-drug-chip", "compound in the DepMap screen"),
            span(class = "help-drug-badge", "pre-trained"),
            span(class = "help-drug-badge help-drug-badge-demo", "demo")
          ),
          div(class = "help-drug-controls",
            textInput(ns("drug_search"), NULL, placeholder = "Search compounds...",
                      width = "260px"),
            radioButtons(ns("drug_filter"), NULL, inline = TRUE,
                         choices = c("Pre-trained (44)" = "pt",
                                     "Demo" = "demo",
                                     "All compounds" = "all"),
                         selected = "pt")
          ),
          uiOutput(ns("drug_results"))
        ),

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
              p("Downloaded and loaded via the ", code("Download & Load"), " button on the Data tab (or a pre-downloaded .RDS upload). Includes bulk expression, single-cell expression, and drug response (AUC) data from DepMap."),
              tags$span(class = "status-badge unloaded", "Manual load")
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
              p("File with columns: ", code("cell_id"), " and ", code("patient_id"), " (case-insensitive). Maps each cell to its patient. Clones are auto-detected via Seurat clustering. Accepts CSV / TSV / TXT / Excel / RDS."),
              tags$span(class = "status-badge unloaded", "User upload")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--primary-light);", icon("stethoscope")),
              strong("Clinical Response"),
              p("(Optional) File with columns: ", code("patient"), ", ", code("response"), " (Responder/Non-responder, case-insensitive). Required for ROC curves and boxplots. Accepts CSV / TSV / TXT / Excel / RDS."),
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
              p(code("glmnet"), " (elastic net, recommended) with 3-fold CV, or ", code("rf"), " (random forest).")
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
              p("Weighted maximum across clones, emphasizing the most resistant clones (recommended).")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--primary);", icon("balance-scale")),
              strong("weighted_average"),
              p("Weighted average across all clones, weighted by clone proportion.")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--error);", icon("arrow-down")),
              strong("min"),
              p("Takes the most sensitive clone's prediction (lowest viability, pessimistic).")
            ),
            div(class = "help-data-card",
              div(class = "help-data-icon", style = "color: var(--success-light);", icon("arrow-up")),
              strong("max"),
              p("Takes the most resistant clone's prediction (highest viability, optimistic).")
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
              div(class = "help-viz-icon", icon("chart-pie")),
              div(
                strong("Clone Distribution"),
                p("Stacked bar chart showing clone proportions per patient."),
                tags$span(class = "status-badge unloaded", "Needs: Clones")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("chart-column")),
              div(
                strong("Clone Viability Lollipop"),
                p("Lollipop chart of predicted viability per clone."),
                tags$span(class = "status-badge unloaded", "Needs: Predictions + Clones")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("chart-line")),
              div(
                strong("ROC Curve"),
                p("Receiver Operating Characteristic curve for response prediction."),
                tags$span(class = "status-badge unloaded", "Needs: Patient Pred + Response")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("chart-gantt")),
              div(
                strong("Response Boxplot"),
                p("Boxplot comparing predicted viability between responders and non-responders."),
                tags$span(class = "status-badge unloaded", "Needs: Patient Pred + Response")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("map")),
              div(
                strong("Clone Identity"),
                p("2D embedding (UMAP/t-SNE) colored by clone membership."),
                tags$span(class = "status-badge unloaded", "Needs: Clones + Embedding")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("flask")),
              div(
                strong("Drug Viability"),
                p("2D embedding (UMAP/t-SNE) colored by predicted drug viability."),
                tags$span(class = "status-badge unloaded", "Needs: Predictions + Embedding")
              )
            ),
            div(class = "help-viz-card",
              div(class = "help-viz-icon", icon("dna")),
              div(
                strong("Gene Expression"),
                p("2D embedding (UMAP/t-SNE) colored by selected gene expression."),
                tags$span(class = "status-badge unloaded", "Needs: Predictions + Expression + Embedding")
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
                "Currently, PERCEPTIONx trains on DepMap data. You can use the trained model to predict on your own single-cell data after rank normalization."
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
            ),
            div(class = "help-faq-item",
              div(class = "help-faq-q", icon("question-circle"), " Why doesn't the interface freeze during training / clustering / prediction?"),
              div(class = "help-faq-a",
                "All heavy computation (model training, Seurat clustering, prediction, plot math) runs in ", strong("background worker processes"), ", not in the interface. The UI polls and shows progress, so other pages stay responsive. If a worker ever stops unexpectedly, the app shows \"Background worker stopped\" instead of spinning forever — just submit again."
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
              p(class = "help-subtitle", "If you use this package, please cite both the package and the original PERCEPTION study:")
            )
          ),
          div(class = "citation-box",
            strong("Jia Ding."),
            " PERCEPTIONx: Personalized Drug Response Prediction from Single-Cell Transcriptomics. R package version ",
            pkg_ver, ". ",
            a("https://github.com/WangLabCSU/PERCEPTIONx",
              href = "https://github.com/WangLabCSU/PERCEPTIONx",
              target = "_blank"),
            br(), br(),
            strong("Sinha, S., Vegesna, R., Mukherjee, S."),
            " et al. PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. ",
            em("Nat Cancer"), " 5, 938-952 (2024). ",
            a("DOI: 10.1038/s43018-024-00756-7",
              href = "https://doi.org/10.1038/s43018-024-00756-7",
              target = "_blank")
          )
        )
      )
    )
  )
}

mod_help_server <- function(id, main_session, shared) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- Supported Drugs: searchable catalogue (44 pre-trained marked) ----
    # Debounce the search box so the (possibly 1,448-entry) grid is only
    # re-rendered after the user pauses typing.
    search_q <- debounce(reactive({
      if (is.null(input$drug_search)) {
        ""
      } else {
        trimws(input$drug_search)
      }
    }), 250)

    output$drug_results <- renderUI({
      pt_all <- tolower(PERCEPTIONx:::perception_all_drugs)
      demo <- c("abemaciclib", "erlotinib")
      md <- shared$depmap_meta
      # The full catalogue is shipped in the package; a loaded DepMap refreshes
      # it to the exact contents of the file the user has.
      live <- !is.null(md) && !is.null(md$drugs) && length(md$drugs) > 0
      all <- if (live) {
        sort(unique(trimws(as.character(md$drugs))), na.last = TRUE)
      } else {
        PERCEPTIONx:::perception_depmap_drugs
      }
      all <- all[order(tolower(all))]

      ftype <- if (is.null(input$drug_filter)) "pt" else input$drug_filter
      q <- search_q()
      sel <- switch(ftype,
        pt   = all[tolower(all) %in% pt_all],
        demo = all[tolower(all) %in% demo],
        all  = all)
      if (nzchar(q)) sel <- sel[grepl(tolower(q), tolower(sel), fixed = TRUE)]

      if (length(sel) == 0) {
        return(div(class = "text-muted", style = "text-align: center; padding: 1.2rem;",
                   icon("magnifying-glass"), " No compound matches your search."))
      }

      mk_badges <- function(d) {
        is_pt <- tolower(d) %in% pt_all
        is_demo <- tolower(d) %in% demo
        # In the Pre-trained (44) view every entry is pre-trained, so a badge
        # per cell would be pure noise; only flag the demo pair there. In the
        # All / search views mark both kinds so the 44 stand out.
        bd <- list()
        if (ftype != "pt" && is_pt) {
          bd <- c(bd, list(tags$span(class = "help-drug-badge", "pre-trained")))
        }
        if (is_demo) {
          bd <- c(bd, list(tags$span(class = "help-drug-badge help-drug-badge-demo", "demo")))
        }
        if (length(bd)) div(class = "help-drug-cell-badges", bd)
      }

      items <- lapply(sel, function(d) {
        div(class = "help-drug-cell", title = d,
          span(class = "help-drug-cell-name", d),
          mk_badges(d))
      })
      count_note <- if (ftype == "all" && !nzchar(q)) {
        sprintf("%s compounds in the DepMap screen (%d pre-trained).",
                prettyNum(length(sel), big.mark = ","), sum(tolower(sel) %in% pt_all))
      } else {
        if (length(sel) == 1) {
          plural <- ""
        } else {
          plural <- "es"
        }
        sprintf("%d match%s.", length(sel), plural)
      }

      tagList(
        div(class = "help-drug-count", icon("filter"), " ", count_note),
        div(class = "help-drug-list",
          div(class = "help-drug-grid", items))
      )
    })

    # ---- Guided walkthroughs (Demo / Advanced) ----
    # One slide: optional screenshot(s), kicker, title, paragraphs, bullets,
    # an amber note and/or a red warning. Inline emphasis helpers:
    #   tw_b()   -> bold key term
    #   tw_r()   -> bold + red (only for the single most important point)
    # Trim stray leading/trailing whitespace from the highlighted text, so a
    # span followed by punctuation can never render "word ." with a gap.
    tw_b <- function(...) {
      ch <- lapply(list(...), function(x) {
        if (is.character(x)) trimws(x) else x
      })
      tags$span(class = "tw-b", ch)
    }
    tw_r <- function(...) {
      ch <- lapply(list(...), function(x) {
        if (is.character(x)) trimws(x) else x
      })
      tags$span(class = "tw-r", ch)
    }
    # Small example table card (mimics the Data-tab previews)
    tw_tbl <- function(title, note, headers, rows) {
      div(class = "tw-tbl-card",
        strong(title),
        tags$table(
          tags$thead(tags$tr(lapply(headers, tags$th))),
          tags$tbody(lapply(rows, function(r) tags$tr(lapply(r, tags$td))))
        ),
        span(class = "tw-tbl-note", note)
      )
    }

    tw_slide <- function(img, kicker, title, paras, points = NULL,
                         note = NULL, warn = NULL, chips = NULL,
                         layout = c("lr", "tb"), big = NA) {
      layout <- match.arg(layout)
      imgs <- if (is.null(img)) character(0) else img
      imgs <- imgs[nzchar(imgs)]
      has_img <- length(imgs) > 0
      figs <- if (has_img) {
        shots <- lapply(imgs, function(src) {
          div(class = "tw-shot", img(src = src, alt = title))
        })
        if (length(shots) == 1) shots[[1]] else div(class = "tw-figs", shots)
      }
      # Header is a fixed zone ABOVE the body so every slide keeps its
      # step caption + title pinned to the same top-left spot (no vertical
      # drift between sparse and dense pages).
      slide_head <- div(class = "tw-slide-head",
        div(class = "tw-slide-kicker", kicker),
        h4(title)
      )
      body_text <- tagList(
        lapply(paras, function(t) if (is.character(t)) p(t) else t),
        if (!is.null(chips)) div(class = "tw-chips",
          lapply(chips, function(cp) {
            if (is.character(cp)) tags$span(class = "tw-chip", cp) else cp
          })),
        if (!is.null(points)) tags$ul(class = "tw-points",
          lapply(points, function(x) if (is.character(x)) tags$li(x) else x)),
        # A quiet bottom line (no boxed "warning" component).
        if (!is.null(warn)) div(class = "tw-warn", icon("circle-info"), warn)
      )
      text_col <- div(class = "tw-slide-text", body_text)
      body_class <- if (!has_img) {
        "tw-slide-body tw-slide-body-text"
      } else if (layout == "tb") {
        "tw-slide-body tw-tb"
      } else {
        "tw-slide-body"
      }
      body_children <- if (has_img && layout == "tb") list(figs, text_col) else list(text_col, figs)
      # Only slides explicitly marked `big` get the enlarged base font (the
      # auto word-count heuristic is disabled to avoid any runtime surprises).
      is_big <- !is.na(big) && isTRUE(big)
      div(class = if (is_big) "tw-slide tw-big" else "tw-slide",
        slide_head,
        div(class = body_class, body_children))
    }

    # Builds and opens the centered slideshow modal.
    tw_walkthrough <- function(root_id, header, pages) {
      slides_html <- lapply(seq_along(pages), function(i) {
        htmltools::tagAppendAttributes(
          pages[[i]],
          class = "tw-slide-inner",
          style = if (i == 1) NULL else "display:none;"
        )
      })
      close_req <- ns("tw_close")
      body <- div(
        id = root_id, class = "tw-root",
        div(class = "tw-topbar",
          div(class = "tw-topbar-title", icon("graduation-cap"), header),
          span(class = "tw-close", title = "Close",
               onclick = sprintf("Shiny.setInputValue('%s', 1, {priority:'event'});", close_req),
               icon("xmark"))
        ),
        div(class = "tw-stages",
          div(class = "tw-count", ""),
          slides_html
        ),
        div(class = "tw-nav",
          tags$button(type = "button", class = "tw-btn tw-prev",
                      icon("chevron-left"), " Back"),
          div(class = "tw-dots",
            lapply(seq_along(pages), function(k) {
              span(class = "tw-dot", `data-i` = k - 1L)
            })
          ),
          tags$button(type = "button", class = "tw-btn tw-next tw-btn-primary",
                      "Next ", icon("chevron-right"))
        ),
        tags$script(HTML(sprintf(
          "(function(){\n
           var root=document.getElementById('%1$s');\n
           if(!root) return;\n
           var slides=root.querySelectorAll('.tw-slide-inner');\n
           var dots=root.querySelectorAll('.tw-dot');\n
           var prev=root.querySelector('.tw-prev');\n
           var next=root.querySelector('.tw-next');\n
           var count=root.querySelector('.tw-count');\n
           var stages=root.querySelector('.tw-stages');\n
           var i=0,n=slides.length;\n
           function upd(){\n
             for(var k=0;k<n;k++){ slides[k].style.display=(k===i)?'flex':'none'; }\n
             for(var k=0;k<dots.length;k++){ dots[k].classList.toggle('active', k===i); }\n
             prev.disabled=(i===0);\n
             next.innerHTML=(i===n-1)?'Finish':'Next';\n
             count.textContent=(i+1)+' / '+n;\n
             stages.scrollTop=0;\n
           }\n
           for(var k=0;k<dots.length;k++){(function(k){ dots[k].onclick=function(){ i=k; upd(); }; })(k);}\n
           prev.onclick=function(){ if(i>0){ i--; upd(); } };\n
           next.onclick=function(){ if(i<n-1){ i++; upd(); } else { Shiny.setInputValue('%2$s', 1, {priority:'event'}); } };\n
           upd();\n
           // htmltools pretty-prints tags with newlines between children; a
           // browser collapses each newline into a space, which would put a
           // gap before punctuation ('word .'). Drop whitespace-only text
           // nodes so inline emphasis sits flush with the following character.\n
           var walker=document.createTreeWalker(stages,NodeFilter.SHOW_TEXT,null);\n
           var drop=[], nn;\n
           while((nn=walker.nextNode())){ if(nn.nodeValue.trim()==='') drop.push(nn); }\n
           for(var d=0;d<drop.length;d++){ drop[d].parentNode.removeChild(drop[d]); }\n
           })();",
          root_id, close_req
        )))
      )
      showModal(modalDialog(
        title = NULL, footer = NULL, size = "xl", easyClose = FALSE, fade = TRUE,
        body
      ))
    }

    demo_pages <- list(
      tw_slide(NULL,
        "Step 1 of 3 \u00b7 Data tab", "Load Demo Data",
        list(p("The demo data is included with the app and needs nothing else to run."),
          tags$ul(class = "tw-points",
            tags$li("49 genes"),
            tags$li("400 cells"),
            tags$li("20 patients"),
            tags$li("Seurat clustering already done"),
            tags$li("Two demo drug models: ", tw_b("abemaciclib"), " and ", tw_b("erlotinib"))),
          p(class = "tw-try", "Try it out:"),
          div(class = "tw-center-btn",
            tags$button(type = "button", class = "btn btn-demo btn-lg",
              onclick = paste0(
                "document.getElementById('demo-overlay').style.display='flex';",
                "Shiny.setInputValue('data-load_demo', 1, {priority:'event'});",
                "var b=this; b.style.opacity='0.6';",
                "setTimeout(function(){b.style.opacity='';},2000);"),
              icon("play"), " Load Demo Data"))),
        big = TRUE,
        warn = "The demo models are SIMULATED to illustrate the workflow. Their predictions must NOT be interpreted as real drug response."),
      tw_slide("figures/predict.png",
        "Step 2 of 3 \u00b7 Predict tab", "Run Prediction: Your Result is Here",
        list(p("On the Predict tab, click ", tw_r("Run Prediction"), ". You get a clone-level table and heatmap, plus patient-level predictions."),
             p("Prediction runs on the two active demo models, ", tw_b("abemaciclib"), " and ", tw_b("erlotinib"), "."),
             p("For the Aggregation Mode, use ", tw_b("Weighted Max (recommended)"), ". The other modes are explained in the Help page."))),
      tw_slide("figures/plot_gallery.png",
        "Step 3 of 3 \u00b7 Visualize tab", "Every plot is now unlocked",
        list(p("Click any plot card to draw it."),
             p("Every plot shows its explanation in the ", tw_b("About This Plot"), " card above it."),
             p("Adjust the figure size, text and DPI before exporting."),
             p("Plots that compare drugs come with a drug switcher, set to ", tw_r("Combination"),
               " by default (the ensemble prediction across all drugs).")))
    )

    adv_pages <- list(
      tw_slide(NULL,
        "Step 1 of 7 \u00b7 Files", "The Three Files You Bring",
        list(p("Bring three tables in the shapes shown below: ", tw_r("expression matrix"),
               " and ", tw_r("cell-patient mapping"), " are required; clinical response is optional (needed for ROC curves and boxplots)."),
          div(class = "tw-tbl-row",
            tw_tbl("1 - Expression matrix",
                   "Rows = genes, columns = cells. Cell-level (per cell) or clone-level (per clone).",
                   c("gene", "CELL_0001", "CELL_0002", "CELL_0003"),
                   list(c("TP53", "2.1", "0.4", "5.6"),
                        c("EGFR", "3.0", "3.2", "1.1"),
                        c("MYC", "0.9", "6.7", "4.2"))),
            tw_tbl("2 - Cell-patient mapping",
                   "One row per cell. cell_id + patient_id required (case-insensitive).",
                   c("cell_id", "patient_id"),
                   list(c("CELL_0001", "PAT_001"),
                        c("CELL_0002", "PAT_001"),
                        c("CELL_0003", "PAT_002"))),
            tw_tbl("3 - Clinical response",
                   "Optional, needed for ROC curves and boxplots. patient + response.",
                   c("patient", "response"),
                   list(c("PAT_001", "Responder"),
                        c("PAT_002", "Non-responder"))))
          )),
      tw_slide("figures/seurat.png",
        "Step 2 of 7 \u00b7 Data tab", "Run Seurat Clustering",
        list(p("We support two input types: cell-level and clone-level:"),
             p(tw_b("Cell-level input: "), " ", tw_r("Run Seurat Clustering"),
               " once, do not forget this step! It detects your clones and builds the 2D map the spatial plots need."),
             p(tw_b("Clone-level input:"), " Seurat is not required, your clones are already defined."),
             p("Check the preview tabs under the tables to confirm your files loaded as expected.")),
        layout = "tb"),
      tw_slide("figures/data_management.png",
        "Step 3 of 7 \u00b7 Models", "Pre-trained Models or Train Your Own",
        list(p("If your drug is one of the ", tw_r("44 pre-trained"),
               " models shipped with the app, it is one click away. Otherwise you train it yourself."),
          div(class = "tw-choice-row",
            div(class = "tw-choice",
              div(class = "tw-choice-title", icon("circle-check"), "In the 44"),
              p("Data tab > Pre-trained Models > select drugs (multi-select) > Download & Load.")),
            div(class = "tw-choice",
              div(class = "tw-choice-title", icon("hammer"), "Not in the 44"),
              p(tw_r("Download & Load DepMap"), " > train it on the Train tab."))),
          p("Every model you download or train lands in ", tw_r("Loaded Models Management"), " on the Data tab, where you can tick which ones stay active."),
          p("Supported Drugs on the Help page lists the full catalogue, with the 44 marked.")),
        layout = "lr"),
      tw_slide("figures/model_roc.png",
        "Step 4 of 7 \u00b7 Train tab (outside the 44)", "Train and Evaluate Your Own Model",
        list(p("If your drug is not in the 44 pre-trained models, check the Supported Drugs list in the Help page. If it is there, download the DepMap data we provide and train the drug on the Train tab."),
             p("Training generates a ", tw_b("Validation ROC"), " and a ", tw_b("Model Performance"), " plot, so you can judge how accurate the model is.")),
        layout = "tb"),
      tw_slide("figures/predict.png",
        "Step 5 of 7 \u00b7 Predict tab", "Run Prediction: Your Result is Here",
        list(p("On the Predict tab, click ", tw_r("Run Prediction"), ". You get a clone-level table and heatmap, plus patient-level predictions."),
             p("Prediction uses the active models and prepared expression from the Loaded Models Management."),
             p("For the Aggregation Mode, use ", tw_b("Weighted Max (recommended)"), ". The other modes are explained in the Help page."))),
      tw_slide("figures/plot_gallery.png",
        "Step 6 of 7 \u00b7 Visualize tab", "Dive Deeper with More Plots",
        list(p("Click any plot card to draw it."),
             p("Every plot shows its explanation in the ", tw_b("About This Plot"), " card above it."),
             p("Adjust the figure size, text and DPI before exporting."),
             p("Plots that compare drugs come with a drug switcher, set to ", tw_r("Combination"),
               " by default (the ensemble prediction across all drugs)."))),
      tw_slide("figures/help.png",
        "Step 7 of 7 \u00b7 Help & resources", "More details live in the Help page",
        list(p("The Help page still has more for you: aggregation modes, data formats, every parameter, the full drug list and the FAQ."),
            p("🎉 Explore at your own pace, and have fun with your data!")),
        big = TRUE
      )
    )

    observeEvent(input$tour_demo, {
      tw_walkthrough(ns("tour_demo_root"), "Demo walkthrough", demo_pages)
    })
    observeEvent(input$tour_advanced, {
      tw_walkthrough(ns("tour_adv_root"), "Your own data walkthrough", adv_pages)
    })
    # Close (X button or "Finish" on the last slide) of either walkthrough.
    observeEvent(input$tw_close, { removeModal() })

    # Section navigation scroll - use namespaced IDs
    observeEvent(input$nav_overview, {
      session$sendCustomMessage("scroll-to", ns("section_overview"))
    })
    observeEvent(input$nav_drugs, {
      session$sendCustomMessage("scroll-to", ns("section_drugs"))
    })
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
