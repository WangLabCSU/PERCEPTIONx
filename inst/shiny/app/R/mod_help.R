# Help Module — Reference Documentation
mod_help_ui <- function(id) {
  ns <- NS(id)
  # Single source for the drug list (identical to the Data-tab Pre-trained
  # Models selector). Demo only uses two of them - they are badged below.
  pt_drugs <- PERCEPTIONx:::perception_all_drugs
  demo_drugs <- c("abemaciclib", "erlotinib")
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
              p(class = "help-subtitle", "New to the app? Pick a tour \u2014 it opens as an interactive slideshow in the center of the screen")
            )
          ),

          div(class = "tw-entry-grid",
            div(class = "tw-entry",
              div(class = "tw-entry-top",
                div(class = "tw-entry-ic tw-entry-ic-demo", icon("wand-magic-sparkles")),
                div(class = "tw-entry-title",
                  strong("Try the demo"),
                  span("~1 minute \u00b7 no files needed")
                )
              ),
              p("Load the built-in demo data, run the predictions, then draw and export the plots - nothing to upload or train."),
              actionButton(ns("tour_demo"), "Start the demo tour", class = "btn-primary btn-sm", width = "100%")
            ),
            div(class = "tw-entry",
              div(class = "tw-entry-top",
                div(class = "tw-entry-ic tw-entry-ic-adv", icon("flask-vial")),
                div(class = "tw-entry-title",
                  strong("Analyze your own data"),
                  span("~2 minutes \u00b7 uploads required")
                )
              ),
              p("Walk through uploading expression / mapping / response files, (optional) training, prediction and how to read the plots."),
              actionButton(ns("tour_advanced"), "Start the data tour", class = "btn-primary btn-sm", width = "100%")
            )
          ),

          p(class = "text-muted", style = "font-size: 0.8rem; margin-top: 0.6rem; margin-bottom: 0;",
            icon("book"), " Detailed reference for every parameter and format sits in the sections below.")
        ),

        # Supported Drugs - the pre-trained 44 at a glance
        div(id = ns("section_drugs"), class = "help-ref-section",
          div(class = "help-section-header",
            div(class = "help-section-icon", icon("pills")),
            div(
              h5("Supported Drugs"),
              p(class = "help-subtitle", "Which drugs ship with ready-made models - and what to do about the rest")
            )
          ),
          p("These 44 drugs have ", strong("pre-trained models"), " and can be loaded directly from the Pre-trained Models card on the Data tab (multi-select supported). The two badged drugs are the ones used by the demo."),
          div(class = "help-drug-legend",
            span(class = "help-drug-chip", "pre-trained model available"),
            span(class = "help-drug-chip help-drug-chip-demo", "also used by the demo")
          ),
          div(class = "help-drug-grid",
            lapply(pt_drugs, function(d) {
              demo <- d %in% demo_drugs
              span(class = if (demo) "help-drug-chip help-drug-chip-demo" else "help-drug-chip",
                   d,
                   if (demo) tags$span(class = "help-drug-badge", "demo"))
            })
          ),
          div(class = "info-box", style = "margin-top: 0.8rem;",
            icon("flask-vial"),
            " Drug not in the list? ", strong("It can still be trained"), " if your DepMap data contains its response screen - load DepMap, then enter the drug name on the Train tab (comma- or newline-separated).")
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
            " PERCEPTIONx: Personalized Drug Response Prediction from Single-Cell Transcriptomics. R package version 0.2.0. ",
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

mod_help_server <- function(id, main_session) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Guided walkthroughs (Demo / Advanced) ----
    # One slide: optional screenshot(s), kicker, title, paragraphs, bullets,
    # an amber note and/or a red warning. Inline emphasis helpers:
    #   tw_b()   -> bold key term
    #   tw_r()   -> bold + red (only for the single most important point)
    tw_b <- function(...) tags$span(class = "tw-b", ...)
    tw_r <- function(...) tags$span(class = "tw-r", ...)
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
                         note = NULL, warn = NULL) {
      imgs <- if (is.null(img)) character(0) else img
      imgs <- imgs[nzchar(imgs)]
      has_img <- length(imgs) > 0
      figs <- if (has_img) {
        shots <- lapply(imgs, function(src) {
          div(class = "tw-shot", img(src = src, alt = title))
        })
        if (length(shots) == 1) shots[[1]] else div(class = "tw-figs", shots)
      }
      div(class = "tw-slide",
        div(class = "tw-slide-kicker", kicker),
        h4(title),
        div(class = if (has_img) "tw-slide-body" else "tw-slide-body tw-slide-body-text",
          div(class = "tw-slide-text",
            lapply(paras, function(t) if (is.character(t)) p(t) else t),
            if (!is.null(note)) div(class = "tw-note", icon("triangle-exclamation"), note),
            if (!is.null(warn)) div(class = "tw-note tw-note-warn", icon("circle-exclamation"), warn),
            if (!is.null(points)) tags$ul(class = "tw-points",
              lapply(points, function(x) if (is.character(x)) tags$li(x) else x))
          ),
          figs
        )
      )
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
           var i=0,n=slides.length;\n
           function upd(){\n
             for(var k=0;k<n;k++){ slides[k].style.display=(k===i)?'block':'none'; }\n
             for(var k=0;k<dots.length;k++){ dots[k].classList.toggle('active', k===i); }\n
             prev.disabled=(i===0);\n
             next.innerHTML=(i===n-1)?'Finish':'Next';\n
             count.textContent=(i+1)+' / '+n;\n
           }\n
           for(var k=0;k<dots.length;k++){(function(k){ dots[k].onclick=function(){ i=k; upd(); }; })(k);}\n
           prev.onclick=function(){ if(i>0){ i--; upd(); } };\n
           next.onclick=function(){ if(i<n-1){ i++; upd(); } else { Shiny.setInputValue('%2$s', 1, {priority:'event'}); } };\n
           upd();\n
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
        "1 - Data tab", "Load the demo data",
        list(p("On the Data tab, click ", tw_b("Load Demo Data"),
               " - the demo is ready to use immediately: Seurat has already been run, so clones and the 2D map are in place.")),
        list("49 genes x 400 cells across 20 patients",
             tags$li("Two demo drug models: ", tw_b("abemaciclib"), " and ", tw_b("erlotinib")),
             "Seurat clustering already done - nothing else to compute"),
        warn = "The demo models are SIMULATED to illustrate the workflow. Their predictions must NOT be interpreted as real drug response."),
      tw_slide("figures/shiny-predict-heatmap.png",
        "2 - Predict tab", "Run Prediction to see the viability scores",
        list(p("Go to the Predict tab and click ", tw_b("Run Prediction"), ". The right panel fills with the clone-level predictions - a table plus a viability heatmap you can hover over."),
             p("These scores are what every plot you draw next will read from.")),
        list("Clone- and patient-level predictions appear as table + heatmap",
             "Uses the two demo models (abemaciclib, erlotinib)")),
      tw_slide("figures/shiny-clone-distribution.png",
        "3 - Visualize tab", "Every plot is now unlocked",
        list(p("All plot types are ready. Click any card to draw it - the bar chart on the right, for example, is a ", tw_b("Clone Distribution"),
               " plot (clone proportions per patient)."),
             p("Each plot carries an explanation above it, and you can adjust the canvas size, text and DPI below the figure before exporting PNG (300 dpi) or PDF.")),
        list("Click a card - the plot appears instantly",
             "Canvas / text / DPI adjustable; PNG (300 dpi) and PDF export",
             "To repeat these steps on your own data, close this tour and start 'Analyze your own data'"))
    )

    adv_pages <- list(
      tw_slide(NULL,
        "Before you start", "The three files you bring",
        list(p("PERCEPTION compares the clones in your tumor against drug-response models. Prepare three tables - the same shapes you will preview in the Data tab. ",
               tw_b("Only the first two are mandatory.")),
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
                   "Optional - needed for ROC curves and boxplots. patient + response.",
                   c("patient", "response"),
                   list(c("PAT_001", "Responder"),
                        c("PAT_002", "Non-responder"))))
          ),
        list("Cell-level input: single cells - clones will be found by Seurat (next slide)",
             "Clone-level input: expression already averaged per clone - Seurat can be skipped unless you want the 2D map")),
      tw_slide("figures/shiny-data-clustering.png",
        "2 - Data tab", "Run Seurat Clustering",
        list(p(tw_r("Cell-level input - do not forget this step: "),
               "run Seurat Clustering once. Only it detects your clones and builds the 2D map that the spatial plots need."),
             p(tw_b("Clone-level input"), " - Seurat is not required; your clones are already defined."),
             p("Before running, preview the loaded tables with the View buttons (Data Review) and fix column names if needed.")),
        list(tags$li(tw_r("Cell-level: "), "run Seurat Clustering and watch the live progress"),
             "Clone-level: skip Seurat unless you also want the 2D map",
             "Use View under the data tables to inspect everything first")),
      tw_slide(NULL,
        "Models", "Pre-trained models - or train your own",
        list(p("If your drug is one of the ", tw_b("44 pre-trained models"), " that ship with the app, you get it in one click:")),
        list(tags$li(tw_b("In the 44: "), "Data tab -> Pre-trained Models -> select drugs (multi-select) -> Download & Load"),
             tags$li(tw_b("Not in the 44: "), "check that the drug exists in DepMap, load DepMap, then train it yourself (next slide)"),
             "The full list, with the 44 marked, is in the Supported Drugs section below"),
        note = "A drug outside the 44 can still be trained if your DepMap data contains its response screen."),
      tw_slide(c("figures/shiny-validation-roc.png", "figures/shiny-model-performance.png"),
        "3 - Optional - Train tab", "Training a model - the key fields",
        list(p("Training is only needed for drugs outside the 44. Load DepMap first, then on the Train tab:"),
             p(tw_b("Drug"), " - one or several, comma- or newline-separated."),
             p(tw_b("Cancer type"), " - ", tw_b("PanCan gives the best results"), " (it trains on all cancer types - the most and most diverse data, which usually generalises best). Picking one type, e.g. Breast, makes sense when the drug mainly matters for that indication, but you learn from less data."),
             p(tw_b("Genes of interest"), " - ", tw_b("leave empty to use ALL genes (best default)"), ". Fill it only to restrict feature selection to genes you care about."),
             p(tw_b("k_features"), " - how many top-ranked genes the model keeps (default 100)."),
             p(tw_b("Model type"), " - ", tw_b("glmnet (elastic net) is the default and what the paper validated"), "; random forest is a good alternative when you expect non-linear effects."),
             p("When training finishes, judge the model from the card on the right - the Validation ROC and the Performance Curve.")),
        list("Trained models appear on the right and are immediately usable in Prediction")),
      tw_slide("figures/shiny-predict-heatmap.png",
        "4 - Predict tab", "Run Prediction - this is your result",
        list(p("On the Predict tab click ", tw_b("Run Prediction"), ". The right panel is the deliverable: a clone-level table and heatmap plus the patient-level predictions."),
             p("Prediction always uses the ", tw_b("active models"), " from the Data tab (Loaded Models Management) and the expression prepared there - load, upload or train models first, then come back here."),
             p("Choose the ", tw_b("Aggregation Mode"), " - how clone scores collapse into one patient score:")),
        list(tags$li(tw_b("Weighted Max (recommended)"), ": maximum clone viability weighted by clone proportion - the approach validated in the paper"),
             tags$li(tw_b("Weighted Average"), ": weighted average across all clones"),
             tags$li(tw_b("Min / Max"), ": only the most sensitive / the most resistant clone"))),
      tw_slide("figures/shiny-lollipop.png",
        "5 - Visualize tab", "Dive deeper with more plots",
        list(p("All seven plot types are now available. Several - Clone Viability Lollipop, ROC, Boxplot and the UMAP views - carry a ", tw_b("drug switcher"),
               " in the parameter panel. ", tw_b("Its default is Combination"), ", the ensemble prediction across all drugs."),
             p("The lollipop on the right is the Clone Viability plot: one dot per clone, taller = predicted more resistant.")),
        list("Switch the drug to redraw the plot for a single drug",
             "Every plot explains itself in the 'About This Plot' card above it")),
      tw_slide(NULL,
        "Next steps", "More details live in the Help page",
        list(p("Terminology (aggregation modes, data formats), every parameter, the full supported-drug list and answers to common issues are documented in the sections below on this page.")),
        list("Supported Drugs: the 44 pre-trained models at a glance",
             "FAQ: downloads, gene-name matching and background workers"))
    )

    observeEvent(input$tour_demo, {
      tw_walkthrough(ns("tour_demo_root"), "Demo walkthrough", demo_pages)
    })
    observeEvent(input$tour_advanced, {
      tw_walkthrough(ns("tour_adv_root"), "Your own data - walkthrough", adv_pages)
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
