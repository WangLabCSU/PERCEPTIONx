# Data Loading Module

# ---- Upload helpers: flexible table reading + column name normalization ----

# Read an uploaded table by file extension. Supported: rds, csv, tsv/txt, xlsx.
# Leading "comment" rows are skipped automatically: any row at the top with
# fewer than 2 non-empty fields (e.g. a title line above the real header) is
# treated as a comment, and the first row with >= 2 fields becomes the header.
# The number of skipped rows is attached as attribute "skipped_rows".
read_uploaded_table <- function(file) {
  ext <- tolower(tools::file_ext(file$name))

  # Count fields in a delimited line, ignoring delimiters inside double quotes.
  count_fields <- function(line, sep) {
    if (!nzchar(trimws(line))) return(0L)
    pat <- paste0(sep, "(?=([^\"]*\"[^\"]*\")*[^\"]*$)")
    length(strsplit(line, pat, perl = TRUE)[[1]])
  }

  # First row with >= 2 non-empty cells = the real header row.
  first_header_row <- function(rows) {
    for (i in seq_along(rows)) {
      if (sum(!is.na(rows[[i]]) & nzchar(trimws(as.character(rows[[i]])))) >= 2) return(i)
    }
    1L
  }

  read_text <- function(path, sep) {
    lines <- readLines(path, warn = FALSE)
    nf <- vapply(lines, count_fields, integer(1), sep = sep)
    header_idx <- which(nf >= 2)[1]
    if (is.na(header_idx)) {
      header_idx <- 1L
      n_skip <- 0L
    } else {
      n_skip <- header_idx - 1L
    }
    con <- textConnection(paste(lines[seq.int(header_idx, length(lines))], collapse = "\n"))
    on.exit(close(con))
    out <- read.table(con, header = TRUE,
                      sep = if (sep == "[ \t]+") "" else sep,
                      stringsAsFactors = FALSE, check.names = FALSE,
                      comment.char = "", quote = "\"", fill = TRUE)
    attr(out, "skipped_rows") <- n_skip
    out
  }

  read_excel_tbl <- function(path) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Excel upload requires the 'readxl' package. Please run: install.packages('readxl')")
    }
    raw <- readxl::read_excel(path, col_names = FALSE, .name_repair = "minimal")
    rows <- lapply(seq_len(nrow(raw)), function(i)
      as.character(unlist(raw[i, ], use.names = FALSE)))
    header_idx <- first_header_row(rows)
    out <- as.data.frame(raw[seq.int(header_idx + 1L, nrow(raw)), , drop = FALSE],
                         stringsAsFactors = FALSE, check.names = FALSE)
    colnames(out) <- rows[[header_idx]]
    attr(out, "skipped_rows") <- header_idx - 1L
    out
  }

  switch(ext,
    rds = {
      obj <- readRDS(file$datapath)
      if (!isS4(obj)) attr(obj, "skipped_rows") <- 0L
      obj
    },
    csv = read_text(file$datapath, ","),
    tsv = read_text(file$datapath, "\t"),
    txt = read_text(file$datapath, "[ \t]+"),
    xlsx = read_excel_tbl(file$datapath),
    xls  = read_excel_tbl(file$datapath),
    stop("Unsupported file type: .", ext, ". Please use .csv, .tsv, .txt, .xlsx or .rds.")
  )
}

# Case-insensitive column matching to standard names (renames matched columns in place).
standardize_columns <- function(df, std_cols) {
  nms <- names(df)
  for (std in std_cols) {
    hit <- which(tolower(trimws(nms)) == tolower(std))
    if (length(hit) > 0) nms[hit[1]] <- std
  }
  names(df) <- nms
  df
}

# Map common response spellings to the canonical two-class labels.
normalize_response_labels <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[y %in% c("responder", "response", "responsive", "r", "sensitive", "sensitivity")] <- "Responder"
  y[y %in% c("non-responder", "nonresponder", "non responder", "non-responsive",
             "nonresponsive", "nr", "resistant", "resistance", "progressor",
             "progression", "non", "non_responder")] <- "Non-responder"
  # Longitudinal treatment time points (Maynard et al. Cell 2020; PERCEPTION
  # Fig. 4 lung cohort): kept as canonical 3-group labels rather than being
  # collapsed into R/NR.
  y[y %in% c("tn", "treatment naive", "treatment-naive", "naive", "untreated")] <- "TN"
  y[y %in% c("rd", "residual disease", "residual")] <- "RD"
  y[y %in% c("pd", "progressive disease", "progressed", "progression disease")] <- "PD"
  y
}

# Human-readable suffix for the number of auto-skipped top comment rows.
skipped_rows_msg <- function(n) {
  if (is.null(n) || n <= 0) return("")
  paste0(" (auto-skipped ", n, " top comment/title row", if (n > 1) "s", ")")
}

# Coerce various R object shapes into a long-format mapping data.frame.
# Accepts: data.frame (cell_id/patient_id), named list (patient -> cell ids),
#          named character vector (cell id -> patient id), or a 2-column matrix.
coerce_mapping_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    # Named list: names = patients, elements = cell-id vectors (empty entries dropped)
    if (!is.null(names(x))) {
      len <- lengths(x)
      if (sum(len) > 0) {
        return(data.frame(
          cell_id    = unlist(x, use.names = FALSE),
          patient_id = rep(names(x), len),
          stringsAsFactors = FALSE
        ))
      }
    }
    stop("Mapping list must be named, with cell-id vectors per patient.")
  }
  if (is.matrix(x)) {
    if (ncol(x) == 2) return(as.data.frame(x, stringsAsFactors = FALSE))
    stop("Mapping matrix must have 2 columns (cell_id, patient_id).")
  }
  if (is.atomic(x) && !is.null(names(x))) {
    # Named vector: names = cell ids, values = patient ids
    return(data.frame(cell_id = names(x), patient_id = as.character(x),
                      stringsAsFactors = FALSE))
  }
  stop("Unrecognized mapping format. Upload a table with 'cell_id' and 'patient_id' columns.")
}

# Coerce response input to a patient/response data.frame.
# Accepts a data.frame or a named vector (names = patients, values = labels).
coerce_response_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.atomic(x) && !is.null(names(x))) {
    return(data.frame(patient = names(x), response = as.character(x),
                      stringsAsFactors = FALSE))
  }
  stop("Unrecognized response format. Upload a table with 'patient' and 'response' columns.")
}

mod_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Upload Your Data — the primary entry point, placed first so users know
    # their own data is the main task (DepMap/models are only needed for
    # training and prediction, so they live in the later "Data Management").
    fluidRow(
      column(12,
        div(class = "section-header",
          icon("upload"),
          h4("Upload Your Data")
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
            "Click the button below to load demo data (49 genes × 400 cells, 20 patients) and explore all features without uploading anything."
          ),
          actionButton(ns("load_demo"), "Load Demo Data", class = "btn-demo btn-sm", icon = icon("play"))
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
              "Single-cell expression matrix (genes × cells). Raw counts or normalized values. CSV / TSV / TXT / Excel / RDS format."),
            p(style = "font-size: 0.8rem; font-weight: 700; color: var(--text); margin: 0.3rem 0 0.5rem;",
              "Header is auto-detected as the first row with ≥2 columns (gene / cell names); single-field title/comment rows above it are skipped. Avoid a multi-column title row above the header."),
            div(class = "data-format-hint",
              icon("table"), " Format: rows = genes, columns = cells",
              tags$pre(class = "data-format-example", 
              "CELL_001  CELL_002  CELL_003
TP53      2.1       0.0       5.3
BRCA1     0.0       1.8       3.2
EGFR      4.7       2.1       0.0")
            ),
            radioButtons(ns("expr_format"), "Data format",
                         choices = c(
                           "Cell-level — run Seurat clustering to detect clones" = "cell",
                           "Clone-level — each column is one clone, skip clustering" = "clone"
                         ),
                         selected = "cell", width = "100%"),
            fileInput(ns("expr_file"), "Upload Expression",
                      accept = c(".csv", ".tsv", ".txt", ".xlsx", ".rds", ".RDS"), width = "100%"),
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
              "File with columns: ", tags$code("cell_id"), " and ", tags$code("patient_id"),
              " (case-insensitive). Maps each cell (or clone) to its patient. In Clone-level mode the data is auto-prepared without clustering. CSV / TSV / TXT / Excel / RDS."),
            p(style = "font-size: 0.8rem; font-weight: 700; color: var(--text); margin: 0.3rem 0 0.5rem;",
              "Header is auto-detected as the first row with ≥2 columns (cell_id / patient_id); single-field title/comment rows above it are skipped."),
            # Cell-level example (default). Clone-level mode swaps to the
            # count-column example via the conditional panel below.
            div(class = "data-format-hint",
              icon("users"), " Format: cell_id + patient_id",
              tags$pre(class = "data-format-example", 
              "cell_id    patient_id
CELL_001   PAT_001
CELL_002   PAT_001
CELL_003   PAT_002")
            ),
            conditionalPanel(
              condition = paste0("input['", ns("expr_format"), "'] == 'clone'"),
              div(class = "data-format-hint", style = "margin-top: 0.4rem;",
                icon("users"), " Format: cell_id + patient_id + count",
                tags$pre(class = "data-format-example", 
                "cell_id       patient_id   count
CLONE_001     PAT_001      450
CLONE_002     PAT_001      50
CLONE_003     PAT_002      300"),
                tags$small(class = "text-muted",
                  "count = real cell number per clone; proportions are then shown accurately. Without it, clone proportions fall back to equal (1/n).")
              )
            ),
            fileInput(ns("mapping_file"), "Upload Mapping",
                      accept = c(".csv", ".tsv", ".txt", ".xlsx", ".rds", ".RDS"), width = "100%"),
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
              "Patient response data for evaluation. File with columns: patient, response (Responder/Non-responder, case-insensitive). CSV / TSV / TXT / Excel / RDS."),
            p(style = "font-size: 0.8rem; font-weight: 700; color: var(--text); margin: 0.3rem 0 0.5rem;",
              "Header is auto-detected as the first row with ≥2 columns (patient / response); single-field title/comment rows above it are skipped."),
            div(class = "data-format-hint",
              icon("heartbeat"), " Format: patient + response",
              tags$pre(class = "data-format-example", "patient    response
PAT_001    Responder
PAT_002    Non-responder
PAT_003    Responder")
            ),
            fileInput(ns("response_file"), "Upload Response",
                      accept = c(".csv", ".tsv", ".txt", ".xlsx", ".rds", ".RDS"), width = "100%"),
            uiOutput(ns("response_status"))
          )
        )
      )
    ),

    # Seurat Clustering & Preprocessing
    # NOTE: deliberately no animate-fade-in-up here — the entrance animation
    # (opacity/translate fill-mode) intermittently left this card's top half
    # invisible until a scroll forced a repaint.
    fluidRow(style = "margin-top: 1rem;",
      column(12,
        div(class = "card",
          div(class = "card-header",
            icon("shapes"), " Clone Detection & Preprocessing"
          ),
          div(class = "card-body seurat-body",
            # Left: description (always visible)
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
                "By default Seurat clusters cells into subclones. ",
                "Choose \"Clone-level\" in the Expression Matrix card to skip clustering — ",
                "the data is then prepared automatically after upload."
              )
            ),
            # Right: clone-level note (vertically centered) + clustering controls
            div(class = "seurat-right",
              div(id = ns("clone_note"), style = "display: none; flex: 1 1 auto; align-items: center; justify-content: center;",
                div(class = "info-box", style = "font-size: 0.8rem; padding: 0.5rem 0.7rem; text-align: left;",
                  icon("bolt"),
                  "Clone-level mode: data is prepared automatically after upload — no clustering needed."
                )
              ),
              div(id = ns("seurat_controls"),
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
                             title = "Load expression matrix and patient-cell mapping first")
              ),
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
Shiny.addCustomMessageHandler('expr-format-state-", ns("expr_format"), "', function(cloneMode) {
  var ctl = document.getElementById('", ns("seurat_controls"), "');
  var note = document.getElementById('", ns("clone_note"), "');
  if (ctl) ctl.style.display = cloneMode ? 'none' : '';
  if (note) note.style.display = cloneMode ? 'flex' : 'none';
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

    # --- Data Preview ---
    fluidRow(style = "margin-top: 1.5rem;",
      column(12,
        div(class = "section-header",
          icon("eye"),
          h4("Data Preview")
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
    ),

    # --- Data Management (DepMap / pre-trained models) ---
    fluidRow(style = "margin-top: 1.5rem;",
      column(12,
        div(class = "section-header",
          icon("database"),
          h4("Data Management")
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
              tags$a(href = "https://github.com/WangLabCSU/PERCEPTIONx/releases/tag/depmap",
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
                fileInput(ns("depmap_file"), NULL, accept = c(".RDS", ".rds"), width = "100%",
                          placeholder = "Select a .RDS file - loads automatically")
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
              "Load pre-trained drug response models from the PERCEPTIONx GitHub Release repository. 44 models are available, each trained on DepMap bulk expression with Elastic Net regression and 5-fold cross-validation. Models are cached locally after first download."),
            tags$small(class = "text-muted", style = "display: block; margin-top: 0.3rem;",
              "To download manually, visit ",
              tags$a(href = "https://github.com/WangLabCSU/PERCEPTIONx/releases/tag/models-v1",
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
                             options = list(placeholder = "Select one or more drugs...", maxOptions = 50,
                                            plugins = list("remove_button"))),
              actionButton(ns("load_model"), "Download & Load",
                           class = "btn-primary", icon = icon("download")),
              checkboxInput(ns("model_mirror"), "Use mirror", value = TRUE)
            ),
            div(style = "margin-top: 0.6rem; border-top: 1px dashed var(--border); padding-top: 0.6rem;",
              tags$small(class = "text-muted", "Or load a pre-downloaded model .RDS:"),
              div(style = "display: flex; gap: 0.5rem; align-items: center; margin-top: 0.3rem;",
                fileInput(ns("model_file"), NULL, accept = c(".RDS", ".rds"), width = "100%",
                          placeholder = "Select a .RDS file - loads automatically")
              )
            ),
            div(style = "margin-top: 0.75rem;",
              uiOutput(ns("model_status"))
            )
          )
        )
      )
    ),

    # --- Loaded Models Management ---
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
    )
  )
}

mod_data_server <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- Auto-prepare in clone-level mode (skip clustering) ---
    # Runs whenever expression + mapping are both available and the user has
    # selected "Clone-level" in the Expression Matrix card. No button needed.
    auto_prepare_if_clone <- function() {
      if (!identical(input$expr_format, "clone")) return(invisible(NULL))
      if (is.null(shared$user_expr) || is.null(shared$user_mapping)) return(invisible(NULL))
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Preparing clone-level data..."),
          p(class = "text-muted", "Mapping cells to patients and rank-normalizing")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      tryCatch({
        prepared <- PERCEPTIONx::prepare_data(
          expression_matrix = shared$user_expr,
          patient_mapping   = shared$user_mapping,
          skip_clustering   = TRUE
        )
        shared$prepared_data <- prepared
        shared$user_clones <- prepared$cell_clone_map
        w$hide()
        showNotification(paste0("Data ready (clone-level): ",
                                ncol(prepared$clone_expression_rnorm),
                                " clones across ", nrow(prepared$clone_counts), " patients"),
                         type = "message", duration = 5)
      }, error = function(e) {
        w$hide()
        showNotification(paste("Auto-prepare error:", e$message), type = "error", duration = 10)
      })
    }

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
      n_patients <- 20
      cell_names <- paste0("CELL_", sprintf("%04d", 1:n_cells))
      patient_names <- paste0("PAT_", sprintf("%03d", 1:n_patients))

      # Clinical response — defined FIRST so expression can be structured around it
      # 10 Responders + 10 Non-responders for meaningful box plots
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

            # y = signed mean of marker gene expression, interpreted as VIABILITY
            #   (high = resistant, low = sensitive):
            #   direction = -1: high marker expr → low viability (responder sensitive)
            #   direction = +1: high marker expr → high viability (responder resistant)
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

          # abemaciclib: markers HIGH in responders → direction -1 maps high expr
          #   to LOW viability (responders sensitive)
          # erlotinib: markers LOW in responders → direction +1 maps low expr
          #   to LOW viability (responders sensitive)
          shared$models <- list(
            abemaciclib = make_drug_model("abemaciclib", 101, abemaciclib_markers, direction = -1),
            erlotinib   = make_drug_model("erlotinib",   202, erlotinib_markers,   direction = +1)
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
        shared$depmap <- PERCEPTIONx::load_depmap(dest = tempdir(), read = TRUE, mirror = use_mirror,
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
    observeEvent(input$depmap_file, {
      req(input$depmap_file)
      file <- input$depmap_file

      # Pre-flight: deserializing an RDS needs roughly 2-4x its compressed
      # size as free RAM. Warn early instead of crashing mid-read (the
      # packaged 567 MB release is the known-good reference point).
      fsize_mb <- file.size(file$datapath) / 1024^2
      if (is.finite(fsize_mb) && fsize_mb > 700) {
        showNotification(
          sprintf("Large file (%.0f MB). Loading needs roughly 2-4x that in free RAM. Close other applications first, or use the built-in 567 MB 'Download & Load' above.",
                  fsize_mb),
          type = "warning", duration = 12)
      }

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
        gc()  # reclaim memory before the big allocation
        DepMap <- readRDS(file$datapath)

        # Validate the expected PERCEPTIONx structure so training does not
        # later fail with a cryptic error about a missing component.
        if (!is.list(DepMap)) {
          stop("The file is not a list. Expected a DepMap.RDS from the PERCEPTIONx release (or an object with the same structure).")
        }
        required_fields <- c("secondary_prism", "secondary_screen_drugAnnotation",
                             "expression_rnorm", "scRNA_complete", "scRNA_subset_rnorm",
                             "CPM_scRNA_CCLE_rnorm", "annotation_20Q4",
                             "metadata_CPM_scRNA", "expression_20Q4")
        missing_fields <- setdiff(required_fields, names(DepMap))
        if (length(missing_fields) > 0) {
          stop("The RDS is missing required components: ",
               paste(missing_fields, collapse = ", "),
               ". Re-save it with these components, or use the built-in 'Download & Load' (567 MB, known-good).")
        }

        depmap_env <- getFromNamespace(".depmap_env", "PERCEPTIONx")
        assign("DepMap", DepMap, envir = depmap_env)
        # Also write to .GlobalEnv so train_models() (which checks the global
        # env, and whose PSOCK workers export from it) can see DepMap.
        assign("DepMap", DepMap, envir = .GlobalEnv)
        shared$depmap <- DepMap
        w$hide()
        showNotification(paste("DepMap loaded from file:", length(DepMap), "datasets"), type = "message")
      }, error = function(e) {
        w$hide()
        msg <- conditionMessage(e)
        if (grepl("cannot allocate|failed to allocate", msg, ignore.case = TRUE)) {
          showNotification("Out of memory: the file needs more free RAM than this machine has right now.", type = "error", duration = 15)
          showNotification("Tips: close other applications; re-save the RDS containing only the components PERCEPTIONx needs; or use the built-in 567 MB 'Download & Load' instead.",
                           type = "warning", duration = 20)
        } else {
          showNotification(paste("Error reading file:", msg), type = "error", duration = 10)
        }
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
            new_model <- PERCEPTIONx::load_model(drug, dest = tempdir(), read = TRUE, mirror = use_mirror,
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
    observeEvent(input$model_file, {
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
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Reading expression matrix..."),
          p(class = "text-muted", "Parsing file and normalizing values")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      res <- tryCatch({
        mat <- read_uploaded_table(file)
        n_skip <- attr(mat, "skipped_rows"); if (is.null(n_skip)) n_skip <- 0L
        # RDS may be a matrix; text/Excel uploads come in as a data.frame
        if (!is.matrix(mat)) {
          mat <- as.matrix(mat)
          # First column holding gene names (e.g. from Excel/csv export) -> rownames
          if (ncol(mat) > 1 && is.character(mat[, 1]) &&
              anyNA(suppressWarnings(as.numeric(mat[, 1])))) {
            rownames(mat) <- mat[, 1]
            mat <- mat[, -1, drop = FALSE]
          }
        }
        # Ensure numeric (Seurat LogMap requires numeric, not integer)
        if (is.integer(mat)) mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat),
                                           dimnames = dimnames(mat))
        storage.mode(mat) <- "numeric"
        shared$user_expr <- mat
        shared$prepared_data <- NULL  # Reset prepared data when expression changes
        list(ok = TRUE,
             msg = paste0("Expression matrix loaded: ", nrow(mat), " genes x ",
                          ncol(mat), " cells", skipped_rows_msg(n_skip)))
      }, error = function(e) {
        list(ok = FALSE, msg = conditionMessage(e))
      })
      w$hide()
      if (isTRUE(res$ok)) {
        showNotification(res$msg, type = "message")
        auto_prepare_if_clone()
      } else {
        showNotification(paste("Error:", res$msg), type = "error")
      }
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
        df <- coerce_mapping_df(read_uploaded_table(file))
        n_skip <- attr(df, "skipped_rows"); if (is.null(n_skip)) n_skip <- 0L
        df <- standardize_columns(df, c("cell_id", "patient_id"))
        if (!all(c("cell_id", "patient_id") %in% names(df))) {
          stop("Mapping must contain columns 'cell_id' and 'patient_id' (case-insensitive). Found: ",
               paste(names(df), collapse = ", "))
        }
        df$cell_id    <- as.character(df$cell_id)
        df$patient_id <- as.character(df$patient_id)

        # Optional per-clone abundance column (count / cells / n_cells /
        # abundance): if present, expand each row to that many cells so clone
        # proportions reflect TRUE cell counts (matching the paper's figures)
        # instead of equal 1/n per clone. The abundance column is consumed here
        # and not passed downstream.
        count_col <- intersect(c("count", "cells", "n_cells", "n_cells_per_clone",
                                 "abundance"), tolower(names(df)))
        shared$mapping_has_count <- length(count_col) > 0
        if (length(count_col) > 0) {
          cnt_col <- names(df)[tolower(names(df)) == count_col[1]][1]
          cnt <- suppressWarnings(as.numeric(df[[cnt_col]]))
          cnt[is.na(cnt) | cnt < 1] <- 1L
          df <- df[rep(seq_len(nrow(df)), cnt), , drop = FALSE]
          df[[cnt_col]] <- NULL
          rownames(df) <- NULL
        }

        shared$user_mapping <- df
        shared$prepared_data <- NULL  # Reset prepared data when mapping changes
        showNotification(paste0("Patient-cell mapping loaded: ", nrow(shared$user_mapping),
                                " cells", skipped_rows_msg(n_skip)), type = "message")
        auto_prepare_if_clone()
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })

    output$mapping_status <- renderUI({
      if (is.null(shared$user_mapping)) return(NULL)
      n_patients <- length(unique(shared$user_mapping$patient_id))
      clone_mode <- identical(input$expr_format, "clone")
      # In clone-level mode without a count column there are no per-clone cell
      # counts, so clone proportions fall back to equal 1/n per clone. Be
      # explicit about this so users know what the plots actually show.
      equal_weight_note <- if (clone_mode && !isTRUE(shared$mapping_has_count)) {
        div(class = "info-box", style = "border-left-color: var(--warning, #d97706); margin-top: 0.4rem; font-size: 0.78rem; padding: 0.4rem 0.6rem;",
          icon("triangle-exclamation"),
          " No count column — clone proportions are shown as equal (1/n per clone). Add a <code>count</code> column with real cell numbers to show true proportions."
        )
      } else {
        NULL
      }
      tagList(
        span(class = "status-badge loaded", span(class = "status-dot green"), "Loaded"),
        tags$small(class = "text-muted",
          paste0(nrow(shared$user_mapping), " cells, ", n_patients, " patients")
        ),
        equal_weight_note
      )
    })

    # --- Expression format (cell-level vs clone-level) ---
    observeEvent(input$expr_format, {
      clone_mode <- identical(input$expr_format, "clone")
      session$sendCustomMessage(paste0("expr-format-state-", ns("expr_format")), clone_mode)
      # Mode switch invalidates any previous preparation; in clone-level mode
      # the data is re-prepared automatically.
      shared$prepared_data <- NULL
      shared$user_clones <- NULL
      if (clone_mode) auto_prepare_if_clone()
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
        prepared <- PERCEPTIONx::prepare_data(
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
        if (identical(input$expr_format, "clone")) {
          tags$span(class = "status-badge unloaded",
            span(class = "status-dot gray"),
            "Clone-level mode — will auto-prepare once expression + mapping are loaded"
          )
        } else {
          tags$span(class = "status-badge unloaded",
            span(class = "status-dot gray"),
            "Not run yet — upload expression + mapping, then click Run Seurat Clustering"
          )
        }
      } else {
        pd <- shared$prepared_data
        tagList(
          span(class = "status-badge loaded", span(class = "status-dot green"), "Complete"),
          tags$small(class = "text-muted",
            if (identical(pd$reduction_method, "none"))
              "(clone-level, clustering skipped) " else "",
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
        df <- coerce_response_df(read_uploaded_table(file))
        n_skip <- attr(df, "skipped_rows"); if (is.null(n_skip)) n_skip <- 0L
        df <- standardize_columns(df, c("patient", "response"))
        if (!all(c("patient", "response") %in% names(df))) {
          stop("Response file must contain columns 'patient' and 'response' (case-insensitive). Found: ",
               paste(names(df), collapse = ", "))
        }
        df$patient  <- as.character(df$patient)
        df$response <- normalize_response_labels(df$response)
        shared$user_response <- df
        # Warn if patient IDs do not overlap with the loaded mapping
        if (!is.null(shared$user_mapping)) {
          mapped_patients <- unique(as.character(shared$user_mapping$patient_id))
          overlap <- intersect(as.character(df$patient), mapped_patients)
          if (length(overlap) == 0) {
            showNotification("Warning: no patient IDs in the response file match the mapping's patient_id. Predictions cannot be validated.", type = "warning", duration = 10)
          }
        }
        showNotification(paste0("Clinical response loaded: ", nrow(shared$user_response),
                                " patients", skipped_rows_msg(n_skip)), type = "message")
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
    # --- DepMap preview: human-readable descriptions + click-to-inspect ---
    # Known members of the DepMap.RDS reference object, with short plain-English
    # descriptions so users know what each table actually contains.
    depmap_help <- c(
      expression_20q4              = "Bulk RNA expression (TPM) of CCLE cell lines (DepMap 20Q4). Rows = genes, columns = cell lines.",
      expression_rnorm             = "Rank-normalized bulk expression (genes x cell lines). The main input used to build drug response models.",
      scrna_complete               = "Single-cell RNA expression of CCLE cell lines (genes x cells). Used to tune models at single-cell resolution.",
      scrna_subset_rnorm           = "Rank-normalized subset of single-cell expression (genes x cells), aligned with the model feature space.",
      cpm_scrna_ccle_rnorm         = "Rank-normalized single-cell CPM expression (cells x genes). The full single-cell reference used for clone-level prediction.",
      metadata_cpm_scrna           = "Cell-level metadata for the single-cell reference (cell barcodes, cell line, quality metrics, etc.).",
      mutations_matrix             = "Binarized mutation matrix (genes x cell lines): 1 = non-silent mutation, 0 = wild type.",
      annotation_20q4              = "Cell line annotations (DepMap 20Q4): lineage, primary disease, growth conditions, etc.",
      drugcategory                 = "Drug classification table: drug name, target pathway, mechanism of action, etc.",
      secondary_prism              = "PRISM secondary-screen drug response values (viability/AUC) across cell lines.",
      secondary_screen_drugannotation = "Annotations for the PRISM secondary-screen drugs."
    )
    # Safe lookup: names not in depmap_help return NULL instead of an error.
    depmap_lookup <- function(n) {
      n <- tolower(n)
      if (n %in% names(depmap_help)) depmap_help[[n]] else NULL
    }

    output$depmap_preview <- renderDT({
      req(shared$depmap)
      nms <- names(shared$depmap)
      summary_df <- data.frame(
        Dataset = nms,
        Rows = sapply(shared$depmap, function(x) nrow(x)),
        Cols = sapply(shared$depmap, function(x) ncol(x)),
        Description = vapply(tolower(nms), function(n) {
          d <- depmap_lookup(n)
          if (is.null(d)) "Click View to inspect contents." else d
        }, character(1)),
        View = vapply(nms, function(n) {
          sprintf('<a href="#" onclick="Shiny.setInputValue(\'%s\', \'%s\', {priority: \'event\'}); return false;">View</a>',
                  ns("depmap_view"), n)
        }, character(1)),
        stringsAsFactors = FALSE
      )
      datatable(summary_df, options = list(pageLength = 10, dom = "tp"),
                rownames = FALSE, class = "display",
                escape = c(TRUE, TRUE, TRUE, TRUE, FALSE))
    })

    # Clicking "View" opens a modal with a description, column list, and a
    # small data preview so users know what the table actually contains.
    observeEvent(input$depmap_view, {
      req(shared$depmap, input$depmap_view)
      nm <- input$depmap_view
      obj <- shared$depmap[[nm]]
      if (is.null(obj)) return(NULL)
      desc <- depmap_lookup(tolower(nm))
      if (is.null(desc)) desc <- paste0("R object '", nm, "'.")
      col_head <- paste(head(colnames(obj), 12), collapse = ", ")
      if (ncol(obj) > 12) col_head <- paste0(col_head, ", ...")
      showModal(modalDialog(
        title = tags$strong(icon("database"), nm),
        size = "l",
        easyClose = TRUE,
        footer = uiOutput(ns("depmap_footer")),
        tags$div(
          p(style = "margin-bottom: 0.6rem;", desc),
          p(class = "text-muted", style = "font-size: 0.82rem; margin-bottom: 0.8rem;",
            sprintf("Dimensions: %s rows x %s columns.", nrow(obj), ncol(obj)),
            " | Columns: ", col_head),
          # tableOutput instead of DTOutput: the preview is a small fixed slice,
          # and a plain Shiny table has no client-side state to go stale across
          # successive clicks (the DataTables "column not found" bug).
          tableOutput(ns("depmap_detail_table"))
        )
      ))
    })

    # Modal footer: export button for small/medium tables; large reference
    # tables (millions of cells, e.g. the 53k-column single-cell matrix)
    # would produce enormous CSVs, so they are excluded with a short note.
    output$depmap_footer <- renderUI({
      req(shared$depmap, input$depmap_view)
      obj <- shared$depmap[[input$depmap_view]]
      if (is.null(obj)) return(modalButton("Close"))
      if (nrow(obj) * ncol(obj) > 3e6) {
        tagList(
          tags$span(class = "text-muted", style = "font-size: 0.8rem; margin-right: 0.5rem;",
            "Table too large for CSV export."),
          modalButton("Close")
        )
      } else {
        tagList(
          downloadButton(ns("depmap_export"), "Export CSV",
                         class = "btn-outline-primary btn-sm", icon = icon("download")),
          modalButton("Close")
        )
      }
    })

    output$depmap_export <- downloadHandler(
      filename = function() {
        paste0(input$depmap_view, "_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        req(shared$depmap, input$depmap_view)
        obj <- shared$depmap[[input$depmap_view]]
        write.csv(as.data.frame(obj), file, row.names = TRUE)
      }
    )

    output$depmap_detail_table <- renderTable({
      req(shared$depmap, input$depmap_view)
      obj <- shared$depmap[[input$depmap_view]]
      prev <- as.data.frame(obj)
      prev <- prev[seq_len(min(8, nrow(prev))),
                   seq_len(min(12, ncol(prev))), drop = FALSE]
      prev
    }, rownames = TRUE, digits = 4, width = "100%",
       caption = "Preview: first 8 rows (first 12 columns)")

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