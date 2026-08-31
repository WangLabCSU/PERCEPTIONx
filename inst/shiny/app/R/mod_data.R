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
          # Static button (NOT uiOutput/renderUI): server-rendered content
          # arrives only after a round trip, so the button would flash in a
          # second late on first load. The server swaps its HTML in place via
          # the 'set-html' handler when the demo loads / is cleared.
          div(id = ns("demo_btn"),
            actionButton(ns("load_demo"), "Load Demo Data",
                         class = "btn-demo btn-sm", icon = icon("play"),
                         # Instant client-side feedback: disable on click so
                         # the ~1s Shiny round trip before the waiter overlay
                         # appears does not look like a dead button.
                         onclick = "document.getElementById('demo-overlay').style.display='flex'; this.style.opacity='0.6'; var b=this; setTimeout(function(){b.style.opacity='';},2000);"))
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
              "Load pre-trained drug response models from the PERCEPTIONx GitHub Release repository. 44 models are available, each trained on DepMap bulk expression with Elastic Net regression and 3-fold cross-validation. Models are cached locally after first download."),
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

    # Demo state: whether the loaded data/models came from the demo button.
    # Drives the "Load Demo Data" <-> "Clear Demo Data" toggle. Stored on
    # shared so the Home page's "Load Demo" (mod_home.R) can mark it too.
    # NOTE: no explicit initialization here — reading an unset reactiveValues
    # field OUTSIDE a reactive consumer errors (app.R's onSessionEnded note),
    # and every read below happens inside a reactive context where
    # isTRUE(NULL) == FALSE is the correct default anyway.

    # Guards against double-submitting a background prepare task (rapid double
    # clicks on Run Seurat / auto-prepare would otherwise queue duplicate jobs
    # whose results race each other).
    prep_busy <- reactiveVal(FALSE)

    # Active DepMap download job: list(jobid, destfile, expected) or NULL when
    # idle. The progress observer below polls the growing file while set.
    download_progress <- reactiveVal(NULL)

    # Demo button state. The button is rendered STATICALLY in the UI (no
    # round-trip flash on first paint); this observer swaps it in place when
    # the demo loads (-> "Clear Demo Data") and back when it is cleared or
    # the user switches to their own data.
    observe({
      html <- if (isTRUE(shared$demo_loaded)) {
        as.character(actionButton(ns("load_demo"), "Clear Demo Data",
                                  class = "btn-danger btn-sm", icon = icon("trash"),
                                  onclick = "document.getElementById('demo-overlay').style.display='flex'; this.style.opacity='0.6'; var b=this; setTimeout(function(){b.style.opacity='';},2000);"))
      } else {
        as.character(actionButton(ns("load_demo"), "Load Demo Data",
                                  class = "btn-demo btn-sm", icon = icon("play"),
                                  onclick = "document.getElementById('demo-overlay').style.display='flex'; this.style.opacity='0.6'; var b=this; setTimeout(function(){b.style.opacity='';},2000);"))
      }
      session$sendCustomMessage("set-html", list(id = ns("demo_btn"), html = html))
    })

    # Auto-prepare in clone-level mode (skip clustering) ---
    # Runs whenever expression + mapping are both available and the user has
    # selected "Clone-level" in the Expression Matrix card. No button needed.
    auto_prepare_if_clone <- function() {
      if (!identical(input$expr_format, "clone")) return(invisible(NULL))
      if (is.null(shared$user_expr) || is.null(shared$user_mapping)) return(invisible(NULL))
      if (prep_busy()) return(invisible(NULL))   # a prepare task is already running
      prep_busy(TRUE)
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Preparing clone-level data..."),
          p(class = "text-muted", "Mapping cells to patients and rank-normalizing")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      jobid <- tryCatch(
        submit_session_task(shared, "prepare", list(
          expression_matrix = shared$user_expr,
          patient_mapping   = shared$user_mapping,
          skip_clustering   = TRUE
        )),
        error = function(e) {
          prep_busy(FALSE)
          w$hide()
          showNotification(paste("Failed to start the background task:", e$message),
                           type = "error", duration = 12)
          NULL
        }
      )
      if (is.null(jobid)) return(invisible(NULL))
      poll_task(shared, session, jobid,
        on_done = function(prepared) {
          prep_busy(FALSE)
          w$hide()
          shared$prepared_data <- prepared
          shared$user_clones <- prepared$cell_clone_map
          showNotification(paste0("Data ready (clone-level): ",
                                  ncol(prepared$clone_expression_rnorm),
                                  " clones across ", nrow(prepared$clone_counts), " patients"),
                           type = "message", duration = 5)
        },
        on_error = function(msg) {
          prep_busy(FALSE)
          w$hide()
          showNotification(paste("Auto-prepare error:", msg), type = "error", duration = 10)
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

    # --- Load / Clear Demo Data ---
    observeEvent(input$load_demo, {
      # Guard against double-submitting: the ~1s round trip before the waiter
      # overlay covers the page would otherwise let a second click queue a
      # duplicate demo task.
      if (isTRUE(shared$demo_busy)) {
        # Overlay was shown by onclick — hide it again since nothing new starts.
        session$sendCustomMessage("hide-demo-overlay", list())
        showNotification("Demo data is already being prepared...", type = "warning", duration = 5)
        return()
      }
      shared$demo_busy <- TRUE
      # Toggle: once the demo is loaded, the button becomes "Clear Demo Data".
      if (isTRUE(shared$demo_loaded)) {
        shared$user_response <- NULL
        shared$user_mapping  <- NULL
        shared$user_expr     <- NULL
        shared$prepared_data <- NULL
        shared$user_clones   <- NULL
        shared$predictions   <- NULL
        shared$patient_pred  <- NULL
        shared$models        <- NULL
        shared$model_cache   <- list()
        shared$model_active  <- list()
        shared$demo_loaded <- FALSE
        shared$demo_busy <- FALSE
        session$sendCustomMessage("hide-demo-overlay", list())
        showNotification("Demo data cleared. Upload your own data or load the demo again.",
                         type = "message", duration = 5)
        return()
      }
      # Shared demo pipeline (also used by the Home "Load Demo" button).
      run_demo_pipeline(shared, session,
        on_success = function() { shared$demo_loaded <- TRUE; shared$demo_busy <- FALSE },
        on_error = function() shared$demo_busy <- FALSE)
    })

    # --- Load DepMap (download) ---
    # Runs in the per-session background worker (submit_session_task) so the
    # main Shiny thread never blocks on the ~567 MB transfer. A progress
    # observer below polls the growing file and drives an overlay progress bar
    # styled like the training overlay (mod_train.R).
    observeEvent(input$load_depmap, {
      if (!is.null(download_progress())) {
        showNotification("A DepMap download is already in progress.", type = "warning", duration = 5)
        return()
      }
      use_mirror <- isTRUE(input$depmap_mirror)
      # Unified cache dir: options(PERCEPTIONx.depmap_cache_dir) / env var win;
      # otherwise derived from options(PERCEPTIONx.cache_root) (or its env var);
      # otherwise the R user data dir. Resolution lives in R/paths.R so the
      # package functions and the Shiny app share one rule.
      cache_dir <- PERCEPTIONx:::perception_depmap_cache_dir()
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      destfile <- file.path(cache_dir, "DepMap.RDS")

      # Cache TTL: drop the cached file if it has not been USED (each cache
      # hit refreshes a dedicated last-used flag file) for more than N hours,
      # so a stale ~567 MB file does not linger forever. Overridable via
      # options(PERCEPTIONx.depmap_cache_ttl_hours) or the
      # PERCEPTIONx_DEPMAP_CACHE_TTL_HOURS env var. Set it to 0 (or any value
      # <= 0) to DISABLE expiry entirely.
      #
      # Default: when a cache_root IS explicitly configured (deployment
      # servers, where the DepMap is pre-downloaded into a persistent dir),
      # expiry is DISABLED by default — a pre-cached file must survive until
      # it is deliberately replaced, no matter how long it sits unused.
      # Without cache_root (personal local use, data in tempdir) the default
      # stays 12 h idle cleanup so a 567 MB file does not linger forever.
      # If the file is locked by a running worker, unlink() fails harmlessly
      # and the next click retries.
      #
      # NB: the "last used" time lives in a SEPARATE flag file, NEVER in
      # DepMap.RDS's own mtime — extract_depmap_meta() only trusts the meta
      # sidecar when it is newer than the source RDS, so touching the source
      # would invalidate the meta cache and force a full 567 MB re-read on
      # every click.
      last_used <- file.path(cache_dir, "DepMap_used.flag")
      cache_ttl_h <- suppressWarnings(as.numeric(getOption(
        "PERCEPTIONx.depmap_cache_ttl_hours",
        Sys.getenv("PERCEPTIONx_DEPMAP_CACHE_TTL_HOURS", ""))))
      if (is.na(cache_ttl_h)) {
        cache_ttl_h <- if (is.null(PERCEPTIONx:::perception_cache_root())) 12 else 0
      }
      ref_file <- if (file.exists(last_used)) last_used else destfile
      # Only expire when a POSITIVE TTL is set; 0 / negative disables expiry
      # (pre-cached deployment files must survive).
      if (file.exists(destfile) && !is.na(cache_ttl_h) && cache_ttl_h > 0) {
        age_h <- as.numeric(difftime(Sys.time(), file.info(ref_file)$mtime, units = "hours"))
        if (is.finite(age_h) && age_h > cache_ttl_h) {
          unlink(destfile)
          unlink(file.path(cache_dir, "DepMap_meta.RDS"))
          unlink(last_used)
          showNotification(sprintf("Cached DepMap expired (unused for %.0f h > %.0f h) — re-downloading.",
                                   age_h, cache_ttl_h), type = "message", duration = 6)
        }
      }

      # Already cached: serve the metadata sidecar directly (kilobytes) when
      # it is fresh. If the sidecar is missing/stale, re-extract it in the
      # background worker — the main process must never deserialize the
      # 567 MB file (blocks the UI + spikes RAM).
      if (file.exists(destfile) && file.size(destfile) > 0) {
        meta_cache <- file.path(cache_dir, "DepMap_meta.RDS")
        meta <- read_cached_meta(destfile, meta_cache)
        if (!is.null(meta)) {
          # Cache hit: refresh the last-used timestamp in the SEPARATE flag
          # file — never touch DepMap.RDS's own mtime (see the TTL note above).
          suppressWarnings(try(file.create(last_used), silent = TRUE))
          shared$depmap_meta <- meta
          shared$depmap_path <- destfile
          # Cached standard download = trusted standard file -> shared worker pool.
          shared$depmap_is_standard <- TRUE
          notify_master_depmap(destfile)
          showNotification("DepMap data loaded from cache", type = "message")
          return()
        }
        # Meta cache unusable -> re-extract in the worker (never the main
        # process: reading 567 MB here would block the UI).
        ew <- Waiter$new(
          html = tagList(
            div(class = "spinner-ring"),
            h4("Reading DepMap data..."),
            p(class = "text-muted", "Extracting metadata (a few seconds)")
          ),
          color = "rgba(255,255,255,0.85)"
        )
        ew$show()
        jobid <- tryCatch(
          submit_session_task(shared, "extract_meta", list(
            depmap_path = destfile,
            cache_file  = meta_cache
          )),
          error = function(e) { ew$hide(); NULL }
        )
        if (is.null(jobid)) return()
        poll_task(shared, session, jobid,
          on_done = function(meta) {
            suppressWarnings(try(file.create(last_used), silent = TRUE))
            shared$depmap_meta <- meta
            shared$depmap_path <- destfile
            shared$depmap_is_standard <- TRUE
            notify_master_depmap(destfile)
            ew$hide()
            showNotification("DepMap data loaded", type = "message")
          },
          on_error = function(msg) {
            ew$hide()
            showNotification(paste("Error:", msg), type = "error", duration = 10)
            showNotification("Tip: Try enabling 'Use mirror' or download DepMap.RDS manually and use the 'Load' button below.",
                             type = "warning", duration = 15)
          })
        return()
      }

      # ONE overlay layer: spinner + stage text + progress bar live together
      # inside the same waiter overlay (same pattern as the training overlay).
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Downloading DepMap data..."),
          p(id = ns("dwn_stage"), class = "text-muted", "Starting download..."),
          p(id = ns("dwn_detail"), class = "text-muted",
            style = "font-size: 0.82rem; opacity: 0.75;",
            if (use_mirror) "567 MB via mirror (mirrors tried in order)"
            else "567 MB from GitHub"),
          # The id sits on the TRACK (see mod_train.R note): set-html replaces
          # the bar inside it with a fresh .train-progress-bar.
          div(id = ns("dwn_bar"), class = "train-progress-track",
            div(class = "train-progress-bar", style = "width: 2%;")
          )
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()

      base_urls <- if (use_mirror) PERCEPTIONx::get_mirrors() else "https://github.com"
      urls <- paste0(base_urls, "/WangLabCSU/PERCEPTIONx/releases/download/depmap/DepMap.RDS")

      jobid <- tryCatch(
        submit_session_task(shared, "download", list(
          urls            = urls,
          destfile        = destfile,
          timeout_seconds = 600,
          retries         = 2
        )),
        error = function(e) { w$hide(); NULL }
      )
      if (is.null(jobid)) return()

      # Expected size of the official DepMap.RDS release asset (bytes) — used
      # only as the progress-bar denominator. The real file is validated by
      # the parsing step below, so a wrong estimate can never corrupt data.
      download_progress(list(jobid = jobid, destfile = destfile, expected = 594589700))

      poll_task(shared, session, jobid,
        on_done = function(path) {
          download_progress(NULL)
          # Download hit 100%: close the download overlay, then immediately
          # show a "parsing" overlay. Reading the 567 MB file + extracting
          # metadata runs as a SECOND background task, so the UI never blocks
          # and the user always knows which stage we are in.
          w$hide()
          pw <- Waiter$new(
            html = tagList(
              div(class = "spinner-ring"),
              h4("Parsing DepMap data..."),
              p(class = "text-muted", "Reading the file and extracting metadata (a few seconds)")
            ),
            color = "rgba(255,255,255,0.85)"
          )
          pw$show()
          jobid2 <- tryCatch(
            submit_session_task(shared, "extract_meta", list(
              depmap_path = path,
              cache_file  = file.path(cache_dir, "DepMap_meta.RDS")
            )),
            error = function(e) { pw$hide(); e }
          )
          if (inherits(jobid2, "error")) {
            showNotification(paste("Download finished but parsing failed to start:", conditionMessage(jobid2)),
                             type = "error", duration = 10)
            return()
          }
          poll_task(shared, session, jobid2,
            on_done = function(meta) {
              suppressWarnings(try(file.create(last_used), silent = TRUE))
              shared$depmap_meta <- meta
              shared$depmap_path <- path
              # Built-in download = trusted standard file -> shared worker pool.
              shared$depmap_is_standard <- TRUE
              notify_master_depmap(path)
              pw$hide()
              showNotification("DepMap data loaded successfully", type = "message")
            },
            on_error = function(msg) {
              pw$hide()
              # A corrupt/truncated file from an interrupted download would
              # otherwise poison EVERY later session. Detect and delete it so
              # the next click re-downloads a clean copy. Memory errors must
              # NOT delete anything — the file may be fine.
              corrupt <- grepl(paste0("error reading from connection|unknown input format|",
                                     "cannot open the connection|not a list|",
                                     "missing required components"), msg, ignore.case = TRUE)
              if (corrupt) {
                unlink(path)
                unlink(file.path(cache_dir, "DepMap_meta.RDS"))
                showNotification("Downloaded DepMap.RDS is corrupt (incomplete download?) — it was deleted. Please click 'Download & Load' again.",
                                 type = "error", duration = 10)
              } else {
                showNotification(paste("Error:", msg), type = "error", duration = 10)
              }
              showNotification("Tip: Try enabling 'Use mirror' or download DepMap.RDS manually and use the 'Load' button below.",
                               type = "warning", duration = 15)
            })
        },
        on_error = function(msg) {
          download_progress(NULL)
          w$hide()
          showNotification(paste("Download failed:", msg), type = "error", duration = 10)
          showNotification("Tip: Try enabling 'Use mirror' or download DepMap.RDS manually and use the 'Load' button below.",
                           type = "warning", duration = 15)
        })
    })

    # Drive the download progress bar: poll the growing file once per second.
    # Stops as soon as the job finishes (download_progress -> NULL).
    observe({
      st <- download_progress()
      if (is.null(st)) return()
      tw <- shared$task_worker
      if (!is.null(tw) && !tw$is_alive()) { download_progress(NULL); return() }
      invalidateLater(1000)
      size <- if (file.exists(st$destfile)) file.size(st$destfile) else 0
      expected <- max(1, st$expected)
      pct <- round(100 * pmin(1, size / expected))
      stage <- sprintf("Downloaded %.0f MB of %.0f MB (%d%%)",
                       size / 1024^2, expected / 1024^2, pct)
      session$sendCustomMessage("set-html", list(id = ns("dwn_stage"), html = stage))
      session$sendCustomMessage("set-html",
        list(id = ns("dwn_bar"),
             html = paste0("<div class='train-progress-bar' style='width: ", pct, "%;'></div>")))
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
      jobid <- tryCatch(
        submit_session_task(shared, "extract_meta", list(
          depmap_path = file$datapath,
          cache_file  = NULL
        )),
        error = function(e) {
          w$hide()
          showNotification(paste("Failed to start the background task:", e$message),
                           type = "error", duration = 12)
          NULL
        }
      )
      if (is.null(jobid)) return()
      poll_task(shared, session, jobid,
        on_done = function(meta) {
          # Record the source file so the background training worker (callr)
          # can load the identical DepMap from disk in its own process.
          shared$depmap_meta <- meta
          shared$depmap_path <- file$datapath
          # User-uploaded file = untrusted -> isolated per-session worker,
          # never the shared pool (a wrong upload must not affect others).
          shared$depmap_is_standard <- FALSE
          w$hide()
          showNotification(paste("DepMap loaded from file:",
                                 length(meta$components), "datasets"),
                           type = "message")
        },
        on_error = function(msg) {
          w$hide()
          if (grepl("cannot allocate|failed to allocate|Error in readRDS",
                    msg, ignore.case = TRUE)) {
            showNotification("Failed to read the file (possible out of memory in the worker).",
                             type = "error", duration = 15)
            showNotification("Tips: close other applications; re-save the RDS containing only the components PERCEPTIONx needs; or use the built-in 567 MB 'Download & Load' instead.",
                             type = "warning", duration = 20)
          } else {
            showNotification(paste("Error reading file:", msg), type = "error", duration = 10)
          }
        })
    })

    output$depmap_status <- renderUI({
      if (is.null(shared$depmap_meta)) {
        tagList(span(class = "status-badge unloaded", span(class = "status-dot gray"), "Not loaded"))
      } else {
        n_datasets <- length(shared$depmap_meta$components)
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
        # Pre-trained models live in a persistent cache dir when
        # options(PERCEPTIONx.cache_root) is set; otherwise tempdir() keeps the
        # historical behaviour. load_model() skips drugs already on disk, so a
        # cached model is not downloaded again on the next app session.
        model_dir <- PERCEPTIONx:::perception_model_dir()
        if (is.null(model_dir)) model_dir <- tempdir()
        loaded_count <- 0
        failed_drugs <- c()
        for (drug in drug_list) {
          tryCatch({
            new_model <- PERCEPTIONx::load_model(drug, dest = model_dir, read = TRUE, mirror = use_mirror,
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
            # NOTE: the status dot must go in `label`, not `icon` — Shiny >= 1.10
            # validates the icon arg (validateIcon) and errors "Invalid icon"
            # for any non-icon() tag.
            actionButton(ns(btn_id), tagList(
              tags$span(class = "status-dot", class = dot_class,
                        style = "width:8px;height:8px;border-radius:50%;display:inline-block;"),
              label
            ),
            class = paste("btn-model-toggle btn-sm", state_class))
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
      shared$demo_loaded <- FALSE  # user switched to their own data
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
      shared$demo_loaded <- FALSE  # user switched to their own data
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

    # --- Run Seurat Clustering (prepare_data) — background worker ---
    observeEvent(input$run_seurat, {
      req(shared$user_expr, shared$user_mapping)
      if (prep_busy()) {
        showNotification("A preparation task is already running", type = "warning", duration = 5)
        return()
      }
      prep_busy(TRUE)
      w <- Waiter$new(
        html = tagList(
          div(class = "spinner-ring"),
          h4("Running Seurat clustering..."),
          p(class = "text-muted", "Detecting subclones and preparing data")
        ),
        color = "rgba(255,255,255,0.85)"
      )
      w$show()
      jobid <- tryCatch(
        submit_session_task(shared, "prepare", list(
          method            = input$seurat_method,
          expression_matrix = shared$user_expr,
          patient_mapping   = shared$user_mapping,
          seurat_resolution = input$seurat_resolution,
          seurat_dims       = input$seurat_dims
        )),
        error = function(e) {
          prep_busy(FALSE)
          w$hide()
          # Make a failed worker spawn VISIBLE — otherwise the overlay just
          # vanishes after ~0.5s with no error and the status never changes.
          showNotification(paste("Failed to start the background task:", e$message),
                           type = "error", duration = 12)
          NULL
        }
      )
      if (is.null(jobid)) return()
      poll_task(shared, session, jobid,
        on_done = function(res) {
          prep_busy(FALSE)
          w$hide()
          shared$prepared_data <- res
          shared$user_clones <- res$cell_clone_map
          # Keep shared$user_expr as cell-level expression for biomarker plots
          showNotification(paste0("Clustering complete: ", ncol(res$clone_expression_rnorm),
                                  " clones detected across ", nrow(res$clone_counts), " patients"),
                           type = "message", duration = 5)
        },
        on_error = function(msg) {
          prep_busy(FALSE)
          w$hide()
          showNotification(paste("Clustering error:", msg), type = "error", duration = 10)
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
      shared$demo_loaded <- FALSE  # user switched to their own data
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
        # NO trailing comma on the last element: renderUI captures the
        # expression with rlang, which turns a trailing comma into an empty
        # 7th argument -> "Error in list: argument 7 is empty".
        list(name = "DepMap Data", loaded = !is.null(shared$depmap_meta))
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
      req(shared$depmap_meta)
      comps <- shared$depmap_meta$components
      nms <- names(comps)
      summary_df <- data.frame(
        Dataset = nms,
        Rows = vapply(comps, function(x) x$nrow, numeric(1)),
        Cols = vapply(comps, function(x) x$ncol, numeric(1)),
        Description = vapply(tolower(nms), function(n) {
          d <- depmap_lookup(n)
          if (is.null(d)) "Click View to inspect details." else d
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

    # Clicking "View" shows a modal with the description, dimensions, and the
    # first column names. The main process keeps ONLY metadata (the full
    # multi-GB object lives in the training worker), so no data preview is
    # rendered here.
    observeEvent(input$depmap_view, {
      req(shared$depmap_meta, input$depmap_view)
      nm <- input$depmap_view
      comp <- shared$depmap_meta$components[[nm]]
      if (is.null(comp)) return(NULL)
      desc <- depmap_lookup(tolower(nm))
      if (is.null(desc)) desc <- paste0("R object '", nm, "'.")
      col_head <- paste(comp$cols_preview, collapse = ", ")
      if (!is.na(comp$ncol) && comp$ncol > 12) col_head <- paste0(col_head, ", ...")
      showModal(modalDialog(
        title = tags$strong(icon("database"), nm),
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$div(
          p(style = "margin-bottom: 0.6rem;", desc),
          p(class = "text-muted", style = "font-size: 0.82rem; margin-bottom: 0.8rem;",
            sprintf("Dimensions: %s rows x %s columns.", comp$nrow, comp$ncol),
            " | Columns: ",
            if (length(col_head) == 0 || !nzchar(col_head)) "—" else col_head),
          p(class = "text-muted", style = "font-size: 0.78rem;",
            "No data preview — the main process keeps only DepMap metadata; the full object lives in the background training worker.")
        )
      ))
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