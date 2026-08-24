# ---------------------------------------------------------------------------
# Background (async) training via callr
#
# The teacher's suggestion: heavy compute (train_models with CV + single-cell
# refinement) must NOT run on the Shiny main thread — on a multi-user server a
# long synchronous call would block every session. Instead each training job
# runs in a SEPARATE R process (callr::r_bg) that keeps DepMap cached, writes
# progress to a file, and hands back the finished model list as an RDS.
#
# Job protocol (all under tempdir()/perception_jobs/<jobid>/):
#   params.rds   -> training parameters (written by the app)
#   progress.txt -> latest "phase|i|n|drug" line (written by the worker)
#   status.txt   -> "running" / "done" / "error|<message>"
#   result.rds   -> finished model list (env-stripped, written by the worker)
# ---------------------------------------------------------------------------

# Job id: timestamp + random salt (teacher: user-keyed + timestamp uuid-ish).
make_job_id <- function(user_key = "") {
  salt <- paste(sample(c(letters, LETTERS, 0:9), 8, replace = TRUE), collapse = "")
  paste0("job_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
         substr(paste0(user_key, salt), 1, 12))
}

# ---- Worker main: runs inside the background R process -------------------
# Loads the live repo code + DepMap once, then polls the jobs dir forever.
train_worker_main <- function(pkg_root, depmap_path, jobs_dir, poll_secs = 1L) {
  suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  message("[worker] loading DepMap from ", depmap_path)
  DepMap <- readRDS(depmap_path)
  # Same dual assignment as the Data-tab loader (mod_data.R): the package
  # cache env AND the global env, so every train_models() code path sees it.
  depmap_env <- tryCatch(getFromNamespace(".depmap_env", "PERCEPTIONx"),
                         error = function(e) .GlobalEnv)
  assign("DepMap", DepMap, envir = depmap_env)
  assign("DepMap", DepMap, envir = .GlobalEnv)
  message("[worker] DepMap ready (", length(DepMap), " datasets), polling ", jobs_dir)

  while (TRUE) {
    jobs <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
    for (job in jobs) {
      status_file <- file.path(job, "status.txt")
      if (file.exists(status_file)) next          # already done/running
      progress_file <- file.path(job, "progress.txt")
      params <- readRDS(file.path(job, "params.rds"))
      writeLines("running", status_file)

      out <- tryCatch(
        PERCEPTIONx::train_models(
          drug_list        = params$drug_list,
          cancer_type      = params$cancer_type,
          exclude_cancer   = params$exclude_cancer,
          GOI              = params$GOI,
          k_features_values = params$k_features_values,
          model_type       = params$model_type,
          ncores           = params$ncores,
          output_dir       = params$output_dir,
          progress_cb = function(phase, i, n, drug) {
            writeLines(sprintf("%s|%d|%d|%s", phase, i, n, drug), progress_file)
          }
        ),
        error = function(e) structure(list(message = conditionMessage(e)),
                                      class = "perception_job_error")
      )

      if (inherits(out, "perception_job_error")) {
        writeLines(sprintf("error|%s", out$message), status_file)
      } else {
        # Strip formula/call environments before serialization (same trick as
        # train.R's strip_model_env) so saveRDS never walks into the session.
        strip_env <- function(m) {
          if (is.list(m) && inherits(m$model, "train")) {
            m$model$call <- NULL
            if (!is.null(m$model$terms)) attr(m$model$terms, ".Environment") <- baseenv()
          }
          m
        }
        saveRDS(lapply(out, strip_env), file.path(job, "result.rds"))
        writeLines("done", status_file)
      }
    }
    Sys.sleep(poll_secs)
  }
}

# ---- App-side job manager -------------------------------------------------

# Spawn the worker once per session; it loads DepMap once and stays resident.
ensure_train_worker <- function(shared) {
  if (isTRUE(shared$worker_alive) && !is.null(shared$train_worker) &&
      shared$train_worker$is_alive()) {
    return(invisible(TRUE))
  }
  if (is.null(shared$depmap_path)) {
    stop("DepMap file path is unknown. Please load DepMap data first (Data tab).")
  }
  jobs_dir <- file.path(tempdir(), "perception_jobs")
  dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)
  shared$jobs_dir <- jobs_dir
  worker <- callr::r_bg(
    func = train_worker_main,
    args = list(
      pkg_root    = getOption("perception.pkg_root", "."),
      depmap_path = shared$depmap_path,
      jobs_dir    = jobs_dir
    ),
    supervise = TRUE,
    poll_connection = FALSE      # avoid blocking the worker on a full stdout pipe
  )
  shared$train_worker <- worker
  shared$worker_alive <- TRUE
  invisible(TRUE)
}

# Submit a training job and return its job id. The worker picks it up as soon
# as it is ready (first job waits for DepMap to load in the worker).
submit_train_job <- function(shared, params) {
  ensure_train_worker(shared)
  jobid <- make_job_id()
  job_dir <- file.path(shared$jobs_dir, jobid)
  dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(params, file.path(job_dir, "params.rds"))
  file.create(file.path(job_dir, "progress.txt"))
  shared$active_job <- jobid
  invisible(jobid)
}

# Read the latest state of a job (queued / running / done / error).
read_job_state <- function(shared, jobid) {
  job_dir <- file.path(shared$jobs_dir, jobid)
  status_file <- file.path(job_dir, "status.txt")
  progress_file <- file.path(job_dir, "progress.txt")

  if (!file.exists(status_file)) {
    return(list(status = "queued", phase = "rank", i = 0L, n = 1L, drug = ""))
  }
  status <- readLines(status_file, warn = FALSE)
  if (length(status) == 0) status <- "running"
  s <- status[1]
  if (s == "done") {
    return(list(status = "done"))
  }
  if (startsWith(s, "error")) {
    return(list(status = "error", message = sub("^error\\|", "", s)))
  }
  # running: latest progress line
  prog <- "rank|0|1|"
  if (file.exists(progress_file)) {
    lines <- readLines(progress_file, warn = FALSE)
    if (length(lines) > 0) prog <- lines[length(lines)]
  }
  parts <- strsplit(prog, "\\|")[[1]]
  list(status = "running",
       phase = if (length(parts) > 0) parts[1] else "rank",
       i     = suppressWarnings(as.integer(parts[2])),
       n     = suppressWarnings(as.integer(parts[3])),
       drug  = if (length(parts) > 3) parts[4] else "")
}
