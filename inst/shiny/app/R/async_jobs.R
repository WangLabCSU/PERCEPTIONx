# ---------------------------------------------------------------------------
# Background (async) training — worker pool
#
# Design (agreed with the user):
#   1. LAZY master: one background R process ("master") is spawned only when
#      the first user loads the STANDARD DepMap. It reads DepMap once (the
#      only full copy in the server) and keeps it.
#   2. FORK PARALLELISM (Linux only): pending training jobs are dispatched
#      with parallel::mclapply — forked children inherit DepMap via
#      copy-on-write, so N concurrent jobs share ONE physical copy.
#      On Windows (no fork) jobs run sequentially inside the master — correct,
#      just not parallel, so self-hosted Windows users keep full support.
#   3. CONTAMINATION ISOLATION: users who UPLOAD their own DepMap file never
#      touch the shared master. They get a dedicated per-session worker that
#      loads THEIR file, so a wrong upload can only hurt that one user.
#   4. IDLE RELEASE: after <idle_minutes> without any training activity the
#      master drops the DepMap (rm + gc) so memory returns to the server.
#      The next job reloads it from the last known path (cold start, ~30-60s,
#      shown as "waiting for worker" in the progress overlay).
#
# IMPORTANT (callr serialization): callr::r_bg() ships the worker function to
# a fresh R process. Only the function itself travels — TOP-LEVEL helpers do
# NOT. Therefore every worker entry function (train_master_main /
# train_custom_main) is fully SELF-CONTAINED: all helpers are defined INSIDE
# its body, so they exist at runtime in the child.
#
# Job protocol (jobs_dir/<jobid>/):
#   params.rds   <- training parameters (written by the app)
#   progress.txt <- latest "phase|i|n|drug" line (written by the worker)
#   status.txt   <- "running" / "done" / "error|<message>"
#   result.rds   <- finished model list (env-stripped, written by the worker)
# ---------------------------------------------------------------------------

# Job id: timestamp + random chars so ids never collide.
make_job_id <- function(user_key = "") {
  salt <- paste(sample(c(letters, LETTERS, 0:9), 8, replace = TRUE), collapse = "")
  paste0("job_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
         substr(paste0(user_key, salt), 1, 12))
}

# ---- Master worker: runs inside the shared background R process ------------
# Fully self-contained (see note above).

train_master_main <- function(pkg_root, jobs_dir, max_parallel, idle_minutes) {
  # Deployment fallback: during development the repo lives locally (load_all);
  # on a server the package is installed (library) and pkg_root is NULL/absent.
  if (!is.null(pkg_root) && nzchar(pkg_root) && dir.exists(pkg_root)) {
    suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  } else {
    library(PERCEPTIONx)
  }

  # --- helpers (local to this process) ---
  write_progress <- function(progress_file, phase, i, n, drug) {
    writeLines(sprintf("%s|%d|%d|%s", phase, i, n, drug), progress_file)
  }
  strip_model_env <- function(m) {
    if (is.list(m) && inherits(m$model, "train")) {
      m$model$call <- NULL
      if (!is.null(m$model$terms)) attr(m$model$terms, ".Environment") <- baseenv()
    }
    m
  }
  assign_depmap <- function(DepMap) {
    depmap_env <- tryCatch(getFromNamespace(".depmap_env", "PERCEPTIONx"),
                           error = function(e) .GlobalEnv)
    assign("DepMap", DepMap, envir = depmap_env)
    assign("DepMap", DepMap, envir = .GlobalEnv)
    invisible(TRUE)
  }
  run_train_job <- function(job_dir) {
    status_file <- file.path(job_dir, "status.txt")
    progress_file <- file.path(job_dir, "progress.txt")
    params <- readRDS(file.path(job_dir, "params.rds"))
    out <- tryCatch(
      PERCEPTIONx::train_models(
        drug_list         = params$drug_list,
        cancer_type       = params$cancer_type,
        exclude_cancer    = params$exclude_cancer,
        GOI               = params$GOI,
        k_features_values = params$k_features_values,
        model_type        = params$model_type,
        ncores            = params$ncores,
        output_dir        = params$output_dir,
        progress_cb = function(phase, i, n, drug) {
          write_progress(progress_file, phase, i, n, drug)
        }
      ),
      error = function(e) structure(list(message = conditionMessage(e)),
                                    class = "perception_job_error")
    )
    if (inherits(out, "perception_job_error")) {
      writeLines(sprintf("error|%s", out$message), status_file)
    } else {
      saveRDS(lapply(out, strip_model_env), file.path(job_dir, "result.rds"))
      writeLines("done", status_file)
    }
    invisible(TRUE)
  }

  # --- main loop ---
  DepMap <- NULL
  depmap_path <- NULL
  request_file <- file.path(jobs_dir, "depmap_request.rds")
  last_activity <- Sys.time()
  message("[master] started, polling ", jobs_dir, " (idle release after ",
          idle_minutes, " min)")

  while (TRUE) {
    # 1. DepMap: honour a new request; otherwise reload the last known path
    #    after an idle release.
    if (file.exists(request_file)) {
      p <- tryCatch(readRDS(request_file), error = function(e) NULL)
      file.remove(request_file)
      if (!is.null(p) && !identical(p, depmap_path)) depmap_path <- p
    }
    if (!is.null(depmap_path) && is.null(DepMap)) {
      message("[master] loading DepMap from ", depmap_path)
      DepMap <- readRDS(depmap_path)
      assign_depmap(DepMap)
      message("[master] DepMap ready (", length(DepMap), " datasets)")
      last_activity <- Sys.time()
    }

    # 2. Claim pending jobs (file.create wins; only this master claims, so no
    #    double-processing).
    pending <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
    pending <- pending[!vapply(pending, function(j)
      file.exists(file.path(j, "status.txt")), logical(1))]
    claimed <- character(0)
    for (j in pending) {
      if (file.create(file.path(j, "status.txt"))) claimed <- c(claimed, j)
    }
    if (length(claimed) > 0) {
      message("[master] dispatching ", length(claimed), " job(s)")
      if (.Platform$OS.type == "unix" && max_parallel > 1L &&
          length(claimed) > 1L) {
        # Linux: fork — children inherit DepMap copy-on-write, one physical copy
        parallel::mclapply(claimed, run_train_job,
                           mc.cores = min(max_parallel, length(claimed)))
      } else {
        # Windows (no fork) or a single job: run in-process, serial
        lapply(claimed, run_train_job)
      }
      last_activity <- Sys.time()
    }

    # 3. Idle release: R's gc() never returns heap pages to the OS, so merely
    #    dropping DepMap would NOT free the ~8GB. The only way to give memory
    #    back to the server is to EXIT the process. The next submitted job then
    #    re-spawns the master (ensure_master detects !is_alive) and cold-loads
    #    the DepMap from the last known path — shown as "waiting for worker".
    if (!is.null(DepMap) &&
        as.numeric(difftime(Sys.time(), last_activity, units = "mins")) > idle_minutes) {
      message("[master] idle for ", idle_minutes, " min - exiting to release memory")
      quit(save = "no")
    }

    Sys.sleep(2)
  }
}

# ---- Dedicated worker for users who UPLOADED their own DepMap --------------
# A wrong upload must never reach the shared master; this per-session worker
# loads the user's own file and dies with the session (supervise).

train_custom_main <- function(pkg_root, depmap_path, jobs_dir, poll_secs = 1L) {
  if (!is.null(pkg_root) && nzchar(pkg_root) && dir.exists(pkg_root)) {
    suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  } else {
    library(PERCEPTIONx)
  }

  write_progress <- function(progress_file, phase, i, n, drug) {
    writeLines(sprintf("%s|%d|%d|%s", phase, i, n, drug), progress_file)
  }
  strip_model_env <- function(m) {
    if (is.list(m) && inherits(m$model, "train")) {
      m$model$call <- NULL
      if (!is.null(m$model$terms)) attr(m$model$terms, ".Environment") <- baseenv()
    }
    m
  }
  assign_depmap <- function(DepMap) {
    depmap_env <- tryCatch(getFromNamespace(".depmap_env", "PERCEPTIONx"),
                           error = function(e) .GlobalEnv)
    assign("DepMap", DepMap, envir = depmap_env)
    assign("DepMap", DepMap, envir = .GlobalEnv)
    invisible(TRUE)
  }
  run_train_job <- function(job_dir) {
    status_file <- file.path(job_dir, "status.txt")
    progress_file <- file.path(job_dir, "progress.txt")
    params <- readRDS(file.path(job_dir, "params.rds"))
    out <- tryCatch(
      PERCEPTIONx::train_models(
        drug_list         = params$drug_list,
        cancer_type       = params$cancer_type,
        exclude_cancer    = params$exclude_cancer,
        GOI               = params$GOI,
        k_features_values = params$k_features_values,
        model_type        = params$model_type,
        ncores            = params$ncores,
        output_dir        = params$output_dir,
        progress_cb = function(phase, i, n, drug) {
          write_progress(progress_file, phase, i, n, drug)
        }
      ),
      error = function(e) structure(list(message = conditionMessage(e)),
                                    class = "perception_job_error")
    )
    if (inherits(out, "perception_job_error")) {
      writeLines(sprintf("error|%s", out$message), status_file)
    } else {
      saveRDS(lapply(out, strip_model_env), file.path(job_dir, "result.rds"))
      writeLines("done", status_file)
    }
    invisible(TRUE)
  }

  message("[custom-worker] loading DepMap from ", depmap_path)
  DepMap <- readRDS(depmap_path)
  assign_depmap(DepMap)
  message("[custom-worker] DepMap ready, polling ", jobs_dir)
  while (TRUE) {
    pending <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
    for (j in pending) {
      if (!file.exists(file.path(j, "status.txt")) &&
          file.create(file.path(j, "status.txt"))) {
        run_train_job(j)
      }
    }
    Sys.sleep(poll_secs)
  }
}

# ---- App-side job manager (runs in the main Shiny process) -----------------

# Server-wide master: spawn once (lazily), reused by every session.
ensure_master <- function() {
  m <- getOption("perception.master")
  if (!is.null(m) && inherits(m, "r_process") && m$is_alive()) {
    return(invisible(TRUE))
  }
  jobs_dir <- file.path(tempdir(), "perception_jobs")
  dir.create(jobs_dir, recursive = TRUE, showWarnings = FALSE)
  options(perception.jobs_dir = jobs_dir)
  max_par <- suppressWarnings(as.integer(Sys.getenv("PERCEPTION_WORKERS", "16")))
  if (is.na(max_par) || max_par < 1L) max_par <- 16L
  idle <- suppressWarnings(as.integer(Sys.getenv("PERCEPTION_WORKER_IDLE_MINUTES", "720")))
  if (is.na(idle) || idle < 1L) idle <- 720L
  worker <- callr::r_bg(
    func = train_master_main,
    # NOTE: default must be NULL (not ".") — a NULL pkg_root makes the worker
    # fall back to library(PERCEPTIONx); "." would trigger devtools::load_all
    # on an unrelated working directory and crash the worker.
    args = list(pkg_root = getOption("perception.pkg_root", NULL),
                jobs_dir = jobs_dir,
                max_parallel = max_par,
                idle_minutes = idle),
    supervise = TRUE,
    poll_connection = FALSE
  )
  options(perception.master = worker)
  invisible(TRUE)
}

# Tell the master which standard DepMap file to load (called when a user loads
# the standard DepMap via the built-in button; idempotent).
notify_master_depmap <- function(depmap_path) {
  ensure_master()
  saveRDS(depmap_path, file.path(getOption("perception.jobs_dir"),
                                 "depmap_request.rds"))
  invisible(TRUE)
}

# Per-session worker for a user who UPLOADED a custom DepMap file.
ensure_custom_worker <- function(shared) {
  if (isTRUE(shared$worker_alive) && !is.null(shared$train_worker) &&
      shared$train_worker$is_alive()) {
    return(invisible(TRUE))
  }
  if (is.null(shared$depmap_path)) {
    stop("DepMap file path is unknown. Please load DepMap data first (Data tab).")
  }
  jobs_dir <- file.path(tempdir(), "perception_jobs_custom", make_job_id())
  dir.create(jobs_dir, recursive = TRUE, showWarnings = FALSE)
  worker <- callr::r_bg(
    func = train_custom_main,
    args = list(pkg_root = getOption("perception.pkg_root", NULL),
                depmap_path = shared$depmap_path,
                jobs_dir = jobs_dir),
    supervise = TRUE,
    poll_connection = FALSE
  )
  shared$jobs_dir <- jobs_dir
  shared$train_worker <- worker
  shared$worker_alive <- TRUE
  invisible(TRUE)
}

# Submit a training job. Standard DepMap -> shared master pool; custom upload
# -> this session's dedicated worker.
submit_train_job <- function(shared, params) {
  if (isTRUE(shared$depmap_is_standard)) {
    ensure_master()
    jobs_dir <- getOption("perception.jobs_dir")
    # Safety net: make sure the master knows the path even if the load-time
    # notification was missed (idempotent — master skips identical paths).
    notify_master_depmap(shared$depmap_path)
  } else {
    ensure_custom_worker(shared)
    jobs_dir <- shared$jobs_dir
  }
  jobid <- make_job_id()
  job_dir <- file.path(jobs_dir, jobid)
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(params, file.path(job_dir, "params.rds"))
  file.create(file.path(job_dir, "progress.txt"))
  shared$active_job <- jobid
  shared$jobs_dir <- jobs_dir
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
  if (s == "done") return(list(status = "done"))
  if (startsWith(s, "error")) {
    return(list(status = "error", message = sub("^error\\|", "", s)))
  }
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
