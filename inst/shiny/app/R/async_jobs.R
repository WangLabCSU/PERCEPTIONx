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

# Spawn a callr background worker with retries. processx_exec can transiently
# fail under resource pressure (first spawn in a busy container, cgroup limits,
# fork contention with other sessions) even though the R install itself is
# fine — so we retry briefly before giving up. On final failure the error
# message carries the diagnostics needed to debug a stubborn deployment
# (Rscript path, PATH, cwd, tempdir) instead of a bare "processx_exec failed".
#
# WHY A LADDER OF CONFIGS: restricted containers reject individual processx
# configurations in deployment-specific ways — TRUE/FALSE stdio got EACCES
# (fixed by /dev/null on some servers), the supervisor's socketpair or a
# stderr file redirect can be rejected on others. Pinning one config makes the
# app depend on which exact failure mode the host hits. Instead we try
# progressively simpler configs and keep the first one the environment
# accepts. Each config gets the full retry budget, and a failed attempt runs
# gc() so half-open connections from the failed spawn do not accumulate.
spawn_r_bg_retry <- function(func, args, jobs_dir, retries = 2L, delay_ms = 500L) {
  configs <- list(
    # 1. Preferred: supervisor on, stdout discarded, stderr to the job log.
    list(supervise = TRUE,  stdout = "/dev/null", stderr = file.path(jobs_dir, "worker.log")),
    # 2. Supervisor and stderr file are the two most container-hostile bits;
    #    drop both, keep stdout discarded. (On Windows "/dev/null" does not
    #    exist, so configs 1-2 fail there and 3 below is what runs locally.)
    list(supervise = FALSE, stdout = "/dev/null", stderr = "/dev/null"),
    # 3. Most basic: inherit stdio (no pipes to fill, nothing to redirect,
    #    no supervisor). Always accepted unless fork itself is blocked.
    list(supervise = FALSE, stdout = NULL, stderr = NULL)
  )
  last_err <- NULL
  config_err <- character(0)
  for (cfg in configs) {
    for (attempt in seq_len(retries + 1L)) {
      tryCatch(
        return(callr::r_bg(
          func = func,
          args = args,
          supervise = cfg$supervise,
          poll_connection = FALSE,
          stdout = cfg$stdout,
          stderr = cfg$stderr
        )),
        error = function(e) {
          last_err <<- e
          if (attempt <= retries) Sys.sleep(delay_ms / 1000)
          NULL
        }
      )
    }
    config_err <- c(config_err, conditionMessage(last_err))
    gc()  # release any half-open connections before the next config
  }
  stop(
    "Failed to start the background worker after ", length(configs), " configs x ",
    retries + 1L, " attempts. ",
    "Last errors: ", paste(config_err, collapse = " | "), ". ",
    "Diagnostics: R.home(bin)=", R.home("bin"),
    "; Rscript exists=", file.exists(file.path(R.home("bin"), "Rscript")),
    "; Sys.which(Rscript)=", Sys.which("Rscript"),
    "; PATH=", Sys.getenv("PATH"),
    "; getwd()=", getwd(),
    "; tempdir=", tempdir(),
    "; jobs_dir exists=", dir.exists(jobs_dir),
    "; jobs_dir writable=", isTRUE(file.access(jobs_dir, 2) == 0),
    "; open connections=", length(connections()),
    "; R version=", paste(R.version$major, R.version$minor, sep = "."),
    "; callr=", tryCatch(as.character(utils::packageVersion("callr")), error = function(e) "?"),
    "; processx=", tryCatch(as.character(utils::packageVersion("processx")), error = function(e) "?"),
    call. = FALSE
  )
}

# ---- Master worker: runs inside the shared background R process ------------
# Fully self-contained (see note above).

train_master_main <- function(pkg_root, jobs_dir, max_parallel, idle_minutes) {
  # Loading strategy: dev mode (repo on disk) -> load_all; otherwise the
  # package is expected to be INSTALLED (e.g. via
  # remotes::install_github("WangLabCSU/PERCEPTIONx")). load_all() would fail
  # on a deployed machine because there is no source tree (no DESCRIPTION);
  # library() is the correct path once the package is installed into .libPaths().
  # A dir counts as the source tree ONLY with DESCRIPTION + R/*.R — an installed
  # package dir has DESCRIPTION but no .R sources, and load_all() on it yields
  # an empty namespace ("X is not an exported object from namespace:PERCEPTIONx").
  pkg_is_src <- !is.null(pkg_root) && nzchar(pkg_root) && dir.exists(pkg_root) &&
    length(list.files(file.path(pkg_root, "R"), pattern = "\\.R$")) > 0
  if (pkg_is_src) {
    suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  } else if (requireNamespace("PERCEPTIONx", quietly = TRUE)) {
    library(PERCEPTIONx)
  } else {
    stop("Package 'PERCEPTIONx' is not installed. Install it with ",
         "remotes::install_github(\"WangLabCSU/PERCEPTIONx\"), then restart the app.")
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
    params <- tryCatch(
      readRDS(file.path(job_dir, "params.rds")),
      error = function(e) structure(list(message = conditionMessage(e)),
                                    class = "perception_job_error")
    )
    if (inherits(params, "perception_job_error")) {
      writeLines(sprintf("error|Failed to read job params: %s", params$message),
                 status_file)
      return(invisible(TRUE))
    }
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
    # Pending jobs are needed both for DepMap-load error reporting and for
    # claiming below, so compute them once per loop iteration.
    # Only dirs with a "ready" marker are claimable: the app writes
    # params.rds FIRST, then the marker. Without this handshake a big,
    # slow saveRDS (e.g. a large expression matrix) could be claimed
    # mid-write and read as a truncated/corrupt file.
    pending <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
    pending <- pending[!vapply(pending, function(j)
      file.exists(file.path(j, "status.txt")), logical(1))]
    pending <- pending[vapply(pending, function(j)
      file.exists(file.path(j, "ready")), logical(1))]

    if (!is.null(depmap_path) && is.null(DepMap)) {
      message("[master] loading DepMap from ", depmap_path)
      loaded <- tryCatch(
        readRDS(depmap_path),
        error = function(e) structure(list(message = conditionMessage(e)),
                                      class = "perception_job_error")
      )
      if (inherits(loaded, "perception_job_error")) {
        # A bad file must never crash the master (that would strand every
        # queued job as "queued" forever). Report the failure to each pending
        # job and drop the path so a fresh notify can retry.
        message("[master] DepMap load FAILED: ", loaded$message)
        for (j in pending) {
          writeLines(sprintf("error|DepMap load failed: %s", loaded$message),
                     file.path(j, "status.txt"))
        }
        DepMap <- NULL
        depmap_path <- NULL
        last_activity <- Sys.time()
      } else {
        DepMap <- loaded
        assign_depmap(DepMap)
        message("[master] DepMap ready (", length(DepMap), " datasets)")
        last_activity <- Sys.time()
      }
    }

    # 2. Claim pending jobs (file.create wins; only this master claims, so no
    #    double-processing).
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
      # SAFE here (unlike in the Shiny session): this runs inside the worker's
      # OWN process, which is about to exit anyway, so closing its leftover
      # connections only avoids "closing unused connection N" noise in the
      # server log. It can NOT corrupt any other process. (In the Shiny app
      # process closeAllConnections() is fatal — it kills the session's data
      # channel, making callr workers fail with "invalid connection" and the
      # page disconnect — so app.R must never do this.)
      suppressWarnings(try(closeAllConnections(), silent = TRUE))
      quit(save = "no")
    }

    Sys.sleep(2)
  }
}

# ---- Dedicated worker for users who UPLOADED their own DepMap --------------
# A wrong upload must never reach the shared master; this per-session worker
# loads the user's own file and dies with the session (supervise).

train_custom_main <- function(pkg_root, depmap_path, jobs_dir, poll_secs = 1L) {
  # Source tree check (DESCRIPTION + R/*.R): an installed package dir has
  # DESCRIPTION but no .R sources, and load_all() on it yields an EMPTY
  # namespace ("X is not an exported object from namespace:PERCEPTIONx").
  pkg_is_src <- !is.null(pkg_root) && nzchar(pkg_root) && dir.exists(pkg_root) &&
    length(list.files(file.path(pkg_root, "R"), pattern = "\\.R$")) > 0
  if (pkg_is_src) {
    suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  } else if (requireNamespace("PERCEPTIONx", quietly = TRUE)) {
    library(PERCEPTIONx)
  } else {
    stop("Package 'PERCEPTIONx' is not installed. Install it with ",
         "remotes::install_github(\"WangLabCSU/PERCEPTIONx\"), then restart the app.")
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
    params <- tryCatch(
      readRDS(file.path(job_dir, "params.rds")),
      error = function(e) structure(list(message = conditionMessage(e)),
                                    class = "perception_job_error")
    )
    if (inherits(params, "perception_job_error")) {
      writeLines(sprintf("error|Failed to read job params: %s", params$message),
                 status_file)
      return(invisible(TRUE))
    }
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
  DepMap <- tryCatch(
    readRDS(depmap_path),
    error = function(e) structure(list(message = conditionMessage(e)),
                                  class = "perception_job_error")
  )
  if (inherits(DepMap, "perception_job_error")) {
    # The uploaded file may already be gone (e.g. its session ended and
    # Shiny cleaned up the upload). Nothing useful left to do: exit, and
    # let the app's worker-death detection surface the failure.
    message("[custom-worker] DepMap load FAILED: ", DepMap$message)
    quit(save = "no")
  }
  assign_depmap(DepMap)
  message("[custom-worker] DepMap ready, polling ", jobs_dir)
  while (TRUE) {
    # Only claim dirs whose params.rds is fully written (see master note).
    pending <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
    pending <- pending[vapply(pending, function(j)
      file.exists(file.path(j, "ready")), logical(1))]
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
  # Shared master jobs dir: under options(PERCEPTIONx.cache_root)/jobs when a
  # cache root is set (deployment-friendly), else tempdir()/perception_jobs.
  jobs_dir <- PERCEPTIONx:::perception_jobs_dir()
  dir.create(jobs_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(jobs_dir) || file.access(jobs_dir, 2) != 0) {
    stop("Cannot use master jobs dir at '", jobs_dir,
         "' (unwritable or not creatable). Fix the permission, or unset ",
         "PERCEPTIONx.cache_root to fall back to tempdir().")
  }
  options(perception.jobs_dir = jobs_dir)
  max_par <- suppressWarnings(as.integer(Sys.getenv("PERCEPTION_WORKERS", "16")))
  if (is.na(max_par) || max_par < 1L) max_par <- 16L
  idle <- suppressWarnings(as.integer(Sys.getenv("PERCEPTION_WORKER_IDLE_MINUTES", "720")))
  if (is.na(idle) || idle < 1L) idle <- 720L
  worker <- spawn_r_bg_retry(
    train_master_main,
    # NOTE: pkg_root default must be NULL (not ".") — NULL makes the worker
    # fall back to library(PERCEPTIONx); "." would trigger devtools::load_all
    # on an unrelated working directory and crash the worker.
    args = list(pkg_root = getOption("perception.pkg_root", NULL),
                jobs_dir = jobs_dir,
                max_parallel = max_par,
                idle_minutes = idle),
    jobs_dir = jobs_dir
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
    # The worker already loaded SOME DepMap into memory. If the user has since
    # uploaded a DIFFERENT file, keep the worker but make it reload: training
    # must use the current upload, never a stale copy.
    if (identical(shared$worker_depmap_path, shared$depmap_path)) {
      return(invisible(TRUE))
    }
    old_dir <- shared$jobs_dir
    shared$train_worker$kill()
    shared$worker_alive <- FALSE
    if (!is.null(old_dir) && dir.exists(old_dir)) unlink(old_dir, recursive = TRUE)
  }
  if (is.null(shared$depmap_path)) {
    stop("DepMap file path is unknown. Please load DepMap data first (Data tab).")
  }
  jobs_dir <- file.path(tempdir(), "perception_jobs_custom", make_job_id())
  dir.create(jobs_dir, recursive = TRUE, showWarnings = FALSE)
  worker <- spawn_r_bg_retry(
    train_custom_main,
    args = list(pkg_root = getOption("perception.pkg_root", NULL),
                depmap_path = shared$depmap_path,
                jobs_dir = jobs_dir),
    jobs_dir = jobs_dir
  )
  shared$jobs_dir <- jobs_dir
  shared$train_worker <- worker
  shared$worker_alive <- TRUE
  shared$worker_depmap_path <- shared$depmap_path
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
  # Handshake: params.rds must be fully written BEFORE the worker is allowed
  # to claim this job (workers only pick up dirs with the ready marker).
  file.create(file.path(job_dir, "ready"))
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

# ---------------------------------------------------------------------------
# Generic per-session task worker (prepare / demo / predict / plot /
# download / extract_meta)
#
# Training goes to the shared master (or a dedicated worker for uploads)
# because it needs the big DepMap. Everything else that computes on
# SESSION-PRIVATE objects (user expression, models, prepared data) runs in a
# light per-session worker so the main Shiny thread never blocks:
#   - "prepare"      -> PERCEPTIONx::prepare_data()          (Seurat clustering)
#   - "demo"         -> full demo pipeline (structured data + clustering + 2 models)
#   - "predict"      -> predict_drugs() + predict_patients()
#   - "plot"         -> plot math for the Visualize tab
#   - "download"     -> DepMap download (mirror fallback)
#   - "extract_meta" -> lightweight metadata from an existing DepMap.RDS
# Inputs travel with the job (params.rds <- list(task, args)) because the
# worker is a separate process and cannot see the parent's memory. Same file
# protocol as training: status.txt / progress.txt / result.rds.
# Self-contained rule (callr): all helpers live INSIDE the entry function.
# ---------------------------------------------------------------------------

session_task_main <- function(pkg_root, jobs_dir, poll_secs = 0.25) {
  # Source tree check (DESCRIPTION + R/*.R): an installed package dir has
  # DESCRIPTION but no .R sources, and load_all() on it yields an EMPTY
  # namespace ("X is not an exported object from namespace:PERCEPTIONx").
  pkg_is_src <- !is.null(pkg_root) && nzchar(pkg_root) && dir.exists(pkg_root) &&
    length(list.files(file.path(pkg_root, "R"), pattern = "\\.R$")) > 0
  if (pkg_is_src) {
    suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  } else if (requireNamespace("PERCEPTIONx", quietly = TRUE)) {
    library(PERCEPTIONx)
  } else {
    stop("Package 'PERCEPTIONx' is not installed. Install it with ",
         "remotes::install_github(\"WangLabCSU/PERCEPTIONx\"), then restart the app.")
  }

  run_prepare <- function(a) {
    PERCEPTIONx::prepare_data(
      method            = a$method,
      expression_matrix = a$expression_matrix,
      patient_mapping   = a$patient_mapping,
      seurat_resolution = a$seurat_resolution,
      seurat_dims       = a$seurat_dims,
      seurat_nfeatures  = if (!is.null(a$seurat_nfeatures)) a$seurat_nfeatures else 2000,
      seurat_seed       = 42,
      skip_clustering   = isTRUE(a$skip_clustering)
    )
  }

  run_demo <- function(a) {
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

    response_labels <- c(rep("Responder", 10), rep("Non-responder", 10))
    clinical_response <- data.frame(patient = patient_names, response = response_labels,
                                    stringsAsFactors = FALSE)

    patient_assignment <- rep(patient_names, each = ceiling(n_cells / n_patients))[1:n_cells]
    patient_mapping <- data.frame(cell_id = cell_names, patient_id = patient_assignment,
                                  stringsAsFactors = FALSE)

    is_responder_cell <- clinical_response$response[
      match(patient_assignment, clinical_response$patient)] == "Responder"
    abemaciclib_markers <- gene_names[1:5]
    erlotinib_markers   <- gene_names[6:10]
    noise_genes         <- gene_names[11:length(gene_names)]

    expr_matrix <- matrix(0.1, nrow = length(gene_names), ncol = n_cells)
    rownames(expr_matrix) <- gene_names
    colnames(expr_matrix) <- cell_names
    for (g in abemaciclib_markers) {
      expr_matrix[g, is_responder_cell] <- pmax(rnorm(sum(is_responder_cell), mean = 8, sd = 3), 0.1)
      expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 3, sd = 2), 0.1)
    }
    for (g in erlotinib_markers) {
      expr_matrix[g, is_responder_cell] <- pmax(rnorm(sum(is_responder_cell), mean = 3, sd = 2), 0.1)
      expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 8, sd = 3), 0.1)
    }
    for (g in noise_genes) expr_matrix[g, ] <- runif(n_cells, 0.5, 8)
    storage.mode(expr_matrix) <- "numeric"

    prepared <- PERCEPTIONx::prepare_data(
      method = "umap", expression_matrix = expr_matrix, patient_mapping = patient_mapping,
      seurat_resolution = 0.8, seurat_dims = 10,
      seurat_nfeatures = min(2000, length(gene_names)))

    make_structured_training <- function(marker_genes, direction, n_train = 160) {
      x_train <- matrix(0, nrow = length(gene_names), ncol = n_train)
      rownames(x_train) <- gene_names
      half <- n_train %/% 2
      responder_like <- seq_len(half)
      nonresponder_like <- seq.int(half + 1, n_train)
      for (g in marker_genes) {
        x_train[g, responder_like] <- pmax(rnorm(half, mean = 8, sd = 3), 0.1)
        x_train[g, nonresponder_like] <- pmax(rnorm(half, mean = 3, sd = 2), 0.1)
      }
      for (g in setdiff(gene_names, marker_genes)) x_train[g, ] <- runif(n_train, 0.5, 8)
      y_train <- direction * colMeans(x_train[marker_genes, , drop = FALSE])
      y_train <- y_train + rnorm(n_train, sd = 0.3)
      as.data.frame(cbind(y = y_train, t(x_train)))
    }
    make_drug_model <- function(drug_name, seed, marker_genes, direction) {
      set.seed(seed)
      train_df <- make_structured_training(marker_genes, direction, n_train = 160)
      caret_model <- suppressWarnings(caret::train(
        y ~ ., data = train_df, method = "glmnet",
        trControl = caret::trainControl(method = "cv", number = 3),
        tuneLength = 3))
      obj <- list(
        model = caret_model,
        performance_in_scRNA = data.frame(estimate.cor = c(0.45, 0.38), p.value = c(0.001, 0.01)),
        performance_in_bulk = data.frame(estimate.cor = c(0.55, 0.48), p.value = c(0.0001, 0.001)),
        performance_in_pseudo_bulk = data.frame(estimate.cor = c(0.50, 0.42), p.value = c(0.0005, 0.005)),
        predVSgroundTruth = list(pred_gt_scRNA = data.frame(Observed = rnorm(20), Test_pred_sc = rnorm(20))),
        single_best = marker_genes[1])
      attr(obj, "drug_name") <- drug_name
      # SIMULATED demo model: trained on random noise, NOT on real drug
      # response data. UI must label it so users never mistake its output for
      # a real prediction.
      attr(obj, "model_source") <- "demo"
      obj
    }

    models <- list(
      abemaciclib = make_drug_model("abemaciclib", 101, abemaciclib_markers, direction = -1),
      erlotinib   = make_drug_model("erlotinib",   202, erlotinib_markers,   direction = +1))

    list(
      user_response = clinical_response,
      user_mapping  = patient_mapping,
      user_expr     = expr_matrix,
      prepared_data = prepared,
      user_clones   = prepared$cell_clone_map,
      models        = models
    )
  }

  run_predict <- function(a) {
    clone_pred <- PERCEPTIONx::predict_drugs(model_list = a$model_list, expr = a$expr)
    patient_pred <- NULL
    if (!is.null(a$patient_input)) {
      # Standard path: prepared_data-derived template + counts (small).
      patient_pred <- PERCEPTIONx::predict_patients(clone_pred, a$patient_input,
                                                    mode = a$mode)
    } else if (!is.null(a$legacy_clones)) {
      # Legacy path (no prepared_data): aggregate from the user clone map.
      clone_data <- a$legacy_clones
      clone_rows <- rownames(clone_pred)
      clone_to_patient <- split(clone_data$patient, clone_data$clone_id)
      clone_to_patient <- lapply(clone_to_patient, unique)
      clone_viability_list <- lapply(clone_rows, function(cl) {
        pat <- clone_to_patient[[cl]]
        if (is.null(pat) || length(pat) == 0) return(NULL)
        if (length(pat) > 1) pat <- pat[1]
        df <- data.frame(patient = pat, clone_id = cl, stringsAsFactors = FALSE)
        for (drug in colnames(clone_pred)) df[[drug]] <- clone_pred[cl, drug]
        df
      })
      clone_viability_df <- do.call(rbind, clone_viability_list)
      if (is.null(clone_viability_df) || nrow(clone_viability_df) == 0) {
        stop("No matching clones between prediction and annotation")
      }
      clone_counts <- as.data.frame.matrix(table(clone_data$patient, clone_data$clone_id))
      clone_counts$patients <- rownames(clone_counts)
      patient_pred <- PERCEPTIONx::predict_patients(clone_viability_df, clone_counts,
                                                    mode = a$mode)
    }
    list(clone_pred = clone_pred, patient_pred = patient_pred)
  }

  run_plot <- function(a) {
    pt <- a$plot_type
    d <- a$data
    p <- a$params

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
    safe_range01 <- function(x) {
      if (length(x) == 0L) stop("empty input to safe_range01")
      r <- tryCatch(PERCEPTIONx::range01(x), error = function(e) NULL)
      if (!is.null(r) && length(r) == length(x)) return(r)
      rng <- range(x, na.rm = TRUE)
      # Constant vector: map to 0 ("no signal"), not the old 0.5 midpoint.
      if (rng[2] == rng[1]) rep(0, length(x)) else (x - rng[1]) / (rng[2] - rng[1])
    }
    get_embedding_xy <- function(coords, cell_ids) {
      if (!is.data.frame(coords)) stop("umap_coords is not a data frame")
      if (!"cell_id" %in% names(coords)) stop("umap_coords has no 'cell_id' column")
      x_col <- intersect(c("dim_1", "umap_1"), names(coords))[1]
      y_col <- intersect(c("dim_2", "umap_2"), names(coords))[1]
      if (is.na(x_col) || is.na(y_col))
        stop("Embedding coordinate columns not found in umap_coords")
      idx <- match(cell_ids, coords$cell_id)
      if (all(is.na(idx))) stop("0 matching cell IDs between query and umap_coords")
      list(X = coords[[x_col]][idx], Y = coords[[y_col]][idx])
    }
    combo_clone_frame <- function() {
      pred_mat <- d$predictions
      if (is.null(pred_mat) || ncol(pred_mat) < 1) stop("No clone-level predictions to combine.")
      pred_mat <- as.matrix(pred_mat)
      z_mat <- pred_mat
      for (j in seq_len(ncol(pred_mat))) {
        col <- pred_mat[, j]
        s <- stats::sd(col, na.rm = TRUE)
        z_mat[, j] <- if (is.na(s) || s == 0) rep(0, nrow(pred_mat)) else
          (col - mean(col, na.rm = TRUE)) / s
      }
      comb <- apply(z_mat, 1, min, na.rm = TRUE)
      tmpl <- d$prepared$clone_viability_template
      if (!is.null(tmpl)) {
        clone_viability_df <- data.frame(patient = tmpl$patient, clone_id = tmpl$clone_id,
                                         comb_viability = unname(comb), stringsAsFactors = FALSE)
      } else {
        parsed <- PERCEPTIONx::parse_clone_keys(names(comb))
        clone_viability_df <- data.frame(patient = parsed$patient, clone_id = parsed$clone_id,
                                         comb_viability = unname(comb), stringsAsFactors = FALSE)
      }
      clone_data <- d$user_clones
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
      if (!is.null(d$user_response)) {
        clone_viability_df$response <- label_resp(d$user_response$response[
          match(clone_viability_df$patient, d$user_response$patient)])
      }
      clone_viability_df
    }
    combo_patient_frame <- function(clone_viability_df) {
      patients <- unique(clone_viability_df$patient)
      scores <- vapply(patients, function(p) {
        sub <- clone_viability_df[clone_viability_df$patient == p, ]
        sub <- sub[!is.na(sub$comb_viability) & !is.na(sub$weights), ]
        if (nrow(sub) == 0) return(NA_real_)
        max(sub$comb_viability * sub$weights)
      }, numeric(1))
      data.frame(patient = patients, combination = scores, stringsAsFactors = FALSE)
    }
    response_groups <- function() {
      cr <- d$user_response
      if (is.null(cr) || is.null(cr$response)) return(character(0))
      grps <- unique(label_resp(cr$response))
      grps[!is.na(grps) & nzchar(trimws(grps)) & grps != "NA"]
    }

    is_combo <- isTRUE(p$combo)
    drug <- p$drug
    plot_obj <- NULL
    msg <- NULL

    if (pt == "clone_dist") {
      clone_data <- d$user_clones
      clone_dist_list <- lapply(unique(clone_data$patient), function(pp) {
        p_cells <- clone_data[clone_data$patient == pp, ]
        p_clones <- unique(p_cells$clone_id)
        n_p_cells <- nrow(p_cells)
        display_clones <- vapply(p_clones, function(cl) {
          if (!startsWith(cl, pp)) return(cl)
          rest <- sub("^[\\._-]", "", substring(cl, nchar(pp) + 1))
          if (nchar(rest) == 0) cl else rest
        }, character(1))
        data.frame(patients = pp, clones = display_clones,
                   weights = sapply(p_clones, function(cl) sum(p_cells$clone_id == cl) / n_p_cells),
                   stringsAsFactors = FALSE)
      })
      clone_distribution <- do.call(rbind, clone_dist_list)
      if (!is.null(d$user_response)) {
        clone_distribution$response <- label_resp(d$user_response$response[
          match(clone_distribution$patients, d$user_response$patient)])
        clone_distribution$response <- label_resp(clone_distribution$response)
        clone_distribution <- clone_distribution[!is.na(clone_distribution$response), ]
      }
      plot_obj <- PERCEPTIONx::plot_clone_distribution(clone_distribution, response_var = "response")

    } else if (pt == "clone_kill") {
      clone_data <- d$user_clones
      pred_mat <- d$predictions
      drug <- if (!is.null(drug) && nchar(drug) > 0) drug else colnames(pred_mat)[1]
      if (is_combo) {
        combo_df <- combo_clone_frame()
        combo_df <- combo_df[!is.na(combo_df$response), ]
        if (nrow(combo_df) == 0) stop("No clones with a valid R/NR response for the combination plot.")
        plot_obj <- PERCEPTIONx::plot_clone_viability(combo_df, viability_var = "comb_viability",
                                                      weights_var = "weights", response_var = "response",
                                                      drug = "Combination")
      } else {
        clone_viability_df <- NULL
        tmpl <- d$prepared$clone_viability_template
        if (!is.null(tmpl)) {
          clone_viability_df <- data.frame(
            patient = tmpl$patient, clone_id = tmpl$clone_id,
            comb_viability = pred_mat[rownames(pred_mat), drug],
            stringsAsFactors = FALSE)
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
          parsed <- PERCEPTIONx::parse_clone_keys(rownames(pred_mat))
          pred_clone_ids <- parsed$clone_id
          pred_patients  <- parsed$patient
          clone_kill_list <- lapply(unique(clone_data$patient), function(pp) {
            p_clones <- unique(clone_data$clone_id[clone_data$patient == pp])
            p_clones <- intersect(p_clones, pred_clone_ids)
            if (length(p_clones) == 0) return(NULL)
            n_p_cells <- sum(clone_data$patient == pp)
            pred_rows <- which(pred_clone_ids %in% p_clones & pred_patients == pp)
            if (length(pred_rows) == 0) pred_rows <- match(p_clones, pred_clone_ids)
            data.frame(patient = pp, clone_id = pred_clone_ids[pred_rows],
                       comb_viability = pred_mat[pred_rows, drug],
                       weights = sapply(pred_clone_ids[pred_rows], function(cl) sum(clone_data$clone_id == cl) / n_p_cells),
                       stringsAsFactors = FALSE)
          })
          clone_viability_df <- do.call(rbind, clone_kill_list)
        }
        if (is.null(clone_viability_df) || nrow(clone_viability_df) == 0) {
          stop("No matching clones between prediction matrix and clone annotation.")
        }
        if (!is.null(d$user_response)) {
          clone_viability_df$response <- label_resp(d$user_response$response[
            match(clone_viability_df$patient, d$user_response$patient)])
          clone_viability_df <- clone_viability_df[!is.na(clone_viability_df$response), ]
        }
        clone_viability_df$comb_viability <- as.numeric(scale(clone_viability_df$comb_viability))
        plot_obj <- PERCEPTIONx::plot_clone_viability(clone_viability_df, viability_var = "comb_viability",
                                                      weights_var = "weights", response_var = "response",
                                                      drug = drug)
      }

    } else if (pt == "roc") {
      cr <- d$user_response
      if (is_combo) {
        combo_df <- combo_clone_frame()
        pat_df <- combo_patient_frame(combo_df)
        rv <- label_resp(cr$response[match(pat_df$patient, cr$patient)])
        predictor_vec <- pat_df$combination
      } else {
        pp <- d$patient_pred
        drug <- if (!is.null(drug) && nchar(drug) > 0) drug else colnames(pp)[1]
        rv <- label_resp(cr$response[match(rownames(pp), cr$patient)])
        predictor_vec <- pp[[drug]]
      }
      keep <- !is.na(rv) & !is.na(predictor_vec)
      rv <- rv[keep]
      predictor_vec <- predictor_vec[keep]
      grps <- response_groups()
      if (length(grps) > 2) {
        a <- p$roc_group_a; b <- p$roc_group_b
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
      smooth <- length(response_vec) >= 10
      title <- if (is_combo) "Combination" else drug
      if (length(grps) > 2) title <- paste0(title, " — ", paste(sel, collapse = " vs "))
      plot_obj <- PERCEPTIONx::plot_roc_curve(response = response_vec, predictor = predictor_vec,
                                              smooth_curve = smooth, title = title)

    } else if (pt == "boxplot") {
      cr <- d$user_response
      if (is_combo) {
        combo_df <- combo_clone_frame()
        pat_df <- combo_patient_frame(combo_df)
        rv <- label_resp(cr$response[match(pat_df$patient, cr$patient)])
        predictor_vec <- pat_df$combination
      } else {
        pp <- d$patient_pred
        drug <- if (!is.null(drug) && nchar(drug) > 0) drug else colnames(pp)[1]
        rv <- label_resp(cr$response[match(rownames(pp), cr$patient)])
        predictor_vec <- pp[[drug]]
      }
      keep <- !is.na(rv) & !is.na(predictor_vec)
      rv <- rv[keep]
      predictor_vec <- predictor_vec[keep]
      predictor_vec <- if (length(predictor_vec) > 1 && is.finite(sd(predictor_vec)) && sd(predictor_vec) > 0) {
        as.numeric(scale(predictor_vec))
      } else predictor_vec
      exp_vs_pred <- data.frame(response = factor(rv, levels = unique(rv)),
                                predicted_viability = predictor_vec, stringsAsFactors = FALSE)
      plot_obj <- PERCEPTIONx::plot_response_boxplot(exp_vs_pred)
      plot_obj <- plot_obj + ggplot2::ggtitle(if (is_combo) "Combination" else drug)

    } else if (pt %in% c("umap_gene", "umap_viability", "umap_clone")) {
      umap_coords <- d$prepared$umap_coords
      clone_data <- d$user_clones
      common_cells <- intersect(clone_data$cell_id, umap_coords$cell_id)
      if (length(common_cells) == 0) {
        stop("No matching cells between embedding coordinates and clone annotation.")
      }
      if (pt == "umap_gene") {
        gene <- p$umap_gene
        # The app pre-extracts only this gene's per-cell expression vector
        # (named by cell id) — the whole matrix never travels to the worker.
        expr_vec <- d$umap_gene_expr
        if (is.null(gene) || nchar(gene) == 0) {
          msg <- "Select a gene first."
        } else if (is.null(expr_vec) || length(expr_vec) == 0) {
          msg <- paste0("Gene '", gene, "' not found in expression matrix.")
        } else {
          cell_expr <- as.numeric(expr_vec[match(common_cells, names(expr_vec))])
          xy <- get_embedding_xy(umap_coords, common_cells)
          n <- length(common_cells)
          stopifnot(length(xy$X) == n, length(xy$Y) == n, length(cell_expr) == n,
                    !anyNA(cell_expr))
          umap_data <- data.frame(X = xy$X, Y = xy$Y, expression = safe_range01(cell_expr),
                                  row.names = common_cells)
          plot_obj <- PERCEPTIONx::plot_tsne_response(umap_data, color_var = "expression",
                                                      title = gene, color_label = "Expression (0-1)",
                                                      colors = c("#e0e0e0", "#c13232"),
                                                      limits = c(0, 1), base_size = 11)
        }
      } else if (pt == "umap_viability") {
        pred_mat <- d$predictions
        if (is.null(pred_mat)) {
          msg <- "No clone-level predictions found. Run predictions first."
        } else {
          drug <- if (!is.null(p$umap_drug) && nchar(p$umap_drug) > 0) p$umap_drug else colnames(pred_mat)[1]
          cell_keys <- PERCEPTIONx::build_clone_key(clone_data$patient, clone_data$clone_id)
          pred_keys <- rownames(pred_mat)
          cell_viability <- setNames(pred_mat[match(cell_keys, pred_keys), drug], clone_data$cell_id)
          kill_common <- intersect(names(cell_viability), umap_coords$cell_id)
          kill_common <- kill_common[!is.na(cell_viability[kill_common])]
          if (length(kill_common) == 0) {
            msg <- "No matching cells between embedding coordinates and prediction data."
          } else {
            raw_vals <- cell_viability[kill_common]
            scaled_vals <- if (length(raw_vals) > 1 && is.finite(sd(raw_vals)) && sd(raw_vals) > 0) {
              as.numeric(scale(raw_vals))
            } else raw_vals
            xy <- get_embedding_xy(umap_coords, kill_common)
            n <- length(kill_common)
            stopifnot(length(xy$X) == n, length(xy$Y) == n, length(scaled_vals) == n)
            umap_data <- data.frame(X = xy$X, Y = xy$Y, viability_scaled = scaled_vals,
                                    row.names = kill_common)
            plot_obj <- PERCEPTIONx::plot_tsne_response(umap_data, color_var = "viability_scaled",
                                                        title = drug, color_label = "Predicted Viability (z-score)",
                                                        palette = "diverging", midpoint = 0,
                                                        base_size = 11)
          }
        }
      } else {  # umap_clone
        idx <- match(common_cells, clone_data$cell_id)
        xy <- get_embedding_xy(umap_coords, common_cells)
        n <- length(common_cells)
        stopifnot(length(xy$X) == n, length(xy$Y) == n, length(idx) == n)
        umap_data <- data.frame(X = xy$X, Y = xy$Y, clone_id = clone_data$clone_id[idx],
                                stringsAsFactors = FALSE)
        plot_obj <- PERCEPTIONx::plot_clone_umap(umap_data, title = "Clone Identity")
      }
    } else {
      stop("Unknown plot type: ", pt)
    }

    list(plot = plot_obj,
         message = if (!is.null(msg)) msg else
           if (is.null(plot_obj)) "Selected plot type is not available with current data" else "")
  }

  # Download with mirror fallback. Downloading AND parsing (reading the 567 MB
  # file) are split into two sequential tasks so the UI can show a proper
  # download progress bar, then switch to a "parsing" overlay — the main
  # process never deserializes the big object either way.
  run_download <- function(a) {
    ok <- PERCEPTIONx:::download_with_mirrors(
      a$urls, a$destfile, quiet = TRUE,
      timeout_seconds = a$timeout_seconds,
      retries = a$retries,
      expected_size = a$expected_size)
    if (!ok) stop("Automatic download failed after trying all mirrors")
    if (!file.exists(a$destfile) || file.size(a$destfile) == 0) {
      stop("Downloaded file is missing or empty")
    }
    a$destfile
  }

  # Rebuild/refresh the metadata sidecar for an EXISTING DepMap.RDS (used when
  # the meta cache is missing or stale). The full read stays in the worker.
  run_extract_meta <- function(a) {
    PERCEPTIONx:::extract_depmap_meta(a$depmap_path, cache_file = a$cache_file)
  }

  # Warm-up task: pre-load the heavy namespaces (Seurat / caret / glmnet) in
  # the background so the FIRST prepare/demo task does not pay the ~3-4 s
  # package load inline (Seurat is only touched on the first Seurat:: call).
  run_warmup <- function(a) {
    invisible(requireNamespace("Seurat", quietly = TRUE))
    invisible(requireNamespace("caret", quietly = TRUE))
    invisible(requireNamespace("glmnet", quietly = TRUE))
    TRUE
  }

  message("[session-worker] started, polling ", jobs_dir)
  while (TRUE) {
    # Only claim dirs whose params.rds is fully written (see master note).
    pending <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
    pending <- pending[vapply(pending, function(j)
      file.exists(file.path(j, "ready")), logical(1))]
    for (j in pending) {
      status_file <- file.path(j, "status.txt")
      if (!file.exists(status_file) && file.create(status_file)) {
        params <- tryCatch(
          readRDS(file.path(j, "params.rds")),
          error = function(e) structure(list(message = conditionMessage(e)),
                                        class = "perception_job_error")
        )
        if (inherits(params, "perception_job_error")) {
          # A corrupt params file must not crash the worker (that would strand
          # every later job in this session's queue).
          writeLines(sprintf("error|Failed to read job params: %s", params$message),
                     status_file)
          next
        }
        out <- tryCatch(
          switch(params$task,
            prepare = run_prepare(params$args),
            demo    = run_demo(params$args),
            predict = run_predict(params$args),
            plot    = run_plot(params$args),
            download = run_download(params$args),
            extract_meta = run_extract_meta(params$args),
            warmup  = run_warmup(params$args),
            stop("Unknown task: ", params$task)),
          error = function(e) structure(list(message = conditionMessage(e)),
                                        class = "perception_job_error")
        )
        if (inherits(out, "perception_job_error")) {
          writeLines(sprintf("error|%s", out$message), status_file)
        } else {
          saveRDS(out, file.path(j, "result.rds"))
          writeLines("done", status_file)
        }
      }
    }
    Sys.sleep(poll_secs)
  }
}

# ---- App-side manager for session tasks --------------------------------

# Lazily spawn ONE light worker per session for prepare/demo/predict jobs.
ensure_session_worker <- function(shared) {
  if (isTRUE(shared$task_worker_alive) && !is.null(shared$task_worker) &&
      shared$task_worker$is_alive()) {
    return(invisible(TRUE))
  }
  jobs_dir <- file.path(tempdir(), "perception_tasks", make_job_id())
  # If a previous worker for this session died, its job dir is abandoned —
  # drop it so stale task files (possibly large result.rds) do not accumulate
  # in tempdir across repeated worker restarts.
  if (!is.null(shared$task_jobs_dir) && nzchar(shared$task_jobs_dir)) {
    unlink(shared$task_jobs_dir, recursive = TRUE)
  }
  dir.create(jobs_dir, recursive = TRUE, showWarnings = FALSE)
  worker <- tryCatch(
    spawn_r_bg_retry(
      session_task_main,
      args = list(pkg_root = getOption("perception.pkg_root", NULL),
                  jobs_dir = jobs_dir),
      jobs_dir = jobs_dir
    ),
    error = function(e) {
      # Spawn failed after retries: drop the abandoned job dir so stale task
      # files cannot accumulate, then surface the diagnostic error.
      unlink(jobs_dir, recursive = TRUE)
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  shared$task_jobs_dir <- jobs_dir
  shared$task_worker <- worker
  shared$task_worker_alive <- TRUE
  invisible(TRUE)
}

# Submit a generic task. args are serialized with the job, so pass only the
# objects the worker needs (they must be serializable).
submit_session_task <- function(shared, task, args) {
  ensure_session_worker(shared)
  jobid <- make_job_id()
  job_dir <- file.path(shared$task_jobs_dir, jobid)
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(task = task, args = args), file.path(job_dir, "params.rds"))
  # Handshake: params.rds must be fully written BEFORE the worker is allowed
  # to claim this job (workers only pick up dirs with the ready marker).
  file.create(file.path(job_dir, "ready"))
  file.create(file.path(job_dir, "progress.txt"))
  shared$active_task <- jobid
  invisible(jobid)
}

# Read the state of a generic task (queued / running / done / error).
# "done" carries the deserialized result object back into the app.
read_task_state <- function(shared, jobid) {
  job_dir <- file.path(shared$task_jobs_dir, jobid)
  status_file <- file.path(job_dir, "status.txt")
  if (!file.exists(status_file)) return(list(status = "queued"))
  s <- readLines(status_file, warn = FALSE)
  if (length(s) == 0) return(list(status = "running"))
  s <- s[1]
  if (s == "done") return(list(status = "done", result = readRDS(file.path(job_dir, "result.rds"))))
  if (startsWith(s, "error")) {
    return(list(status = "error", message = sub("^error\\|", "", s)))
  }
  list(status = "running")
}

# Poll a generic task every second until it finishes, then run on_done(result)
# or on_error(message). The observer destroys itself when the task resolves,
# so one call per submission is enough (no manual cleanup needed). If the
# session worker died mid-task (is_alive == FALSE) the poll gives up with an
# actionable error instead of spinning forever.
poll_task <- function(shared, session, jobid, on_done, on_error = NULL,
                      poll_ms = 250) {
  obs <- NULL
  obs <- observe({
    st <- read_task_state(shared, jobid)
    # Helper: drop the finished job's on-disk files (result is already in
    # memory). Guard against a NULL jobs dir — file.path(NULL, id) would
    # degrade to a bare relative path that unlink could resolve in cwd.
    cleanup_job_dir <- function() {
      tj <- shared$task_jobs_dir
      if (!is.null(tj) && nzchar(tj)) unlink(file.path(tj, jobid), recursive = TRUE)
    }
    if (st$status == "done") {
      obs$destroy()
      # Result is already in memory (st$result): the on-disk job dir is no
      # longer needed. Remove it so tempdir does not fill with big result.rds
      # files over a long-running server (params.rds + result.rds can be
      # many MB per predict/plot/prepare call).
      cleanup_job_dir()
      on_done(st$result)
    } else if (st$status == "error") {
      obs$destroy()
      cleanup_job_dir()
      if (!is.null(on_error)) {
        on_error(st$message)
      } else {
        showNotification(paste("Background task failed:", st$message),
                         type = "error", duration = 10)
      }
    } else {
      w <- shared$task_worker
      if (!is.null(w) && !w$is_alive()) {
        # Worker crashed/exited without finishing: give up with a clear error.
        obs$destroy()
        cleanup_job_dir()
        if (!is.null(on_error)) {
          on_error("Background worker stopped unexpectedly. Please try again.")
        } else {
          showNotification("Background worker stopped unexpectedly. Please try again.",
                           type = "error", duration = 10)
        }
      } else {
        # NOTE: no explicit session arg here. Inside observe() the default
        # reactive domain IS this observer's session; on Shiny >= 1.10 an
        # explicit promise from the outer function frame can evaluate as a
        # missing argument at flush time ("argument 'session' is missing").
        invalidateLater(poll_ms)
      }
    }
  })
  invisible(obs)
}
