# End-to-end test of the ASYNC TRAINING pipeline on a real server.
#
# It drives the app's OWN code path — ensure_master() -> submit_train_job()
# -> read_job_state() — in a separate R session (like the Shiny main process),
# spawning the real callr background worker. No hand-written simulation of the
# job protocol.
#
# Usage:
#   R_LIBS=/data/home/dingjia/R/library Rscript dev/e2e_train.R
#       -> standard user: DepMap from the fixed cache dir (load_depmap flow)
#   R_LIBS=/data/home/dingjia/R/library Rscript dev/e2e_train.R --uploaded /path/to/DepMap.RDS
#       -> upload user: dedicated custom worker with the given DepMap file
#
# Asserts: both jobs finish (done), result.rds exists per job, and the master
# exits itself after idle_minutes (memory release). Exits non-zero on failure.

.libPaths(c('/data/home/dingjia/R/library', .libPaths()))
source('/data/home/dingjia/R/library/PERCEPTIONx/shiny/app/R/async_jobs.R')
options(perception.pkg_root = NULL)

# Test-only: shrink the master's idle-release window so this run also verifies
# the memory release. The app default stays 720 minutes.
Sys.setenv(PERCEPTION_WORKER_IDLE_MINUTES = '1')

args <- commandArgs(trailingOnly = TRUE)
uploaded <- NULL
if (length(args) >= 2 && args[1] == '--uploaded') uploaded <- args[2]

# ---- build the app-side `shared` state ------------------------------------
# An environment (reference semantics, like reactiveValues) — a plain list is
# pass-by-value, so the worker layer's `shared$jobs_dir <- ...` writes would
# be lost. reactiveValues itself errors outside a reactive context on Shiny
# 1.10, so an environment is the faithful non-Shiny stand-in.
cache_dir <- getOption(
  'PERCEPTIONx.depmap_cache_dir',
  Sys.getenv('PERCEPTIONx_DEPMAP_CACHE_DIR', tools::R_user_dir('PERCEPTIONx', 'data'))
)
shared <- new.env(parent = emptyenv())
if (!is.null(uploaded)) {
  shared$depmap_is_standard <- FALSE
  shared$depmap_path <- normalizePath(uploaded)
  cat('[e2e] SCENARIO: uploaded-DepMap user, file =', uploaded, '\n')
} else {
  shared$depmap_is_standard <- TRUE
  shared$depmap_path <- file.path(cache_dir, 'DepMap.RDS')
  cat('[e2e] SCENARIO: standard user (load_depmap cache), file =', shared$depmap_path, '\n')
}
if (!file.exists(shared$depmap_path)) {
  stop('[e2e] DepMap file not found: ', shared$depmap_path)
}
cat('[e2e] DepMap size:', round(file.info(shared$depmap_path)$size / 1e6), 'MB\n')

# ---- start worker + submit jobs via the app's own functions ---------------
t0 <- Sys.time()
if (isTRUE(shared$depmap_is_standard)) {
  ensure_master()
  worker <- getOption('perception.master')
  jobs_dir <- getOption('perception.jobs_dir')
  cat('[e2e] master spawned, jobs_dir =', jobs_dir, '\n')
} else {
  ensure_custom_worker(shared)
  worker <- shared$train_worker
  jobs_dir <- shared$jobs_dir
  cat('[e2e] custom worker spawned, jobs_dir =', jobs_dir, '\n')
}

mk_params <- function(drug, tag) {
  # PanCan for both include/exclude = the app's UI defaults. Do NOT pass NULL
  # for exclude_cancer: feature_ranking_bulk does grep(exclude_cancer, lineage)
  # and grep(NULL, ...) errors, which err_handle swallows into a silent
  # "feature ranking failed" (all models skipped).
  list(drug_list         = drug,
       cancer_type       = 'PanCan',
       exclude_cancer    = 'PanCan',
       GOI               = c('TP53', 'BRCA1'),
       k_features_values = c(50),
       model_type        = 'glmnet',
       ncores            = 1,
       output_dir        = file.path(jobs_dir, paste0('out_', tag)))
}

# erlotinib/gemcitabine verified to train cleanly on the clean release asset
# (abemaciclib hits a train_models data-level bug: "attempt to select less
# than one element in get1index" — a package robustness issue, not an async one).
jid1 <- submit_train_job(shared, mk_params('erlotinib', 'a'))
jid2 <- submit_train_job(shared, mk_params('gemcitabine', 'e'))
# shared is an environment: submit_train_job mutated it in place, re-read the
# jobs_dir it recorded (covers the uploaded scenario where it differs).
jobs_dir <- if (isTRUE(shared$depmap_is_standard)) getOption('perception.jobs_dir') else shared$jobs_dir
cat('[e2e] submitted jobs:', jid1, '|', jid2, '\n')

# ---- poll ----------------------------------------------------------------
deadline <- Sys.time() + 1500
done <- character(0)
while (Sys.time() < deadline) {
  if (!worker$is_alive()) {
    cat('[e2e] WORKER DIED unexpectedly; stderr tail:\n')
    cat(tail(readLines(file.path(jobs_dir, 'worker.log'), warn = FALSE), 30), sep = '\n')
    quit(status = 1)
  }
  for (jid in c(jid1, jid2)) {
    st <- read_job_state(shared, jid)
    if (st$status %in% c('done', 'error') && !(jid %in% done)) {
      done <- c(done, jid)
      msg <- if (st$status == 'error') paste0(' ERROR: ', st$message) else ''
      res <- file.exists(file.path(jobs_dir, jid, 'result.rds'))
      cat(sprintf('[e2e] %s -> %s (result.rds=%s) after %.0fs%s\n',
                  jid, st$status, res,
                  as.numeric(difftime(Sys.time(), t0, units = 'secs')), msg))
    } else if (st$status == 'running') {
      cat(sprintf('[e2e] %s running: phase=%s %d/%d drug=%s\n',
                  jid, st$phase, st$i, st$n, st$drug))
    }
  }
  if (length(done) == 2) break
  Sys.sleep(5)
}

if (length(done) != 2) {
  cat('[e2e] TIMEOUT; final states:\n')
  for (jid in c(jid1, jid2)) {
    st <- read_job_state(shared, jid)
    cat(jid, ': status =', st$status,
        if (!is.null(st$message)) paste0(' (', st$message, ')') else '', '\n')
  }
  quit(status = 1)
}

# result content sanity check
for (jid in c(jid1, jid2)) {
  res <- tryCatch(readRDS(file.path(jobs_dir, jid, 'result.rds')), error = function(e) NULL)
  if (is.null(res)) {
    cat('[e2e] FAIL: result.rds missing for', jid, '\n')
    quit(status = 1)
  }
  n <- if (is.list(res)) length(res) else NA
  cat(sprintf('[e2e] %s result: len=%s names=%s\n', jid, n,
              paste(if (is.list(res)) names(res) else '', collapse = ', ')))
  if (is.list(res) && length(res) > 0) {
    nm <- names(res)[1]
    m <- res[[1]]
    cat(sprintf('[e2e]   %s: class=%s fields=%s\n', nm,
                paste(class(m), collapse = '/'),
                paste(names(m), collapse = ', ')))
  }
}

# ---- worker must release itself after idle ---------------------------------
if (isTRUE(shared$depmap_is_standard)) {
  cat('[e2e] jobs done, waiting for master idle exit (idle_minutes=1)...\n')
  deadline2 <- Sys.time() + 180
  while (Sys.time() < deadline2) {
    if (!worker$is_alive()) {
      cat('[e2e] master exited after',
          round(as.numeric(difftime(Sys.time(), t0, units = 'secs'))), 's\n')
      break
    }
    Sys.sleep(2)
  }
  if (worker$is_alive()) {
    cat('[e2e] FAIL: master still alive after 180s (idle exit broken)\n')
    worker$kill()
    quit(status = 1)
  }
}
cat('[e2e] ALL PASS\n')
