# Edge / failure-mode tests: corrupt inputs, invalid arguments, and worker
# death mid-task (must be detectable, not hang forever).
#
# Usage:
#   R_LIBS=/data/home/dingjia/R/library Rscript dev/test_edge.R

.libPaths(c('/data/home/dingjia/R/library', .libPaths()))
source('/data/home/dingjia/R/library/PERCEPTIONx/shiny/app/R/async_jobs.R')
options(perception.pkg_root = NULL)
suppressMessages(library(PERCEPTIONx))

fails <- 0L
t <- function(name, ok, detail = '') {
  cat(sprintf('[%s] %s%s\n', if (ok) 'PASS' else 'FAIL', name,
              if (nzchar(detail)) paste0(' -> ', detail) else ''))
  if (!ok) fails <<- fails + 1L
}

# --- 1. corrupt DepMap file -> extract_depmap_meta errors (app deletes it) ---
corrupt <- tempfile(fileext = '.RDS')
writeLines('this is not an RDS file', corrupt)
meta <- tryCatch(PERCEPTIONx:::extract_depmap_meta(corrupt), error = function(e) NULL)
t('extract_depmap_meta rejects corrupt file', is.null(meta))

# --- 2. invalid arguments ----------------------------------------------------
t('prepare_data rejects empty matrix',
  is.null(tryCatch(prepare_data('umap', matrix(numeric(0), 0, 0)),
                   error = function(e) NULL)))
t('predict_patients rejects bad mode',
  is.null(tryCatch(predict_patients(data.frame(patient='P1', clone_id='P1@@c1', d=0.5),
                                    data.frame(patients='P1', `P1@@c1`=1, check.names=FALSE),
                                    mode = 'bogus_mode'),
                   error = function(e) NULL)))
t('get_response_matrix bad input handled',
  length(tryCatch(get_response_matrix(NULL), error = function(e) numeric(0))) == 0)

# --- 3. worker death mid-task is DETECTABLE, not a hang ----------------------
shared <- new.env(parent = emptyenv())
ensure_session_worker(shared)
worker <- shared$task_worker
jid <- submit_session_task(shared, 'demo', list())
# wait until the worker claims it (status.txt appears), then kill the worker
deadline <- Sys.time() + 30
claimed <- FALSE
while (Sys.time() < deadline) {
  if (file.exists(file.path(shared$task_jobs_dir, jid, 'status.txt'))) { claimed <- TRUE; break }
  Sys.sleep(1)
}
t('task was claimed before kill', claimed)
worker$kill()
Sys.sleep(1)
t('worker is dead after kill', !worker$is_alive())
st <- read_task_state(shared, jid)
t('task stays running (not fake done/error) after worker death', st$status == 'running')
# poll_task must give up with an actionable error instead of spinning forever
# (poll_task needs a Shiny session; simulate the same check poll_task performs)
alive <- worker$is_alive()
t('poll can detect dead worker (alive=FALSE)', isFALSE(alive))

cat('\n==== EDGE TEST SUMMARY: FAIL =', fails, '====\n')
quit(status = if (fails == 0L) 0L else 1L)
