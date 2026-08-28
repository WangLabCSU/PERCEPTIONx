# End-to-end test of the per-session task worker (demo / prepare / predict /
# error propagation) via the app's own code path: ensure_session_worker ->
# submit_session_task -> read_task_state. Runs on the server.
#
# Usage:
#   R_LIBS=/data/home/dingjia/R/library Rscript dev/e2e_session.R

.libPaths(c('/data/home/dingjia/R/library', .libPaths()))
source('/data/home/dingjia/R/library/PERCEPTIONx/shiny/app/R/async_jobs.R')
options(perception.pkg_root = NULL)

# environment = reference semantics (reactiveValues stand-in, see e2e_train.R)
shared <- new.env(parent = emptyenv())

fails <- 0L
t <- function(name, ok, detail = '') {
  cat(sprintf('[%s] %s%s\n', if (ok) 'PASS' else 'FAIL', name,
              if (nzchar(detail)) paste0(' -> ', detail) else ''))
  if (!ok) fails <<- fails + 1L
}

ensure_session_worker(shared)
worker <- shared$task_worker
jobs_dir <- shared$task_jobs_dir
cat('[e2e-session] worker spawned, jobs_dir =', jobs_dir, '\n')
t('session worker alive', !is.null(worker) && worker$is_alive())

wait_done <- function(jobid, timeout = 600) {
  deadline <- Sys.time() + timeout
  while (Sys.time() < deadline) {
    if (!worker$is_alive()) return(list(status = 'worker_died'))
    st <- read_task_state(shared, jobid)
    if (st$status %in% c('done', 'error')) return(st)
    Sys.sleep(2)
  }
  list(status = 'timeout')
}

# --- 1. demo task (full pipeline: synthetic data + Seurat + 2 models) --------
jid <- submit_session_task(shared, 'demo', list())
cat('[e2e-session] demo submitted:', jid, '\n')
r <- wait_done(jid)
t('demo task done', r$status == 'done', if (r$status == 'error') r$message else '')
if (r$status == 'done') {
  res <- r$result
  t('demo result has models', !is.null(res$models) && length(res$models) >= 2)
  t('demo result has prepared_data', !is.null(res$prepared_data))
  t('demo result has expr/mapping/response',
    !is.null(res$user_expr) && !is.null(res$user_mapping) && !is.null(res$user_response))
}

# --- 2. prepare task (user data through Seurat in worker) --------------------
set.seed(1)
gene_names <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "RB1",
                "APC", "BRAF", "CDH1", "CDKN2A", "ERBB2", "FGFR1", "ALK",
                "MET", "RET", "ROS1", "NRAS", "HRAS", "MAP2K1", "MAPK1",
                "JAK2", "STAT3", "MTOR", "AKT1", "AKT2", "CTNNB1", "SMAD4",
                "VHL", "NF1", "NF2", "STK11", "FBXW7", "ARID1A", "KDM5C",
                "KMT2D", "SETD2", "BAP1", "PBRM1", "NOTCH1", "NOTCH2",
                "DICER1", "TET2", "IDH1", "IDH2", "DNMT3A", "FLT3", "NPM1")
cell_names <- paste0('C', seq_len(400))
expr49 <- matrix(rpois(49 * 400, lambda = 5), 49, 400,
                 dimnames = list(gene_names, cell_names))
pm49 <- data.frame(cell_id = cell_names,
                   patient_id = sample(paste0('PT', 1:20), 400, replace = TRUE),
                   stringsAsFactors = FALSE)
jid <- submit_session_task(shared, 'prepare',
  list(method = 'umap', expression_matrix = expr49, patient_mapping = pm49,
       seurat_resolution = 0.8, seurat_dims = 10, seurat_nfeatures = 50))
cat('[e2e-session] prepare submitted:', jid, '\n')
r <- wait_done(jid)
t('prepare task done', r$status == 'done', if (r$status == 'error') r$message else '')
prepared <- if (r$status == 'done') r$result else NULL
if (!is.null(prepared)) {
  t('prepare returns clone map (400 cells)',
    nrow(prepared$cell_clone_map) == 400)
}

# --- 3. predict task (clone-level + patient-level) ---------------------------
# Same shape as mod_predict.R: expr = clone_expression_rnorm, patient_input =
# list(clone_viability_template, clone_counts).
jid <- submit_session_task(shared, 'predict',
  list(model_list = list(erlotinib = res$models$erlotinib),
       expr = prepared$clone_expression_rnorm,
       patient_input = list(
         clone_viability_template = prepared$clone_viability_template,
         clone_counts = prepared$clone_counts),
       mode = 'weighted_max'))
cat('[e2e-session] predict submitted:', jid, '\n')
r <- wait_done(jid)
t('predict task done', r$status == 'done', if (r$status == 'error') r$message else '')
if (r$status == 'done') {
  pr <- r$result
  t('predict returns clone_pred matrix', !is.null(pr$clone_pred) && is.matrix(pr$clone_pred))
  t('predict returns patient_pred', !is.null(pr$patient_pred))
}

# --- 4. unknown task type -> clean error propagation -------------------------
jid <- submit_session_task(shared, 'no_such_task', list())
cat('[e2e-session] bad task submitted:', jid, '\n')
r <- wait_done(jid)
t('unknown task errors cleanly', r$status == 'error' && nzchar(r$message), r$message)

cat('\n==== SESSION WORKER TEST SUMMARY: FAIL =', fails, '====\n')
quit(status = if (fails == 0L) 0L else 1L)
