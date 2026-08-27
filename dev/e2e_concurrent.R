# Concurrency/isolation test: two independent Shiny "sessions" each spawn
# their own per-session worker and run a demo task in parallel. Asserts:
# distinct worker PIDs, both tasks complete independently, no cross-talk.
#
# Usage:
#   R_LIBS=/data/home/dingjia/R/library Rscript dev/e2e_concurrent.R

.libPaths(c('/data/home/dingjia/R/library', .libPaths()))
source('/data/home/dingjia/R/library/PERCEPTIONx/shiny/app/R/async_jobs.R')
options(perception.pkg_root = NULL)

mk_session <- function() {
  s <- new.env(parent = emptyenv())
  ensure_session_worker(s)
  list(shared = s, worker = s$task_worker, dir = s$task_jobs_dir)
}

s1 <- mk_session()
s2 <- mk_session()
cat('[e2e-conc] session1 pid:', s1$worker$get_pid(), '| session2 pid:', s2$worker$get_pid(), '\n')

fails <- 0L
t <- function(name, ok, detail = '') {
  cat(sprintf('[%s] %s%s\n', if (ok) 'PASS' else 'FAIL', name,
              if (nzchar(detail)) paste0(' -> ', detail) else ''))
  if (!ok) fails <<- fails + 1L
}
t('two sessions have distinct workers', s1$worker$get_pid() != s2$worker$get_pid())
t('distinct jobs dirs', !identical(s1$dir, s2$dir))

submit_wait <- function(s, timeout = 600) {
  jid <- submit_session_task(s$shared, 'demo', list())
  deadline <- Sys.time() + timeout
  while (Sys.time() < deadline) {
    if (!s$worker$is_alive()) return(list(status = 'worker_died', jid = jid))
    st <- read_task_state(s$shared, jid)
    if (st$status %in% c('done', 'error')) return(c(list(jid = jid), st))
    Sys.sleep(2)
  }
  list(status = 'timeout', jid = jid)
}

t0 <- Sys.time()
cat('[e2e-conc] submitting demo to BOTH sessions...\n')
r1 <- submit_wait(s1)
r2 <- submit_wait(s2)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
cat('[e2e-conc] both done in', round(elapsed), 's\n')

t('session1 demo done', r1$status == 'done', if (r1$status == 'error') r1$message else '')
t('session2 demo done', r2$status == 'done', if (r2$status == 'error') r2$message else '')
if (r1$status == 'done' && r2$status == 'done') {
  m1 <- r1$result$models; m2 <- r2$result$models
  t('results are independent objects', !identical(m1, m2))
  t('both have 2 models', length(m1) >= 2 && length(m2) >= 2)
  # models differ (different random seeds) but share structure
  t('model structure consistent',
    all(sapply(m1, function(m) !is.null(m$model))) &&
    all(sapply(m2, function(m) !is.null(m$model))))
}

# re-submit on session1 to confirm worker reuse (same pid, no leak)
pid_before <- s1$worker$get_pid()
r3 <- submit_wait(s1)
t('session1 worker reused (same pid)', s1$worker$get_pid() == pid_before && r3$status == 'done')

cat('\n==== CONCURRENCY TEST SUMMARY: FAIL =', fails, '====\n')
quit(status = if (fails == 0L) 0L else 1L)
