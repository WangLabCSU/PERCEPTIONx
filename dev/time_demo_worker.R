# time_demo_worker.R — measures the REAL wall time of a "demo" task through
# the exact async machinery, INCLUDING the new warm-up task: spawn the worker
# (source mode -> devtools::load_all), submit "warmup" (pre-loads Seurat/
# caret/glmnet), then submit the demo job and time it. Demonstrates that the
# ~3-4 s package load no longer sits inside the demo task.
# Usage: Rscript dev/time_demo_worker.R
suppressMessages({
  source("inst/shiny/app/R/async_jobs.R", local = TRUE)
  library(callr)
})

jobs_dir <- file.path(tempdir(), "perception_tasks", make_job_id())
dir.create(jobs_dir, recursive = TRUE, showWarnings = FALSE)

t_spawn <- proc.time()[3]
worker <- callr::r_bg(
  func = session_task_main,
  args = list(pkg_root = "c:/Users/Lenovo/Desktop/PERCEPTION", jobs_dir = jobs_dir),
  supervise = TRUE, poll_connection = FALSE,
  stdout = FALSE, stderr = file.path(jobs_dir, "worker.log"))
cat(sprintf("worker spawned (returned): %.1fs\n", proc.time()[3] - t_spawn))

deadline <- proc.time()[3] + 120
while (proc.time()[3] < deadline) {
  if (file.exists(file.path(jobs_dir, "worker.log"))) {
    log <- tryCatch(readLines(file.path(jobs_dir, "worker.log"), warn = FALSE), error = function(e) "")
    if (any(grepl("started, polling", log))) break
  }
  Sys.sleep(0.5)
}
cat(sprintf("worker ready (PERCEPTIONx loaded, polling): %.1fs\n", proc.time()[3] - t_spawn))

# Helper: submit a task and wait for done (mirrors app behaviour)
submit_wait <- function(task, args) {
  jobid <- make_job_id()
  job_dir <- file.path(jobs_dir, jobid)
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(task = task, args = args), file.path(job_dir, "params.rds"))
  file.create(file.path(job_dir, "ready"))
  file.create(file.path(job_dir, "progress.txt"))
  t0 <- proc.time()[3]
  repeat {
    sf <- file.path(job_dir, "status.txt")
    if (file.exists(sf)) {
      s <- readLines(sf, warn = FALSE)
      if (length(s) > 0 && s[1] == "done") break
      if (length(s) > 0 && startsWith(s[1], "error")) stop("task error: ", s[1])
    }
    if (proc.time()[3] - t0 > 300) stop("timeout")
    Sys.sleep(0.2)
  }
  proc.time()[3] - t0
}

t_w <- submit_wait("warmup", list())
cat(sprintf("warmup task (Seurat/caret load): %.1fs\n", t_w))
t_d <- submit_wait("demo", list())
cat(sprintf("demo task AFTER warmup: %.1fs\n", t_d))
cat(sprintf("TOTAL (spawn -> demo done): %.1fs\n", proc.time()[3] - t_spawn))
worker$kill()
unlink(jobs_dir, recursive = TRUE)
