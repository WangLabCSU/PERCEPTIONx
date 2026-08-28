# time_demo_worker.R — measures the REAL wall time of a "demo" task through
# the exact async machinery: spawn session_task_main via callr::r_bg (source
# mode -> devtools::load_all), submit the demo job via the file protocol,
# poll until done. This captures process spawn + package load + Seurat/caret.
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
cat(sprintf("worker spawn+startup (returned): %.1fs (process loading in background)\n",
            proc.time()[3] - t_spawn))

# Wait for the worker to be polling (log line written after package load)
deadline <- proc.time()[3] + 120
while (proc.time()[3] < deadline) {
  if (file.exists(file.path(jobs_dir, "worker.log"))) {
    log <- tryCatch(readLines(file.path(jobs_dir, "worker.log"), warn = FALSE), error = function(e) "")
    if (any(grepl("started, polling", log))) break
  }
  Sys.sleep(0.5)
}
cat(sprintf("worker ready (packages loaded, polling): %.1fs\n", proc.time()[3] - t_spawn))

# Submit the demo job (same as submit_session_task)
jobid <- make_job_id()
job_dir <- file.path(jobs_dir, jobid)
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(task = "demo", args = list()), file.path(job_dir, "params.rds"))
file.create(file.path(job_dir, "ready"))
file.create(file.path(job_dir, "progress.txt"))

t_job <- proc.time()[3]
while (TRUE) {
  sf <- file.path(job_dir, "status.txt")
  if (file.exists(sf)) {
    s <- readLines(sf, warn = FALSE)
    if (length(s) > 0 && s[1] %in% c("done", "error", "running")) {
      if (s[1] == "done") break
      if (startsWith(s[1], "error")) { cat("JOB ERROR:", s, "\n"); break }
    }
  }
  if (proc.time()[3] - t_job > 300) { cat("TIMEOUT after 300s\n"); break }
  Sys.sleep(0.5)
}
cat(sprintf("demo job wall time: %.1fs\n", proc.time()[3] - t_job))
cat(sprintf("TOTAL (spawn->demo done): %.1fs\n", proc.time()[3] - t_spawn))
if (file.exists(file.path(jobs_dir, "worker.log"))) {
  cat("--- worker.log tail ---\n")
  cat(tail(readLines(file.path(jobs_dir, "worker.log"), warn = FALSE), 5), sep = "\n")
}
worker$kill()
unlink(jobs_dir, recursive = TRUE)
