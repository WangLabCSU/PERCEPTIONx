# Download the official DepMap release asset into a FIXED cache directory so
# that the Shiny app (and e2e tests) can reuse it across restarts instead of
# re-downloading 567MB each time.
#
# Usage (server):
#   R_LIBS=/data/home/dingjia/R/library nohup Rscript dev/download_depmap.R \
#     /data/home/dingjia/DepMap_cache > download_depmap.log 2>&1 &
#
# The app reads this directory via the same load_depmap() call — the function
# skips the download when the file already exists.

.libPaths(c('/data/home/dingjia/R/library', .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
dest <- if (length(args) >= 1) args[1] else file.path(Sys.getenv('HOME'), 'DepMap_cache')
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(PERCEPTIONx))
cat('R:', R.version.string, '\n')
cat('PERCEPTIONx version:', as.character(packageVersion('PERCEPTIONx')), '\n')
cat('dest:', dest, '\n')

# mirror = TRUE: direct GitHub is ~0 on this server, the gh-proxy mirrors are
# usable (~50KB/s). The app's "Download & Load" button passes the same flag
# when the user enables 'Use mirror'.
PERCEPTIONx::load_depmap(dest = dest, read = FALSE, mirror = TRUE,
                         timeout_seconds = 600, retries = 2)

f <- file.path(dest, 'DepMap.RDS')
if (file.exists(f)) {
  cat('DONE:', f, round(file.info(f)$size / 1e6), 'MB\n')
} else {
  cat('FAILED: no DepMap.RDS produced\n')
  quit(status = 1)
}
