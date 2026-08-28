# Server environment self-check for PERCEPTION-shiny.
# Run BEFORE starting the app or e2e tests:
#   R_LIBS=/data/home/dingjia/R/library Rscript dev/check_env.R
# Exits non-zero if anything critical fails. Prints a summary either way.

.libPaths(c('/data/home/dingjia/R/library', .libPaths()))

status <- function(label, ok, detail = '') {
  cat(sprintf('[%s] %s%s\n', if (ok) ' OK ' else 'FAIL', label,
              if (nzchar(detail)) paste0(' -> ', detail) else ''))
  ok
}

fails <- 0L

# --- R / library paths ------------------------------------------------------
cat('R version      :', R.version.string, '\n')
cat('.libPaths      :', paste(.libPaths(), collapse = ' | '), '\n')

# --- Installed packages -----------------------------------------------------
pkgs <- c('PERCEPTIONx', 'shiny', 'bslib', 'DT', 'callr', 'waiter', 'shinycssloaders')
for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf('%-14s: %s\n', p, as.character(packageVersion(p))))
  } else {
    cat(sprintf('%-14s: NOT INSTALLED\n', p))
  }
}

# --- PERCEPTIONx exports required by the app --------------------------------
# extract_depmap_meta() is an internal (@noRd) package function used by the
# app's background "extract_meta" task; check it exists in the namespace so a
# stale installed package is caught here instead of failing in the worker.
required_exports <- c('load_depmap', 'train_models', 'prepare_data',
                      'predict_drugs', 'predict_patients')
if (requireNamespace('PERCEPTIONx', quietly = TRUE)) {
  exp <- getNamespaceExports('PERCEPTIONx')
  for (f in required_exports) {
    ok <- f %in% exp
    if (!ok) fails <- fails + 1L
    status(sprintf('export %s()', f), ok)
  }
  ok_internal <- exists('extract_depmap_meta', envir = asNamespace('PERCEPTIONx'))
  if (!ok_internal) fails <- fails + 1L
  status('internal extract_depmap_meta()', ok_internal)
} else {
  fails <- fails + length(required_exports)
  status('PERCEPTIONx installed', FALSE)
}

# App-internal helpers the app sources at startup.
app_r <- '/data/home/dingjia/R/library/PERCEPTIONx/shiny/app/R'
for (f in c('shiny_helpers.R', 'async_jobs.R', 'mod_data.R', 'mod_train.R',
            'mod_predict.R', 'mod_visualize.R', 'mod_home.R', 'mod_help.R')) {
  ok <- file.exists(file.path(app_r, f))
  if (!ok) fails <- fails + 1L
  status(sprintf('app source %s', f), ok)
}

# --- DepMap cache -----------------------------------------------------------
cache_dir <- getOption(
  'PERCEPTIONx.depmap_cache_dir',
  # Accept both env-var spellings (docs use PERCEPTIONX..., early
  # deployments set PERCEPTIONx...); Sys.getenv is case-sensitive on Linux.
  Sys.getenv('PERCEPTIONX_DEPMAP_CACHE_DIR',
             Sys.getenv('PERCEPTIONx_DEPMAP_CACHE_DIR',
                        tools::R_user_dir('PERCEPTIONx', 'data')))
)
cat('DepMap cache dir :', cache_dir, '\n')
f <- file.path(cache_dir, 'DepMap.RDS')
if (file.exists(f)) {
  sz <- round(file.info(f)$size / 1e6)
  ok <- file.info(f)$size > 5e8  # > ~500MB is a plausible complete file
  if (!ok) fails <- fails + 1L
  status(sprintf('DepMap.RDS cache (%d MB)', sz), ok,
         if (!ok) 'suspiciously small — likely incomplete download' else '')
} else {
  fails <- fails + 1L
  status('DepMap.RDS cache', FALSE,
         'not downloaded yet — run dev/download_depmap.R <dir> first')
}

# --- Mirror reachability / speed (quick probe: first 8MB, max 20s) ----------
cat('\nMirror download speed probe (first 8MB):\n')
urls <- c(
  'https://gh-proxy.com/https://github.com',
  'https://v6.gh-proxy.org/https://github.com',
  'https://github.akams.cn/https://github.com'
)
for (u in urls) {
  sp <- tryCatch(system2('curl', c('-L', '-s', '--max-time', '20',
                                   '-r', '0-8000000', '-o', '/dev/null',
                                   '-w', '%{speed_download}',
                                   file.path(u, 'WangLabCSU/PERCEPTIONx/releases/download/depmap/DepMap.RDS')),
                        stdout = TRUE, stderr = TRUE),
                 error = function(e) 'ERR')
  kb <- suppressWarnings(as.numeric(sp[1]))
  cat(sprintf('  %-55s: %s KB/s\n', u,
              if (is.finite(kb)) round(kb / 1024) else 'unreachable'))
}

cat('\nResult:', if (fails == 0L) 'ALL CHECKS PASSED' else paste(fails, 'FAILURE(S)'), '\n')
quit(status = if (fails == 0L) 0L else 1L)
