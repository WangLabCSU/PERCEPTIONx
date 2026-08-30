# ---------------------------------------------------------------------------
# Unified cache / data directory resolution
#
# All on-disk caches (DepMap, pre-trained models, training job dirs) can be
# pinned under ONE root directory, so a deployment can keep every file in a
# chosen location instead of the per-user defaults:
#
#   options(PERCEPTIONx.cache_root = "/data/perceptionx")
#   # or the environment variable (both spellings accepted):
#   #   PERCEPTIONX_CACHE_ROOT=/data/perceptionx
#
# Resulting layout:
#   <root>/depmap/  DepMap.RDS + DepMap_meta.RDS + DepMap_used.flag
#   <root>/models/  pre-trained model .RDS files
#   <root>/jobs/    background training job dirs (shared master pool)
#
# Precedence for the DepMap directory specifically (kept from the original
# app-only option): the dedicated option/env var wins; otherwise it is derived
# from cache_root; otherwise the R user data dir is used.
# When no cache_root is set, package-level functions keep their historical
# defaults ("./" for load_depmap, "./models" for load_model) and the Shiny app
# keeps using tempdir() for models / jobs — i.e. behaviour is unchanged.
# ---------------------------------------------------------------------------

# Resolve the unified cache root (options > env var). NULL when unset.
perception_cache_root <- function() {
  root <- getOption("PERCEPTIONx.cache_root")
  if (is.null(root) || !nzchar(root)) {
    root <- Sys.getenv("PERCEPTIONX_CACHE_ROOT",
                       Sys.getenv("PERCEPTIONx_CACHE_ROOT", ""))
  }
  if (!is.null(root) && nzchar(root)) root else NULL
}

# DepMap cache directory used by the Shiny app.
perception_depmap_cache_dir <- function() {
  d <- getOption("PERCEPTIONx.depmap_cache_dir")
  if (is.null(d) || !nzchar(d)) {
    d <- Sys.getenv("PERCEPTIONX_DEPMAP_CACHE_DIR",
                    Sys.getenv("PERCEPTIONx_DEPMAP_CACHE_DIR", ""))
  }
  if (!is.null(d) && nzchar(d)) return(d)
  root <- perception_cache_root()
  if (!is.null(root)) return(file.path(root, "depmap"))
  tools::R_user_dir("PERCEPTIONx", "data")
}

# Pre-trained model cache directory when a cache_root is set (NULL otherwise).
perception_model_dir <- function() {
  root <- perception_cache_root()
  if (is.null(root)) NULL else file.path(root, "models")
}

# Default dest for load_depmap(): cache_root-derived, else historical ".".
perception_default_depmap_dir <- function() {
  root <- perception_cache_root()
  if (is.null(root)) "." else file.path(root, "depmap")
}

# Default dest for load_model(): cache_root-derived, else historical "./models".
perception_default_model_dir <- function() {
  root <- perception_cache_root()
  if (is.null(root)) "./models" else file.path(root, "models")
}

# Shared master jobs directory: cache_root-derived, else tempdir()/perception_jobs.
perception_jobs_dir <- function() {
  root <- perception_cache_root()
  if (is.null(root)) file.path(tempdir(), "perception_jobs") else file.path(root, "jobs")
}
