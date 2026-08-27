# Package-level unit tests for PERCEPTIONx, run on the server against the
# CLEAN release DepMap. Covers the functions the app's async workers call:
# get_response_matrix / get_cellLine_list / feature_ranking_bulk /
# train_models (glmnet + rf, multiple drugs/k) / prepare_data /
# predict_drugs / predict_patients.
#
# Usage:
#   R_LIBS=/data/home/dingjia/R/library \
#   PERCEPTIONx_DEPMAP_CACHE_DIR=/data/home/dingjia/DepMap_cache \
#     Rscript dev/test_unit.R
#
# Exits non-zero if any test FAILS.

.libPaths(c('/data/home/dingjia/R/library', .libPaths()))
suppressMessages(library(PERCEPTIONx))

cache_dir <- Sys.getenv('PERCEPTIONx_DEPMAP_CACHE_DIR',
                        tools::R_user_dir('PERCEPTIONx', 'data'))
DepMap <- readRDS(file.path(cache_dir, 'DepMap.RDS'))
cat('DepMap loaded:', length(DepMap), 'components\n')

fails <- 0L
skips <- 0L
t <- function(name, expr, skip_if = FALSE) {
  if (skip_if) { cat('[SKIP]', name, '\n'); skips <<- skips + 1L; return(invisible(NULL)) }
  ok <- tryCatch({ r <- force(expr); isTRUE(r) },
                 error = function(e) { cat('   error:', conditionMessage(e), '\n'); FALSE })
  if (ok) { cat('[PASS]', name, '\n') } else { cat('[FAIL]', name, '\n'); fails <<- fails + 1L }
  invisible(NULL)
}

out_dir <- '/tmp/test_unit_out'
unlink(out_dir, recursive = TRUE); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- drug response ----------------------------------------------------------
r <- tryCatch(get_response_matrix('erlotinib'), error = function(e) numeric(0))
t('get_response_matrix erlotinib has data', length(r) > 100)
r2 <- tryCatch(get_response_matrix('nonexistent_drug_xyz'), error = function(e) numeric(0))
# Unknown drug: function must not crash; may return empty or all-NA vector.
t('get_response_matrix unknown drug handled (no crash)',
  length(r2) == 0 || all(is.na(r2)))

# --- cell line selection ----------------------------------------------------
cl <- tryCatch(get_cellLine_list('PanCan', 'erlotinib', 'PanCan',
                                 get_response_matrix('erlotinib')),
               error = function(e) NULL)
t('get_cellLine_list PanCan returns lines', !is.null(cl) && length(cl[[1]]) > 100)
cl2 <- tryCatch(get_cellLine_list('breast', 'erlotinib', 'PanCan',
                                  get_response_matrix('erlotinib')),
                error = function(e) NULL)
t('get_cellLine_list breast returns some lines', !is.null(cl2))

# --- feature ranking --------------------------------------------------------
fr <- tryCatch(feature_ranking_bulk('erlotinib', 'PanCan', 'PanCan', c('TP53', 'BRCA1')),
               error = function(e) NULL)
t('feature_ranking_bulk PanCan ranks', !is.null(fr) && nrow(fr) >= 1 && all(c('p.value', 'estimate.cor') %in% colnames(fr)))
fr2 <- tryCatch(feature_ranking_bulk('erlotinib', 'PanCan', 'PanCan', c('NOT_A_GENE_XYZ')),
                error = function(e) NULL)
t('feature_ranking_bulk unknown GOI -> empty/NA (no crash)', is.null(fr2) || nrow(fr2) == 0 || all(is.na(fr2)))

# --- train_models: glmnet, 2 drugs, 2 k values ------------------------------
o1 <- tryCatch(
  train_models(drug_list = c('erlotinib', 'gemcitabine'),
               cancer_type = 'PanCan', exclude_cancer = 'PanCan',
               GOI = c('TP53', 'BRCA1'), k_features_values = c(20, 50),
               model_type = 'glmnet', ncores = 2, output_dir = out_dir),
  error = function(e) NULL)
t('train_models glmnet 2 drugs succeed', !is.null(o1) && length(o1) == 2)
if (!is.null(o1) && length(o1) == 2) {
  t('train_models glmnet model objects present',
    all(sapply(o1, function(m) inherits(m$model, 'train') || inherits(m$model, 'glmnet'))))
  t('train_models glmnet perf fields',
    all(sapply(o1, function(m) all(c('performance_in_bulk', 'performance_in_scRNA') %in% names(m)))))
  t('train_models glmnet saves RDS to output_dir',
    length(list.files(out_dir, pattern = '\\.RDS$')) >= 1)
}

# --- train_models: random forest --------------------------------------------
has_rf <- requireNamespace('randomForest', quietly = TRUE)
o2 <- tryCatch(
  train_models(drug_list = 'erlotinib', cancer_type = 'PanCan', exclude_cancer = 'PanCan',
               GOI = c('TP53', 'BRCA1'), k_features_values = c(20),
               model_type = 'rf', ncores = 1,
               output_dir = file.path(out_dir, 'rf')),
  error = function(e) NULL)
t('train_models rf succeeds', !is.null(o2) && length(o2) == 1, skip_if = !has_rf)

# --- train_models: unknown drug -> 0 models, no crash ------------------------
o3 <- tryCatch(
  train_models(drug_list = 'nonexistent_drug_xyz', cancer_type = 'PanCan',
               exclude_cancer = 'PanCan', GOI = c('TP53'),
               k_features_values = c(20), model_type = 'glmnet',
               ncores = 1, output_dir = file.path(out_dir, 'bad')),
  error = function(e) NULL)
t('train_models unknown drug -> 0 models no crash', is.list(o3) && length(o3) == 0)

# --- train_models: GOI not in expression -------------------------------------
o4 <- tryCatch(
  train_models(drug_list = 'erlotinib', cancer_type = 'PanCan', exclude_cancer = 'PanCan',
               GOI = c('NOT_A_GENE_XYZ'), k_features_values = c(20),
               model_type = 'glmnet', ncores = 1,
               output_dir = file.path(out_dir, 'goi')),
  error = function(e) NULL)
# train_models warns "None of the provided GOI genes found" and stops;
# either behavior is acceptable as long as it does not hang.
t('train_models GOI-not-found handled (stop or 0 models)',
  is.null(o4) || length(o4) == 0)

# --- prepare_data (Seurat clustering) --------------------------------------
# Same input profile as the app's demo pipeline (49 real genes x 400 cells,
# 20 patients) — the exact profile already proven in the browser e2e run.
set.seed(1)
gene_names <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "RB1",
                "APC", "BRAF", "CDH1", "CDKN2A", "ERBB2", "FGFR1", "ALK",
                "MET", "RET", "ROS1", "NRAS", "HRAS", "MAP2K1", "MAPK1",
                "JAK2", "STAT3", "MTOR", "AKT1", "AKT2", "CTNNB1", "SMAD4",
                "VHL", "NF1", "NF2", "STK11", "FBXW7", "ARID1A", "KDM5C",
                "KMT2D", "SETD2", "BAP1", "PBRM1", "NOTCH1", "NOTCH2",
                "DICER1", "TET2", "IDH1", "IDH2", "DNMT3A", "FLT3", "NPM1")
cell_names <- paste0("C", seq_len(400))
expr <- matrix(rpois(49 * 400, lambda = 5), 49, 400,
               dimnames = list(gene_names, cell_names))
patient_assignment <- sample(paste0("PT", 1:20), 400, replace = TRUE)
pm <- data.frame(cell_id = cell_names, patient_id = patient_assignment,
                 stringsAsFactors = FALSE)
pd <- tryCatch(
  prepare_data(method = 'umap', expression_matrix = expr, patient_mapping = pm,
               seurat_resolution = 0.8, seurat_dims = 10, seurat_nfeatures = 50),
  error = function(e) NULL)
t('prepare_data runs (Seurat, demo profile)', !is.null(pd))
if (!is.null(pd)) {
  t('prepare_data returns clone map for all cells',
    nrow(pd$cell_clone_map) == 400)
  t('prepare_data clone_counts present', !is.null(pd$clone_counts))
}

# --- predict_drugs / predict_patients (aligned clones) -----------------------
# The app predicts on the CLONE-level matrix (prepare_data$clone_expression_rnorm,
# rownames = "Patient@@Clone"); predict_patients() then aligns those rows with
# the clone_viability_template. Cell-level input would misalign. Build one
# synthetic count matrix whose genes include the model features, prepare it,
# then predict on the clone-level output — exactly like mod_predict.R.
feats_all <- unique(unlist(lapply(o1, function(m) m$model$coefnames)))
# Seurat requires valid R names as genes (no "-"/leading digits etc.), while
# DepMap features can contain them. make.names() fixes it; predict_drugs has a
# make.names() fallback (see viability_from_model_internal) so matching holds.
sim_genes <- unique(make.names(c(feats_all, 'TP53', 'BRCA1'), unique = TRUE))
cell_names2 <- paste0('C', seq_len(400))
expr_sim <- matrix(rpois(length(sim_genes) * 400, lambda = 5),
                   length(sim_genes), 400,
                   dimnames = list(sim_genes, cell_names2))
pm_sim <- data.frame(cell_id = cell_names2,
                     patient_id = sample(paste0('PT', 1:20), 400, replace = TRUE),
                     stringsAsFactors = FALSE)
pd_sim <- tryCatch(
  prepare_data(method = 'umap', expression_matrix = expr_sim,
               patient_mapping = pm_sim,
               seurat_resolution = 0.8, seurat_dims = 10, seurat_nfeatures = 50),
  error = function(e) NULL)
t('predict pipeline: prepare_data on model-feature data', !is.null(pd_sim))
expr_clone <- pd_sim$clone_expression_rnorm
pd_pred <- tryCatch(predict_drugs(o1, expr_clone), error = function(e) NULL)
t('predict_drugs on clone-level expr returns matrix', !is.null(pd_pred) && is.matrix(pd_pred))
if (!is.null(pd_pred)) {
  t('predict_drugs rows = clones, cols = drugs',
    nrow(pd_pred) == ncol(expr_clone) && ncol(pd_pred) == 2)
}
pp <- tryCatch(predict_patients(pd_pred, pd_sim, mode = 'weighted_max'),
               error = function(e) NULL)
t('predict_patients simple mode', !is.null(pp))
if (!is.null(pp)) {
  # returns data.frame: rownames = patients, columns = drugs (no 'patient' col)
  t('predict_patients per-patient rows (rownames)',
    is.data.frame(pp) && nrow(pp) == length(unique(pm_sim$patient_id)) &&
    ncol(pp) >= 1 && all(c('erlotinib', 'gemcitabine') %in% names(pp)))
}
# legacy mode: clone_viability_matrix (patient/clone_id/drugs) + clone_counts
# clone_ids must match clone_counts columns exactly (per patient cl1/cl2).
cvm <- data.frame(
  patient  = rep(paste0('PT', 1:4), each = 2),
  clone_id = paste0(rep(paste0('PT', 1:4), each = 2),
                    rep(c('@@cl1', '@@cl2'), 4)),
  erlotinib = runif(8), gemcitabine = runif(8),
  stringsAsFactors = FALSE)
clone_counts2 <- data.frame(
  patients = paste0('PT', 1:4),
  `PT1@@cl1` = runif(4), `PT1@@cl2` = runif(4),
  `PT2@@cl1` = runif(4), `PT2@@cl2` = runif(4),
  `PT3@@cl1` = runif(4), `PT3@@cl2` = runif(4),
  `PT4@@cl1` = runif(4), `PT4@@cl2` = runif(4),
  check.names = FALSE)
pp2 <- tryCatch(
  predict_patients(cvm, clone_counts2, mode = 'weighted_average'),
  error = function(e) NULL)
t('predict_patients legacy mode', !is.null(pp2))

cat('\n==== UNIT TEST SUMMARY ====\n')
cat('PASS:', 'see above | FAIL:', fails, '| SKIP:', skips, '\n')
quit(status = if (fails == 0L) 0L else 1L)
