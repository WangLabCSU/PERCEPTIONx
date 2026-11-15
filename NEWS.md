# PERCEPTIONx 0.4.0

## New features

- `preprocess_sc()` implements the full raw-counts to clone-label workflow:
  MAD-based per-patient quality control, normalisation, variable-feature
  selection, a shared PCA embedding, and clustering. `cluster_scope = "global"`
  (default) clusters all cells together, as in the original PERCEPTION method,
  so a clone id can occur in several patients; `cluster_scope = "patient"`
  clusters inside each patient only, so a clone can never span two patients.
  The two are identical for a single-patient dataset. With the patient scope,
  harmony is applied inside a patient that has more than one sample
  (e.g. longitudinal biopsies) so a batch effect cannot split one clone in two;
  `harmony` is optional and the step is skipped with a warning when absent.
- `preprocess_sc()` gains `skip_qc = FALSE`. With `skip_qc = TRUE` the MAD
  fences are not applied and every cell is kept, which reproduces the original
  PERCEPTION pipeline on an already-curated matrix; the QC table still reports
  the input but marks the fences as `off`.
- `plot_qc_report()` draws the per-patient QC fences, so the automatic
  filtering can be audited without knowing the MAD rule.
- `prepare_data()` gains two alternative inputs. `data = <Seurat object>` reads
  the matrix from `assay`/`layer` and labels from `clone_col` or `Idents()`;
  `clone_labels = <data.frame>` takes a pre-computed cell-to-clone table. Both
  skip clustering.
- `prepare_data()` now records where the clone labels came from in
  `$label_source`: `auto_cluster`, `user_labels`, `seurat_idents`,
  `seurat_metadata`, or `provided_matrix_columns`.
- `preprocess_sc()` gains `remove_doublets = FALSE`: when enabled, cells
  classified as doublets by `scDblFinder` are removed after QC, one patient at
  a time. Doublets pass every QC fence (normal gene counts and mitochondrial
  fraction) yet can form a spurious clone made of two fused clones.
  `scDblFinder` and `SingleCellExperiment` are optional (`Suggests`); without
  them the step warns and is skipped. Removed cells are returned in
  `$doublet_cells` and counted in `qc$n_doublet_removed`.
- `preprocess_sc()` gains `exclude_markers` / `exclude_min_expr`: cells that
  strongly express a user-supplied marker set are removed before clustering, so
  immune and stromal populations do not form clusters that later read as
  resistance clones. This is an exclusion rule, not a tumour-cell classifier.
  Default `NULL` keeps every cell. Removed cells are returned in
  `$excluded_cells`.
- `confidence` in `preprocess_sc()` output is renamed `cluster_size_class`
  (`large` / `small`), because it only describes cluster size and says nothing
  about how credible a clone is.

## Improvements

- The MAD fence default is now `nmads = 5` (was 3). Benchmarked against
  CellRanger's own cell calls on real unfiltered matrices: PBMC 3k Jaccard
  0.914 -> 0.969, lung squamous carcinoma 0.711 -> 0.806, precision unchanged
  in both. On the labelled Maynard cohort the patient-level AUC returns to the
  no-QC value (0.821) instead of dropping to 0.750, i.e. QC stops deleting
  real tumour cells.
- Mitochondrial genes are now detected when row names are Ensembl IDs
  (10x GRCh38 `features.tsv` column 1, many GEO uploads). Previously the
  `MT-` prefix never matched, so `pct_mito` was silently 0 and the
  mitochondrial fence did nothing.
- Internal code migrated to `data.table`; user-facing results remain
  `data.frame`.
- Messages and errors use `cli`; `run_seurat_pipeline()` and its
  `lifecycle` dependency were removed in favour of `prepare_data()` and the
  shared `cluster_embedding()` core.
- Drug catalogues are exported as documented objects:
  `perception_depmap_drugs`, `perception_fda_approved_drugs`,
  `perception_all_drugs`.
- `glmnet` and `randomForest` remain in `Imports`: `caret::train(method =
  "glmnet")` is the default path and loads them at runtime.
- `Matrix` added to `Imports`; `harmony` added to `Suggests`.
- Added `air` formatting, `jarl` linting, spelling checks, and repository
  templates (CONTRIBUTING / SUPPORT / CODE_OF_CONDUCT).

# PERCEPTIONx 0.3.0

## New features
- Help tab gained a "Supported Drugs" catalogue: the full ~1,450-compound
  DepMap screen with the 44 pre-trained models marked, filterable by
  Pre-trained / Demo / All compounds with debounced search.
- New guided walkthroughs on the Help tab (Demo and Advanced tours) with
  retaken, correctly-sized screenshots.
- Reworked the Advanced tour steps: three input files, Seurat clustering,
  pre-trained vs. trained models, training output, prediction and
  visualisation reading order.

## Improvements
- README, vignettes and tutorial sections rearranged and synced with the app.

# PERCEPTIONx 0.2.0

## New features
- Redesigned the visualization tab as a "plot gallery": all 7 plot types
  visible as clickable cards, a data-readiness bar showing what is loaded,
  and click-to-draw (no Generate button). Parameter changes (drug / gene /
  ROC pair) redraw the plot automatically; hovering a card previews its
  description before you commit; generated-at timestamps are shown in UTC.
- Seurat clustering now reports live stage progress (normalizing / variable
  features / PCA / clustering / UMAP / ...) in the Load Demo overlay.
- Demo data is more realistic and reproducible: tumor-heterogeneity mixing
  (responder/non-responder cells mixed 65/35 per patient) brings the ROC AUC
  from a perfect 1.0 to ~0.86, tighter noise-gene range sharpens the UMAP
  clone layout, and a fixed seed makes every Load Demo identical.
- Layout widths unified across Train / Predict / Visualize (sidebar ratio
  4:8); Data readiness order now groups the three user-uploaded inputs first.
- `annotate_clones()` / `prepare_data()` gained an optional `progress_cb`
  and a `variable_selection` option (vst default; dispersion/mvp for speed).

## Bug fixes
- Background workers could not start in container deployments
  ("Failed to start the background worker ... invalid connection").
  Spawning now tries a ladder of progressively simpler callr configs
  (supervise + stderr file -> /dev/null -> inherited stdio) and reports
  rich diagnostics on final failure.
- `stdout = FALSE` was rejected with EACCES in the rocker/shiny-verse
  container; stdout now redirects to `/dev/null` instead.
- Removed `closeAllConnections()` from the Shiny session-end callback; it
  corrupted the process connection table and killed later spawns.
- DepMap load no longer disconnects the session; demo load failures fixed.
- Clearing demo data now also drops stale predictions; invalid expression
  files are rejected up front.
- ROC curve corners (0,0) / (1,1) no longer render as a broken line
  (smoothed coords are densified onto a uniform FPR grid).
- ROC / boxplot defaulted to the first drug on first click instead of the
  selected "Combination" (selectize init race); Gene Expression no longer
  shows "Select a gene first" on first draw.
- DepMap cache is never auto-expired on disk (removed the TTL cleanup that
  could delete a 567 MB pre-downloaded cache); `cache_ttl_hours <= 0` no
  longer deletes the file immediately.
- Removed a stray `inst/shiny/app/FALSE` binary file from the repository.

## Changes
- Removed the on-disk TTL mechanism for the DepMap cache; memory is instead
  released by the master worker's idle timeout.

# PERCEPTIONx 0.1.0

Initial release.
