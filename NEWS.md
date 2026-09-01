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
