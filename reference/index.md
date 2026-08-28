# Package index

## Data Loading

Download and load DepMap reference data and pre-trained models.

- [`load_depmap()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/load_depmap.md)
  : Download filtered DepMap data
- [`load_model()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/load_model.md)
  : Load pre-built model of provided drugs
- [`get_depmap()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_depmap.md)
  : Get the DepMap dataset
- [`get_mirrors()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_mirrors.md)
  : Get current download mirrors
- [`add_mirrors()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/add_mirrors.md)
  : Add custom download mirrors
- [`list_mirrors()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/list_mirrors.md)
  : List current download mirrors
- [`reset_mirrors()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/reset_mirrors.md)
  : Reset mirrors to default

## Preprocessing

- [`rank_normalization_mat()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/rank_normalization_mat.md)
  : Rank-normalize each column of a matrix
- [`range01()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/range01.md)
  : Change range to 0-1
- [`zscore_viability()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/zscore_viability.md)
  : Z-score scale viability values across patients
- [`prepare_data()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/prepare_data.md)
  : Prepare patient data for PERCEPTIONx prediction

## Model Training

- [`train_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/train_models.md)
  : Train PERCEPTIONx models for multiple drugs
- [`get_response_matrix()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_response_matrix.md)
  : Get specific drug response data for cell-lines
- [`get_cellLine_list()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_cellLine_list.md)
  : Determine training and test cell-lines for a given drug
- [`feature_ranking_bulk()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/feature_ranking_bulk.md)
  : Feature ranking for a single drug using bulk expression
- [`run_parallel_feature_ranking_bulk()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/run_parallel_feature_ranking_bulk.md)
  : Parallel feature ranking for multiple drugs
- [`build_on_BULK_v2()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/build_on_BULK_v2.md)
  : Build PERCEPTIONx model on bulk expression data

## Prediction

- [`predict_drugs()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_drugs.md)
  : Predict drug response for cells or clones
- [`predict_patients()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/predict_patients.md)
  : Predict drug response at patient level

## Evaluation

- [`compare_performance()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/compare_performance.md)
  : Compare performance of multiple trained models
- [`get_significant_models()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_significant_models.md)
  : Get best performing models
- [`get_performance()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/get_performance.md)
  : Load performance metrics from a saved model file
- [`each_patient_pseudo_bulk()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/each_patient_pseudo_bulk.md)
  : Compute pseudo-bulk expression for a patient

## Clone & Seurat Pipeline

- [`run_seurat_pipeline()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/run_seurat_pipeline.md)
  : PERCEPTIONx Patient Data Annotation Functions
- [`annotate_clones()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/annotate_clones.md)
  : Annotate cells with clone IDs via Seurat clustering
- [`build_clone_counts()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/build_clone_counts.md)
  : Build clone abundance table from cell-clone mapping
- [`build_clone_key()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/build_clone_key.md)
  : Build a clone key from patient and clone id
- [`parse_clone_keys()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/parse_clone_keys.md)
  : Parse clone keys in "Patient@Clone" format
- [`clone_mean_expression()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/clone_mean_expression.md)
  : Compute clone-level mean expression from single-cell data

## Visualization

- [`plot_tsne_response()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_tsne_response.md)
  : Plot UMAP with drug response overlay
- [`plot_tsne_biomarker_viability()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_tsne_biomarker_viability.md)
  : Plot UMAP side-by-side for biomarker and viability
- [`plot_clone_distribution()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_clone_distribution.md)
  : Plot clone distribution as stacked bar
- [`plot_clone_viability()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_clone_viability.md)
  : Plot clone-level viability (lollipop plot)
- [`plot_clone_umap()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_clone_umap.md)
  : Plot UMAP colored by clone identity
- [`plot_model_performance()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_model_performance.md)
  : Plot model performance across datasets
- [`plot_model_roc()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_model_roc.md)
  : Plot validation ROC curves for trained models
- [`plot_roc_curve()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_roc_curve.md)
  : Plot ROC curve with AUC annotation
- [`plot_response_boxplot()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_response_boxplot.md)
  : Plot predicted vs observed response boxplot
- [`plot_seurat_clustering()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_seurat_clustering.md)
  : Run Seurat clustering and plot 2D embedding
- [`plot_patient_response_panel()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/plot_patient_response_panel.md)
  : Complete patient response visualization pipeline

## Shiny Application

- [`run_perception_app()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/run_perception_app.md)
  : Launch PERCEPTIONx Shiny Dashboard

## Utilities

- [`hypergeometric_test_for_twolists()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/hypergeometric_test_for_twolists.md)
  : Hypergeometric test for gene list overlap
- [`fdrcorr()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/fdrcorr.md)
  : FDR correction for multiple testing
- [`cor.test_trimmed_v0()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/cor.test_trimmed_v0.md)
  : Fast correlation test (generic)
- [`export_plot_cairo()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/export_plot_cairo.md)
  : Export plot to file via Cairo device
