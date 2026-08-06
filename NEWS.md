# PERCEPTIONx 0.1.0

* Initial release of the R implementation of the PERCEPTION approach
  (predicting personalized drug response from single-cell transcriptomics).
* Train drug-response models on large-scale cell-line screens (`train_models`)
  or load pre-trained models (`load_model`).
* Predict response at clone level (`predict_drugs`) and patient level
  (`predict_patients`) with multiple aggregation strategies.
* Interactive visualizations with hover tooltips (`ggiraph`) and
  publication-quality export (600 dpi PNG, PDF, SVG).
* Shiny web application (`run_perception_app`) covering the whole pipeline.
* Vignette "The PERCEPTIONx Pipeline" documents the complete workflow.
