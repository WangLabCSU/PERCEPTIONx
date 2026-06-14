#' PERCEPTION Visualization Functions
#'
#' This module provides visualization functions for PERCEPTION model results,
#' including t-SNE/UMAP plots, clone distribution plots, ROC curves, and more.
#' These functions are designed to work seamlessly with trained PERCEPTION models
#' and patient/single-cell expression data.
#'
#' @name plot_perception
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_boxplot geom_bar geom_hline geom_line geom_vline coord_cartesian theme_bw labs theme element_text element_rect facet_grid scale_colour_gradientn scale_fill_brewer scale_size margin unit guide_colourbar guides guide_legend ggtitle annotate rel
#' @importFrom rlang .data
#' @importFrom pROC ggroc roc auc smooth
NULL

# Column names used in aes() - declare as global variables to suppress R CMD check notes
utils::globalVariables(c("X", "Y", "clones", "weights", "patients",
                         "clone_id", "Predictibility", "drugsCount", "dataused",
                         "pred_viab"))

#' Plot t-SNE/UMAP with drug response overlay
#'
#' Visualizes single cells in t-SNE/UMAP space with color overlay representing
#' either biomarker expression or predicted drug sensitivity.
#'
#' @param tsne_data Data frame with columns: X, Y (coordinates), and optional
#'        biomarker/killing columns.
#' @param color_var Character. Name of the column to use for color mapping.
#'        Default = "killing_scaled".
#' @param title Character. Plot title. Default = "Drug Killing".
#' @param color_label Character. Legend label for color. Default = "Predicted Killing".
#' @param point_size Numeric. Size of points. Default = 0.5.
#' @param colors Character vector. Gradient colors (low, mid, high).
#'        Default = c("#F8766D", "lightgrey", "#00BFC4").
#' @param base_size Numeric. Base font size for theme. Default = 8.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   # After predicting killing for single cells
#'   tsne_data <- data.frame(
#'     X = lung_tSNE$X,
#'     Y = lung_tSNE$Y,
#'     killing_scaled = range01(rank(-viability_pred))
#'   )
#'   plot_tsne_response(tsne_data, color_var = "killing_scaled")
#' }
#'
#' @export
plot_tsne_response <- function(tsne_data,
                               color_var = "killing_scaled",
                               title = NULL,
                               color_label = "Predicted Killing",
                               point_size = 0.5,
                               colors = c("#F8766D", "lightgrey", "#00BFC4"),
                               base_size = 8) {

  if (!all(c("X", "Y", color_var) %in% colnames(tsne_data))) {
    stop("tsne_data must contain columns: X, Y, and ", color_var)
  }

  p <- ggplot(tsne_data, aes(x = X, y = Y, color = .data[[color_var]])) +
    geom_point(size = point_size) +
    theme_bw(base_size = base_size) +
    labs(color = color_label, x = "", y = "") +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradientn(colours = colors)

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  return(p)
}

#' Plot clone distribution as stacked bar
#'
#' Visualizes the proportion of each clone across patients as a stacked bar plot.
#' Useful for understanding tumor heterogeneity and clonal architecture.
#'
#' @param clone_distribution Data frame with columns: patients, clones, weights.
#' @param response_var Character. Optional column name for response annotation.
#'        If provided, facets by response. Default = NULL.
#' @param base_size Numeric. Base font size. Default = 15.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   # After computing clone weights
#'   clone_dist <- data.frame(
#'     patients = c("P1", "P1", "P1", "P2", "P2", "P2"),
#'     clones = c("c1", "c2", "c3", "c1", "c2", "c3"),
#'     weights = c(0.3, 0.5, 0.2, 0.6, 0.3, 0.1)
#'   )
#'   plot_clone_distribution(clone_dist)
#' }
#'
#' @export
plot_clone_distribution <- function(clone_distribution,
                                    response_var = NULL,
                                    base_size = 15) {

  required_cols <- c("patients", "clones", "weights")
  if (!all(required_cols %in% colnames(clone_distribution))) {
    stop("clone_distribution must contain columns: ", paste(required_cols, collapse = ", "))
  }

  p <- ggplot(clone_distribution, aes(fill = clones, y = weights, x = patients)) +
    geom_bar(position = "stack", stat = "identity") +
    theme_bw(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "top") +
    labs(y = "Clone Proportion", x = "Patients")

  if (!is.null(response_var) && response_var %in% colnames(clone_distribution)) {
    p <- p + facet_grid(. ~ .data[[response_var]], shrink = TRUE,
                        scales = "free", space = "free_x")
  }

  return(p)
}

#' Plot clone-level killing (lollipop plot)
#'
#' Visualizes predicted drug sensitivity for each clone within patients.
#' Each clone is represented as a point with a stem (lollipop style).
#' Useful for identifying resistant clones within heterogeneous tumors.
#'
#' @param clone_killing Data frame with columns: patient, clone_id, killing (or drug-specific).
#' @param killing_var Character. Column name for killing values. Default = "comb_killing".
#' @param weights_var Character. Optional column name for clone weights (point size).
#'        Default = NULL.
#' @param response_var Character. Optional column for response annotation.
#'        Default = NULL.
#' @param base_size Numeric. Base font size. Default = 15.
#' @param y_limits Numeric vector. Y-axis limits. Default = c(-3, 1.2).
#' @param viridis_scale Logical. Use viridis color scale. Default = TRUE.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   clone_kill <- data.frame(
#'     patient = c("P1", "P1", "P2", "P2"),
#'     clone_id = c("c1", "c2", "c1", "c2"),
#'     comb_killing = c(-0.5, 0.8, -1.2, 0.3)
#'   )
#'   plot_clone_killing(clone_kill, killing_var = "comb_killing")
#' }
#'
#' @export
plot_clone_killing <- function(clone_killing,
                               killing_var = "comb_killing",
                               weights_var = NULL,
                               response_var = NULL,
                               base_size = 15,
                               y_limits = c(-3, 1.2),
                               viridis_scale = TRUE) {

  if (!all(c("patient", "clone_id", killing_var) %in% colnames(clone_killing))) {
    stop("clone_killing must contain columns: patient, clone_id, and ", killing_var)
  }

  # Build aes mapping
  aes_mapping <- aes(y = .data[[killing_var]], x = clone_id)

  if (!is.null(weights_var) && weights_var %in% colnames(clone_killing)) {
    aes_mapping <- aes(y = .data[[killing_var]], x = clone_id,
                       color = .data[[killing_var]], size = .data[[weights_var]])
  } else {
    aes_mapping <- aes(y = .data[[killing_var]], x = clone_id,
                       color = .data[[killing_var]])
  }

  y_data <- clone_killing[[killing_var]]
  y_min <- min(y_data, na.rm = TRUE)
  y_max <- max(y_data, na.rm = TRUE)
  y_bottom <- min(0, y_min)
  y_top <- y_max + (y_max - y_bottom) * 0.1

  p <- ggplot(clone_killing, aes_mapping) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_segment(aes(x = clone_id, xend = clone_id, y = 0, yend = .data[[killing_var]]),
                 color = "black", linewidth = 0.4) +
    geom_point() +
    coord_cartesian(ylim = c(y_bottom, y_top)) +
    theme_bw(base_size = base_size) +
    labs(x = "Clones", y = "Predicted Viability (z-score)",
         color = "Predicted Viability", size = "Proportion in Tumor") +
    theme(legend.position = "top",
          strip.placement = "outside",
          strip.background = element_rect(fill = "white", linewidth = 1, color = "white"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.box = "horizontal",
          legend.box.spacing = unit(8, "pt"),
          legend.key.size = unit(14, "pt"),
          legend.text = element_text(size = rel(0.7)),
          legend.title = element_text(size = rel(0.8), vjust = 0.5),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

  if (viridis_scale) {
    p <- p + scale_colour_gradientn(
        colours = c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"),
        breaks = function(limits) pretty(limits, n = 4)
    ) +
      guides(color = guide_colourbar(
        barwidth = unit(120, "pt"),
        barheight = unit(8, "pt"),
        label.hjust = 0.5,
        ticks.colour = "grey50",
        title.position = "top"
      ))
  }

  if (!is.null(weights_var) && weights_var %in% colnames(clone_killing)) {
    p <- p + scale_size(range = c(1, 8))
  } else {
    p <- p + guides(size = "none")
  }

  if (!is.null(response_var) && response_var %in% colnames(clone_killing)) {
    p <- p + facet_grid(reformulate(c(response_var, "patient"), "."),
                        scales = "free_x", shrink = TRUE,
                        drop = TRUE, space = "free_x", switch = "x",
                        as.table = TRUE)
  } else {
    p <- p + facet_grid(reformulate("patient", "."),
                        scales = "free_x", shrink = TRUE,
                        drop = TRUE, space = "free_x", switch = "x",
                        as.table = TRUE)
  }

  p
}

#' Plot ROC curve with AUC annotation
#'
#' Generates a ROC curve from predicted vs observed response with AUC value annotation.
#'
#' @param response Factor or numeric. True response labels (e.g., "R"/"NR" or 0/1).
#' @param predictor Numeric. Predicted values (e.g., killing scores).
#' @param smooth_curve Logical. Whether to smooth the ROC curve. Default = TRUE.
#' @param base_size Numeric. Base font size. Default = 15.
#' @param auc_digits Integer. Number of digits for AUC display. Default = 3.
#' @param title Character. Plot title. Default = NULL.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   response <- factor(c("R", "NR", "R", "NR", "R"))
#'   predictor <- c(0.8, 0.2, 0.7, 0.3, 0.9)
#'   plot_roc_curve(response, predictor)
#' }
#'
#' @export
plot_roc_curve <- function(response,
                           predictor,
                           smooth_curve = TRUE,
                           base_size = 15,
                           auc_digits = 3,
                           title = NULL) {

  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required for ROC curve plotting. Install with: install.packages('pROC')")
  }

  rocobj <- pROC::roc(response = response, predictor = predictor)

  if (smooth_curve) {
    rocobj_smooth <- pROC::smooth(rocobj)
    p <- ggroc(rocobj_smooth)
  } else {
    p <- ggroc(rocobj)
  }

  p <- p +
    theme_bw(base_size = base_size) +
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,
             color = "grey", linetype = "dashed") +
    annotate("text", x = 0.2, y = 0.1,
             label = paste0("AUC=", round(rocobj$auc, auc_digits)),
             size = 5, color = "black")

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  return(p)
}

#' Plot predicted vs observed response boxplot
#'
#' Creates a boxplot comparing predicted viability between response groups
#' (e.g., Responders vs Non-Responders) with statistical annotation.
#'
#' @param exp_vs_pred Data frame with columns: response, predicted_killing.
#' @param response_var Character. Column name for response labels. Default = "response".
#' @param predicted_var Character. Column name for predicted values. Default = "predicted_killing".
#' @param y_label Character. Y-axis label. Default = "Predicted Viability (z-score)".
#' @param base_size Numeric. Base font size. Default = 15.
#' @param compare_method Character. Statistical test method. Default = "wilcox.test".
#' @param alternative Character. Alternative hypothesis direction. Default = "greater".
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   exp_pred <- data.frame(
#'     response = factor(c("R", "NR", "R", "NR")),
#'     predicted_killing = c(0.8, 0.2, 0.7, 0.3)
#'   )
#'   plot_response_boxplot(exp_pred)
#' }
#'
#' @export
plot_response_boxplot <- function(exp_vs_pred,
                                  response_var = "response",
                                  predicted_var = "predicted_killing",
                                  y_label = "Predicted Viability (z-score)",
                                  base_size = 15,
                                  compare_method = "wilcox.test",
                                  alternative = "greater") {

  if (!all(c(response_var, predicted_var) %in% colnames(exp_vs_pred))) {
    stop("exp_vs_pred must contain columns: ", response_var, " and ", predicted_var)
  }

  # Ensure response is factor with correct order
  exp_vs_pred[[response_var]] <- factor(exp_vs_pred[[response_var]])
  if (length(levels(exp_vs_pred[[response_var]])) == 2) {
    exp_vs_pred[[response_var]] <- factor(exp_vs_pred[[response_var]],
                                          rev(levels(exp_vs_pred[[response_var]])))
  }

  p <- ggplot(exp_vs_pred, aes(y = .data[[predicted_var]], x = .data[[response_var]],
                               color = .data[[response_var]])) +
    geom_boxplot() +
    geom_point(size = 1, alpha = 0.5) +
    ggpubr::stat_compare_means(method.args = list(alternative = alternative),
                       size = 5, label = "p",
                       label.x.npc = 0.95, label.y.npc = 0.95) +
    theme_bw(base_size = base_size) +
    labs(y = y_label, x = "Patients") +
    theme(legend.position = "top")

  return(p)
}

#' Plot model performance across datasets
#'
#' Visualizes the number of drugs achieving different correlation thresholds
#' across bulk, pseudo-bulk, and single-cell datasets.
#'
#' @param performance_list Named list of model performance objects (from train_perception_models).
#' @param threshold_range Numeric vector. Correlation thresholds to evaluate.
#'        Default = seq(0.1, 0.6, 0.01).
#' @param base_size Numeric. Base font size. Default = 20.
#' @param highlight_threshold Numeric. Threshold to highlight with vertical line.
#'        Default = 0.3.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   models <- train_perception_models(c("abemaciclib", "erlotinib"), ...)
#'   plot_model_performance(models)
#' }
#'
#' @export
plot_model_performance <- function(performance_list,
                                   threshold_range = seq(0.1, 0.6, 0.01),
                                   base_size = 20,
                                   highlight_threshold = 0.3) {

  # Extract performance metrics
  performance_in_scRNA <- data.frame(do.call(rbind,
    lapply(performance_list, function(x) x$performance_in_scRNA)))
  performance_in_bulk <- data.frame(do.call(rbind,
    lapply(performance_list, function(x) x$performance_in_bulk)))
  performance_in_pseudo_bulk <- data.frame(do.call(rbind,
    lapply(performance_list, function(x) x$performance_in_pseudo_bulk)))

  # Build summary data frame
  df2plot <- rbind(
    data.frame(
      drugsCount = sapply(threshold_range, function(x) sum(performance_in_scRNA$estimate.cor > x)),
      dataused = "scRNA-seq",
      Predictibility = threshold_range
    ),
    data.frame(
      drugsCount = sapply(threshold_range, function(x) sum(performance_in_bulk$estimate.cor > x)),
      dataused = "bulk",
      Predictibility = threshold_range
    ),
    data.frame(
      drugsCount = sapply(threshold_range, function(x) sum(performance_in_pseudo_bulk$estimate.cor > x)),
      dataused = "pseudo-bulk",
      Predictibility = threshold_range
    )
  )

  p <- ggplot(df2plot, aes(x = Predictibility, y = drugsCount, color = dataused)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = highlight_threshold, linetype = "dashed") +
    theme_bw(base_size = base_size) +
    labs(y = "Number of Drugs", color = "Validation Dataset",
         x = "Predictibility (Pearson Correlation)") +
    theme(legend.position = "top")

  return(p)
}

#' Run Seurat clustering and plot UMAP
#'
#' Performs Seurat clustering on an expression matrix and generates UMAP visualization.
#' Useful for identifying subclones within patient tumor samples.
#'
#' @param expression_matrix Matrix. Gene expression matrix (genes as rows, cells as columns).
#' @param min_cells Integer. Minimum cells per feature. Default = 3.
#' @param min_features Integer. Minimum features per cell. Default = 200.
#' @param nfeatures Integer. Number of variable features. Default = 2000.
#' @param dims Integer. Number of PCA dimensions for clustering. Default = 10.
#' @param resolution Numeric. Clustering resolution. Default = 0.8.
#' @param seed Integer. Random seed. Default = 1.
#'
#' @return A list containing:
#'   \item{seurat_object}{Seurat object with clustering results}
#'   \item{umap_plot}{ggplot UMAP visualization}
#'   \item{cluster_ids}{Named vector of cluster IDs per cell}
#'
#' @examples
#' \dontrun{
#'   result <- plot_seurat_clustering(patient_expression)
#'   result$umap_plot
#'   result$cluster_ids
#' }
#'
#' @export
plot_seurat_clustering <- function(expression_matrix,
                                  min_cells = 3,
                                  min_features = 200,
                                  nfeatures = 2000,
                                  dims = 10,
                                  resolution = 0.8,
                                  seed = 1) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Install with: install.packages('Seurat')")
  }

  set.seed(seed)

  # Create Seurat object
  so <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                   project = "PERCEPTION",
                                   min.cells = min_cells,
                                   min.features = min_features)

  # Standard workflow
  so <- Seurat::NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
  so <- Seurat::FindVariableFeatures(so, selection.method = "vst", nfeatures = nfeatures)
  so <- Seurat::ScaleData(so)
  so <- Seurat::RunPCA(so, features = Seurat::VariableFeatures(object = so))
  so <- Seurat::FindNeighbors(so, dims = 1:dims)
  so <- Seurat::FindClusters(so, resolution = resolution)
  so <- Seurat::RunUMAP(so, dims = 1:dims)

  # Generate plot
  umap_plot <- Seurat::DimPlot(so, reduction = "umap")

  # Extract cluster IDs
  cluster_ids <- Seurat::Idents(so)

  return(list(
    seurat_object = so,
    umap_plot = umap_plot,
    cluster_ids = cluster_ids
  ))
}

#' Complete patient response visualization pipeline
#'
#' Generates a comprehensive visualization panel for patient drug response prediction,
#' including clone distribution, clone-level killing, response boxplot, and ROC curve.
#' This is a convenience function that combines multiple plot functions.
#'
#' @param clone_distribution Data frame. Clone weights per patient.
#' @param clone_killing Data frame. Killing scores per clone.
#' @param exp_vs_pred Data frame. Predicted vs observed response.
#' @param response_col Character. Response column name. Default = "response".
#' @param killing_col Character. Killing column name. Default = "comb_killing".
#' @param predicted_col Character. Predicted values column name. Default = "predicted_killing".
#' @param weights_col Character. Weights column name. Default = "weights".
#' @param layout_matrix Matrix. Layout for grid.arrange. Default = NULL (auto).
#'
#' @return A gtable object from grid.arrange.
#'
#' @examples
#' \dontrun{
#'   # After running prediction pipeline
#'   panel <- plot_patient_response_panel(
#'     clone_distribution = clone_dist_df,
#'     clone_killing = clone_kill_df,
#'     exp_vs_pred = response_df
#'   )
#'   ggsave(panel, filename = "patient_response.pdf", height = 15, width = 10)
#' }
#'
#' @export
plot_patient_response_panel <- function(clone_distribution,
                                        clone_killing,
                                        exp_vs_pred,
                                        response_col = "response",
                                        killing_col = "comb_killing",
                                        predicted_col = "predicted_killing",
                                        weights_col = "weights",
                                        layout_matrix = NULL) {

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required. Install with: install.packages('gridExtra')")
  }

  # Panel 1: Clone distribution
  p1 <- plot_clone_distribution(clone_distribution, response_var = response_col)

  # Panel 2: Clone-level killing (uses killing_col from clone_killing data)
  p2 <- plot_clone_killing(clone_killing, killing_var = killing_col,
                           weights_var = weights_col, response_var = response_col)

  # Panel 3: Response boxplot (uses predicted_col from exp_vs_pred data)
  p3 <- plot_response_boxplot(exp_vs_pred, response_var = response_col,
                              predicted_var = predicted_col)

  # Panel 4: ROC curve (uses predicted_col from exp_vs_pred data)
  p4 <- plot_roc_curve(response = exp_vs_pred[[response_col]],
                       predictor = exp_vs_pred[[predicted_col]])

  # Arrange panels
  if (is.null(layout_matrix)) {
    layout_matrix <- rbind(
      c(1, 1),
      c(2, 2),
      c(2, 2),
      c(3, 4)
    )
  }

  panel <- gridExtra::grid.arrange(p1, p2, p3, p4, layout_matrix = layout_matrix)

  return(panel)
}

#' Plot t-SNE/UMAP side-by-side for biomarker and killing
#'
#' Creates a side-by-side comparison of biomarker expression and predicted killing
#' in t-SNE/UMAP space. Useful for visualizing correlation between marker and response.
#'
#' @param tsne_data Data frame with X, Y coordinates and both biomarker/killing columns.
#' @param biomarker_var Character. Column name for biomarker expression. Default = "biomarker_scaled".
#' @param killing_var Character. Column name for killing values. Default = "killing_scaled".
#' @param biomarker_label Character. Legend label for biomarker. Default = "Biomarker Exp".
#' @param killing_label Character. Legend label for killing. Default = "Drug Killing".
#' @param nrow Integer. Number of rows in arrangement. Default = 1.
#' @param base_size Numeric. Base font size. Default = 8.
#'
#' @return A gtable object from grid.arrange.
#'
#' @examples
#' \dontrun{
#'   tsne_data <- data.frame(
#'     X = lung_tSNE$X,
#'     Y = lung_tSNE$Y,
#'     biomarker_scaled = range01(rank(MDM2_expression)),
#'     killing_scaled = range01(rank(-viability_pred))
#'   )
#'   plot_tsne_biomarker_killing(tsne_data)
#' }
#'
#' @export
plot_tsne_biomarker_killing <- function(tsne_data,
                                        biomarker_var = "biomarker_scaled",
                                        killing_var = "killing_scaled",
                                        biomarker_label = "Biomarker Exp",
                                        killing_label = "Drug Killing",
                                        nrow = 1,
                                        base_size = 8) {

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required. Install with: install.packages('gridExtra')")
  }

  p1 <- plot_tsne_response(tsne_data, color_var = biomarker_var,
                           color_label = biomarker_label, base_size = base_size)
  p2 <- plot_tsne_response(tsne_data, color_var = killing_var,
                           color_label = killing_label, base_size = base_size)

  combined <- gridExtra::grid.arrange(p1, p2, nrow = nrow)

  return(combined)
}