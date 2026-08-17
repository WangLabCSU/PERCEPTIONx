#' PERCEPTIONx Visualization Functions
#'
#' This module provides visualization functions for PERCEPTIONx model results,
#' including UMAP plots, clone distribution plots, ROC curves, and more.
#' These functions are designed to work seamlessly with trained PERCEPTIONx models
#' and patient/single-cell expression data.
#'
#' @name plot_perception
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_boxplot geom_violin geom_bar geom_hline geom_line geom_vline geom_jitter geom_text coord_cartesian coord_equal theme theme_bw theme_void labs element_text element_rect element_line element_blank facet_grid facet_wrap vars scale_colour_gradientn scale_colour_gradient2 scale_fill_manual scale_color_manual scale_x_discrete scale_size margin unit guide_colourbar guides guide_legend ggtitle annotate rel
#' @importFrom rlang .data
#' @importFrom pROC ggroc roc auc smooth
#' @importFrom stats t.test wilcox.test
NULL

# Column names used in aes() - declare as global variables to suppress R CMD check notes
utils::globalVariables(c("X", "Y", "clones", "weights", "patients",
                         "clone_id", "Predictability", "drugsCount", "dataused",
                         "pred_viab", "expression",
                         "tooltip_text", "data_id", "fpr", "tpr",
                         "x", "y", "label", "yintercept"))

# ---------------------------------------------------------------------------
# Design system
# ---------------------------------------------------------------------------

# Unified theme for all PERCEPTIONx plots: white background, subtle grid,
# consistent font hierarchy and legend styling.
theme_perception <- function(base_size = 11, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.background    = element_rect(fill = "white", color = NA),
      plot.background     = element_rect(fill = "white", color = NA),
      panel.border        = element_rect(color = "#dfe3ee", linewidth = 0.5),
      panel.grid.major    = element_line(color = "#eef0f6", linewidth = 0.4),
      panel.grid.minor    = element_blank(),
      axis.ticks          = element_line(color = "#b6becd"),
      axis.title          = element_text(face = "bold", color = "#1e2a4a"),
      axis.text           = element_text(color = "#3d4a63"),
      plot.title          = element_text(face = "bold", color = "#1e2a4a", hjust = 0,
                                         margin = margin(b = 4)),
      plot.subtitle       = element_text(color = "#5a6a8a", hjust = 0,
                                         margin = margin(b = 8)),
      legend.background   = element_rect(fill = "white", color = NA),
      legend.key          = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA),
      legend.title        = element_text(face = "bold", color = "#1e2a4a", size = rel(0.85)),
      legend.text         = element_text(color = "#3d4a63", size = rel(0.8)),
      strip.background    = element_rect(fill = "#f6f8fc", color = "#dfe3ee", linewidth = 0.5),
      strip.text          = element_text(face = "bold", color = "#1e2a4a", size = rel(0.85)),
      plot.margin         = margin(t = 6, r = 8, b = 6, l = 6)
    )
}

# Qualitative palette for many clones, zero extra dependencies.
# Base of 15 curated colors (blue -> purple -> pink -> coral -> orange ->
# yellow -> green -> teal -> indigo -> slate), kept at moderate saturation
# so neighboring clones stay distinguishable without clashing.
clone_palette_base <- c(
  "#c1dc55ff", "#20bb3dff", "#28a5dfff", "#2c6a9dff", "#5C6BC0",
  "#764BA2", "#AB47BC", "#c977d2ff", "#910013ff", "#d42e2bff",
  "#FF7043", "#FFCA28", "#18861dff", "#23acbeff", "#7e645bff"
)
clone_palette <- function(n) {
  if (n <= length(clone_palette_base)) {
    return(clone_palette_base[seq_len(n)])
  }
  # More clones: evenly space hues around the wheel, anchored on the deep
  # blue of the first curated color (hue ~214 deg), with moderate chroma and
  # high-ish lightness for a soft, harmonious look.
  hues <- (214 + seq(0, n - 1) * 360 / n) %% 360
  unname(grDevices::hcl(h = hues, c = 65, l = 70))
}

# Diverging scale for z-score style values: blue = low, white = neutral, red = high
diverging_colors <- c("#2166AC", "#F7F7F7", "#B2182B")
# Sequential scale for percentile/rank values (viridis family)
sequential_colors <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")

# ---------------------------------------------------------------------------
# Plot: UMAP / t-SNE response
# ---------------------------------------------------------------------------

#' Plot UMAP with drug response overlay
#'
#' Visualizes single cells in UMAP space with color overlay representing
#' either biomarker expression or predicted drug sensitivity.
#'
#' @param tsne_data Data frame with columns: X, Y (coordinates), and optional
#'        biomarker/viability columns.
#' @param color_var Character. Name of the column to use for color mapping.
#'        Default = "viability_scaled".
#' @param title Character. Plot title. Default = NULL.
#' @param color_label Character. Legend label for color. Default = "Predicted Viability".
#' @param point_size Numeric. Point size. If NULL (default), auto-adapts to
#'        the number of cells to avoid overplotting.
#' @param colors Character vector. Custom gradient colors (low, mid, high)
#'        for the sequential palette. Default = NULL (uses built-in viridis).
#' @param palette Character. One of \code{"viridis"} (sequential, default) or
#'        \code{"diverging"} (blue-white-red centered at \code{midpoint}).
#' @param midpoint Numeric. Center value for diverging palette. Default = 0.
#' @param base_size Numeric. Base font size for theme. Default = 11.
#' @param tooltip Logical. If TRUE (default) and \pkg{ggiraph} is installed,
#'        points get hover tooltips (cell id if present, else the colored value).
#' @param tooltip_col Character. Optional existing column used as the tooltip
#'        text (e.g. "cell_id"). Default = NULL (auto-builds from color_var).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   # After predicting viability for single cells
#'   tsne_data <- data.frame(
#'     X = lung_tSNE$X,
#'     Y = lung_tSNE$Y,
#'     viability_scaled = range01(rank(-viability_pred))
#'   )
#'   plot_tsne_response(tsne_data, color_var = "viability_scaled")
#' }
#'
#' @export
plot_tsne_response <- function(tsne_data,
                               color_var = "viability_scaled",
                               title = NULL,
                               color_label = "Predicted Viability",
                               point_size = NULL,
                               colors = NULL,
                               palette = c("viridis", "diverging"),
                               midpoint = 0,
                               base_size = 11,
                               tooltip = TRUE,
                               tooltip_col = NULL) {

  if (!all(c("X", "Y", color_var) %in% colnames(tsne_data))) {
    stop("tsne_data must contain columns: X, Y, and ", color_var)
  }

  if (!is.numeric(tsne_data[[color_var]])) {
    stop(color_var, " must be a numeric column for continuous color mapping")
  }

  palette <- match.arg(palette)
  n_points <- nrow(tsne_data)

  # Adaptive point size / alpha to prevent overplotting on large scRNA data.
  if (is.null(point_size)) {
    point_size <- if (n_points > 50000) 0.4
                  else if (n_points > 20000) 0.5
                  else if (n_points > 8000) 0.6
                  else 0.8
  }
  point_alpha <- if (n_points > 50000) 0.35
                 else if (n_points > 20000) 0.5
                 else if (n_points > 8000) 0.65
                 else 0.8

  # Downsample extremely large datasets for smooth interactive rendering.
  if (n_points > 50000) {
    set.seed(42)
    keep <- sample.int(n_points, 50000)
    tsne_data <- tsne_data[keep, , drop = FALSE]
    message("plot_tsne_response: downsampled to 50,000 points for rendering (n = ", n_points, ")")
  }

  # Tooltip for ggiraph interactivity.
  tt_ok <- tooltip && requireNamespace("ggiraph", quietly = TRUE)
  if (!is.null(tooltip_col) && tooltip_col %in% colnames(tsne_data)) {
    tsne_data$tooltip_text <- sprintf("%s\n%s: %.3f",
                                      tsne_data[[tooltip_col]], color_label,
                                      tsne_data[[color_var]])
  } else {
    tsne_data$tooltip_text <- sprintf("%s: %.3f", color_label,
                                      tsne_data[[color_var]])
  }
  tsne_data$data_id <- paste0("cell_", seq_len(nrow(tsne_data)))

  p <- ggplot(tsne_data, aes(x = X, y = Y, color = .data[[color_var]])) +
    { if (tt_ok) ggiraph::geom_point_interactive(
        aes(tooltip = tooltip_text, data_id = data_id),
        size = point_size, alpha = point_alpha)
      else geom_point(size = point_size, alpha = point_alpha) } +
    theme_perception(base_size = base_size) +
    labs(color = color_label, x = "", y = "") +
    theme(legend.position = "top",
          legend.key.height = unit(0.9, "lines"),
          legend.key.width = unit(2.6, "lines"))

  if (palette == "diverging") {
    p <- p + scale_colour_gradient2(low = diverging_colors[1],
                                    mid = diverging_colors[2],
                                    high = diverging_colors[3],
                                    midpoint = midpoint)
  } else {
    cols <- if (is.null(colors)) sequential_colors else colors
    p <- p + scale_colour_gradientn(colours = cols)
  }

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  return(p)
}

# ---------------------------------------------------------------------------
# Plot: Clone distribution
# ---------------------------------------------------------------------------

#' Plot clone distribution as stacked bar
#'
#' Visualizes the proportion of each clone across patients as a stacked bar plot.
#' Useful for understanding tumor heterogeneity and clonal architecture.
#'
#' @param clone_distribution Data frame with columns: patients, clones, weights.
#' @param response_var Character. Optional column name for response annotation.
#'        If provided, facets by response. Default = NULL.
#' @param base_size Numeric. Base font size. Default = 15.
#' @param tooltip Logical. If TRUE (default) and \pkg{ggiraph} is installed,
#'        bar segments get hover tooltips (clone + proportion).
#' @param tooltip_col Character. Optional existing column used as the tooltip
#'        text. Default = NULL (auto-builds a rich tooltip).
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
                                    base_size = 15,
                                    tooltip = TRUE,
                                    tooltip_col = NULL) {

  required_cols <- c("patients", "clones", "weights")
  if (!all(required_cols %in% colnames(clone_distribution))) {
    stop("clone_distribution must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Pre-convert the response column to an R/NR factor BEFORE building the ggplot
  # object (ggplot captures data at build time, so the factor must exist first;
  # this avoids the %+% data-update operator entirely).
  if (!is.null(response_var) && response_var %in% colnames(clone_distribution)) {
    resp_vec <- clone_distribution[[response_var]]
    if (all(toupper(resp_vec[!is.na(resp_vec)]) %in% c("R", "NR"))) {
      clone_distribution[[response_var]] <- factor(resp_vec, levels = c("R", "NR"))
    }
  }

  # Order patients by appearance (stable across facets)
  clone_distribution$patients <- factor(clone_distribution$patients,
                                        levels = unique(clone_distribution$patients))

  n_clones <- length(unique(clone_distribution$clones))
  n_patients <- length(unique(clone_distribution$patients))

  # Adaptive x-axis labels: rotate 45 deg only when many patients would
  # otherwise overlap. The FACET STRIP labels (R/NR on top) stay horizontal.
  x_angle <- if (n_patients <= 8) 0 else 45
  x_hjust <- if (n_patients <= 8) 0.5 else 1
  x_vjust <- if (n_patients <= 8) 0.5 else 1
  x_breaks <- NULL
  if (n_patients > 24) {
    keep_idx <- seq(1, n_patients, by = 2)
    x_breaks <- levels(clone_distribution$patients)[keep_idx]
  }

  # Tooltip for ggiraph interactivity (hover shows clone + proportion).
  tt_ok <- tooltip && requireNamespace("ggiraph", quietly = TRUE)
  if (is.null(tooltip_col) || !(tooltip_col %in% colnames(clone_distribution))) {
    clone_distribution$tooltip_text <- sprintf(
      "Patient: %s\nClone: %s\nProportion: %.2f",
      clone_distribution$patients,
      clone_distribution$clones,
      clone_distribution$weights
    )
  }
  tt_col <- if (!is.null(tooltip_col) && tooltip_col %in% colnames(clone_distribution)) {
    tooltip_col
  } else {
    "tooltip_text"
  }

  p <- ggplot(clone_distribution, aes(fill = clones, y = weights, x = patients)) +
    { if (tt_ok) ggiraph::geom_bar_interactive(
        aes(tooltip = .data[[tt_col]], data_id = .data[[tt_col]]),
        position = "stack", stat = "identity", width = 0.75,
        color = "white", linewidth = 0.15)
      else geom_bar(position = "stack", stat = "identity", width = 0.75,
                    color = "white", linewidth = 0.15) } +
    scale_fill_manual(values = clone_palette(n_clones)) +
    scale_x_discrete(breaks = x_breaks) +
    theme_perception(base_size = base_size) +
    theme(axis.text.x = element_text(angle = x_angle, hjust = x_hjust, vjust = x_vjust,
                                     size = rel(0.8), margin = margin(t = 4)),
          legend.position = "bottom",
          legend.text = element_text(size = rel(0.7)),
          legend.key.height = unit(0.7, "lines"),
          strip.placement = "outside",
          strip.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                    size = rel(0.8))) +
    labs(y = "Clone Proportion", x = "Patients", fill = "Clone")

  # Facet by response group (R/NR) — Y axis kept fixed (0-1) so proportions
  # are directly comparable across facets; only X is free.
  if (!is.null(response_var) && response_var %in% colnames(clone_distribution)) {
    p <- p +
      facet_grid(cols = vars(.data[[response_var]]), shrink = TRUE,
                 scales = "free_x", space = "free_x")
  }

  return(p)
}

# ---------------------------------------------------------------------------
# Plot: Clone-level viability (lollipop)
# ---------------------------------------------------------------------------

#' Plot clone-level viability (lollipop plot)
#'
#' Visualizes predicted drug sensitivity for each clone within patients.
#' Each clone is represented as a point with a stem (lollipop style).
#' Useful for identifying resistant clones within heterogeneous tumors.
#'
#' Facet strategy is adaptive:
#' \itemize{
#'   \item <= 12 patients: single-row grid, strips on the bottom (45 deg).
#'   \item > 12 patients: \code{facet_wrap} grid, one compact panel per patient.
#' }
#'
#' @param clone_viability Data frame with columns: patient, clone_id, viability (or drug-specific).
#' @param viability_var Character. Column name for viability values. Default = "comb_viability".
#' @param weights_var Character. Optional column name for clone weights (point size).
#'        Default = NULL.
#' @param response_var Character. Optional column for response annotation.
#'        Default = NULL.
#' @param drug Character. Drug name, used as plot title. Default = NULL.
#' @param base_size Numeric. Base font size. Default = 11.
#' @param y_limits Numeric vector. Y-axis limits. Default = c(-3, 1.2).
#' @param viridis_scale Logical. If TRUE (default), uses the standard PERCEPTIONx
#'        diverging scale (blue = sensitive, red = resistant) centred at 0.
#'        If FALSE, uses a plain viridis sequential scale.
#' @param tooltip Logical. If TRUE (default) and \pkg{ggiraph} is installed,
#'        points get hover tooltips (clone + viability + proportion).
#' @param tooltip_col Character. Optional existing column used as the tooltip
#'        text. Default = NULL (auto-builds a rich tooltip).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   clone_kill <- data.frame(
#'     patient = c("P1", "P1", "P2", "P2"),
#'     clone_id = c("c1", "c2", "c1", "c2"),
#'     comb_viability = c(-0.5, 0.8, -1.2, 0.3)
#'   )
#'   plot_clone_viability(clone_kill, viability_var = "comb_viability")
#' }
#'
#' @export
plot_clone_viability <- function(clone_viability,
                               viability_var = "comb_viability",
                               weights_var = NULL,
                               response_var = NULL,
                               drug = NULL,
                               base_size = 11,
                               y_limits = c(-3, 1.2),
                               viridis_scale = TRUE,
                               tooltip = TRUE,
                               tooltip_col = NULL) {

  if (!all(c("patient", "clone_id", viability_var) %in% colnames(clone_viability))) {
    stop("clone_viability must contain columns: patient, clone_id, and ", viability_var)
  }

  # --- Prepare facet labels ---
  facet_col <- "patient"
  if (!is.null(response_var) && response_var %in% colnames(clone_viability)) {
    resp_map <- clone_viability[!duplicated(clone_viability$patient), ]
    rownames(resp_map) <- resp_map$patient
    # Order: Responders first
    patient_order <- unique(clone_viability$patient[
      order(clone_viability[[response_var]], decreasing = TRUE)
    ])
    clone_viability$patient <- factor(clone_viability$patient, levels = patient_order)
    clone_viability$facet_label <- vapply(as.character(clone_viability$patient), function(pat) {
      rv <- as.character(resp_map[pat, response_var])
      tag <- if (toupper(rv) %in% c("R", "RESPONDER")) "R" else "NR"
      paste0(tag, " ", pat)
    }, character(1))
    # Two-line variant used by facet_wrap (top strip, no rotation)
    clone_viability$facet_label_wrap <- vapply(as.character(clone_viability$patient), function(pat) {
      rv <- as.character(resp_map[pat, response_var])
      tag <- if (toupper(rv) %in% c("R", "RESPONDER")) "R" else "NR"
      paste0(tag, "\n", pat)
    }, character(1))
  } else {
    clone_viability$facet_label <- as.character(clone_viability$patient)
    clone_viability$facet_label_wrap <- as.character(clone_viability$patient)
  }

  # Preserve patient ordering via factor levels on facet_label
  if (!is.null(response_var) && response_var %in% colnames(clone_viability)) {
    clone_viability$facet_label <- factor(clone_viability$facet_label,
                                        levels = unique(clone_viability$facet_label[
                                          order(clone_viability$patient)
                                        ]))
    clone_viability$facet_label_wrap <- factor(clone_viability$facet_label_wrap,
                                             levels = unique(clone_viability$facet_label_wrap[
                                               order(clone_viability$patient)
                                             ]))
  }

  # --- Sort clones within each patient by proportion (descending) ---
  if (!is.null(weights_var) && weights_var %in% colnames(clone_viability)) {
    clone_viability <- clone_viability[order(clone_viability$facet_label,
                                         -clone_viability[[weights_var]]), ]
  }

  # Build aes mapping. NOTE: `size` (clone proportion) is mapped ONLY on
  # geom_point — a global size mapping would be inherited by the line layers
  # (geom_segment/geom_hline), which triggers the ggplot2 "size for lines is
  # deprecated" warning.
  has_weights <- !is.null(weights_var) && weights_var %in% colnames(clone_viability)
  aes_mapping <- aes(y = .data[[viability_var]], x = clone_id,
                     color = .data[[viability_var]])

  # Tooltip for ggiraph interactivity (hover shows clone, viability, proportion).
  tt_ok <- tooltip && requireNamespace("ggiraph", quietly = TRUE)
  if (is.null(tooltip_col) || !(tooltip_col %in% colnames(clone_viability))) {
    clone_viability$tooltip_text <- if (has_weights) {
      sprintf("Clone: %s\nPatient: %s\nViability: %.2f\nProportion: %.2f",
              clone_viability$clone_id, clone_viability$facet_label,
              clone_viability[[viability_var]], clone_viability[[weights_var]])
    } else {
      sprintf("Clone: %s\nPatient: %s\nViability: %.2f",
              clone_viability$clone_id, clone_viability$facet_label,
              clone_viability[[viability_var]])
    }
  }
  tt_col <- if (!is.null(tooltip_col) && tooltip_col %in% colnames(clone_viability)) {
    tooltip_col
  } else {
    "tooltip_text"
  }

  # Adaptive Y limits
  y_data <- clone_viability[[viability_var]]
  y_min <- min(y_data, na.rm = TRUE)
  y_max <- max(y_data, na.rm = TRUE)
  y_bottom <- min(0, y_min)
  y_top <- y_max + (y_max - y_bottom) * 0.1

  # Adaptive facet layout: wrap for many patients, single-row grid otherwise
  n_patients <- length(unique(clone_viability$facet_label))
  wrap_mode <- n_patients > 12

  # hline data (one per facet) so ggplotly can locate the facet column
  facet_var <- if (wrap_mode) "facet_label_wrap" else "facet_label"
  hline_data <- clone_viability[!duplicated(clone_viability[[facet_var]]), facet_var, drop = FALSE]
  hline_data$yintercept <- 0

  p <- ggplot(clone_viability, aes_mapping) +
    geom_hline(data = hline_data, mapping = aes(yintercept = yintercept),
               color = "grey60", linewidth = 0.3) +
    geom_segment(aes(x = clone_id, xend = clone_id, y = 0, yend = .data[[viability_var]]),
                 color = "grey40", linewidth = 0.35) +
    { if (tt_ok && has_weights) ggiraph::geom_point_interactive(
        aes(size = .data[[weights_var]], tooltip = .data[[tt_col]],
            data_id = .data[[tt_col]]), alpha = 0.85)
      else if (tt_ok) ggiraph::geom_point_interactive(
        aes(tooltip = .data[[tt_col]], data_id = .data[[tt_col]]), alpha = 0.85)
      else if (has_weights) geom_point(aes(size = .data[[weights_var]]), alpha = 0.85)
      else geom_point(alpha = 0.85) } +
    coord_cartesian(ylim = c(y_bottom, y_top)) +
    theme_perception(base_size = base_size) +
    labs(x = "Clones", y = "Predicted Viability (z-score)",
         color = "Predicted Viability",
         size = if (has_weights) "Proportion in Tumor" else NULL) +
    theme(legend.position = "top",
          legend.box = "horizontal",
          legend.box.spacing = unit(8, "pt"),
          legend.key.size = unit(13, "pt"),
          axis.title.y = element_text(margin = margin(r = 6)))

  # Color scale: diverging blue-white-red (sensitive -> neutral -> resistant)
  if (viridis_scale) {
    p <- p + scale_colour_gradient2(low = diverging_colors[1],
                                    mid = diverging_colors[2],
                                    high = diverging_colors[3],
                                    midpoint = 0,
                                    breaks = function(limits) pretty(limits, n = 4))
  } else {
    p <- p + scale_colour_gradientn(
      colours = sequential_colors,
      breaks = function(limits) pretty(limits, n = 4)
    )
  }

  if (has_weights) {
    p <- p + scale_size(range = c(1, 7))
  } else {
    p <- p + guides(size = "none")
  }

  if (wrap_mode) {
    ncol <- min(6, max(2, ceiling(sqrt(n_patients) * 1.4)))
    p <- p + facet_wrap(reformulate(facet_var), scales = "free_x",
                        ncol = ncol, drop = TRUE, strip.position = "top") +
      theme(strip.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                      size = rel(0.62), lineheight = 0.9),
            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                       size = rel(0.55)),
            panel.spacing = unit(0.3, "lines"))
  } else {
    p <- p + facet_grid(reformulate(facet_var, "."),
                        scales = "free_x", shrink = TRUE,
                        drop = TRUE, space = "free_x", switch = "x",
                        as.table = TRUE) +
      theme(strip.text = element_text(angle = 45, hjust = 0, vjust = 0.5,
                                      size = rel(0.75), lineheight = 0.9,
                                      margin = margin(t = 2, b = 2, l = 4, r = 4, unit = "pt")),
            panel.spacing = unit(0.15, "lines"))
  }

  if (!is.null(drug) && nchar(drug) > 0) {
    p <- p + ggtitle(drug, subtitle = sprintf("N = %d patients", n_patients))
  }

  p
}

# ---------------------------------------------------------------------------
# Plot: ROC curve
# ---------------------------------------------------------------------------

#' Plot ROC curve with AUC annotation
#'
#' Generates a ROC curve from predicted vs observed response with AUC value annotation.
#'
#' @param response Factor or numeric. True response labels (e.g., "R"/"NR" or 0/1).
#' @param predictor Numeric. Predicted values (e.g., viability scores).
#' @param smooth_curve Logical. Whether to smooth the ROC curve. Default = TRUE.
#' @param base_size Numeric. Base font size. Default = 15.
#' @param auc_digits Integer. Number of digits for AUC display. Default = 3.
#' @param title Character. Plot title. Default = NULL.
#' @param tooltip Logical. If TRUE (default) and \pkg{ggiraph} is installed,
#'        curve points get hover tooltips (FPR / TPR).
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
                           title = NULL,
                           tooltip = TRUE) {

  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required for ROC curve plotting. Install with: install.packages('pROC')")
  }

  # Single-class response cannot build a ROC curve — return an informative placeholder
  if (length(unique(response[!is.na(response)])) < 2) {
    stop("ROC curve requires at least two response classes (e.g. both R and NR).")
  }

  rocobj <- pROC::roc(response = response, predictor = predictor)

  if (smooth_curve) {
    rocobj_smooth <- tryCatch(
      pROC::smooth(rocobj),
      error = function(e) {
        warning("ROC curve not smoothable: ", e$message, ". Using unsmoothed curve.")
        NULL
      }
    )
    if (!is.null(rocobj_smooth)) {
      rocobj <- rocobj_smooth
    }
  }

  # Extract curve coordinates for both static & interactive rendering.
  curve_df <- as.data.frame(
    pROC::coords(rocobj, x = "all", ret = c("specificity", "sensitivity"),
                 transpose = FALSE)
  )
  curve_df$fpr <- 1 - curve_df$specificity
  curve_df$tpr <- curve_df$sensitivity
  curve_df$tooltip_text <- sprintf("FPR: %.3f\nTPR: %.3f", curve_df$fpr, curve_df$tpr)
  curve_df$data_id <- paste0("roc_", seq_len(nrow(curve_df)))

  tt_ok <- tooltip && requireNamespace("ggiraph", quietly = TRUE)

  if (tt_ok) {
    p <- ggplot(curve_df, aes(x = fpr, y = tpr)) +
      ggiraph::geom_line_interactive(aes(tooltip = tooltip_text, data_id = data_id),
                                     linewidth = 1.1, color = "#2E86AB") +
      ggiraph::geom_point_interactive(aes(tooltip = tooltip_text, data_id = data_id),
                                      size = 1.4, color = "#2E86AB", alpha = 0.7)
  } else {
    p <- ggplot(curve_df, aes(x = fpr, y = tpr)) +
      geom_line(linewidth = 1.1, color = "#2E86AB") +
      geom_point(size = 1.4, color = "#2E86AB", alpha = 0.7)
  }

  p <- p +
    theme_perception(base_size = base_size) +
    labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") +
    coord_equal() +
    # Random-chance diagonal from (0,0) to (1,1)
    annotate("segment", x = 0, xend = 1, y = 0, yend = 1,
             color = "grey60", linewidth = 0.5, linetype = "dashed") +
    # AUC badge in a rounded, semi-transparent box (bottom-right corner,
    # inside the region below the curve).
    annotate("rect", xmin = 0.60, xmax = 0.98, ymin = 0.05, ymax = 0.20,
             fill = "white", alpha = 0.85, color = "#dfe3ee", linewidth = 0.3) +
    annotate("text", x = 0.79, y = 0.125,
             label = paste0("AUC = ", round(rocobj$auc, auc_digits)),
             size = base_size * 0.28, color = "#1e2a4a", fontface = "bold")

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  return(p)
}

# ---------------------------------------------------------------------------
# Plot: Response boxplot
# ---------------------------------------------------------------------------

#' Plot predicted vs observed response boxplot
#'
#' Creates a combined violin + box + jitter plot comparing predicted viability
#' between response groups (e.g., Responders vs Non-Responders) with a
#' statistical significance bracket.
#'
#' @param exp_vs_pred Data frame with columns: response, predicted_viability.
#' @param response_var Character. Column name for response labels. Default = "response".
#' @param predicted_var Character. Column name for predicted values. Default = "predicted_viability".
#' @param y_label Character. Y-axis label. Default = "Predicted Viability (z-score)".
#' @param base_size Numeric. Base font size. Default = 15.
#' @param compare_method Character. Statistical test method. One of
#'        \code{"wilcox.test"} (default) or \code{"t.test"}.
#' @param alternative Character. Alternative hypothesis direction. Default = "greater".
#' @param tooltip Logical. If TRUE (default) and \pkg{ggiraph} is installed,
#'        jittered points get hover tooltips.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   exp_pred <- data.frame(
#'     response = factor(c("R", "NR", "R", "NR")),
#'     predicted_viability = c(0.8, 0.2, 0.7, 0.3)
#'   )
#'   plot_response_boxplot(exp_pred)
#' }
#'
#' @export
plot_response_boxplot <- function(exp_vs_pred,
                                  response_var = "response",
                                  predicted_var = "predicted_viability",
                                  y_label = "Predicted Viability (z-score)",
                                  base_size = 15,
                                  compare_method = "wilcox.test",
                                  alternative = "greater",
                                  tooltip = TRUE) {

  if (!all(c(response_var, predicted_var) %in% colnames(exp_vs_pred))) {
    stop("exp_vs_pred must contain columns: ", response_var, " and ", predicted_var)
  }

  # Ensure response is a factor preserving R/NR order (R first, as passed in).
  # Do NOT reverse levels — reversing would flip the x-order AND invert the
  # direction of the one-sided test below (e.g. test "NR > R" instead of "R > NR").
  resp_vec <- exp_vs_pred[[response_var]]
  if (is.factor(resp_vec)) {
    resp_vec <- factor(resp_vec)  # drop unused levels, keep order
  } else if (all(toupper(as.character(resp_vec[!is.na(resp_vec)])) %in% c("R", "NR"))) {
    resp_vec <- factor(resp_vec, levels = c("R", "NR"))
  } else {
    resp_vec <- factor(resp_vec)
  }
  exp_vs_pred[[response_var]] <- resp_vec

  # Group colors: responder = blue, non-responder = red
  grp_colors <- c("#2E86AB", "#C13232")
  if (length(levels(exp_vs_pred[[response_var]])) == 2) {
    names(grp_colors) <- levels(exp_vs_pred[[response_var]])
  }

  # Tooltip for ggiraph interactivity (hover shows group + predicted value).
  tt_ok <- tooltip && requireNamespace("ggiraph", quietly = TRUE)
  exp_vs_pred$tooltip_text <- sprintf("%s: %.3f",
                                      exp_vs_pred[[response_var]],
                                      exp_vs_pred[[predicted_var]])
  exp_vs_pred$data_id <- paste0("box_", seq_len(nrow(exp_vs_pred)))

  # Layer order & styling are deliberate (simplified, non-overlapping):
  #   1) violin = density body (semi-transparent, fills the background),
  #   2) box    = outline ONLY (transparent fill) so it never covers the
  #               violin's silhouette,
  #   3) jitter = raw points, spread across the violin width so the dots
  #               visually match the density shape.
  p <- ggplot(exp_vs_pred, aes(x = .data[[response_var]], y = .data[[predicted_var]])) +
    geom_violin(aes(fill = .data[[response_var]]), alpha = 0.22,
                width = 0.80, color = "grey55", linewidth = 0.30) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = NA,
                 color = "grey30", linewidth = 0.45) +
    { if (tt_ok) ggiraph::geom_jitter_interactive(
        aes(color = .data[[response_var]], tooltip = tooltip_text, data_id = data_id),
        width = 0.18, size = 1.5, alpha = 0.55)
      else geom_jitter(aes(color = .data[[response_var]]),
                       width = 0.18, size = 1.5, alpha = 0.55) } +
    theme_perception(base_size = base_size) +
    labs(y = y_label, x = "") +
    theme(legend.position = "none") +
    scale_fill_manual(values = grp_colors) +
    scale_color_manual(values = grp_colors)

  # Add significance bracket + p-value annotation
  if (length(levels(exp_vs_pred[[response_var]])) == 2) {
    lvls <- levels(exp_vs_pred[[response_var]])
    grp1 <- exp_vs_pred[[predicted_var]][exp_vs_pred[[response_var]] == lvls[1]]
    grp2 <- exp_vs_pred[[predicted_var]][exp_vs_pred[[response_var]] == lvls[2]]
    if (length(grp1) > 1 && length(grp2) > 1) {
      test_fun <- if (tolower(compare_method) == "t.test") t.test else wilcox.test
      wt <- tryCatch(test_fun(grp1, grp2, alternative = alternative),
                     error = function(e) NULL)
      if (!is.null(wt)) {
        p_label <- if (wt$p.value < 0.001) {
          "p < 0.001"
        } else {
          sprintf("p = %.3f", wt$p.value)
        }
        y_max <- max(exp_vs_pred[[predicted_var]], na.rm = TRUE)
        if (is.finite(y_max)) {
          y_pos <- y_max + (y_max - min(exp_vs_pred[[predicted_var]], na.rm = TRUE)) * 0.10
          p <- p +
            annotate("segment", x = 1, xend = 1, y = y_pos * 0.985, yend = y_pos,
                     color = "#1e2a4a", linewidth = 0.4) +
            annotate("segment", x = 2, xend = 2, y = y_pos * 0.985, yend = y_pos,
                     color = "#1e2a4a", linewidth = 0.4) +
            annotate("segment", x = 1, xend = 2, y = y_pos, yend = y_pos,
                     color = "#1e2a4a", linewidth = 0.4) +
            annotate("text", x = 1.5, y = y_pos * 1.03, label = p_label,
                     size = base_size * 0.3, hjust = 0.5, fontface = "italic",
                     color = "#1e2a4a") +
            coord_cartesian(ylim = c(NA, y_pos * 1.06))
        }
      }
    }
  }

  return(p)
}

# ---------------------------------------------------------------------------
# Plot: Model performance
# ---------------------------------------------------------------------------

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
#' @param tooltip Logical. If TRUE (default) and \pkg{ggiraph} is installed,
#'        points get hover tooltips (dataset, threshold, drug count).
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
                                   highlight_threshold = 0.3,
                                   tooltip = TRUE) {

  # Extract performance metrics with safety: skip entries lacking the expected list
  extract_col <- function(lst, what, col) {
    vals <- lapply(lst, function(x) {
      el <- x[[what]]
      if (is.null(el)) return(NULL)
      # Training output is a named vector (unlist(cor.test(...))), while the
      # built-in demo models use a one-row data.frame — support both.
      if (is.data.frame(el)) {
        if (!(col %in% colnames(el))) return(NULL)
        el[[col]]
      } else {
        if (is.null(names(el)) || !(col %in% names(el))) return(NULL)
        el[[col]]
      }
    })
    vals <- vals[!sapply(vals, is.null)]
    if (length(vals) == 0) stop("No performance metrics found in performance_list (missing ", what, ").")
    unlist(vals)
  }

  cor_scRNA <- extract_col(performance_list, "performance_in_scRNA", "estimate.cor")
  cor_bulk  <- extract_col(performance_list, "performance_in_bulk", "estimate.cor")
  cor_pseudo <- extract_col(performance_list, "performance_in_pseudo_bulk", "estimate.cor")

  # Build summary data frame
  df2plot <- rbind(
    data.frame(
      drugsCount = sapply(threshold_range, function(x) sum(cor_scRNA > x)),
      dataused = "scRNA-seq",
      Predictability = threshold_range
    ),
    data.frame(
      drugsCount = sapply(threshold_range, function(x) sum(cor_bulk > x)),
      dataused = "bulk",
      Predictability = threshold_range
    ),
    data.frame(
      drugsCount = sapply(threshold_range, function(x) sum(cor_pseudo > x)),
      dataused = "pseudo-bulk",
      Predictability = threshold_range
    )
  )
  df2plot$dataused <- factor(df2plot$dataused,
                             levels = c("scRNA-seq", "pseudo-bulk", "bulk"))

  dataset_colors <- c("scRNA-seq" = "#2E86AB", "pseudo-bulk" = "#27AE60",
                      "bulk" = "#E67E22")

  # Tooltip for ggiraph interactivity.
  tt_ok <- tooltip && requireNamespace("ggiraph", quietly = TRUE)
  df2plot$tooltip_text <- sprintf("Dataset: %s\nThreshold: %.2f\nDrugs: %d",
                                  df2plot$dataused, df2plot$Predictability,
                                  df2plot$drugsCount)
  df2plot$data_id <- paste0("perf_", seq_len(nrow(df2plot)))

  p <- ggplot(df2plot, aes(x = Predictability, y = drugsCount, color = dataused)) +
    geom_vline(xintercept = highlight_threshold, linetype = "dashed",
               color = "#C13232", alpha = 0.5, linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    { if (tt_ok) ggiraph::geom_point_interactive(
        aes(tooltip = tooltip_text, data_id = data_id), size = 1.8, alpha = 0.85)
      else geom_point(size = 1.8, alpha = 0.85) } +
    scale_color_manual(values = dataset_colors) +
    theme_perception(base_size = base_size) +
    labs(y = "Number of Drugs", color = "Validation Dataset",
         x = "Predictability (Pearson Correlation)") +
    theme(legend.position = "top",
          legend.box = "horizontal",
          legend.box.spacing = unit(8, "pt"))

  return(p)
}

# ---------------------------------------------------------------------------
# Plot: Seurat clustering
# ---------------------------------------------------------------------------

#' Run Seurat clustering and plot 2D embedding
#'
#' Performs Seurat clustering on an expression matrix and generates a 2D
#' embedding visualization (UMAP or t-SNE).
#' Useful for identifying subclones within patient tumor samples.
#'
#' @param method Character. Dimensionality reduction method. One of \code{"umap"}
#'        (default) or \code{"tsne"}.
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
#'   \item{embedding_plot}{ggplot 2D embedding visualization}
#'   \item{cluster_ids}{Named vector of cluster IDs per cell}
#'
#' @examples
#' \dontrun{
#'   result <- plot_seurat_clustering(patient_expression)
#'   result$embedding_plot
#'   result$cluster_ids
#' }
#'
#' @export
plot_seurat_clustering <- function(method = c("umap", "tsne"),
                                    expression_matrix,
                                    min_cells = 3,
                                    min_features = 200,
                                    nfeatures = 2000,
                                    dims = 10,
                                    resolution = 0.8,
                                    seed = 1) {

  method <- match.arg(method)

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Install with: install.packages('Seurat')")
  }

  set.seed(seed)

  # Create Seurat object
  so <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                   project = "PERCEPTIONx",
                                   min.cells = min_cells,
                                   min.features = min_features)

  # Standard workflow
  so <- Seurat::NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
  so <- Seurat::FindVariableFeatures(so, selection.method = "vst", nfeatures = nfeatures)
  so <- Seurat::ScaleData(so)
  so <- Seurat::RunPCA(so, features = Seurat::VariableFeatures(object = so))
  so <- Seurat::FindNeighbors(so, dims = 1:dims)
  so <- Seurat::FindClusters(so, resolution = resolution)

  # Dimensionality reduction
  if (method == "umap") {
    so <- Seurat::RunUMAP(so, dims = 1:dims)
  } else {
    so <- Seurat::RunTSNE(so, dims = 1:dims, check_duplicates = FALSE)
  }

  # Generate plot
  embedding_plot <- Seurat::DimPlot(so, reduction = method)

  # Extract cluster IDs
  cluster_ids <- Seurat::Idents(so)

  return(list(
    seurat_object = so,
    embedding_plot = embedding_plot,
    cluster_ids = cluster_ids
  ))
}

# ---------------------------------------------------------------------------
# Plot: Complete patient response panel
# ---------------------------------------------------------------------------

#' Complete patient response visualization pipeline
#'
#' Generates a comprehensive visualization panel for patient drug response prediction,
#' including clone distribution, clone-level viability, response boxplot, and ROC curve.
#' This is a convenience function that combines multiple plot functions.
#'
#' @param clone_distribution Data frame. Clone weights per patient.
#' @param clone_viability Data frame. Viability scores per clone.
#' @param exp_vs_pred Data frame. Predicted vs observed response.
#' @param response_col Character. Response column name. Default = "response".
#' @param viability_col Character. Viability column name. Default = "comb_viability".
#' @param predicted_col Character. Predicted values column name. Default = "predicted_viability".
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
#'     clone_viability = clone_kill_df,
#'     exp_vs_pred = response_df
#'   )
#'   ggsave(panel, filename = "patient_response.pdf", height = 15, width = 10)
#' }
#'
#' @export
plot_patient_response_panel <- function(clone_distribution,
                                        clone_viability,
                                        exp_vs_pred,
                                        response_col = "response",
                                        viability_col = "comb_viability",
                                        predicted_col = "predicted_viability",
                                        weights_col = "weights",
                                        layout_matrix = NULL) {

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required. Install with: install.packages('gridExtra')")
  }

  # Panel 1: Clone distribution
  p1 <- plot_clone_distribution(clone_distribution, response_var = response_col)

  # Panel 2: Clone-level viability (uses viability_col from clone_viability data)
  p2 <- plot_clone_viability(clone_viability, viability_var = viability_col,
                           weights_var = weights_col, response_var = response_col)

  # Panel 3: Response boxplot (uses predicted_col from exp_vs_pred data)
  p3 <- plot_response_boxplot(exp_vs_pred, response_var = response_col,
                              predicted_var = predicted_col)

  # Panel 4: ROC curve — auto-disable smoothing for small samples + safety net
  n_pts <- sum(!is.na(exp_vs_pred[[response_col]]) & !is.na(exp_vs_pred[[predicted_col]]))
  p4 <- tryCatch({
    plot_roc_curve(response = exp_vs_pred[[response_col]],
                   predictor = exp_vs_pred[[predicted_col]],
                   smooth_curve = n_pts >= 10)
  }, error = function(e) {
    # Last-resort fallback: render an empty placeholder ggplot with error message
    ggplot(data.frame(x = 0, y = 0, label = paste("ROC unavailable:", e$message))) +
      geom_text(aes(x, y, label = label), hjust = 0.5, vjust = 0.5, size = 3) +
      theme_void() +
      labs(title = "ROC Curve (unavailable)")
  })

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

# ---------------------------------------------------------------------------
# Plot: UMAP biomarker vs viability side-by-side
# ---------------------------------------------------------------------------

#' Plot UMAP side-by-side for biomarker and viability
#'
#' Creates a side-by-side comparison of biomarker expression and predicted viability
#' in UMAP space. Useful for visualizing correlation between marker and response.
#'
#' @param tsne_data Data frame with X, Y coordinates and both biomarker/viability columns.
#' @param biomarker_var Character. Column name for biomarker expression. Default = "biomarker_scaled".
#' @param viability_var Character. Column name for viability values. Default = "viability_scaled".
#' @param biomarker_label Character. Legend label for biomarker. Default = "Biomarker Exp".
#' @param viability_label Character. Legend label for viability. Default = "Drug Viability".
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
#'     viability_scaled = range01(rank(-viability_pred))
#'   )
#'   plot_tsne_biomarker_viability(tsne_data)
#' }
#'
#' @export
plot_tsne_biomarker_viability <- function(tsne_data,
                                        biomarker_var = "biomarker_scaled",
                                        viability_var = "viability_scaled",
                                        biomarker_label = "Biomarker Exp",
                                        viability_label = "Drug Viability",
                                        nrow = 1,
                                        base_size = 8) {

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required. Install with: install.packages('gridExtra')")
  }

  p1 <- plot_tsne_response(tsne_data, color_var = biomarker_var,
                           color_label = biomarker_label, base_size = base_size,
                           palette = "diverging")
  p2 <- plot_tsne_response(tsne_data, color_var = viability_var,
                           color_label = viability_label, base_size = base_size,
                           palette = "viridis")

  combined <- gridExtra::grid.arrange(p1, p2, nrow = nrow)

  return(combined)
}

# ---------------------------------------------------------------------------
# Plot: Cairo export
# ---------------------------------------------------------------------------

#' Export plot to file via Cairo device
#'
#' High-resolution plot export using the Cairo rendering engine for
#' superior anti-aliasing and cross-platform font rendering.
#'
#' @param file Character. Output file path.
#' @param plot A ggplot or grid object to export.
#' @param format Character. One of \code{"png"}, \code{"svg"}, or \code{"pdf"}.
#'        Default = \code{"png"}.
#' @param width Numeric. Plot width in inches. Default = 7.
#' @param height Numeric. Plot height in inches. Default = 5.
#' @param res Numeric. Output resolution (DPI) for PNG. Default = 600, minimum 96.
#' @param draw_fun Function. Optional custom draw function (e.g. \code{gridExtra::grid.arrange}
#'        for multi-panel layouts). If \code{NULL}, \code{print(plot)} is used.
#'
#' @return Invisibly returns the file path. Called for its side effect of creating the file.
#'
#' @examples
#' \dontrun{
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
#' export_plot_cairo("plot.png", p)
#' export_plot_cairo("plot.pdf", p, format = "pdf")
#' }
#'
#' @export
export_plot_cairo <- function(file, plot, format = "png",
                               width = 7, height = 5, res = 600,
                               draw_fun = NULL) {
  format <- tolower(format)
  has_cairo <- capabilities("cairo")
  if (format == "png") {
    if (has_cairo) {
      grDevices::png(file, width = width, height = height,
        units = "in", res = max(res, 96), type = "cairo")
    } else {
      grDevices::png(file, width = width, height = height,
        units = "in", res = max(res, 96))
    }
  } else if (format == "svg") {
    grDevices::svg(file, width = width, height = height, antialias = "default")
  } else {
    if (has_cairo) {
      grDevices::cairo_pdf(file, width = width, height = height)
    } else {
      grDevices::pdf(file, width = width, height = height)
    }
  }
  if (is.null(draw_fun)) {
    print(plot)
  } else {
    draw_fun(plot)
  }
  grDevices::dev.off()
  invisible(file)
}
