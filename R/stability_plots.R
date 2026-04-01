# stability_plots.R
# ggplot2 plotting helpers for stability analysis report integration.

#' Plot DEG count distribution across perturbations.
plot_stability_deg_counts <- function(iteration_metrics) {
  ggplot2::ggplot(iteration_metrics, ggplot2::aes(x = .data$n_deg_iter, fill = .data$mode)) +
    ggplot2::geom_histogram(bins = 20, alpha = 0.7, color = "white") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = "DEG Count Across Perturbation Runs",
      x = "Number of strong DEGs per iteration",
      y = "Frequency",
      fill = "Mode"
    )
}

#' Plot DEG recovery-rate distribution.
plot_stability_recovery <- function(iteration_metrics) {
  ggplot2::ggplot(iteration_metrics, ggplot2::aes(x = .data$recovery_rate, fill = .data$mode)) +
    ggplot2::geom_histogram(bins = 20, alpha = 0.7, color = "white") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = "Reference DEG Recovery Across Perturbations",
      x = "Recovery rate",
      y = "Frequency",
      fill = "Mode"
    )
}

#' Plot top-N overlap distribution.
plot_stability_top_overlap <- function(iteration_metrics) {
  ggplot2::ggplot(iteration_metrics, ggplot2::aes(x = .data$top_n_overlap, fill = .data$mode)) +
    ggplot2::geom_histogram(bins = 20, alpha = 0.7, color = "white") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = "Top-N DEG Overlap Across Perturbations",
      x = "Top-N overlap",
      y = "Frequency",
      fill = "Mode"
    )
}

#' Plot per-sample influence scores from LOO runs.
plot_sample_influence <- function(sample_influence) {
  ggplot2::ggplot(sample_influence, ggplot2::aes(
    x = stats::reorder(.data$removed_sample, .data$influence_score),
    y = .data$influence_score,
    fill = .data$collapse_flag
  )) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#4C78A8", "TRUE" = "#E45756")) +
    ggplot2::labs(
      title = "Single-Sample Influence (Leave-One-Out)",
      x = "Removed sample",
      y = "Influence score (higher = more influential)",
      fill = "Signal collapse"
    )
}

#' Plot gene-level stability scatter.
plot_gene_stability <- function(gene_stability, highlight_n = 20) {
  to_label <- gene_stability |>
    dplyr::arrange(dplyr::desc(.data$significance_frequency), .data$log2FC_mad) |>
    dplyr::slice_head(n = highlight_n)

  ggplot2::ggplot(gene_stability, ggplot2::aes(
    x = .data$significance_frequency,
    y = .data$log2FC_mad,
    color = .data$stability_label
  )) +
    ggplot2::geom_point(alpha = 0.7) +
    ggrepel::geom_text_repel(
      data = to_label,
      ggplot2::aes(label = .data$gene_id),
      size = 3,
      max.overlaps = 50
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = "Gene-Level Stability Landscape",
      x = "Significance frequency across perturbations",
      y = "log2FC MAD across perturbations",
      color = "Stability label"
    )
}

#' Optional heatmap for top genes across iterations (log2FC values).
plot_stability_heatmap <- function(gene_delta_long,
                                   gene_ids,
                                   value_col = "iter_log2FC") {
  sub <- gene_delta_long |>
    dplyr::filter(.data$gene_id %in% gene_ids) |>
    dplyr::select(.data$iteration_index, .data$gene_id, value = dplyr::all_of(value_col))

  wide <- tidyr::pivot_wider(sub, names_from = .data$iteration_index, values_from = .data$value)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$gene_id

  pheatmap::pheatmap(
    mat,
    color = colorRampPalette(c("#313695", "#FFFFFF", "#A50026"))(100),
    main = "Top-Gene log2FC Stability Heatmap",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    silent = TRUE
  )
}
