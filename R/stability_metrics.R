# stability_metrics.R
# Metrics for comparing perturbation DESeq2 results against reference analysis.

#' Assign DEG significance class for each gene.
#' @param res_tbl DESeq2 results table with gene_id, padj, log2FoldChange.
#' @param padj_threshold Adjusted p-value threshold.
#' @param log2fc_threshold Absolute log2FC threshold for strong DEG.
#' @return tibble with ref_class/significant flags.
classify_significance <- function(res_tbl,
                                  padj_threshold = 0.05,
                                  log2fc_threshold = 1.0) {
  res_tbl |>
    dplyr::mutate(
      padj = suppressWarnings(as.numeric(.data$padj)),
      log2FoldChange = suppressWarnings(as.numeric(.data$log2FoldChange)),
      ref_class = dplyr::case_when(
        !is.na(.data$padj) & !is.na(.data$log2FoldChange) &
          .data$padj < padj_threshold & abs(.data$log2FoldChange) >= log2fc_threshold ~ "strong_DEG",
        !is.na(.data$padj) & !is.na(.data$log2FoldChange) &
          .data$padj < padj_threshold & abs(.data$log2FoldChange) < log2fc_threshold ~ "weak_effect",
        TRUE ~ "not_significant"
      ),
      is_strong_deg = .data$ref_class == "strong_DEG",
      is_weak_effect = .data$ref_class == "weak_effect",
      lfc_sign = dplyr::case_when(
        is.na(.data$log2FoldChange) ~ NA_integer_,
        .data$log2FoldChange > 0 ~ 1L,
        .data$log2FoldChange < 0 ~ -1L,
        TRUE ~ 0L
      )
    )
}

#' Build ranked top-gene vectors from a DE table.
#' @param res_tbl Classified DE table.
#' @param top_n Number of top genes.
#' @return list with top_all/top_up/top_down vectors.
rank_top_genes <- function(res_tbl, top_n = 50) {
  ordered <- res_tbl |>
    dplyr::mutate(
      padj_rank = dplyr::if_else(is.na(.data$padj), Inf, .data$padj),
      abs_lfc = abs(dplyr::coalesce(.data$log2FoldChange, 0))
    ) |>
    dplyr::arrange(.data$padj_rank, dplyr::desc(.data$abs_lfc))

  top_all <- head(ordered$gene_id, top_n)

  up_tbl <- ordered |>
    dplyr::filter(!is.na(.data$log2FoldChange), .data$log2FoldChange > 0)
  down_tbl <- ordered |>
    dplyr::filter(!is.na(.data$log2FoldChange), .data$log2FoldChange < 0)

  list(
    top_all = top_all,
    top_up = head(up_tbl$gene_id, top_n),
    top_down = head(down_tbl$gene_id, top_n),
    top_combined = top_all
  )
}

.safe_cor <- function(x, y) {
  idx <- which(!is.na(x) & !is.na(y))
  if (length(idx) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx], method = "pearson"))
}

.safe_rmse <- function(x, y) {
  idx <- which(!is.na(x) & !is.na(y))
  if (length(idx) < 1) return(NA_real_)
  sqrt(mean((x[idx] - y[idx])^2))
}

#' Compare one perturbation iteration to the full-data reference.
#' @param iter_tbl DE table from one iteration.
#' @param ref_obj Reference object from run_reference_deseq.
#' @param top_n Number of top genes used in overlap metrics.
#' @param padj_threshold padj threshold.
#' @param log2fc_threshold log2FC threshold.
#' @return list(iter_metrics, gene_delta)
compare_iteration_to_reference <- function(iter_tbl,
                                           ref_obj,
                                           top_n = 50,
                                           padj_threshold = 0.05,
                                           log2fc_threshold = 1.0) {
  iter_class <- classify_significance(iter_tbl, padj_threshold, log2fc_threshold)

  iter_deg <- iter_class |>
    dplyr::filter(.data$is_strong_deg) |>
    dplyr::pull(.data$gene_id)

  ref_deg <- ref_obj$reference_deg_set
  n_ref_deg <- length(ref_deg)
  deg_metrics_applicable <- n_ref_deg > 0
  signal_detectable <- deg_metrics_applicable
  lack_of_detectable_signal <- !signal_detectable

  overlap <- intersect(ref_deg, iter_deg)
  union_set <- union(ref_deg, iter_deg)

  iter_top <- rank_top_genes(iter_class, top_n = top_n)
  ref_top <- ref_obj$reference_top

  merged <- ref_obj$reference_table |>
    dplyr::select(
      .data$gene_id,
      ref_log2FC = .data$log2FoldChange,
      ref_padj = .data$padj,
      ref_class = .data$ref_class,
      ref_sign = .data$lfc_sign,
      ref_is_deg = .data$is_strong_deg
    ) |>
    dplyr::left_join(
      iter_class |>
        dplyr::select(
          .data$gene_id,
          iter_log2FC = .data$log2FoldChange,
          iter_padj = .data$padj,
          iter_class = .data$ref_class,
          iter_sign = .data$lfc_sign,
          iter_is_deg = .data$is_strong_deg
        ),
      by = "gene_id"
    ) |>
    dplyr::mutate(
      iter_tested = !is.na(.data$iter_log2FC) | !is.na(.data$iter_padj),
      lfc_delta = .data$iter_log2FC - .data$ref_log2FC,
      abs_lfc_delta = abs(.data$lfc_delta),
      sign_flip = dplyr::if_else(
        .data$ref_is_deg & .data$iter_tested & !is.na(.data$ref_sign) & !is.na(.data$iter_sign),
        .data$ref_sign != .data$iter_sign,
        FALSE,
        FALSE
      )
    )

  n_deg_iter <- length(iter_deg)
  n_ref_recovered <- if (deg_metrics_applicable) length(overlap) else NA_integer_
  collapse_cutoff <- max(0L, floor(0.2 * n_ref_deg))
  signal_collapse_flag <- deg_metrics_applicable && (n_deg_iter <= collapse_cutoff)

  iter_metrics <- tibble::tibble(
    n_deg_iter = n_deg_iter,
    n_ref_deg = n_ref_deg,
    deg_metrics_applicable = deg_metrics_applicable,
    signal_detectable = signal_detectable,
    lack_of_detectable_signal = lack_of_detectable_signal,
    n_ref_deg_recovered = n_ref_recovered,
    recovery_rate = ifelse(deg_metrics_applicable, n_ref_recovered / n_ref_deg, NA_real_),
    jaccard_deg = ifelse(deg_metrics_applicable && length(union_set) > 0, length(overlap) / length(union_set), NA_real_),
    top_n_overlap = length(intersect(ref_top$top_combined, iter_top$top_combined)) / max(1, top_n),
    up_top_n_overlap = length(intersect(ref_top$top_up, iter_top$top_up)) / max(1, top_n),
    down_top_n_overlap = length(intersect(ref_top$top_down, iter_top$top_down)) / max(1, top_n),
    log2fc_correlation = .safe_cor(merged$ref_log2FC, merged$iter_log2FC),
    log2fc_rmse = .safe_rmse(merged$ref_log2FC, merged$iter_log2FC),
    median_abs_delta_log2fc = stats::median(merged$abs_lfc_delta, na.rm = TRUE),
    sign_flip_count = sum(merged$sign_flip, na.rm = TRUE),
    signal_collapse = signal_collapse_flag
  )

  list(
    iter_metrics = iter_metrics,
    gene_delta = merged,
    iter_classified = iter_class,
    iter_top = iter_top
  )
}

#' Aggregate gene-level stability metrics across perturbation iterations.
#' @param reference_table Classified full-data DE table.
#' @param gene_delta_list List of gene-level merged tables from iterations.
#' @param top_gene_hits_list Character vectors of top genes from each iteration.
#' @param padj_threshold padj threshold.
#' @param log2fc_threshold log2FC threshold.
#' @return tibble with per-gene stability metrics and labels.
summarize_gene_stability <- function(reference_table,
                                     gene_delta_list,
                                     top_gene_hits_list,
                                     padj_threshold = 0.05,
                                     log2fc_threshold = 1.0) {
  all_deltas <- dplyr::bind_rows(gene_delta_list, .id = "iteration_index")

  top_hits <- tibble::tibble(
    iteration_index = as.character(seq_along(top_gene_hits_list)),
    gene_id = top_gene_hits_list
  ) |>
    tidyr::unnest(cols = c(.data$gene_id)) |>
    dplyr::mutate(in_top_n = TRUE)

  by_gene <- all_deltas |>
    dplyr::left_join(top_hits, by = c("iteration_index", "gene_id")) |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::summarise(
      ref_log2FC = dplyr::first(.data$ref_log2FC),
      ref_padj = dplyr::first(.data$ref_padj),
      ref_class = dplyr::first(.data$ref_class),
      times_tested = sum(.data$iter_tested, na.rm = TRUE),
      times_significant = sum(.data$iter_is_deg, na.rm = TRUE),
      significance_frequency = dplyr::if_else(.data$times_tested > 0, .data$times_significant / .data$times_tested, NA_real_),
      mean_log2FC = mean(.data$iter_log2FC, na.rm = TRUE),
      sd_log2FC = stats::sd(.data$iter_log2FC, na.rm = TRUE),
      log2FC_mad = stats::mad(.data$iter_log2FC, center = stats::median(.data$iter_log2FC, na.rm = TRUE), na.rm = TRUE),
      cv_like_log2FC = dplyr::if_else(abs(.data$mean_log2FC) > 1e-8, .data$sd_log2FC / abs(.data$mean_log2FC), NA_real_),
      sign_consistency = mean(.data$iter_sign == .data$ref_sign, na.rm = TRUE),
      mean_padj = mean(.data$iter_padj, na.rm = TRUE),
      recovered_in_top_n_frequency = mean(dplyr::coalesce(.data$in_top_n, FALSE), na.rm = TRUE),
      perturb_only_significant = dplyr::first(.data$ref_class) == "not_significant" & sum(.data$iter_is_deg, na.rm = TRUE) > 0,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      significance_frequency = dplyr::coalesce(.data$significance_frequency, 0),
      sign_consistency = dplyr::coalesce(.data$sign_consistency, 0),
      log2FC_mad = dplyr::coalesce(.data$log2FC_mad, Inf),
      stability_label = dplyr::case_when(
        .data$sign_consistency < 0.8 ~ "effect_direction_unstable",
        .data$significance_frequency >= 0.8 & .data$sign_consistency >= 0.9 & .data$log2FC_mad <= 0.35 ~ "highly_stable",
        .data$significance_frequency >= 0.5 ~ "moderately_stable",
        TRUE ~ "unstable"
      )
    )

  # Keep reference genes first, but keep perturbation-only candidates too.
  ref_ids <- reference_table$gene_id
  by_gene |>
    dplyr::mutate(ref_priority = dplyr::if_else(.data$gene_id %in% ref_ids, 0L, 1L)) |>
    dplyr::arrange(.data$ref_priority, dplyr::desc(.data$significance_frequency), .data$gene_id) |>
    dplyr::select(-.data$ref_priority)
}

#' Build per-sample influence table from leave-one-out iterations.
#' @param iteration_metrics Table with per-iteration metrics and metadata.
#' @return tibble per removed sample with influence score.
compute_sample_influence <- function(iteration_metrics, deg_metrics_applicable = TRUE) {
  loo <- iteration_metrics |>
    dplyr::filter(.data$mode == "loo", !is.na(.data$removed_sample))

  if (nrow(loo) == 0) {
    return(tibble::tibble(
      removed_sample = character(),
      iter_deg_count = integer(),
      deg_recovery_rate = numeric(),
      top_n_overlap = numeric(),
      log2fc_correlation = numeric(),
      collapse_flag = logical(),
      influence_mode = character(),
      influence_score = numeric(),
      influence_note = character()
    ))
  }

  deg_norm <- loo$n_deg_iter / max(1, max(loo$n_deg_iter, na.rm = TRUE))
  corr_norm <- pmax(-1, pmin(1, loo$log2fc_correlation))
  corr_norm <- (corr_norm + 1) / 2

  rec_weight <- if (isTRUE(deg_metrics_applicable)) 0.35 else 0.00
  top_weight <- if (isTRUE(deg_metrics_applicable)) 0.30 else 0.00
  corr_weight <- if (isTRUE(deg_metrics_applicable)) 0.20 else 0.00
  deg_weight <- if (isTRUE(deg_metrics_applicable)) 0.15 else 0.00
  rec_component <- if (isTRUE(deg_metrics_applicable)) dplyr::coalesce(loo$recovery_rate, 0) else 0

  influence <- loo |>
    dplyr::mutate(
      collapse_flag = .data$signal_collapse,
      influence_mode = ifelse(isTRUE(deg_metrics_applicable), "deg_based", "effect_or_rank_based"),
      influence_score = ifelse(
        isTRUE(deg_metrics_applicable),
        1 - (
          rec_weight * rec_component +
            top_weight * dplyr::coalesce(.data$top_n_overlap, 0) +
            corr_weight * dplyr::coalesce(corr_norm, 0) +
            deg_weight * dplyr::coalesce(deg_norm, 0)
        ),
        compute_rank_based_sample_influence(loo)
      ),
      influence_score = pmax(0, pmin(1, .data$influence_score)),
      influence_note = dplyr::case_when(
        !isTRUE(deg_metrics_applicable) ~ "Reference DEG signal is absent; influence is estimated from ranking and effect-size consistency only.",
        .data$collapse_flag ~ "Removing this sample collapses detected DEG signal.",
        .data$influence_score >= 0.7 ~ "Removing this sample substantially changes DEG signal.",
        .data$influence_score >= 0.45 ~ "Removing this sample moderately affects DEG signal.",
        TRUE ~ "Removing this sample has limited impact on DEG signal."
      )
    ) |>
    dplyr::transmute(
      removed_sample = .data$removed_sample,
      iter_deg_count = as.integer(.data$n_deg_iter),
      deg_recovery_rate = .data$recovery_rate,
      top_n_overlap = .data$top_n_overlap,
      log2fc_correlation = .data$log2fc_correlation,
      collapse_flag = .data$collapse_flag,
      influence_mode = .data$influence_mode,
      influence_score = .data$influence_score,
      influence_note = .data$influence_note
    ) |>
    dplyr::arrange(dplyr::desc(.data$influence_score), .data$removed_sample)

  influence
}

#' Compute sample influence score from ranking/effect metrics only.
#' @param tbl Leave-one-out iteration table.
#' @return Numeric vector in [0,1] where higher means more influential.
compute_rank_based_sample_influence <- function(tbl) {
  corr_norm_local <- pmax(-1, pmin(1, tbl$log2fc_correlation))
  corr_norm_local <- (corr_norm_local + 1) / 2
  top_component <- dplyr::coalesce(tbl$top_n_overlap, 0)
  corr_component <- dplyr::coalesce(corr_norm_local, 0)
  # In low/no-signal settings, focus on ranking and effect-size consistency only.
  score_local <- 1 - (0.60 * top_component + 0.40 * corr_component)
  pmax(0, pmin(1, score_local))
}
