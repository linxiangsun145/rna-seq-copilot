# stability_scoring.R
# Dataset-level stability aggregation, scalar scoring, and confidence penalty mapping.

.cap01 <- function(x) pmax(0, pmin(1, x))

#' Detect dataset signal state from strong-DEG and weak-effect counts.
#' @param n_ref_deg Number of strong DEGs in reference analysis.
#' @param weak_effect_gene_count Number of weak-effect genes in reference analysis.
#' @param thresholds List of optional thresholds.
#' @return Character state: no_detectable_signal | weak_signal | strong_signal
detect_signal_state <- function(n_ref_deg,
                                weak_effect_gene_count = 0L,
                                thresholds = list()) {
  n_ref <- as.integer(n_ref_deg)
  n_weak <- as.integer(weak_effect_gene_count)

  if (is.na(n_ref) || n_ref <= 0) {
    if (!is.na(n_weak) && n_weak > 0) return("weak_signal")
    return("no_detectable_signal")
  }

  "strong_signal"
}

#' Determine whether effect-size diagnostics are mathematically computable.
#' @param dataset_metrics Dataset-level metric list.
#' @return Logical.
is_effect_metrics_computable <- function(dataset_metrics) {
  !is.na(dplyr::coalesce(dataset_metrics$mean_log2fc_correlation, NA_real_)) ||
    !is.na(dplyr::coalesce(dataset_metrics$median_log2fc_mad, NA_real_))
}

#' Determine stability mode from signal state and effect-metric availability.
#' @param signal_state Output of detect_signal_state().
#' @param effect_metrics_computable Whether effect metrics can be computed.
#' @return Character mode: deg_based | effect_only | not_applicable
determine_stability_mode <- function(signal_state, effect_metrics_computable = FALSE) {
  if (identical(signal_state, "strong_signal")) return("deg_based")
  if (identical(signal_state, "weak_signal") || identical(signal_state, "no_detectable_signal")) {
    if (isTRUE(effect_metrics_computable)) return("effect_only")
    return("not_applicable")
  }
  "not_applicable"
}

# Backward-compatible alias.
compute_stability_mode <- determine_stability_mode

#' Determine run status (completed | limited | failed).
#' @param signal_state Stability signal state.
#' @param stability_mode Stability mode.
#' @param engine_failed Whether perturbation engine failed.
#' @return Character run status.
determine_stability_run_status <- function(signal_state,
                                           stability_mode,
                                           engine_failed = FALSE) {
  if (isTRUE(engine_failed)) return("failed")
  if (identical(signal_state, "strong_signal")) return("completed")
  "limited"
}

#' Compute effect-size-only stability score for no-signal datasets.
#' @param dataset_metrics Dataset-level metric list.
#' @return List(score, components)
compute_effect_size_only_score <- function(dataset_metrics) {
  has_corr <- !is.na(dplyr::coalesce(dataset_metrics$mean_log2fc_correlation, NA_real_))
  has_var <- !is.na(dplyr::coalesce(dataset_metrics$median_log2fc_mad, NA_real_))
  if (!has_corr && !has_var) {
    return(list(
      score = NA_real_,
      components = list(
        lfc_corr_component = NA_real_,
        lfc_variability_component = NA_real_
      )
    ))
  }

  corr_scaled <- ifelse(
    is.na(dataset_metrics$mean_log2fc_correlation),
    0,
    .cap01((dataset_metrics$mean_log2fc_correlation + 1) / 2)
  )

  # Use MAD of log2FC perturbation as variability proxy.
  mad_val <- dplyr::coalesce(dataset_metrics$median_log2fc_mad, NA_real_)
  var_scaled <- ifelse(is.na(mad_val), 1, .cap01(mad_val / 1.5))

  score <- .cap01(0.6 * corr_scaled + 0.4 * (1 - var_scaled))

  list(
    score = score,
    components = list(
      lfc_corr_component = corr_scaled,
      lfc_variability_component = 1 - var_scaled
    )
  )
}

#' Compute DEG-set stability score (valid only when strong signal exists).
#' @param dataset_metrics Dataset-level metric list.
#' @return Numeric score in [0,1] or NA.
compute_deg_stability_score <- function(dataset_metrics) {
  if (!isTRUE(dataset_metrics$deg_metrics_applicable)) return(NA_real_)

  rec <- .cap01(dplyr::coalesce(dataset_metrics$mean_deg_recovery_rate, 0))
  jacc <- .cap01(dplyr::coalesce(dataset_metrics$mean_jaccard_deg, 0))
  top <- .cap01(dplyr::coalesce(dataset_metrics$mean_top_n_overlap, 0))
  collapse <- dplyr::coalesce(dataset_metrics$fraction_of_iterations_with_signal_collapse, 1)
  collapse_component <- 1 - .cap01(collapse)

  .cap01(0.45 * rec + 0.25 * jacc + 0.20 * top + 0.10 * collapse_component)
}

#' Compute effect-size consistency score (always defined when metrics available).
#' @param dataset_metrics Dataset-level metric list.
#' @return Numeric score in [0,1] or NA.
compute_effect_stability_score <- function(dataset_metrics) {
  eff <- compute_effect_size_only_score(dataset_metrics)
  dplyr::coalesce(eff$score, NA_real_)
}

#' Compute public final stability score.
#' @param signal_state Signal state string.
#' @param deg_stability_score DEG stability score.
#' @param effect_stability_score Effect stability score.
#' @return Numeric final score or NA when not publicly applicable.
compute_final_stability_score <- function(signal_state,
                                          deg_stability_score,
                                          effect_stability_score) {
  if (!identical(signal_state, "strong_signal")) {
    return(NA_real_)
  }
  .cap01(0.65 * dplyr::coalesce(deg_stability_score, 0) + 0.35 * dplyr::coalesce(effect_stability_score, 0))
}

#' Assign public-facing stability badge from signal and final score.
#' @param signal_state Signal state string.
#' @param final_stability_score Public final stability score.
#' @return Badge label.
assign_stability_badge <- function(signal_state,
                                   final_stability_score,
                                   stability_run_status = "limited") {
  if (identical(stability_run_status, "failed")) return("unknown")
  if (identical(signal_state, "no_detectable_signal")) return("low_signal")
  if (identical(signal_state, "weak_signal")) return("low_signal")

  score <- dplyr::coalesce(final_stability_score, 0)
  dplyr::case_when(
    score >= 0.75 ~ "high",
    score >= 0.50 ~ "moderate",
    TRUE ~ "low"
  )
}

#' Summarize dataset-level stability metrics from iteration metrics.
#' @param iteration_metrics Per-iteration metrics table.
#' @param reference_deg_count Number of strong DEGs in full-data reference.
#' @param weak_effect_gene_count Number of weak-effect genes in reference analysis.
#' @return List with aggregated metrics and interpretations.
summarize_dataset_stability <- function(iteration_metrics,
                                        reference_deg_count,
                                        weak_effect_gene_count = 0L) {
  if (!is.data.frame(iteration_metrics) || nrow(iteration_metrics) == 0) {
    stop("iteration_metrics must be a non-empty data.frame")
  }

  safe_mean <- function(x) ifelse(all(is.na(x)), NA_real_, mean(x, na.rm = TRUE))
  safe_median <- function(x) ifelse(all(is.na(x)), NA_real_, stats::median(x, na.rm = TRUE))
  safe_sd <- function(x) ifelse(sum(!is.na(x)) >= 2, stats::sd(x, na.rm = TRUE), NA_real_)

  mean_iter_deg <- safe_mean(iteration_metrics$n_deg_iter)
  sd_iter_deg <- safe_sd(iteration_metrics$n_deg_iter)
  effect_metrics_computable <- is_effect_metrics_computable(list(
    mean_log2fc_correlation = safe_mean(iteration_metrics$log2fc_correlation),
    median_log2fc_mad = safe_median(iteration_metrics$median_abs_delta_log2fc)
  ))
  signal_state <- detect_signal_state(
    n_ref_deg = reference_deg_count,
    weak_effect_gene_count = weak_effect_gene_count
  )
  stability_mode <- determine_stability_mode(signal_state, effect_metrics_computable)
  stability_run_status <- determine_stability_run_status(signal_state, stability_mode, engine_failed = FALSE)
  deg_metrics_applicable <- identical(signal_state, "strong_signal")
  signal_detectable <- deg_metrics_applicable

  out <- list(
    reference_deg_count = as.integer(reference_deg_count),
    weak_effect_gene_count = as.integer(weak_effect_gene_count),
    signal_state = signal_state,
    stability_run_status = stability_run_status,
    stability_mode = stability_mode,
    effect_metrics_computable = effect_metrics_computable,
    deg_metrics_applicable = deg_metrics_applicable,
    signal_detectable = signal_detectable,
    lack_of_detectable_signal = !signal_detectable,
    mean_iter_deg_count = mean_iter_deg,
    sd_iter_deg_count = sd_iter_deg,
    cv_iter_deg_count = ifelse(!is.na(mean_iter_deg) && mean_iter_deg > 0, sd_iter_deg / mean_iter_deg, NA_real_),
    mean_deg_recovery_rate = ifelse(deg_metrics_applicable, safe_mean(iteration_metrics$recovery_rate), NA_real_),
    median_deg_recovery_rate = ifelse(deg_metrics_applicable, safe_median(iteration_metrics$recovery_rate), NA_real_),
    mean_jaccard_deg = ifelse(deg_metrics_applicable, safe_mean(iteration_metrics$jaccard_deg), NA_real_),
    mean_top_n_overlap = safe_mean(iteration_metrics$top_n_overlap),
    mean_up_top_n_overlap = safe_mean(iteration_metrics$up_top_n_overlap),
    mean_down_top_n_overlap = safe_mean(iteration_metrics$down_top_n_overlap),
    mean_log2fc_correlation = safe_mean(iteration_metrics$log2fc_correlation),
    mean_log2fc_rmse = safe_mean(iteration_metrics$log2fc_rmse),
    median_log2fc_mad = safe_median(iteration_metrics$median_abs_delta_log2fc),
    fraction_of_iterations_with_signal_collapse = ifelse(deg_metrics_applicable, safe_mean(as.numeric(iteration_metrics$signal_collapse)), NA_real_),
    worst_case_top_n_overlap = ifelse(all(is.na(iteration_metrics$top_n_overlap)), NA_real_, min(iteration_metrics$top_n_overlap, na.rm = TRUE)),
    worst_case_deg_recovery = ifelse(deg_metrics_applicable && !all(is.na(iteration_metrics$recovery_rate)), min(iteration_metrics$recovery_rate, na.rm = TRUE), NA_real_)
  )

  out$metric_interpretation <- list(
    deg_recovery = ifelse(
      !deg_metrics_applicable || is.na(out$mean_deg_recovery_rate),
      "Reference DEG recovery is not defined because no strong reference DEG set is available.",
      sprintf("On average, %.1f%% of reference strong DEGs are recovered under perturbation.", 100 * out$mean_deg_recovery_rate)
    ),
    top_overlap = sprintf("Top-ranked genes (statistical ranking, not necessarily significant DEGs) overlap by %.1f%% across perturbations.", 100 * dplyr::coalesce(out$mean_top_n_overlap, 0)),
    effect_consistency = sprintf(
      "Mean log2FC correlation is %.2f and median absolute log2FC shift is %.3f.",
      dplyr::coalesce(out$mean_log2fc_correlation, NA_real_),
      dplyr::coalesce(out$median_log2fc_mad, NA_real_)
    ),
    collapse = ifelse(
      deg_metrics_applicable,
      sprintf("Signal collapse occurs in %.1f%% of perturbation runs.", 100 * dplyr::coalesce(out$fraction_of_iterations_with_signal_collapse, 0)),
      "Signal collapse is not applicable because no detectable reference DEG signal was present."
    )
  )

  out
}

#' Compute transparent scalar stability score in [0,1].
#' @param dataset_metrics Output list from summarize_dataset_stability.
#' @return List with score and component breakdown.
compute_stability_score <- function(dataset_metrics) {
  signal_state <- dataset_metrics$signal_state
  stability_mode <- determine_stability_mode(
    signal_state,
    effect_metrics_computable = isTRUE(dataset_metrics$effect_metrics_computable)
  )
  deg_stability_score <- compute_deg_stability_score(dataset_metrics)
  effect_stability_score <- compute_effect_stability_score(dataset_metrics)
  final_stability_score <- compute_final_stability_score(
    signal_state = signal_state,
    deg_stability_score = deg_stability_score,
    effect_stability_score = effect_stability_score
  )
  run_status <- determine_stability_run_status(signal_state, stability_mode, engine_failed = FALSE)
  badge <- assign_stability_badge(signal_state, final_stability_score, run_status)

  level <- dplyr::case_when(
    badge == "high" ~ "high",
    badge == "moderate" ~ "moderate",
    badge == "low" ~ "low",
    badge == "low_signal" ~ "low_signal",
    TRUE ~ "not_applicable"
  )

  components <- list(
    deg_stability_score = deg_stability_score,
    effect_stability_score = effect_stability_score,
    final_stability_score = final_stability_score
  )

  formula <- if (identical(stability_mode, "deg_based")) {
    "final = 0.65*deg_stability_score + 0.35*effect_stability_score; deg_stability_score uses DEG recovery/Jaccard/top-overlap/collapse"
  } else {
    "effect_stability_score = 0.6*LFC_correlation + 0.4*(1-normalized_log2FC_variability); final_stability_score is not publicly applicable"
  }

  list(
    # Public-facing overall score (N/A when no strong DEG signal).
    stability_score = final_stability_score,
    deg_stability_score = deg_stability_score,
    effect_stability_score = effect_stability_score,
    final_stability_score = final_stability_score,
    stability_level = level,
    stability_badge = badge,
    stability_run_status = run_status,
    stability_mode = stability_mode,
    formula = formula,
    components = components
  )
}

#' Compute confidence-penalty output from stability results.
#' @param stability_results Full output list from run_stability_analysis.
#' @return List with penalty value and explanation.
compute_stability_penalty <- function(stability_results) {
  score <- stability_results$dataset_stability$final_stability_score
  signal_state <- as.character(stability_results$dataset_stability$signal_state)
  if (is.null(signal_state) || is.na(signal_state)) signal_state <- ""

  if (identical(signal_state, "no_detectable_signal")) {
    return(list(
      stability_penalty = 20,
      penalty_level = "moderate",
      penalty_reason = "No robust differential signal was detectable; stability is evaluated using effect-size consistency only."
    ))
  }

  if (identical(signal_state, "weak_signal")) {
    return(list(
      stability_penalty = 14,
      penalty_level = "moderate",
      penalty_reason = "Only weak differential signal was detected; DEG-based stability is limited under current thresholds."
    ))
  }

  score <- dplyr::coalesce(score, 0)

  penalty <- dplyr::case_when(
    score >= 0.75 ~ 0,
    score >= 0.50 ~ 8,
    score >= 0.30 ~ 18,
    TRUE ~ 30
  )

  level <- dplyr::case_when(
    penalty == 0 ~ "none",
    penalty <= 10 ~ "mild",
    penalty <= 20 ~ "moderate",
    TRUE ~ "severe"
  )

  reason <- dplyr::case_when(
    level == "none" ~ "Stability analysis indicates robust DEG signal under perturbation.",
    level == "mild" ~ "Stability analysis indicates moderate robustness with some perturbation sensitivity.",
    level == "moderate" ~ "Stability analysis shows notable perturbation sensitivity in DEG recovery or ranking.",
    TRUE ~ "Stability analysis indicates substantial instability and possible sample-driven signal."
  )

  list(
    stability_penalty = penalty,
    penalty_level = level,
    penalty_reason = reason
  )
}
