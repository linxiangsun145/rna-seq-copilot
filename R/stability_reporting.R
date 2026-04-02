# stability_reporting.R
# Warning generation and publication-style narrative reporting.

#' Generate structured warning codes from stability metrics.
#' @param dataset_stability Dataset-level stability list.
#' @param gene_stability Gene-level stability table.
#' @param sample_influence Per-sample influence table.
#' @return Character vector of warning codes.
generate_stability_warnings <- function(dataset_stability,
                                        gene_stability,
                                        sample_influence,
                                        pca_separation = "unknown",
                                        qc_context = list()) {
  warnings <- character()

  get_scalar_numeric <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    suppressWarnings(as.numeric(x[[1]]))
  }

  get_scalar_text <- function(x, default = "unknown") {
    if (is.null(x) || length(x) == 0) return(default)
    val <- as.character(x[[1]])
    if (is.na(val) || !nzchar(val)) return(default)
    val
  }

  signal_state <- as.character(dataset_stability$signal_state)
  stability_run_status <- as.character(dataset_stability$stability_run_status)
  deg_metrics_applicable <- isTRUE(dataset_stability$deg_metrics_applicable)
  final_score <- dataset_stability$final_stability_score
  ref_deg_count <- suppressWarnings(as.integer(dataset_stability$reference_deg_count))

  # Technical failure warnings are reserved for true engine/runtime failures.
  if (identical(stability_run_status, "failed")) {
    return("STABILITY_ANALYSIS_FAILED")
  }

  if (identical(signal_state, "no_detectable_signal")) {
    warnings <- c(warnings, "LACK_OF_DETECTABLE_SIGNAL")
  }
  if (!is.na(ref_deg_count) && ref_deg_count == 0) {
    warnings <- c(warnings, "LACK_OF_DETECTABLE_STRONG_SIGNAL")
  }

  score <- dplyr::coalesce(final_score, 0)
  if (score < 0.4 && identical(signal_state, "strong_signal")) warnings <- c(warnings, "STABILITY_LOW")

  top_overlap <- dplyr::coalesce(dataset_stability$mean_top_n_overlap, 0)
  if (top_overlap < 0.4) warnings <- c(warnings, "TOP_GENES_UNSTABLE")

  ref_deg <- gene_stability[gene_stability$ref_class == "strong_DEG", , drop = FALSE]

  median_sign_consistency <- if (nrow(ref_deg) > 0) {
    stats::median(ref_deg$sign_consistency, na.rm = TRUE)
  } else {
    NA_real_
  }

  if (!is.na(median_sign_consistency) && median_sign_consistency < 0.8) {
    warnings <- c(warnings, "EFFECT_SIZE_UNSTABLE")
  }

  collapse_frac <- dplyr::coalesce(dataset_stability$fraction_of_iterations_with_signal_collapse, 0)
  if (deg_metrics_applicable && collapse_frac >= 0.25) warnings <- c(warnings, "SIGNAL_COLLAPSE_RISK")

  mean_corr <- get_scalar_numeric(qc_context$mean_within_group_correlation)
  if (is.na(mean_corr)) mean_corr <- get_scalar_numeric(qc_context$mean_correlation)
  group_correlation_low <- !is.na(mean_corr) && mean_corr < 0.75
  group_correlation_very_low <- !is.na(mean_corr) && mean_corr <= 0.30
  group_correlation_negative <- !is.na(mean_corr) && mean_corr < 0
  if (tolower(as.character(pca_separation)) == "clear" && group_correlation_low) {
    warnings <- c(warnings, "PCA_QC_CONFLICT")
  }
  if (group_correlation_very_low || group_correlation_negative) {
    warnings <- c(warnings, "WITHIN_GROUP_INCONSISTENCY")
  }

  qc_status <- tolower(get_scalar_text(qc_context$qc_status, default = "unknown"))
  qc_high_risk <- qc_status %in% c("critical", "high", "high_risk")
  if ((qc_high_risk && group_correlation_low) ||
      ((group_correlation_very_low || group_correlation_negative) && !identical(signal_state, "strong_signal")) ||
      (!identical(signal_state, "strong_signal") && group_correlation_low)) {
    warnings <- c(warnings, "QC_DOMINATES_SIGNAL")
  }

  if (nrow(sample_influence) > 0) {
    if (max(sample_influence$influence_score, na.rm = TRUE) >= 0.7) {
      warnings <- c(warnings, "SINGLE_SAMPLE_DEPENDENCY")
    }
  }

  unique(warnings)
}

#' Build compact headline and cautious scientific interpretation.
#' @param dataset_stability Dataset-level stability list.
#' @param warnings Warning code vector.
#' @param sample_influence Per-sample influence table.
#' @return List with headline, badge, findings and narrative paragraphs.
interpret_stability_results <- function(dataset_stability,
                                        warnings,
                                        sample_influence,
                                        pca_separation = "unknown",
                                        qc_context = list()) {
  get_scalar_numeric <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    suppressWarnings(as.numeric(x[[1]]))
  }

  signal_state <- as.character(dataset_stability$signal_state)
  stability_run_status <- as.character(dataset_stability$stability_run_status)
  if ((identical(signal_state, "no_detectable_signal") || identical(signal_state, "weak_signal")) &&
      !("STABILITY_ANALYSIS_FAILED" %in% warnings)) {
    stability_run_status <- "limited"
  }
  stability_mode <- as.character(dataset_stability$stability_mode)
  deg_metrics_applicable <- isTRUE(dataset_stability$deg_metrics_applicable)
  effect_metrics_computable <- isTRUE(dataset_stability$effect_metrics_computable)
  badge <- if (identical(stability_run_status, "failed")) {
    "unknown"
  } else if (identical(signal_state, "no_detectable_signal")) {
    "low_signal"
  } else if (identical(signal_state, "weak_signal")) {
    "low_signal"
  } else {
    dplyr::case_when(
      dplyr::coalesce(dataset_stability$final_stability_score, 0) >= 0.75 ~ "high",
      dplyr::coalesce(dataset_stability$final_stability_score, 0) >= 0.50 ~ "moderate",
      TRUE ~ "low"
    )
  }

  headline <- if (identical(stability_run_status, "failed")) {
    "Stability analysis could not be completed due to a technical failure in the perturbation assessment module."
  } else if (identical(signal_state, "no_detectable_signal")) {
    "No robust differential expression signal was detectable under the configured thresholds."
  } else if (identical(signal_state, "weak_signal")) {
    "No robust differential expression (DEG) signal was detected under the configured thresholds; only weak effect-level signal was observed."
  } else {
    dplyr::case_when(
      badge == "high" ~ "DEG results appear stable under sample perturbation.",
      badge == "moderate" ~ "DEG results are moderately stable, with measurable perturbation sensitivity.",
      TRUE ~ "DEG results appear sensitive to perturbation and should be interpreted with caution."
    )
  }

  key_findings <- c()
  if (deg_metrics_applicable) {
    key_findings <- c(
      key_findings,
      sprintf(
        "Mean reference DEG recovery rate across perturbations: %.1f%%.",
        100 * dplyr::coalesce(dataset_stability$mean_deg_recovery_rate, 0)
      )
    )
  } else {
    key_findings <- c(
      key_findings,
      "Reference strong-DEG set is empty; DEG recovery and DEG Jaccard metrics are not applicable."
    )
  }
  key_findings <- c(
    key_findings,
    sprintf(
      "Top-ranked gene overlap across perturbations: %.1f%% (ranking-based, not threshold-based DEG overlap).",
      100 * dplyr::coalesce(dataset_stability$mean_top_n_overlap, 0)
    ),
    sprintf(
      "Mean log2 fold-change correlation: %.2f.",
      dplyr::coalesce(dataset_stability$mean_log2fc_correlation, NA_real_)
    )
  )

  if (!is.na(dplyr::coalesce(dataset_stability$effect_stability_score, NA_real_)) &&
      !is.na(dplyr::coalesce(dataset_stability$mean_log2fc_correlation, NA_real_))) {
    key_findings <- c(
      key_findings,
      "Composite effect-size consistency is derived from multiple perturbation-based measures (including log2 fold-change correlation and variability) and is not identical to the mean log2 fold-change correlation."
    )
  }

  if ("SIGNAL_COLLAPSE_RISK" %in% warnings) {
    key_findings <- c(
      key_findings,
      sprintf(
        "Signal collapse was observed in %.1f%% of perturbation runs.",
        100 * dplyr::coalesce(dataset_stability$fraction_of_iterations_with_signal_collapse, 0)
      )
    )
  }

  influence_mode <- if (nrow(sample_influence) == 0) {
    "not_applicable"
  } else {
    as.character(dplyr::coalesce(sample_influence$influence_mode[1], ifelse(deg_metrics_applicable, "deg_based", "effect_or_rank_based")))
  }

  single_sample_text <- if (deg_metrics_applicable) {
    "No single sample appears to dominate the DEG signal under leave-one-out analysis."
  } else {
    "Sample influence is evaluated using effect-size and ranking sensitivity because DEG-based influence is not applicable in the current signal state."
  }
  mean_corr <- get_scalar_numeric(qc_context$mean_within_group_correlation)
  if (is.na(mean_corr)) mean_corr <- get_scalar_numeric(qc_context$mean_correlation)
  qc_poor <- !is.na(mean_corr) && mean_corr < 0.75

  if (nrow(sample_influence) > 0) {
    top_inf <- sample_influence[1, ]
    if (!is.na(top_inf$influence_score) && top_inf$influence_score >= 0.7) {
      single_sample_text <- if (deg_metrics_applicable) {
        sprintf(
          "Sample %s appears highly influential. Removing this sample substantially alters the DEG signal, suggesting potential over-reliance on individual samples.",
          top_inf$removed_sample
        )
      } else {
        sprintf(
          "Sample %s appears highly influential in ranking/effect patterns under perturbation. In the current low-signal context, this supports cautious interpretation rather than robust DEG-level conclusions.",
          top_inf$removed_sample
        )
      }
      if (qc_poor) {
        single_sample_text <- paste0(single_sample_text, " This may contribute to the observed group inconsistency and indicate technical variability.")
      }
    } else if (!is.na(top_inf$influence_score) && top_inf$influence_score >= 0.45) {
      single_sample_text <- if (deg_metrics_applicable) {
        sprintf(
          "Sample %s shows moderate influence on DEG ranking and recovery.",
          top_inf$removed_sample
        )
      } else {
        sprintf(
          "Sample %s may contribute to group inconsistency and instability in differential signal detection.",
          top_inf$removed_sample
        )
      }
    }
  }

  summary_text <- if (identical(stability_run_status, "failed")) {
    "Stability analysis could not be completed due to a technical failure in the perturbation assessment module."
  } else if (identical(signal_state, "no_detectable_signal") && effect_metrics_computable) {
    sprintf(
      paste0(
        "No robust differential expression (DEG) signal was detected under the configured thresholds; only weak effect-level signal was observed. ",
        "Stability evaluation is limited, as DEG-based stability metrics are not applicable. ",
        "Importantly, substantial within-group inconsistency is present (mean within-group correlation = %.2f), indicating that technical variability may dominate the observed signal. ",
        "Effect-size direction shows partial consistency across perturbations (mean log2 fold-change correlation = %.2f), but this does not provide strong evidence for biologically robust differential expression under current data conditions. ",
        "Composite effect-size consistency is derived from multiple perturbation-based measures (including log2 fold-change correlation and variability) and is not identical to the mean log2 fold-change correlation. ",
        "Top-ranked genes (based on statistical ranking) show %.1f%% overlap and may vary across perturbations. ",
        "Overall, the dataset does not provide stable and reproducible evidence of differential expression."
      ),
      dplyr::coalesce(mean_corr, NA_real_),
      dplyr::coalesce(dataset_stability$mean_log2fc_correlation, NA_real_),
      100 * dplyr::coalesce(dataset_stability$mean_top_n_overlap, 0)
    )
  } else if (identical(signal_state, "no_detectable_signal") && !effect_metrics_computable) {
    sprintf(
      paste0(
        "No robust differential expression (DEG) signal was detected under the configured thresholds; only weak effect-level signal was observed. ",
        "Stability evaluation is limited, as DEG-based stability metrics are not applicable. ",
        "Importantly, substantial within-group inconsistency is present (mean within-group correlation = %.2f), indicating that technical variability may dominate the observed signal. ",
        "Perturbation-based effect-size consistency could not be estimated reliably for this run, and no robust biological interpretation is supported under current data conditions."
      ),
      dplyr::coalesce(mean_corr, NA_real_)
    )
  } else if (identical(signal_state, "weak_signal")) {
    sprintf(
      paste0(
        "No robust differential expression (DEG) signal was detected under the configured thresholds; only weak effect-level signal was observed. ",
        "Stability evaluation is limited, as DEG-based stability metrics are not applicable. ",
        "Importantly, substantial within-group inconsistency is present (mean within-group correlation = %.2f), indicating that technical variability may dominate the observed signal. ",
        "Effect-size direction shows partial consistency across perturbations (mean log2 fold-change correlation = %.2f). However, this consistency should not be interpreted as evidence of robust biological signal under current data conditions. ",
        "These patterns may reflect weak or subtle expression differences; however, substantial within-group inconsistency and lack of robust DEGs significantly limit biological interpretation. ",
        "Top-ranked genes are defined based on statistical ranking and may vary across perturbations. ",
        "Overall, the dataset does not provide stable and reproducible evidence of differential expression."
      ),
      dplyr::coalesce(mean_corr, NA_real_),
      dplyr::coalesce(dataset_stability$mean_log2fc_correlation, NA_real_)
    )
  } else {
    dplyr::case_when(
      badge == "high" ~
        "Bootstrap/leave-one-out stability analysis suggests that the differential expression signal is robust to moderate sample perturbation.",
      badge == "moderate" ~
        "Stability analysis indicates that the differential expression signal is moderately robust, although top-ranked genes show noticeable sensitivity to sample perturbation.",
      TRUE ~
        "Stability analysis indicates that the differential expression signal is unstable across perturbations and should be interpreted with caution."
    )
  }

  top_rank_definition_text <- "Top-ranked genes are defined based on statistical ranking and may not meet differential expression thresholds."

  direction_text <- if (identical(stability_run_status, "failed")) {
    "Effect-size direction consistency was not produced because the perturbation analysis failed technically."
  } else if (!effect_metrics_computable) {
    "Effect-size direction consistency could not be estimated reliably for this run."
  } else if ("EFFECT_SIZE_UNSTABLE" %in% warnings) {
    "Effect-size direction is not consistently preserved across perturbations for a subset of reference DEGs."
  } else {
    "Effect-size direction remains broadly consistent across perturbation runs."
  }

  pca_qc_conflict_text <- if ("PCA_QC_CONFLICT" %in% warnings) {
    "Despite apparent PCA separation, low within-group correlation suggests that the observed clustering may not reflect biologically coherent groups."
  } else {
    ""
  }

  list(
    stability_headline = headline,
    stability_badge = badge,
    key_stability_findings = key_findings,
    signal_state = signal_state,
    stability_run_status = stability_run_status,
    deg_metrics_applicable = deg_metrics_applicable,
    stability_mode = stability_mode,
    influence_mode = influence_mode,
    deg_stability_score = dataset_stability$deg_stability_score,
    effect_stability_score = dataset_stability$effect_stability_score,
    final_stability_score = dataset_stability$final_stability_score,
    summary_text = summary_text,
    directionality_text = direction_text,
    top_rank_definition_text = top_rank_definition_text,
    pca_qc_conflict_text = pca_qc_conflict_text,
    sample_influence_text = single_sample_text,
    warnings = warnings
  )
}

#' Build stability narrative output (alias for report integration).
#' @param dataset_stability Dataset-level stability list.
#' @param warnings Warning code vector.
#' @param sample_influence Per-sample influence table.
#' @return List with headline, badge, findings and narrative paragraphs.
generate_stability_narrative <- function(dataset_stability,
                                         warnings,
                                         sample_influence,
                                         pca_separation = "unknown",
                                         qc_context = list()) {
  interpret_stability_results(
    dataset_stability = dataset_stability,
    warnings = warnings,
    sample_influence = sample_influence,
    pca_separation = pca_separation,
    qc_context = qc_context
  )
}
