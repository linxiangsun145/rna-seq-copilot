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

  signal_state <- as.character(dataset_stability$signal_state)
  deg_metrics_applicable <- isTRUE(dataset_stability$deg_metrics_applicable)
  final_score <- dataset_stability$final_stability_score

  if (identical(signal_state, "no_detectable_signal")) {
    warnings <- c(warnings, "LACK_OF_DETECTABLE_SIGNAL")
  }

  score <- dplyr::coalesce(final_score, 0)
  if (score < 0.4 && identical(signal_state, "strong_signal")) warnings <- c(warnings, "STABILITY_LOW")

  top_overlap <- dplyr::coalesce(dataset_stability$mean_top_n_overlap, 0)
  if (top_overlap < 0.4) warnings <- c(warnings, "TOP_GENES_UNSTABLE")

  ref_deg <- gene_stability |>
    dplyr::filter(.data$ref_class == "strong_DEG")

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

  mean_corr <- suppressWarnings(as.numeric(qc_context$mean_within_group_correlation))
  if (is.na(mean_corr)) mean_corr <- suppressWarnings(as.numeric(qc_context$mean_correlation))
  group_correlation_low <- !is.na(mean_corr) && mean_corr < 0.75
  if (tolower(as.character(pca_separation)) == "clear" && group_correlation_low) {
    warnings <- c(warnings, "PCA_QC_CONFLICT")
  }

  qc_status <- tolower(as.character(qc_context$qc_status))
  qc_high_risk <- qc_status %in% c("critical", "high", "high_risk")
  if (qc_high_risk && !identical(signal_state, "strong_signal") && group_correlation_low) {
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
  signal_state <- as.character(dataset_stability$signal_state)
  stability_mode <- as.character(dataset_stability$stability_mode)
  deg_metrics_applicable <- isTRUE(dataset_stability$deg_metrics_applicable)
  score <- dplyr::coalesce(dataset_stability$final_stability_score, 0)
  badge <- if (identical(signal_state, "no_detectable_signal")) {
    "not_applicable"
  } else if (identical(signal_state, "weak_signal")) {
    "low_signal"
  } else {
    dplyr::case_when(
      score >= 0.75 ~ "high",
      score >= 0.50 ~ "moderate",
      TRUE ~ "low"
    )
  }

  headline <- if (identical(signal_state, "no_detectable_signal")) {
    "No robust differential expression signal was detectable under the configured thresholds."
  } else if (identical(signal_state, "weak_signal")) {
    "Only weak differential signal was detected; DEG-based stability metrics are limited under current thresholds."
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

  if ("SIGNAL_COLLAPSE_RISK" %in% warnings) {
    key_findings <- c(
      key_findings,
      sprintf(
        "Signal collapse was observed in %.1f%% of perturbation runs.",
        100 * dplyr::coalesce(dataset_stability$fraction_of_iterations_with_signal_collapse, 0)
      )
    )
  }

  single_sample_text <- if (deg_metrics_applicable) {
    "No single sample appears to dominate the DEG signal under leave-one-out analysis."
  } else {
    "Single-sample influence on robust DEG recovery is not directly assessable because no strong reference DEG set was detected."
  }
  mean_corr <- suppressWarnings(as.numeric(qc_context$mean_within_group_correlation))
  if (is.na(mean_corr)) mean_corr <- suppressWarnings(as.numeric(qc_context$mean_correlation))
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
      single_sample_text <- sprintf(
        "Sample %s shows moderate influence on DEG ranking and recovery.",
        top_inf$removed_sample
      )
    }
  }

  summary_text <- if (identical(signal_state, "no_detectable_signal")) {
    sprintf(
      paste0(
        "No robust differential expression signal was detectable under the configured thresholds. ",
        "Stability evaluation is therefore limited, as DEG-based stability metrics are not applicable. ",
        "Effect-size direction shows moderate consistency (mean log2FC correlation = %.2f); however, given the low within-group correlation and absence of significant DEGs, this consistency should not be interpreted as evidence of a robust biological signal. ",
        "Top-ranked genes (based on statistical ranking) show limited overlap (%.1f%%), indicating sensitivity to sample perturbation. ",
        "%s ",
        "%s ",
        "Overall, the dataset does not provide stable and reproducible evidence of differential expression under the current thresholds."
      ),
      dplyr::coalesce(dataset_stability$mean_log2fc_correlation, NA_real_),
      100 * dplyr::coalesce(dataset_stability$mean_top_n_overlap, 0),
      ifelse("PCA_QC_CONFLICT" %in% warnings, "Despite apparent PCA separation, low within-group correlation suggests that the observed clustering may not reflect biologically coherent groups.", ""),
      ifelse("QC_DOMINATES_SIGNAL" %in% warnings, "Technical variability appears to dominate the observed expression pattern, limiting confidence in biological interpretation.", "")
    )
  } else if (identical(signal_state, "weak_signal")) {
    "The analysis identified only weak differential signal under the configured thresholds. DEG-based stability assessment is limited, but effect-size direction shows partial consistency across perturbations. These patterns may reflect subtle biological differences; however, the signal remains sensitive to sample perturbation and should be interpreted cautiously."
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

  direction_text <- if ("EFFECT_SIZE_UNSTABLE" %in% warnings) {
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
    deg_metrics_applicable = deg_metrics_applicable,
    stability_mode = stability_mode,
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
