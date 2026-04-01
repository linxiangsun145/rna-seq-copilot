#!/usr/bin/env Rscript

# Minimal regression tests for stability state-machine semantics.
# Run from repo root:
#   Rscript r_scripts/test_stability_state_machine.R

suppressPackageStartupMessages({
  library(dplyr)
})

source("R/stability_scoring.R")
source("R/stability_reporting.R")

assert_equal <- function(actual, expected, msg) {
  if (!identical(actual, expected)) {
    stop(sprintf("[FAIL] %s | expected=%s actual=%s", msg, expected, actual), call. = FALSE)
  }
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) {
    stop(sprintf("[FAIL] %s", msg), call. = FALSE)
  }
}

cat("Running stability state-machine regression checks...\n")

# Case 1: strong signal -> completed, not failed.
signal_1 <- detect_signal_state(n_ref_deg = 12L, weak_effect_gene_count = 0L)
mode_1 <- determine_stability_mode(signal_1, effect_metrics_computable = TRUE)
run_1 <- determine_stability_run_status(signal_1, mode_1, engine_failed = FALSE)
badge_1 <- assign_stability_badge(signal_1, final_stability_score = 0.81, stability_run_status = run_1)
assert_equal(signal_1, "strong_signal", "n_ref_deg > 0 should be strong_signal")
assert_equal(run_1, "completed", "strong_signal should map to completed")
assert_equal(badge_1, "high", "strong signal with high score should be high")

# Case 2: no detectable signal -> limited + low_signal + no failure warning.
signal_2 <- detect_signal_state(n_ref_deg = 0L, weak_effect_gene_count = 0L)
mode_2 <- determine_stability_mode(signal_2, effect_metrics_computable = TRUE)
run_2 <- determine_stability_run_status(signal_2, mode_2, engine_failed = FALSE)
badge_2 <- assign_stability_badge(signal_2, final_stability_score = NA_real_, stability_run_status = run_2)
assert_equal(signal_2, "no_detectable_signal", "n_ref_deg == 0 and no weak effects should be no_detectable_signal")
assert_equal(run_2, "limited", "no_detectable_signal should be limited")
assert_equal(badge_2, "low_signal", "no_detectable_signal should map to low_signal")

warnings_2 <- generate_stability_warnings(
  dataset_stability = list(
    signal_state = signal_2,
    stability_run_status = run_2,
    deg_metrics_applicable = FALSE,
    final_stability_score = NA_real_,
    mean_top_n_overlap = 0.35,
    fraction_of_iterations_with_signal_collapse = NA_real_
  ),
  gene_stability = data.frame(ref_class = character(), sign_consistency = numeric()),
  sample_influence = data.frame(),
  pca_separation = "unknown",
  qc_context = list(
    mean_within_group_correlation = NA_real_,
    mean_correlation = NA_real_,
    qc_status = "unknown"
  )
)
assert_true("LACK_OF_DETECTABLE_SIGNAL" %in% warnings_2, "no-detectable-signal should carry LACK_OF_DETECTABLE_SIGNAL")
assert_true(!("STABILITY_ANALYSIS_FAILED" %in% warnings_2), "no-detectable-signal limited run must not carry technical-failure warning")

# Case 3: weak signal -> limited + low_signal.
signal_3 <- detect_signal_state(n_ref_deg = 0L, weak_effect_gene_count = 7L)
mode_3 <- determine_stability_mode(signal_3, effect_metrics_computable = TRUE)
run_3 <- determine_stability_run_status(signal_3, mode_3, engine_failed = FALSE)
badge_3 <- assign_stability_badge(signal_3, final_stability_score = NA_real_, stability_run_status = run_3)
assert_equal(signal_3, "weak_signal", "n_ref_deg == 0 with weak effects should be weak_signal")
assert_equal(run_3, "limited", "weak_signal should map to limited")
assert_equal(badge_3, "low_signal", "weak_signal should map to low_signal")

# Case 4: true technical failure -> failed + unknown + failure warning.
run_4 <- determine_stability_run_status("no_detectable_signal", "effect_only", engine_failed = TRUE)
badge_4 <- assign_stability_badge("no_detectable_signal", final_stability_score = NA_real_, stability_run_status = run_4)
assert_equal(run_4, "failed", "engine failure must map to failed")
assert_equal(badge_4, "unknown", "failed run must map to unknown badge")

warnings_4 <- generate_stability_warnings(
  dataset_stability = list(
    signal_state = "no_detectable_signal",
    stability_run_status = run_4,
    deg_metrics_applicable = FALSE,
    final_stability_score = NA_real_,
    mean_top_n_overlap = NA_real_,
    fraction_of_iterations_with_signal_collapse = NA_real_
  ),
  gene_stability = data.frame(ref_class = character(), sign_consistency = numeric()),
  sample_influence = data.frame(),
  pca_separation = "unknown",
  qc_context = list(
    mean_within_group_correlation = NA_real_,
    mean_correlation = NA_real_,
    qc_status = "unknown"
  )
)
assert_equal(length(warnings_4), 1L, "failed status should return only STABILITY_ANALYSIS_FAILED")
assert_equal(warnings_4[[1]], "STABILITY_ANALYSIS_FAILED", "failed status should emit technical failure warning")

cat("All stability state-machine checks passed.\n")
