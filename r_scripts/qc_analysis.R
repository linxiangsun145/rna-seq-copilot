# qc_analysis.R — Strict Rule-Based QC + Data Validation Module
#
# Public function:
#   run_qc(dds, vsd, counts_mat, meta, contrast_factor, res_df, plots_dir, results_dir)
#
# Outputs
#   plots_dir/sample_distance_heatmap.png
#   plots_dir/sample_correlation_heatmap.png
#   plots_dir/library_size.png
#   plots_dir/count_distribution.png
#   plots_dir/zero_fraction.png
#   results_dir/qc_report.json

# ── Plot helpers ──────────────────────────────────────────────────────────────
.qc_plot_sample_distance <- function(vsd, meta, outfile) {
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)

  ann_col <- NULL
  if (ncol(meta) > 0) {
    ann_col <- data.frame(
      row.names = rownames(meta),
      Group = as.character(meta[, 1]),
      check.names = FALSE
    )
  }

  colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

  png(outfile, width = 900, height = 760, res = 120)
  on.exit(dev.off(), add = TRUE)
  pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    col = colours,
    annotation_col = ann_col,
    main = "Sample-to-sample distance (Euclidean on VST)",
    fontsize = 10
  )
}

.qc_plot_sample_correlation <- function(vsd, meta, outfile) {
  corr_mat <- cor(assay(vsd), method = "pearson")

  ann_col <- NULL
  if (ncol(meta) > 0) {
    ann_col <- data.frame(
      row.names = rownames(meta),
      Group = as.character(meta[, 1]),
      check.names = FALSE
    )
  }

  colours <- colorRampPalette(c("#313695", "#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#d73027"))(255)

  png(outfile, width = 900, height = 760, res = 120)
  on.exit(dev.off(), add = TRUE)
  pheatmap(
    corr_mat,
    col = colours,
    breaks = seq(-1, 1, length.out = 256),
    annotation_col = ann_col,
    main = "Sample-to-sample Pearson correlation (VST)",
    fontsize = 10
  )
}

.qc_plot_library_size <- function(counts_mat, meta, group_col, outfile) {
  lib_sizes <- colSums(counts_mat)
  df <- data.frame(
    sample = names(lib_sizes),
    lib_size = as.numeric(lib_sizes),
    stringsAsFactors = FALSE
  )

  if (group_col %in% colnames(meta)) {
    df$group <- as.character(meta[df$sample, group_col])
  } else {
    df$group <- "all"
  }

  med <- median(lib_sizes)

  p <- ggplot(df, aes(x = reorder(sample, lib_size), y = lib_size, fill = group)) +
    geom_bar(stat = "identity", alpha = 0.9) +
    geom_hline(yintercept = med, linetype = "dashed", colour = "black") +
    geom_hline(yintercept = 0.5 * med, linetype = "dotted", colour = "#f59e0b") +
    geom_hline(yintercept = 0.3 * med, linetype = "dotted", colour = "#dc2626") +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Library size per sample",
      subtitle = "Dashed: median; dotted orange: 50% median; dotted red: 30% median",
      x = "Sample",
      y = "Total counts",
      fill = group_col
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

  ggsave(outfile, plot = p, width = max(7, 0.55 * nrow(df) + 2), height = 5.5, dpi = 150)
}

.qc_plot_count_distribution <- function(counts_mat, outfile) {
  log_counts <- log2(counts_mat + 1)
  df <- data.frame(
    log2count = as.vector(log_counts),
    sample = rep(colnames(counts_mat), each = nrow(counts_mat)),
    stringsAsFactors = FALSE
  )

  p <- ggplot(df, aes(x = log2count, colour = sample)) +
    geom_density(linewidth = 0.7) +
    labs(
      title = "Count distribution by sample",
      x = expression(log[2](count + 1)),
      y = "Density",
      colour = "Sample"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.text = element_text(size = 7))

  ggsave(outfile, plot = p, width = 8.5, height = 5, dpi = 150)
}

.qc_plot_zero_fraction <- function(counts_mat, outfile) {
  zero_frac <- colMeans(counts_mat == 0)
  df <- data.frame(sample = names(zero_frac), zero_fraction = zero_frac, stringsAsFactors = FALSE)

  p <- ggplot(df, aes(x = reorder(sample, zero_fraction), y = zero_fraction)) +
    geom_bar(stat = "identity", fill = "#2563eb", alpha = 0.9) +
    geom_hline(yintercept = 0.2, linetype = "dotted", colour = "#f59e0b") +
    geom_hline(yintercept = 0.4, linetype = "dotted", colour = "#dc2626") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Zero-count gene fraction per sample",
      subtitle = "Dotted orange: 20%; dotted red: 40%",
      x = "Sample",
      y = "Zero fraction"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

  ggsave(outfile, plot = p, width = max(7, 0.55 * nrow(df) + 2), height = 5, dpi = 150)
}

# ── Stats helpers ─────────────────────────────────────────────────────────────
.qc_add_issue <- function(message, level, qc_warnings, qc_critical) {
  if (level == "critical") {
    qc_critical <- c(qc_critical, message)
  } else {
    qc_warnings <- c(qc_warnings, message)
  }
  list(qc_warnings = qc_warnings, qc_critical = qc_critical)
}

.qc_safe_numeric <- function(x) {
  ifelse(is.finite(x), x, NA_real_)
}

.qc_compute_anova_r2 <- function(score, grp) {
  if (length(unique(grp)) < 2) return(NA_real_)
  fit <- lm(score ~ grp)
  ss_total <- sum((score - mean(score))^2)
  if (ss_total == 0) return(0)
  ss_model <- sum((fitted(fit) - mean(score))^2)
  ss_model / ss_total
}

.qc_is_confounded <- function(batch, condition) {
  tab <- table(batch, condition)
  row_pure <- all(rowSums(tab > 0) == 1)
  col_pure <- all(colSums(tab > 0) == 1)
  row_pure && col_pure
}

.qc_align_group_labels <- function(meta, sample_names, contrast_factor) {
  if (!(contrast_factor %in% colnames(meta))) {
    return(rep("all", length(sample_names)))
  }

  grp <- NULL
  if (!is.null(rownames(meta)) && all(sample_names %in% rownames(meta))) {
    grp <- as.character(meta[sample_names, contrast_factor, drop = TRUE])
  } else if ("sample" %in% colnames(meta) && all(sample_names %in% as.character(meta$sample))) {
    idx <- match(sample_names, as.character(meta$sample))
    grp <- as.character(meta[idx, contrast_factor, drop = TRUE])
  } else if (nrow(meta) == length(sample_names)) {
    cat("[QC][WARN] metadata rownames/sample column missing or mismatched; falling back to row-order alignment\n")
    grp <- as.character(meta[[contrast_factor]])
  } else {
    cat("[QC][WARN] cannot align metadata to sample names; assigning all samples to one group\n")
    grp <- rep("all", length(sample_names))
  }

  grp[is.na(grp) | grp == ""] <- "unknown"
  grp
}

.qc_prepare_vst_for_sample_correlation <- function(vst_mat, sample_names) {
  # We need a genes x samples matrix before calling cor(), so that cor() returns sample x sample.
  if (!is.null(colnames(vst_mat)) && all(sample_names %in% colnames(vst_mat))) {
    aligned <- vst_mat[, sample_names, drop = FALSE]
    return(list(expr = aligned, orientation = "genes_by_samples"))
  }

  if (!is.null(rownames(vst_mat)) && all(sample_names %in% rownames(vst_mat))) {
    aligned <- t(vst_mat[sample_names, , drop = FALSE])
    return(list(expr = aligned, orientation = "samples_by_genes_transposed"))
  }

  if (ncol(vst_mat) == length(sample_names)) {
    colnames(vst_mat) <- sample_names
    return(list(expr = vst_mat, orientation = "genes_by_samples_assumed"))
  }

  if (nrow(vst_mat) == length(sample_names)) {
    aligned <- t(vst_mat)
    colnames(aligned) <- sample_names
    return(list(expr = aligned, orientation = "samples_by_genes_transposed_assumed"))
  }

  stop("Unable to align VST matrix to sample names for correlation computation")
}

# ── Main runner ───────────────────────────────────────────────────────────────
run_qc <- function(dds, vsd, counts_mat, meta, contrast_factor, res_df,
                   plots_dir, results_dir) {
  cat("[QC] Running strict QC + validation checks\n")

  warning_items <- list()
  qc_warnings <- character(0)
  qc_critical <- character(0)

  .add_warning <- function(severity, code, message, sample = NULL, metric = NULL) {
    item <- list(
      type = "qc",
      severity = severity,
      code = code,
      message = message,
      sample = sample,
      metric = metric
    )
    warning_items[[length(warning_items) + 1]] <<- item
    if (severity == "critical") {
      qc_critical <<- c(qc_critical, message)
    } else {
      qc_warnings <<- c(qc_warnings, message)
    }
  }

  vst_mat <- assay(vsd)
  sample_names <- colnames(counts_mat)

  if (is.null(sample_names) || length(sample_names) == 0) {
    stop("counts_mat must have sample columns")
  }

  prep <- .qc_prepare_vst_for_sample_correlation(vst_mat, sample_names)
  vst_for_corr <- prep$expr
  cat(sprintf("[QC] Correlation matrix orientation: %s (genes=%d, samples=%d)\n", prep$orientation, nrow(vst_for_corr), ncol(vst_for_corr)))

  # Required plots
  .qc_plot_sample_distance(vsd, meta, file.path(plots_dir, "sample_distance_heatmap.png"))
  .qc_plot_sample_correlation(vsd, meta, file.path(plots_dir, "sample_correlation_heatmap.png"))
  .qc_plot_library_size(counts_mat, meta, contrast_factor, file.path(plots_dir, "library_size.png"))
  .qc_plot_count_distribution(counts_mat, file.path(plots_dir, "count_distribution.png"))
  .qc_plot_zero_fraction(counts_mat, file.path(plots_dir, "zero_fraction.png"))

  pca_fit <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
  pca_var <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  pca_variance <- list(
    PC1 = round(.qc_safe_numeric(pca_var[1]), 4),
    PC2 = round(.qc_safe_numeric(pca_var[2]), 4)
  )

  corr_mat <- cor(vst_for_corr, method = "pearson", use = "pairwise.complete.obs")
  lib_sizes <- colSums(counts_mat)
  zero_fraction <- colMeans(counts_mat == 0)

  grp <- .qc_align_group_labels(meta, sample_names, contrast_factor)
  grp_tbl <- table(grp)
  replicates_per_group <- as.list(setNames(as.integer(grp_tbl), names(grp_tbl)))

  ratio <- if (length(grp_tbl) > 0) as.numeric(max(grp_tbl) / min(grp_tbl)) else 1
  group_status <- "balanced"
  if (ratio > 2) {
    group_status <- "critical"
    .add_warning("critical", "group_imbalance_critical", sprintf("Group imbalance detected (ratio %.2f > 2.00).", ratio), metric = "group_balance")
  } else if (ratio >= 1.5) {
    group_status <- "warning"
    .add_warning("warning", "group_imbalance_warning", sprintf("Group imbalance detected (ratio %.2f >= 1.50).", ratio), metric = "group_balance")
  }
  group_balance <- list(status = group_status, ratio = round(ratio, 4))

  global_median_lib <- median(lib_sizes)
  library_size_flags <- list()
  zero_fraction_flags <- list()
  correlation_flags <- list()
  batch_flags <- list()
  outlier_samples <- character(0)

  sample_warning_count <- setNames(rep(0L, length(sample_names)), sample_names)
  sample_critical_lib <- setNames(rep(FALSE, length(sample_names)), sample_names)
  sample_critical_corr <- setNames(rep(FALSE, length(sample_names)), sample_names)
  within_group_corr <- setNames(rep(NA_real_, length(sample_names)), sample_names)
  peer_corr_debug <- setNames(vector("list", length(sample_names)), sample_names)
  sample_qc_reasons <- setNames(vector("list", length(sample_names)), sample_names)

  # Library size + zero-fraction flags
  for (s in sample_names) {
    val <- as.numeric(lib_sizes[[s]])
    if (val < 0.3 * global_median_lib) {
      library_size_flags[[length(library_size_flags) + 1]] <- list(
        sample = s, level = "critical", value = val,
        threshold = round(0.3 * global_median_lib, 2),
        rule = "library_size < 30% median"
      )
      sample_critical_lib[[s]] <- TRUE
      .add_warning("critical", "low_library_size_critical", sprintf("Low library size detected for %s (%.0f vs median %.0f).", s, val, global_median_lib), sample = s, metric = "library_size")
    } else if (val < 0.5 * global_median_lib) {
      library_size_flags[[length(library_size_flags) + 1]] <- list(
        sample = s, level = "warning", value = val,
        threshold = round(0.5 * global_median_lib, 2),
        rule = "library_size < 50% median"
      )
      sample_warning_count[[s]] <- sample_warning_count[[s]] + 1L
      .add_warning("warning", "low_library_size_warning", sprintf("Low library size detected for %s (%.0f vs median %.0f).", s, val, global_median_lib), sample = s, metric = "library_size")
    }

    z <- as.numeric(zero_fraction[[s]])
    if (z > 0.4) {
      zero_fraction_flags[[length(zero_fraction_flags) + 1]] <- list(sample = s, level = "critical", value = round(z, 4), threshold = 0.4, rule = "zero_fraction > 40%")
      sample_warning_count[[s]] <- sample_warning_count[[s]] + 1L
      .add_warning("critical", "high_zero_fraction_critical", sprintf("High zero-count fraction detected for %s (%.1f%% > 40%%).", s, 100 * z), sample = s, metric = "zero_fraction")
    } else if (z > 0.2) {
      zero_fraction_flags[[length(zero_fraction_flags) + 1]] <- list(sample = s, level = "warning", value = round(z, 4), threshold = 0.2, rule = "zero_fraction > 20%")
      sample_warning_count[[s]] <- sample_warning_count[[s]] + 1L
      .add_warning("warning", "high_zero_fraction_warning", sprintf("Elevated zero-count fraction detected for %s (%.1f%% > 20%%).", s, 100 * z), sample = s, metric = "zero_fraction")
    }
  }

  # Within-group sample correlation (exclude self, safe for small groups)
  for (g in names(grp_tbl)) {
    group_samples <- sample_names[grp == g]
    if (length(group_samples) < 2) {
      for (s in group_samples) {
        peer_corr_debug[[s]] <- list()
        sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], sprintf("no peer in group '%s'; correlation flag skipped", g))
      }
      next
    }

    for (s in group_samples) {
      peers <- setdiff(group_samples, s)
      peer_vals <- as.numeric(corr_mat[s, peers, drop = TRUE])
      names(peer_vals) <- peers
      peer_vals <- peer_vals[is.finite(peer_vals)]
      peer_corr_debug[[s]] <- as.list(round(peer_vals, 4))

      m <- if (length(peer_vals) > 0) mean(peer_vals, na.rm = TRUE) else NA_real_
      within_group_corr[[s]] <- as.numeric(m)

      if (is.finite(m) && m < 0.75) {
        sample_critical_corr[[s]] <- TRUE
        correlation_flags[[length(correlation_flags) + 1]] <- list(sample = s, level = "critical", mean_correlation = round(m, 4), threshold = 0.75, rule = "mean_within_group_correlation < 0.75")
        .add_warning("critical", "low_within_group_correlation_critical", sprintf("Low within-group correlation detected for %s (%.3f < 0.75).", s, m), sample = s, metric = "mean_within_group_correlation")
        sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], sprintf("critical correlation: %.3f < 0.75", m))
      } else if (is.finite(m) && m < 0.85) {
        sample_warning_count[[s]] <- sample_warning_count[[s]] + 1L
        correlation_flags[[length(correlation_flags) + 1]] <- list(sample = s, level = "warning", mean_correlation = round(m, 4), threshold = 0.85, rule = "mean_within_group_correlation < 0.85")
        .add_warning("warning", "low_within_group_correlation_warning", sprintf("Low within-group correlation detected for %s (%.3f < 0.85).", s, m), sample = s, metric = "mean_within_group_correlation")
        sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], sprintf("warning correlation: %.3f < 0.85", m))
      } else if (is.finite(m)) {
        sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], sprintf("correlation OK: %.3f", m))
      } else {
        sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], "correlation unavailable after finite-value filtering")
      }
    }
  }

  corr_flagged_samples <- unique(c(
    vapply(correlation_flags, function(x) if (!is.null(x$sample)) as.character(x$sample) else NA_character_, character(1), USE.NAMES = FALSE)
  ))
  corr_flagged_samples <- corr_flagged_samples[!is.na(corr_flagged_samples) & corr_flagged_samples != ""]
  corr_flag_ratio <- if (length(sample_names) > 0) length(unique(corr_flagged_samples)) / length(sample_names) else 0
  if (corr_flag_ratio >= 0.8 && length(sample_names) >= 4) {
    .add_warning(
      "warning",
      "correlation_sanity_check_warning",
      sprintf("Sanity check: %.0f%% of samples were correlation-flagged; verify group labels and correlation orientation.", 100 * corr_flag_ratio),
      metric = "mean_within_group_correlation"
    )
    cat(sprintf("[QC][WARN] Sanity check triggered: %d/%d samples flagged by correlation\n", length(unique(corr_flagged_samples)), length(sample_names)))
  }

  # Optional outlier warning (centroid distance) for groups with n>=3
  dist_mat <- as.matrix(dist(t(vst_mat), method = "euclidean"))
  for (g in names(grp_tbl)) {
    group_samples <- sample_names[grp == g]
    if (length(group_samples) < 3) next
    group_dist <- dist_mat[group_samples, group_samples, drop = FALSE]
    sample_to_centroid <- sapply(group_samples, function(s) mean(group_dist[s, setdiff(group_samples, s), drop = TRUE], na.rm = TRUE))
    med <- median(sample_to_centroid, na.rm = TRUE)
    madv <- mad(sample_to_centroid, center = med, constant = 1, na.rm = TRUE)
    thr <- med + 3 * madv
    bad <- names(sample_to_centroid)[sample_to_centroid > thr]
    if (length(bad) > 0) {
      outlier_samples <- unique(c(outlier_samples, bad))
      for (s in bad) {
        sample_warning_count[[s]] <- sample_warning_count[[s]] + 1L
        .add_warning("warning", "distance_outlier_warning", sprintf("Distance outlier detected for %s in group %s (%.3f > %.3f).", s, g, sample_to_centroid[[s]], thr), sample = s, metric = "sample_distance")
      }
    }
  }

  # Batch checks stay in QC type
  if ("batch" %in% colnames(meta)) {
    batch <- as.character(meta$batch)
    condition <- grp

    if (.qc_is_confounded(batch, condition)) {
      batch_flags[[length(batch_flags) + 1]] <- list(level = "critical", rule = "batch_condition_confounded", detail = "Batch and condition are one-to-one confounded")
      .add_warning("critical", "batch_condition_confounded", "Batch is confounded with condition.", metric = "batch")
    } else {
      pc1 <- pca_fit$x[, 1]
      r2_batch <- .qc_compute_anova_r2(pc1, batch)
      r2_cond <- .qc_compute_anova_r2(pc1, condition)
      if (is.finite(r2_batch) && is.finite(r2_cond) && (r2_batch - r2_cond) >= 0.1 && r2_batch >= 0.3) {
        batch_flags[[length(batch_flags) + 1]] <- list(level = "warning", rule = "batch_dominates_pc1", detail = sprintf("PC1 R2 batch=%.3f, condition=%.3f", r2_batch, r2_cond))
        .add_warning("warning", "batch_dominates_pc1", sprintf("Batch effect may dominate PC1 (R2 %.3f vs %.3f).", r2_batch, r2_cond), metric = "batch")
      }
    }
  }

  # Low-quality sample definition:
  # critical library size OR critical correlation OR multiple warnings
  low_quality_samples <- character(0)
  per_sample_qc_metrics <- list()
  for (s in sample_names) {
    qc_flags <- character(0)
    if (isTRUE(sample_critical_lib[[s]])) {
      qc_flags <- c(qc_flags, "low_library_size_critical")
      sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], "critical library size")
    }
    if (isTRUE(sample_critical_corr[[s]])) {
      qc_flags <- c(qc_flags, "low_within_group_correlation_critical")
    }
    if (sample_warning_count[[s]] >= 2) {
      qc_flags <- c(qc_flags, "multiple_qc_warnings")
      sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], sprintf("multiple warnings (%d)", sample_warning_count[[s]]))
    }

    condition <- grp[which(sample_names == s)][1]
    if (is.na(condition) || is.null(condition) || condition == "") condition <- "unknown"

    if (length(qc_flags) > 0) {
      low_quality_samples <- c(low_quality_samples, s)
      sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], "included in low_quality_samples")
    } else {
      sample_qc_reasons[[s]] <- c(sample_qc_reasons[[s]], "not low-quality by composite rule")
    }

    decision <- if (length(qc_flags) > 0) "low_quality" else "retain"
    peer_debug_text <- "none"
    if (length(peer_corr_debug[[s]]) > 0) {
      peer_pairs <- unlist(peer_corr_debug[[s]], use.names = TRUE)
      peer_debug_text <- paste(sprintf("%s=%.4f", names(peer_pairs), as.numeric(peer_pairs)), collapse = ",")
    }

    cat(sprintf("[QC][DEBUG] sample=%s condition=%s lib=%.0f mean_corr=%s peers=[%s] decision=%s reasons=%s\\n",
      s,
      condition,
      as.numeric(lib_sizes[[s]]),
      ifelse(is.finite(within_group_corr[[s]]), sprintf("%.4f", within_group_corr[[s]]), "NA"),
      peer_debug_text,
      decision,
      paste(unique(sample_qc_reasons[[s]]), collapse = "; ")
    ))

    per_sample_qc_metrics[[length(per_sample_qc_metrics) + 1]] <- list(
      sample = s,
      condition = condition,
      library_size = as.numeric(lib_sizes[[s]]),
      zero_fraction = round(as.numeric(zero_fraction[[s]]), 4),
      mean_within_group_correlation = ifelse(is.finite(within_group_corr[[s]]), round(within_group_corr[[s]], 4), NA_real_),
      peer_correlations = peer_corr_debug[[s]],
      qc_decision = decision,
      qc_reasons = as.list(unique(sample_qc_reasons[[s]])),
      qc_flags = as.list(unique(qc_flags))
    )
  }

  qc_metrics <- list(
    library_size = list(
      min = as.numeric(min(lib_sizes)),
      median = as.numeric(median(lib_sizes)),
      min_median_ratio = as.numeric(min(lib_sizes) / median(lib_sizes))
    ),
    correlation = list(
      mean_within_group = as.numeric(mean(unlist(within_group_corr), na.rm = TRUE))
    ),
    zero_fraction = list(
      mean = as.numeric(mean(zero_fraction, na.rm = TRUE))
    )
  )

  # Deduplicate warning items by (type,severity,code,sample,metric)
  if (length(warning_items) > 0) {
    keys <- vapply(warning_items, function(x) paste(x$type, x$severity, x$code, x$sample, x$metric, sep = "|"), character(1))
    warning_items <- warning_items[!duplicated(keys)]
  }

  qc_report <- list(
    outliers = as.list(unique(outlier_samples)),
    low_quality_samples = as.list(unique(low_quality_samples)),
    group_balance = group_balance,
    replicates_per_group = replicates_per_group,
    library_size_flags = library_size_flags,
    zero_fraction_flags = zero_fraction_flags,
    correlation_flags = correlation_flags,
    batch_flags = batch_flags,
    realism_flags = list(),
    qc_warnings = as.list(unique(qc_warnings)),
    qc_critical = as.list(unique(qc_critical)),
    pca_variance = pca_variance,
    qc_metrics = qc_metrics,
    per_sample_qc_metrics = per_sample_qc_metrics,
    warning_items = warning_items
  )

  qc_json_path <- file.path(results_dir, "qc_report.json")
  write(jsonlite::toJSON(qc_report, auto_unbox = TRUE, pretty = TRUE), qc_json_path)
  cat(sprintf("[QC] qc_report.json written: %s\n", qc_json_path))

  invisible(qc_report)
}
