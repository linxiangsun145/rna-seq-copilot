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

# ── Main runner ───────────────────────────────────────────────────────────────
run_qc <- function(dds, vsd, counts_mat, meta, contrast_factor, res_df,
                   plots_dir, results_dir) {
  cat("[QC] Running strict QC + validation checks\n")

  qc_warnings <- character(0)
  qc_critical <- character(0)

  # 0) Base computed objects
  vst_mat <- assay(vsd)
  sample_names <- colnames(counts_mat)

  # 1) Required plots
  .qc_plot_sample_distance(vsd, meta, file.path(plots_dir, "sample_distance_heatmap.png"))
  .qc_plot_sample_correlation(vsd, meta, file.path(plots_dir, "sample_correlation_heatmap.png"))
  .qc_plot_library_size(counts_mat, meta, contrast_factor, file.path(plots_dir, "library_size.png"))
  .qc_plot_count_distribution(counts_mat, file.path(plots_dir, "count_distribution.png"))
  .qc_plot_zero_fraction(counts_mat, file.path(plots_dir, "zero_fraction.png"))

  # 2) PCA + variance explained
  pca_fit <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
  pca_var <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  pca_variance <- list(
    PC1 = round(.qc_safe_numeric(pca_var[1]), 4),
    PC2 = round(.qc_safe_numeric(pca_var[2]), 4)
  )

  # 3) Distances, correlations, library size, zero fraction
  dist_mat <- as.matrix(dist(t(vst_mat), method = "euclidean"))
  corr_mat <- cor(vst_mat, method = "pearson")
  lib_sizes <- colSums(counts_mat)
  zero_fraction <- colMeans(counts_mat == 0)

  # 4) Replicates and group balance
  grp <- if (contrast_factor %in% colnames(meta)) as.character(meta[[contrast_factor]]) else rep("all", ncol(counts_mat))
  grp_tbl <- table(grp)
  replicates_per_group <- as.list(setNames(as.integer(grp_tbl), names(grp_tbl)))

  ratio <- if (length(grp_tbl) > 0) as.numeric(max(grp_tbl) / min(grp_tbl)) else 1
  if (ratio < 1.5) {
    group_status <- "balanced"
  } else if (ratio <= 2) {
    group_status <- "warning"
    msg <- sprintf("Group imbalance warning: ratio=%.3f (1.5 <= ratio <= 2).", ratio)
    out <- .qc_add_issue(msg, "warning", qc_warnings, qc_critical)
    qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
  } else {
    group_status <- "critical"
    msg <- sprintf("Group imbalance critical: ratio=%.3f (>2).", ratio)
    out <- .qc_add_issue(msg, "critical", qc_warnings, qc_critical)
    qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
  }
  group_balance <- list(status = group_status, ratio = round(ratio, 4))

  # 5) Replicate count rules
  for (g in names(grp_tbl)) {
    n <- as.integer(grp_tbl[[g]])
    if (n < 2) {
      msg <- sprintf("Group '%s' has n=%d (<2): critical replicate deficiency.", g, n)
      out <- .qc_add_issue(msg, "critical", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    } else if (n == 2) {
      msg <- sprintf("Group '%s' has n=2: warning for low replicate count.", g)
      out <- .qc_add_issue(msg, "warning", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    }
  }

  # 6) Library size rules
  global_median_lib <- median(lib_sizes)
  library_size_flags <- list()
  for (s in sample_names) {
    val <- as.numeric(lib_sizes[[s]])
    if (val < 0.3 * global_median_lib) {
      library_size_flags[[length(library_size_flags) + 1]] <- list(sample = s, level = "critical", value = val, threshold = round(0.3 * global_median_lib, 2), rule = "library_size < 30% median")
      msg <- sprintf("Sample '%s' critical library size: %.0f < 30%% median (%.0f).", s, val, 0.3 * global_median_lib)
      out <- .qc_add_issue(msg, "critical", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    } else if (val < 0.5 * global_median_lib) {
      library_size_flags[[length(library_size_flags) + 1]] <- list(sample = s, level = "warning", value = val, threshold = round(0.5 * global_median_lib, 2), rule = "library_size < 50% median")
      msg <- sprintf("Sample '%s' warning library size: %.0f < 50%% median (%.0f).", s, val, 0.5 * global_median_lib)
      out <- .qc_add_issue(msg, "warning", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    }
  }

  # 7) Zero-fraction rules
  zero_fraction_flags <- list()
  for (s in sample_names) {
    z <- as.numeric(zero_fraction[[s]])
    if (z > 0.4) {
      zero_fraction_flags[[length(zero_fraction_flags) + 1]] <- list(sample = s, level = "critical", value = round(z, 4), threshold = 0.4, rule = "zero_fraction > 40%")
      msg <- sprintf("Sample '%s' critical zero fraction: %.2f%% > 40%%.", s, 100 * z)
      out <- .qc_add_issue(msg, "critical", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    } else if (z > 0.2) {
      zero_fraction_flags[[length(zero_fraction_flags) + 1]] <- list(sample = s, level = "warning", value = round(z, 4), threshold = 0.2, rule = "zero_fraction > 20%")
      msg <- sprintf("Sample '%s' warning zero fraction: %.2f%% > 20%%.", s, 100 * z)
      out <- .qc_add_issue(msg, "warning", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    }
  }

  # 8) Correlation rules + outlier rules
  correlation_flags <- list()
  outlier_samples <- character(0)

  # global mean correlation (excluding self)
  for (s in sample_names) {
    row_vals <- corr_mat[s, setdiff(sample_names, s), drop = TRUE]
    m <- mean(row_vals, na.rm = TRUE)
    if (m < 0.75) {
      correlation_flags[[length(correlation_flags) + 1]] <- list(sample = s, level = "critical", mean_correlation = round(m, 4), threshold = 0.75, rule = "mean_sample_correlation < 0.75")
      msg <- sprintf("Sample '%s' critical mean correlation: %.3f < 0.75.", s, m)
      out <- .qc_add_issue(msg, "critical", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    } else if (m < 0.85) {
      correlation_flags[[length(correlation_flags) + 1]] <- list(sample = s, level = "warning", mean_correlation = round(m, 4), threshold = 0.85, rule = "mean_sample_correlation < 0.85")
      msg <- sprintf("Sample '%s' warning mean correlation: %.3f < 0.85.", s, m)
      out <- .qc_add_issue(msg, "warning", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    }
  }

  # outlier condition A: distance to same-group centroid > median + 3*MAD
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
        msg <- sprintf("Sample '%s' outlier by centroid distance in group '%s': %.4f > %.4f (median + 3*MAD).", s, g, sample_to_centroid[[s]], thr)
        out <- .qc_add_issue(msg, "warning", qc_warnings, qc_critical)
        qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
      }
    }
  }

  # outlier condition B: within-group mean correlation < 0.85, critical if < 0.75
  for (g in names(grp_tbl)) {
    group_samples <- sample_names[grp == g]
    if (length(group_samples) < 2) next
    for (s in group_samples) {
      other <- setdiff(group_samples, s)
      m <- mean(corr_mat[s, other, drop = TRUE], na.rm = TRUE)
      if (m < 0.75) {
        outlier_samples <- unique(c(outlier_samples, s))
        msg <- sprintf("Sample '%s' outlier (critical) by within-group correlation in group '%s': %.3f < 0.75.", s, g, m)
        out <- .qc_add_issue(msg, "critical", qc_warnings, qc_critical)
        qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
      } else if (m < 0.85) {
        outlier_samples <- unique(c(outlier_samples, s))
        msg <- sprintf("Sample '%s' outlier (warning) by within-group correlation in group '%s': %.3f < 0.85.", s, g, m)
        out <- .qc_add_issue(msg, "warning", qc_warnings, qc_critical)
        qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
      }
    }
  }

  # 9) Batch effect rules
  batch_flags <- list()
  if ("batch" %in% colnames(meta)) {
    batch <- as.character(meta$batch)
    condition <- grp

    if (.qc_is_confounded(batch, condition)) {
      batch_flags[[length(batch_flags) + 1]] <- list(level = "critical", rule = "batch_condition_confounded", detail = "Batch and condition are one-to-one confounded")
      out <- .qc_add_issue("Critical batch effect: batch is confounded with condition.", "critical", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    } else {
      pc1 <- pca_fit$x[, 1]
      r2_batch <- .qc_compute_anova_r2(pc1, batch)
      r2_cond <- .qc_compute_anova_r2(pc1, condition)
      if (is.finite(r2_batch) && is.finite(r2_cond) && (r2_batch - r2_cond) >= 0.1 && r2_batch >= 0.3) {
        batch_flags[[length(batch_flags) + 1]] <- list(level = "warning", rule = "batch_dominates_pc1", detail = sprintf("PC1 R2 batch=%.3f, condition=%.3f", r2_batch, r2_cond))
        out <- .qc_add_issue(sprintf("Batch effect warning: batch explains more PC1 variance than condition (R2 %.3f vs %.3f).", r2_batch, r2_cond), "warning", qc_warnings, qc_critical)
        qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
      }
    }
  }

  # 10) Realism / anti-fake checks
  realism_flags <- list()
  sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, , drop = FALSE]
  deg_count <- nrow(sig)

  canonical <- c("TP53", "MYC", "EGFR", "VEGFA", "IL6", "CXCL8", "GAPDH", "ACTB")
  housekeeping <- c("GAPDH", "ACTB", "RPL13A")

  top_n <- min(20, nrow(res_df))
  top_genes <- toupper(res_df$gene_id[seq_len(top_n)])
  canonical_hits <- intersect(top_genes, canonical)

  if (length(canonical_hits) >= 5 || (top_n > 0 && length(canonical_hits) / top_n >= 0.5)) {
    realism_flags[[length(realism_flags) + 1]] <- list(level = "critical", rule = "canonical_overrepresented_top_genes", genes = as.list(canonical_hits))
    out <- .qc_add_issue(sprintf("Critical realism flag: canonical genes dominate top DEGs (%s).", paste(canonical_hits, collapse = ", ")), "critical", qc_warnings, qc_critical)
    qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
  } else if (length(canonical_hits) >= 3) {
    realism_flags[[length(realism_flags) + 1]] <- list(level = "warning", rule = "canonical_genes_enriched_top_genes", genes = as.list(canonical_hits))
    out <- .qc_add_issue(sprintf("Warning realism flag: canonical genes enriched in top DEGs (%s).", paste(canonical_hits, collapse = ", ")), "warning", qc_warnings, qc_critical)
    qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
  }

  sig_upper <- toupper(sig$gene_id)
  hk_hits <- intersect(sig_upper, housekeeping)
  if (length(hk_hits) > 0) {
    hk_critical <- sig[toupper(sig$gene_id) %in% hk_hits & !is.na(sig$log2FoldChange) & !is.na(sig$padj) & abs(sig$log2FoldChange) >= 2 & sig$padj < 0.01, , drop = FALSE]
    if (nrow(hk_critical) > 0) {
      realism_flags[[length(realism_flags) + 1]] <- list(level = "critical", rule = "housekeeping_strong_deg", genes = as.list(unique(toupper(hk_critical$gene_id))))
      out <- .qc_add_issue(sprintf("Critical realism flag: strong housekeeping DEGs detected (%s).", paste(unique(toupper(hk_critical$gene_id)), collapse = ", ")), "critical", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    } else {
      realism_flags[[length(realism_flags) + 1]] <- list(level = "warning", rule = "housekeeping_in_deg", genes = as.list(hk_hits))
      out <- .qc_add_issue(sprintf("Warning realism flag: housekeeping genes present in DEGs (%s).", paste(hk_hits, collapse = ", ")), "warning", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    }
  }

  if (deg_count < 10) {
    realism_flags[[length(realism_flags) + 1]] <- list(level = "warning", rule = "deg_count_too_low", value = deg_count)
    out <- .qc_add_issue(sprintf("DEG sanity warning: significant DEG count is very low (%d < 10).", deg_count), "warning", qc_warnings, qc_critical)
    qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
  }
  if (deg_count > 5000) {
    realism_flags[[length(realism_flags) + 1]] <- list(level = "warning", rule = "deg_count_too_high", value = deg_count)
    out <- .qc_add_issue(sprintf("DEG sanity warning: significant DEG count is very high (%d > 5000).", deg_count), "warning", qc_warnings, qc_critical)
    qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
  }

  pvals <- res_df$pvalue[!is.na(res_df$pvalue)]
  if (length(pvals) >= 200) {
    pvals_clipped <- pvals[pvals >= 0 & pvals <= 1]
    if (length(pvals_clipped) >= 200) {
      ks <- suppressWarnings(ks.test(pvals_clipped, "punif", 0, 1))
      p_mean <- mean(pvals_clipped)
      if (is.finite(ks$p.value) && ks$p.value > 0.2 && p_mean >= 0.45 && p_mean <= 0.55) {
        realism_flags[[length(realism_flags) + 1]] <- list(level = "warning", rule = "pvalue_distribution_too_uniform", ks_pvalue = round(ks$p.value, 6), mean_pvalue = round(p_mean, 4))
        out <- .qc_add_issue("P-value distribution warning: near-uniform distribution detected.", "warning", qc_warnings, qc_critical)
        qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
      }

      frac_tiny <- mean(pvals_clipped < 1e-10)
      if (frac_tiny > 0.2) {
        realism_flags[[length(realism_flags) + 1]] <- list(level = "warning", rule = "too_many_extremely_small_pvalues", fraction = round(frac_tiny, 4))
        out <- .qc_add_issue(sprintf("P-value distribution warning: %.2f%% p-values < 1e-10.", 100 * frac_tiny), "warning", qc_warnings, qc_critical)
        qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
      }
    }
  }

  log_mat <- log2(counts_mat + 1)
  med_by_sample <- apply(log_mat, 2, median, na.rm = TRUE)
  med_z <- as.numeric(scale(med_by_sample))
  iqr_by_sample <- apply(log_mat, 2, IQR, na.rm = TRUE)
  iqr_med <- median(iqr_by_sample, na.rm = TRUE)

  for (i in seq_along(sample_names)) {
    s <- sample_names[i]
    if (is.finite(med_z[i]) && abs(med_z[i]) > 3) {
      realism_flags[[length(realism_flags) + 1]] <- list(level = "warning", rule = "global_distribution_shift", sample = s, z_median = round(med_z[i], 4))
      out <- .qc_add_issue(sprintf("Expression consistency warning: sample '%s' median expression shifted (z=%.2f).", s, med_z[i]), "warning", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    }
    if (is.finite(iqr_by_sample[i]) && iqr_by_sample[i] < 0.5 * iqr_med) {
      realism_flags[[length(realism_flags) + 1]] <- list(level = "critical", rule = "compressed_dynamic_range", sample = s, iqr = round(iqr_by_sample[i], 4), threshold = round(0.5 * iqr_med, 4))
      out <- .qc_add_issue(sprintf("Expression consistency critical: sample '%s' compressed dynamic range (IQR %.3f < %.3f).", s, iqr_by_sample[i], 0.5 * iqr_med), "critical", qc_warnings, qc_critical)
      qc_warnings <- out$qc_warnings; qc_critical <- out$qc_critical
    }
  }

  # 11) Aggregate sample-level poor quality labels
  flagged_samples <- unique(c(
    outlier_samples,
    vapply(library_size_flags, function(x) x$sample, character(1), USE.NAMES = FALSE),
    vapply(zero_fraction_flags, function(x) x$sample, character(1), USE.NAMES = FALSE),
    vapply(correlation_flags, function(x) x$sample, character(1), USE.NAMES = FALSE)
  ))

  # 12) Assemble strict JSON
  qc_report <- list(
    outliers = as.list(unique(outlier_samples)),
    low_quality_samples = as.list(unique(flagged_samples)),
    group_balance = group_balance,
    replicates_per_group = replicates_per_group,
    library_size_flags = library_size_flags,
    zero_fraction_flags = zero_fraction_flags,
    correlation_flags = correlation_flags,
    batch_flags = batch_flags,
    realism_flags = realism_flags,
    qc_warnings = as.list(unique(qc_warnings)),
    qc_critical = as.list(unique(qc_critical)),
    pca_variance = pca_variance
  )

  qc_json_path <- file.path(results_dir, "qc_report.json")
  write(jsonlite::toJSON(qc_report, auto_unbox = TRUE, pretty = TRUE), qc_json_path)
  cat(sprintf("[QC] qc_report.json written: %s\n", qc_json_path))

  invisible(qc_report)
}
