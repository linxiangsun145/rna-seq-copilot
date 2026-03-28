# summarize_results.R — Build the structured summary JSON consumed by LLM

build_summary <- function(dds, res_df, vsd, meta, contrast, group_col) {

  # Groups
  if (group_col %in% colnames(colData(dds))) {
    group_levels <- levels(colData(dds)[[group_col]])
  } else {
    group_levels <- unique(as.character(meta[[group_col]]))
  }

  # DEG counts
  sig    <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
  deg_up   <- sum(sig$log2FoldChange >  1, na.rm = TRUE)
  deg_down <- sum(sig$log2FoldChange < -1, na.rm = TRUE)
  top_genes <- head(sig$gene_id[order(sig$padj)], 20)

  # PCA separation heuristic — guard against <2 samples per group
  pca_sep <- "unknown"
  tryCatch({
    pca_data <- plotPCA(vsd, intgroup = group_col, returnData = TRUE)
    pct1     <- round(100 * attr(pca_data, "percentVar")[1])
    pca_sep  <- if (pct1 >= 30) "clear" else if (pct1 >= 15) "weak" else "none"
  }, error = function(e) invisible(NULL))

  # Outlier detection via Cook's distance (only if assay is available post-DESeq)
  outliers <- character(0)
  tryCatch({
    if ("cooks" %in% assayNames(dds)) {
      cooks_mat <- assays(dds)[["cooks"]]
      cooks_max <- apply(cooks_mat, 2, max, na.rm = TRUE)
      p <- ncol(model.matrix(design(dds), colData(dds)))
      n <- ncol(dds)
      threshold <- qf(0.99, p, n - p)
      outliers  <- colnames(dds)[cooks_max > threshold]
    }
  }, error = function(e) invisible(NULL))

  # Warnings
  warnings <- character(0)
  if (deg_up + deg_down == 0) warnings <- c(warnings, "No significant DEGs found. Consider relaxing thresholds or checking experimental design.")
  if (pca_sep == "none")  warnings <- c(warnings, "PC1 explains <15% variance — groups do not separate clearly on PCA.")
  if (length(outliers) > 0) warnings <- c(warnings, paste("Potential outlier samples:", paste(outliers, collapse = ", ")))

  # Data issues — compare pre-vs-post low-count filter using nrow(res_df) as pre-filter approximation
  data_issues <- character(0)
  n_total_before_filter <- nrow(res_df)  # res_df was built from filtered dds so this is post-filter
  # Check if many genes were filtered (DESeq2 internal independent filtering removes NA padj rows)
  n_na_padj <- sum(is.na(res_df$padj))
  if (n_na_padj > 0.3 * nrow(res_df)) {
    data_issues <- c(data_issues, sprintf(
      "%d/%d genes (%.0f%%) have NA adjusted p-values due to low counts or independent filtering.",
      n_na_padj, nrow(res_df), 100 * n_na_padj / nrow(res_df)
    ))
  }

  warning_items <- list()
  .add_stat_warning <- function(severity, code, message, sample = NULL, metric = NULL) {
    warning_items[[length(warning_items) + 1]] <<- list(
      type = "statistical",
      severity = severity,
      code = code,
      message = message,
      sample = sample,
      metric = metric
    )
  }

  if (deg_up + deg_down == 0) {
    .add_stat_warning("warning", "no_significant_deg", "No significant differentially expressed genes detected under current thresholds.", metric = "deg_count")
  }
  if (pca_sep == "none") {
    .add_stat_warning("warning", "weak_pca_separation", "PC1 explains <15% variance and group separation is weak.", metric = "pca_separation")
  }
  if (length(outliers) > 0) {
    .add_stat_warning("warning", "summary_outlier_hint", paste("Potential outlier samples:", paste(outliers, collapse = ", ")), metric = "outliers")
  }
  if (n_na_padj > 0.3 * nrow(res_df)) {
    .add_stat_warning("warning", "high_na_padj_fraction", sprintf("High NA padj fraction detected (%d/%d).", n_na_padj, nrow(res_df)), metric = "padj")
  }

  list(
    n_samples       = ncol(dds),
    groups          = as.list(group_levels),
    contrast        = contrast,
    outliers        = as.list(outliers),
    pca_separation  = pca_sep,
    deg_up          = deg_up,
    deg_down        = deg_down,
    top_genes       = as.list(top_genes),
    warnings        = as.list(warnings),
    data_issues     = as.list(data_issues),
    warning_items   = warning_items
  )
}
