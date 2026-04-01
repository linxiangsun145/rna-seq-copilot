# stability_main.R
# Main orchestration for RNA-seq perturbation stability analysis.

#' Ensure required packages are available.
#' @return invisible TRUE
check_stability_dependencies <- function() {
  req <- c("DESeq2", "dplyr", "tibble", "purrr", "tidyr", "jsonlite")
  missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required packages for stability module: %s",
      paste(missing, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#' Run DESeq2 once and return standardized results table.
#' @param counts Integer matrix genes x samples.
#' @param metadata data.frame with aligned sample rows.
#' @param sample_col Metadata sample ID column.
#' @param condition_col Metadata condition column.
#' @param design_formula Formula object or string.
#' @param contrast Character vector c(factor, numerator, denominator).
#' @param min_total_count Gene filter threshold.
#' @param use_lfc_shrink Whether to attempt lfcShrink.
#' @return tibble DE results with stable schema.
run_deseq_once <- function(counts,
                           metadata,
                           sample_col = "sample_id",
                           condition_col = "condition",
                           design_formula = "~ condition",
                           contrast,
                           min_total_count = 10,
                           use_lfc_shrink = TRUE) {
  check_stability_dependencies()

  design_formula <- stats::as.formula(design_formula)
  design_vars <- all.vars(design_formula)

  coldata <- as.data.frame(metadata)
  rownames(coldata) <- coldata[[sample_col]]

  counts <- as.matrix(counts)
  counts <- counts[, rownames(coldata), drop = FALSE]

  for (v in design_vars) {
    if (!v %in% colnames(coldata)) {
      stop(sprintf("Design variable '%s' not found in metadata", v))
    }
    coldata[[v]] <- as.factor(coldata[[v]])
  }

  if (length(contrast) != 3) {
    stop("contrast must be c(factor, numerator, denominator)")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = design_formula
  )

  keep <- rowSums(DESeq2::counts(dds)) >= min_total_count
  dds <- dds[keep, ]
  if (nrow(dds) < 2) {
    stop("Too few genes after filtering for DESeq2 model fitting")
  }

  fit_types <- c("parametric", "local", "mean")
  dds_fit <- NULL
  last_err <- NULL
  fit_used <- NA_character_

  for (ft in fit_types) {
    dds_try <- tryCatch(
      DESeq2::DESeq(dds, fitType = ft),
      error = function(e) {
        last_err <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(dds_try)) {
      dds_fit <- dds_try
      fit_used <- ft
      break
    }
  }

  if (is.null(dds_fit)) {
    stop(sprintf("DESeq2 failed for all fitType options. Last error: %s", last_err))
  }

  res <- DESeq2::results(
    dds_fit,
    contrast = contrast,
    alpha = 0.05
  )

  if (isTRUE(use_lfc_shrink)) {
    shrink_type <- if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal"
    res <- tryCatch(
      DESeq2::lfcShrink(
        dds_fit,
        contrast = contrast,
        res = res,
        type = shrink_type
      ),
      error = function(e) {
        message(sprintf("lfcShrink failed (%s); using unshrunk results.", conditionMessage(e)))
        res
      }
    )
  }

  out <- as.data.frame(res)
  out$gene_id <- rownames(out)

  # Ensure stable schema for downstream JSON/report integration.
  for (nm in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
    if (!nm %in% colnames(out)) out[[nm]] <- NA_real_
  }

  tibble::as_tibble(out[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]) |>
    dplyr::mutate(fit_type = fit_used)
}

#' Run full-data reference DESeq2 and derive baseline objects.
#' @return List reference object for iteration comparisons.
run_reference_deseq <- function(counts,
                                metadata,
                                sample_col = "sample_id",
                                condition_col = "condition",
                                design_formula = "~ condition",
                                contrast,
                                padj_threshold = 0.05,
                                log2fc_threshold = 1.0,
                                top_n = 50) {
  ref_raw <- run_deseq_once(
    counts = counts,
    metadata = metadata,
    sample_col = sample_col,
    condition_col = condition_col,
    design_formula = design_formula,
    contrast = contrast
  )

  ref_tbl <- classify_significance(
    ref_raw,
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold
  )

  ref_deg_set <- ref_tbl |>
    dplyr::filter(.data$is_strong_deg) |>
    dplyr::pull(.data$gene_id)

  ref_top <- rank_top_genes(ref_tbl, top_n = top_n)

  list(
    reference_table = ref_tbl,
    reference_deg_set = ref_deg_set,
    reference_top = ref_top,
    reference_summary = list(
      n_genes_tested = nrow(ref_tbl),
      n_strong_deg = length(ref_deg_set),
      n_weak_effect = sum(ref_tbl$ref_class == "weak_effect", na.rm = TRUE),
      top_n = top_n,
      padj_threshold = padj_threshold,
      log2fc_threshold = log2fc_threshold
    )
  )
}

#' Run one perturbation iteration and compute comparison metrics.
#' @param plan_row Single-row perturbation plan (from build_resample_plan).
#' @param ref_obj Reference object from run_reference_deseq.
#' @return List with iteration_metrics and gene_delta.
run_single_iteration <- function(plan_row,
                                 counts,
                                 metadata,
                                 ref_obj,
                                 sample_col = "sample_id",
                                 condition_col = "condition",
                                 design_formula = "~ condition",
                                 contrast,
                                 padj_threshold = 0.05,
                                 log2fc_threshold = 1.0,
                                 top_n = 50) {
  included <- unlist(plan_row$included_samples[[1]], use.names = FALSE)
  if (length(included) < 4) {
    stop(sprintf("Iteration %s has < 4 included samples", plan_row$iteration_id))
  }

  metadata_iter <- metadata |>
    dplyr::filter(.data[[sample_col]] %in% included)

  # Guard against degenerate per-group sample counts.
  min_group <- metadata_iter |>
    dplyr::count(.data[[condition_col]], name = "n") |>
    dplyr::pull(.data$n) |>
    min()

  if (is.infinite(min_group) || min_group < 2) {
    stop(sprintf(
      "Iteration %s has insufficient group size after perturbation (min group n=%s)",
      plan_row$iteration_id,
      as.character(min_group)
    ))
  }

  counts_iter <- counts[, metadata_iter[[sample_col]], drop = FALSE]

  iter_tbl <- run_deseq_once(
    counts = counts_iter,
    metadata = metadata_iter,
    sample_col = sample_col,
    condition_col = condition_col,
    design_formula = design_formula,
    contrast = contrast
  )

  cmp <- compare_iteration_to_reference(
    iter_tbl = iter_tbl,
    ref_obj = ref_obj,
    top_n = top_n,
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold
  )

  iter_metrics <- cmp$iter_metrics |>
    dplyr::mutate(
      iteration_id = plan_row$iteration_id,
      mode = plan_row$mode,
      removed_sample = dplyr::coalesce(plan_row$removed_sample, NA_character_),
      included_samples = paste(included, collapse = ";"),
      excluded_samples = paste(unlist(plan_row$excluded_samples[[1]], use.names = FALSE), collapse = ";")
    ) |>
    dplyr::select(
      .data$iteration_id,
      .data$mode,
      .data$removed_sample,
      .data$included_samples,
      .data$excluded_samples,
      dplyr::everything()
    )

  list(
    iteration_metrics = iter_metrics,
    gene_delta = cmp$gene_delta,
    iter_top = cmp$iter_top$top_combined
  )
}

#' Run end-to-end stability analysis.
#' @return Nested list containing machine-readable outputs, warnings, and narratives.
run_stability_analysis <- function(counts,
                                   metadata,
                                   sample_col = "sample_id",
                                   condition_col = "condition",
                                   design_formula = "~ condition",
                                   contrast,
                                   padj_threshold = 0.05,
                                   log2fc_threshold = 1.0,
                                   top_n = 50,
                                   modes = NULL,
                                   subsample_fraction = 0.8,
                                   subsample_iterations = 50,
                                   seed = 123L,
                                   pca_separation = "unknown",
                                   qc_context = list(),
                                   parallelizable_hint = TRUE,
                                   verbose = TRUE) {
  check_stability_dependencies()

  aligned <- validate_stability_inputs(
    counts = counts,
    metadata = metadata,
    sample_col = sample_col,
    condition_col = condition_col
  )

  counts2 <- aligned$counts
  metadata2 <- aligned$metadata

  default_modes <- choose_stability_modes(
    metadata2,
    condition_col = condition_col,
    default_subsample_iter = subsample_iterations
  )

  if (is.null(modes)) {
    modes <- default_modes$modes
    subsample_iterations <- default_modes$subsample_iterations
  }

  if (isTRUE(verbose)) {
    message(default_modes$rationale)
    if (isTRUE(parallelizable_hint)) {
      message("Iteration loop is currently sequential but designed for future parallel execution.")
    }
  }

  ref_obj <- run_reference_deseq(
    counts = counts2,
    metadata = metadata2,
    sample_col = sample_col,
    condition_col = condition_col,
    design_formula = design_formula,
    contrast = contrast,
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold,
    top_n = top_n
  )

  plan <- build_resample_plan(
    metadata = metadata2,
    sample_col = sample_col,
    condition_col = condition_col,
    modes = modes,
    subsample_fraction = subsample_fraction,
    subsample_iterations = subsample_iterations,
    seed = seed
  )

  iter_out <- vector("list", nrow(plan))
  failures <- list()

  for (i in seq_len(nrow(plan))) {
    row_i <- plan[i, , drop = FALSE]
    out_i <- tryCatch(
      run_single_iteration(
        plan_row = row_i,
        counts = counts2,
        metadata = metadata2,
        ref_obj = ref_obj,
        sample_col = sample_col,
        condition_col = condition_col,
        design_formula = design_formula,
        contrast = contrast,
        padj_threshold = padj_threshold,
        log2fc_threshold = log2fc_threshold,
        top_n = top_n
      ),
      error = function(e) {
        failures[[length(failures) + 1L]] <<- list(
          iteration_id = row_i$iteration_id,
          mode = row_i$mode,
          error = conditionMessage(e)
        )
        NULL
      }
    )
    iter_out[[i]] <- out_i
  }

  iter_ok <- iter_out[!vapply(iter_out, is.null, logical(1))]
  if (length(iter_ok) == 0) {
    stop("All perturbation iterations failed; no stability metrics produced")
  }

  iteration_metrics <- dplyr::bind_rows(lapply(iter_ok, `[[`, "iteration_metrics"))
  gene_delta_list <- lapply(iter_ok, `[[`, "gene_delta")
  top_hits_list <- lapply(iter_ok, `[[`, "iter_top")

  gene_stability <- summarize_gene_stability(
    reference_table = ref_obj$reference_table,
    gene_delta_list = gene_delta_list,
    top_gene_hits_list = top_hits_list,
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold
  )

  dataset_base <- summarize_dataset_stability(
    iteration_metrics = iteration_metrics,
    reference_deg_count = ref_obj$reference_summary$n_strong_deg
  )

  sample_influence <- compute_sample_influence(
    iteration_metrics,
    deg_metrics_applicable = isTRUE(dataset_base$deg_metrics_applicable)
  )

  score_obj <- compute_stability_score(dataset_base)
  dataset_stability <- utils::modifyList(dataset_base, score_obj)

  warnings <- generate_stability_warnings(
    dataset_stability = dataset_stability,
    gene_stability = gene_stability,
    sample_influence = sample_influence,
    pca_separation = pca_separation,
    qc_context = qc_context
  )

  narrative <- interpret_stability_results(
    dataset_stability = dataset_stability,
    warnings = warnings,
    sample_influence = sample_influence,
    pca_separation = pca_separation,
    qc_context = qc_context
  )

  penalty <- compute_stability_penalty(list(dataset_stability = dataset_stability))

  list(
    config = list(
      sample_col = sample_col,
      condition_col = condition_col,
      design_formula = as.character(design_formula),
      contrast = contrast,
      padj_threshold = padj_threshold,
      log2fc_threshold = log2fc_threshold,
      top_n = top_n,
      modes = modes,
      subsample_fraction = subsample_fraction,
      subsample_iterations = subsample_iterations,
      seed = seed
    ),
    reference_summary = ref_obj$reference_summary,
    iteration_metrics = iteration_metrics,
    gene_stability = gene_stability,
    sample_influence = sample_influence,
    dataset_stability = dataset_stability,
    warnings = warnings,
    narrative = narrative,
    stability_headline = narrative$stability_headline,
    stability_badge = narrative$stability_badge,
    signal_state = dataset_stability$signal_state,
    deg_metrics_applicable = dataset_stability$deg_metrics_applicable,
    stability_mode = dataset_stability$stability_mode,
    stability_penalty = penalty,
    failures = failures
  )
}

#' Export stability outputs to CSV/JSON/RDS for pipeline integration.
#' @param stability_results Output object from run_stability_analysis.
#' @param outdir Output directory path.
#' @param prefix Filename prefix.
#' @return Invisible list of written file paths.
export_stability_outputs <- function(stability_results,
                                     outdir,
                                     prefix = "stability") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  f_iter <- file.path(outdir, sprintf("%s_iteration_metrics.csv", prefix))
  f_gene <- file.path(outdir, sprintf("%s_gene_stability.csv", prefix))
  f_inf <- file.path(outdir, sprintf("%s_sample_influence.csv", prefix))
  f_json <- file.path(outdir, sprintf("%s_summary.json", prefix))
  f_rds <- file.path(outdir, sprintf("%s_full.rds", prefix))

  utils::write.csv(stability_results$iteration_metrics, f_iter, row.names = FALSE)
  utils::write.csv(stability_results$gene_stability, f_gene, row.names = FALSE)
  utils::write.csv(stability_results$sample_influence, f_inf, row.names = FALSE)

  json_payload <- list(
    config = stability_results$config,
    reference_summary = stability_results$reference_summary,
    dataset_stability = stability_results$dataset_stability,
    warnings = stability_results$warnings,
    narrative = stability_results$narrative,
    stability_headline = stability_results$stability_headline,
    stability_badge = stability_results$stability_badge,
    stability_penalty = stability_results$stability_penalty,
    failures = stability_results$failures
  )

  jsonlite::write_json(json_payload, f_json, pretty = TRUE, auto_unbox = TRUE, na = "null")
  saveRDS(stability_results, f_rds)

  invisible(list(
    iteration_metrics_csv = f_iter,
    gene_stability_csv = f_gene,
    sample_influence_csv = f_inf,
    summary_json = f_json,
    full_rds = f_rds
  ))
}
