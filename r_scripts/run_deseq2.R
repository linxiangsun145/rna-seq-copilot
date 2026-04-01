#!/usr/bin/env Rscript
# =============================================================================
# run_deseq2.R
# Main DESeq2 differential expression pipeline.
# Called by the Python backend via subprocess.
#
# Usage:
#   Rscript run_deseq2.R \
#     --counts  /path/to/counts.csv \
#     --metadata /path/to/metadata.csv \
#     --formula "~ condition" \
#     --contrast "condition,treated,control" \
#     --outdir  /path/to/job/dir
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(jsonlite)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Parse arguments ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--counts",   type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--formula",  type = "character", default = "~ condition"),
  make_option("--contrast", type = "character"),   # comma-separated: factor,num,denom
  make_option("--outdir",   type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$counts), !is.null(opt$metadata),
          !is.null(opt$contrast), !is.null(opt$outdir))

outdir     <- opt$outdir
plots_dir  <- file.path(outdir, "plots")
results_dir <- file.path(outdir, "results")
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("[DESeq2] Starting analysis\n")

# ── Load data ─────────────────────────────────────────────────────────────────
sep_counts <- if (grepl("\\.tsv$|\\.txt$", opt$counts, ignore.case = TRUE)) "\t" else ","
sep_meta   <- if (grepl("\\.tsv$|\\.txt$", opt$metadata, ignore.case = TRUE)) "\t" else ","

counts_raw <- read.csv(opt$counts,   sep = sep_counts, row.names = 1, check.names = FALSE)
meta_raw   <- read.csv(opt$metadata, sep = sep_meta,   row.names = 1, check.names = FALSE)

# Keep only common samples, align order
common_samples <- intersect(colnames(counts_raw), rownames(meta_raw))
if (length(common_samples) == 0) stop("No common samples between counts and metadata.")
counts_mat <- counts_raw[, common_samples, drop = FALSE]
meta       <- meta_raw[common_samples, , drop = FALSE]

# Coerce to integer matrix
counts_mat <- as.matrix(round(counts_mat))
mode(counts_mat) <- "integer"

cat(sprintf("[DESeq2] %d genes × %d samples\n", nrow(counts_mat), ncol(counts_mat)))

# ── Parse design formula ──────────────────────────────────────────────────────
design_formula <- as.formula(opt$formula)
design_vars    <- all.vars(design_formula)

for (v in design_vars) {
  if (!v %in% colnames(meta)) {
    stop(sprintf("Design variable '%s' not found in metadata columns: %s",
                 v, paste(colnames(meta), collapse = ", ")))
  }
  meta[[v]] <- factor(meta[[v]])
}

# ── Parse contrast ─────────────────────────────────────────────────────────────
contrast_parts <- trimws(strsplit(opt$contrast, ",")[[1]])
if (length(contrast_parts) != 3) stop("--contrast must be 'factor,numerator,denominator'")
contrast_factor <- contrast_parts[1]
contrast_num    <- contrast_parts[2]
contrast_denom  <- contrast_parts[3]

# Validate that contrast levels exist in the (now-factored) metadata
if (!contrast_factor %in% colnames(meta)) {
  stop(sprintf("Contrast factor '%s' not found in metadata. Available columns: %s",
               contrast_factor, paste(colnames(meta), collapse = ", ")))
}
if (!contrast_num %in% levels(meta[[contrast_factor]])) {
  stop(sprintf("Contrast numerator '%s' is not a level of '%s'. Available levels: %s",
               contrast_num, contrast_factor,
               paste(levels(meta[[contrast_factor]]), collapse = ", ")))
}
if (!contrast_denom %in% levels(meta[[contrast_factor]])) {
  stop(sprintf("Contrast denominator '%s' is not a level of '%s'. Available levels: %s",
               contrast_denom, contrast_factor,
               paste(levels(meta[[contrast_factor]]), collapse = ", ")))
}

# ── Build DESeqDataSet ────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = meta,
  design    = design_formula
)

# Low-count filter: keep genes with total count >= 10
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat(sprintf("[DESeq2] After low-count filter: %d genes remain\n", sum(keep)))

# ── Run DESeq2 (robust fitType fallback) ─────────────────────────────────────
fit_candidates <- c("parametric", "local", "mean")
dds_ok <- NULL
fit_used <- NULL
last_err <- NULL

for (ft in fit_candidates) {
  cat(sprintf("[DESeq2] Trying fitType='%s'...\n", ft))
  attempt <- tryCatch(
    {
      DESeq(dds, fitType = ft)
    },
    error = function(e) {
      last_err <<- e$message
      NULL
    }
  )

  if (!is.null(attempt)) {
    dds_ok <- attempt
    fit_used <- ft
    cat(sprintf("[DESeq2] fitType='%s' succeeded\n", ft))
    break
  }

  cat(sprintf("[DESeq2] fitType='%s' failed: %s\n", ft, last_err))
}

if (is.null(dds_ok)) {
  stop(sprintf(
    "DESeq2 failed for all fitType candidates (%s). Last error: %s",
    paste(fit_candidates, collapse = ", "),
    ifelse(is.null(last_err), "unknown", last_err)
  ))
}

dds <- dds_ok
res <- results(dds,
               contrast = c(contrast_factor, contrast_num, contrast_denom),
               alpha    = 0.05)
res  <- lfcShrink(dds,
                  contrast = c(contrast_factor, contrast_num, contrast_denom),
                  res      = res,
                  type     = if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal")

cat(sprintf("[DESeq2] %d significant DEGs (padj < 0.05), fitType=%s\n",
            sum(!is.na(res$padj) & res$padj < 0.05),
            fit_used))

# ── Save DEG results CSV ───────────────────────────────────────────────────────
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
# lfcShrink drops 'stat'; select only columns that reliably exist
keep_cols <- intersect(c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"),
                       colnames(res_df))
res_df <- res_df[, keep_cols]
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.csv(res_df, file.path(results_dir, "deg_results.csv"), row.names = FALSE, quote = FALSE)

# ── VST for visualisations ────────────────────────────────────────────────────
vsd <- tryCatch(
  vst(dds, blind = FALSE),
  error = function(e) varianceStabilizingTransformation(dds, blind = FALSE)
)

# ── Source helper plot scripts ─────────────────────────────────────────────────
# Robust script-dir detection that works with Rscript on all platforms
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .args[grepl("^--file=", .args)]
if (length(.file_arg) > 0) {
  script_dir <- normalizePath(dirname(sub("^--file=", "", .file_arg[1])))
} else {
  script_dir <- getwd()  # interactive / source() fallback
}

source(file.path(script_dir, "plot_pca.R"))
source(file.path(script_dir, "plot_volcano.R"))
source(file.path(script_dir, "plot_ma.R"))
source(file.path(script_dir, "plot_heatmap.R"))
source(file.path(script_dir, "summarize_results.R"))
source(file.path(script_dir, "qc_analysis.R"))

# Stability module (optional runtime layer; non-fatal if unavailable/fails)
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
stability_dir <- file.path(project_root, "R")
stability_scripts <- c(
  "stability_resampling.R",
  "stability_metrics.R",
  "stability_scoring.R",
  "stability_reporting.R",
  "stability_plots.R",
  "stability_main.R"
)
stability_available <- all(file.exists(file.path(stability_dir, stability_scripts)))
if (stability_available) {
  for (f in stability_scripts) {
    source(file.path(stability_dir, f))
  }
} else {
  cat("[Stability] Module scripts not fully available; skipping stability analysis.\n")
}

# ── Generate DEG plots ────────────────────────────────────────────────────────
plot_pca(vsd,     meta, contrast_factor, file.path(plots_dir, "pca.png"))
plot_volcano(res_df,                     file.path(plots_dir, "volcano.png"),
             contrast_num, contrast_denom)
plot_ma(res_df,                          file.path(plots_dir, "ma.png"))  # res_df: avoids lfcShrink stat-column drop
plot_heatmap(vsd,  meta,                 file.path(plots_dir, "heatmap.png"))

# ── QC module ─────────────────────────────────────────────────────────────────
run_qc(
  dds             = dds,
  vsd             = vsd,
  counts_mat      = counts_mat,
  meta            = meta,
  contrast_factor = contrast_factor,
  res_df          = res_df,
  plots_dir       = plots_dir,
  results_dir     = results_dir
)

# ── Stability / perturbation analysis (LOO + subsampling defaults) ──────────
stability_assessment <- NULL
if (stability_available && exists("run_stability_analysis", mode = "function")) {
  cat("[Stability] Running perturbation stability analysis...\n")
  stability_assessment <- tryCatch({
    metadata_stab <- as.data.frame(meta)
    metadata_stab$sample_id <- rownames(meta)

    qc_context <- list()
    pca_sep_for_stability <- "unknown"
    qc_json_path <- file.path(results_dir, "qc_report.json")
    if (file.exists(qc_json_path)) {
      qc_obj <- tryCatch(jsonlite::fromJSON(qc_json_path, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(qc_obj) && is.list(qc_obj)) {
        group_qc <- qc_obj$group_qc
        group_corr <- c()
        if (!is.null(group_qc) && length(group_qc) > 0) {
          for (nm in names(group_qc)) {
            m <- suppressWarnings(as.numeric(group_qc[[nm]]$mean_correlation))
            if (!is.na(m)) group_corr <- c(group_corr, m)
          }
        }
        if (length(group_corr) > 0) {
          qc_context$mean_within_group_correlation <- mean(group_corr, na.rm = TRUE)
        }

        pc1 <- suppressWarnings(as.numeric(qc_obj$pca_variance$PC1))
        if (!is.na(pc1)) {
          pca_sep_for_stability <- if (pc1 >= 0.30) "clear" else if (pc1 >= 0.15) "weak" else "none"
        }
      }
    }

    stab_results <- run_stability_analysis(
      counts = counts_mat,
      metadata = metadata_stab,
      sample_col = "sample_id",
      condition_col = contrast_factor,
      design_formula = opt$formula,
      contrast = c(contrast_factor, contrast_num, contrast_denom),
      padj_threshold = 0.05,
      log2fc_threshold = 1.0,
      top_n = 50,
      modes = NULL,
      subsample_fraction = 0.8,
      subsample_iterations = 50,
      seed = 123,
      pca_separation = pca_sep_for_stability,
      qc_context = qc_context,
      verbose = FALSE
    )

    stability_outdir <- file.path(results_dir, "stability")
    export_stability_outputs(stab_results, stability_outdir, prefix = "stability")

    ds <- stab_results$dataset_stability
    list(
      stability_score = as.numeric(ds$final_stability_score),
      deg_stability_score = as.numeric(ds$deg_stability_score),
      effect_stability_score = as.numeric(ds$effect_stability_score),
      final_stability_score = as.numeric(ds$final_stability_score),
      stability_level = as.character(ds$stability_level),
      stability_badge = as.character(stab_results$stability_badge),
      signal_state = as.character(stab_results$signal_state),
      stability_run_status = as.character(stab_results$stability_run_status),
      deg_metrics_applicable = as.logical(stab_results$deg_metrics_applicable),
      stability_mode = as.character(stab_results$stability_mode),
      top_rank_overlap = as.numeric(ds$mean_top_n_overlap),
      influence_mode = as.character(stab_results$influence_mode),
      sample_influence = if (nrow(stab_results$sample_influence) > 0) {
        list(
          removed_sample = as.character(stab_results$sample_influence$removed_sample[1]),
          influence_mode = as.character(stab_results$sample_influence$influence_mode[1]),
          influence_score = as.numeric(stab_results$sample_influence$influence_score[1])
        )
      } else {
        NULL
      },
      stability_headline = as.character(stab_results$stability_headline),
      key_stability_findings = as.list(as.character(stab_results$narrative$key_stability_findings)),
      warnings = as.list(as.character(stab_results$warnings)),
      summary_text = as.character(stab_results$narrative$summary_text),
      directionality_text = as.character(stab_results$narrative$directionality_text),
      top_rank_definition_text = as.character(stab_results$narrative$top_rank_definition_text),
      pca_qc_conflict_text = as.character(stab_results$narrative$pca_qc_conflict_text),
      sample_influence_text = as.character(stab_results$narrative$sample_influence_text),
      reference_deg_count = as.integer(ds$reference_deg_count),
      mean_deg_recovery_rate = as.numeric(ds$mean_deg_recovery_rate),
      mean_top_n_overlap = as.numeric(ds$mean_top_n_overlap),
      mean_log2fc_correlation = as.numeric(ds$mean_log2fc_correlation),
      mean_log2fc_rmse = as.numeric(ds$mean_log2fc_rmse),
      signal_collapse_fraction = as.numeric(ds$fraction_of_iterations_with_signal_collapse),
      stability_penalty = list(
        stability_penalty = as.numeric(stab_results$stability_penalty$stability_penalty),
        penalty_level = as.character(stab_results$stability_penalty$penalty_level),
        penalty_reason = as.character(stab_results$stability_penalty$penalty_reason)
      )
    )
  }, error = function(e) {
    cat(sprintf("[Stability] Warning: stability analysis failed: %s\n", e$message))
    list(
      stability_score = NA_real_,
      deg_stability_score = NA_real_,
      effect_stability_score = NA_real_,
      final_stability_score = NA_real_,
      stability_level = "unknown",
      stability_badge = "not_applicable",
      signal_state = "unknown",
      stability_run_status = "failed",
      deg_metrics_applicable = FALSE,
      stability_mode = "not_applicable",
      top_rank_overlap = NA_real_,
      influence_mode = "not_applicable",
      sample_influence = NULL,
      stability_headline = "Stability analysis could not be completed due to a technical failure in the perturbation assessment module.",
      key_stability_findings = as.list(c("Perturbation stability module failed technically; no reliable stability layer was produced.")),
      warnings = as.list(c("STABILITY_ANALYSIS_FAILED")),
      summary_text = "Stability analysis could not be completed due to a technical failure in the perturbation assessment module.",
      directionality_text = "Effect-size direction consistency was not produced because the perturbation analysis failed technically.",
      top_rank_definition_text = "Top-ranked genes are defined based on statistical ranking and may not meet differential expression thresholds.",
      pca_qc_conflict_text = "",
      sample_influence_text = "Single-sample influence was not produced because the perturbation analysis failed technically.",
      stability_penalty = list(
        stability_penalty = 0,
        penalty_level = "none",
        penalty_reason = "No stability penalty applied because stability analysis failed."
      )
    )
  })
}

# ── Write summary JSON ────────────────────────────────────────────────────────
summary_json <- build_summary(
  dds        = dds,
  res_df     = res_df,
  vsd        = vsd,
  meta       = meta,
  contrast   = paste(contrast_parts, collapse = "_vs_"),
  group_col  = contrast_factor,
  stability_assessment = stability_assessment
)
write(toJSON(summary_json, auto_unbox = TRUE, pretty = TRUE, na = "null"),
      file.path(results_dir, "summary.json"))

cat("[DESeq2] Done. Outputs written to:", outdir, "\n")
quit(status = 0, save = "no")
