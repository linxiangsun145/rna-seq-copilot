#!/usr/bin/env Rscript

# Example: run_stability_analysis.R
# Demonstrates end-to-end usage of the stability module in rna-seq-copilot.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

project_root <- normalizePath(file.path(getwd(), "rna-seq-copilot"), mustWork = FALSE)
if (!dir.exists(project_root)) {
  project_root <- normalizePath(getwd())
}

source(file.path(project_root, "R", "stability_resampling.R"))
source(file.path(project_root, "R", "stability_metrics.R"))
source(file.path(project_root, "R", "stability_scoring.R"))
source(file.path(project_root, "R", "stability_reporting.R"))
source(file.path(project_root, "R", "stability_plots.R"))
source(file.path(project_root, "R", "stability_main.R"))

counts_path <- file.path(project_root, "sample_data", "counts.csv")
meta_path <- file.path(project_root, "sample_data", "metadata.csv")

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
metadata_raw <- read.csv(meta_path, check.names = FALSE)

# Align metadata schema to module defaults.
if (!"sample_id" %in% names(metadata_raw)) {
  if ("sample" %in% names(metadata_raw)) {
    metadata_raw$sample_id <- metadata_raw$sample
  } else {
    metadata_raw$sample_id <- metadata_raw[[1]]
  }
}
if (!"condition" %in% names(metadata_raw)) {
  stop("metadata.csv must contain a 'condition' column for this example")
}

metadata <- metadata_raw |>
  as_tibble() |>
  select(sample_id, condition, everything())

# Determine contrast levels from metadata (treated vs control style if available).
conds <- unique(as.character(metadata$condition))
if (length(conds) < 2) {
  stop("Need at least two conditions for DESeq2 contrast")
}

numerator <- if ("treated" %in% conds) "treated" else conds[2]
denominator <- if ("control" %in% conds) "control" else conds[1]

res <- run_stability_analysis(
  counts = counts,
  metadata = metadata,
  sample_col = "sample_id",
  condition_col = "condition",
  design_formula = "~ condition",
  contrast = c("condition", numerator, denominator),
  padj_threshold = 0.05,
  log2fc_threshold = 1.0,
  top_n = 50,
  modes = NULL,
  subsample_fraction = 0.8,
  subsample_iterations = 50,
  seed = 123,
  verbose = TRUE
)

out_dir <- file.path(project_root, "sample_data", "stability_example")
exported <- export_stability_outputs(res, out_dir, prefix = "example")

message(sprintf("Stability score: %.3f (%s)", res$dataset_stability$stability_score, res$stability_badge))
message(sprintf("Headline: %s", res$stability_headline))

p1 <- plot_stability_deg_counts(res$iteration_metrics)
p2 <- plot_stability_recovery(res$iteration_metrics)
p3 <- plot_stability_top_overlap(res$iteration_metrics)
p4 <- if (nrow(res$sample_influence) > 0) plot_sample_influence(res$sample_influence) else NULL
p5 <- plot_gene_stability(res$gene_stability, highlight_n = 20)

ggplot2::ggsave(file.path(out_dir, "example_deg_count_distribution.png"), p1, width = 7, height = 4, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "example_recovery_distribution.png"), p2, width = 7, height = 4, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "example_top_overlap_distribution.png"), p3, width = 7, height = 4, dpi = 150)
if (!is.null(p4)) {
  ggplot2::ggsave(file.path(out_dir, "example_sample_influence.png"), p4, width = 7, height = 5, dpi = 150)
}
ggplot2::ggsave(file.path(out_dir, "example_gene_stability_scatter.png"), p5, width = 7, height = 5, dpi = 150)

message("Exported files:")
print(exported)
