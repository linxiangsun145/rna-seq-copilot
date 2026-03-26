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

# ── Run DESeq2 ────────────────────────────────────────────────────────────────
dds  <- DESeq(dds)
res  <- results(dds,
                contrast = c(contrast_factor, contrast_num, contrast_denom),
                alpha    = 0.05)
res  <- lfcShrink(dds,
                  contrast = c(contrast_factor, contrast_num, contrast_denom),
                  res      = res,
                  type     = if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal")

cat(sprintf("[DESeq2] %d significant DEGs (padj < 0.05)\n",
            sum(!is.na(res$padj) & res$padj < 0.05)))

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

# ── Generate plots ────────────────────────────────────────────────────────────
plot_pca(vsd,     meta, contrast_factor, file.path(plots_dir, "pca.png"))
plot_volcano(res_df,                     file.path(plots_dir, "volcano.png"),
             contrast_num, contrast_denom)
plot_ma(res_df,                          file.path(plots_dir, "ma.png"))  # res_df: avoids lfcShrink stat-column drop
plot_heatmap(vsd,  meta,                 file.path(plots_dir, "heatmap.png"))

# ── Write summary JSON ────────────────────────────────────────────────────────
summary_json <- build_summary(
  dds        = dds,
  res_df     = res_df,
  vsd        = vsd,
  meta       = meta,
  contrast   = paste(contrast_parts, collapse = "_vs_"),
  group_col  = contrast_factor
)
write(toJSON(summary_json, auto_unbox = TRUE, pretty = TRUE),
      file.path(results_dir, "summary.json"))

cat("[DESeq2] Done. Outputs written to:", outdir, "\n")
quit(status = 0, save = "no")
