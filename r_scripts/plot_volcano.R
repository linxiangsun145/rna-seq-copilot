# plot_volcano.R — Volcano plot (log2FC vs -log10 padj)

plot_volcano <- function(res_df, outfile, numerator = "treated", denominator = "control") {
  df <- res_df[!is.na(res_df$padj), ]
  df$neg_log10_padj <- -log10(pmax(df$padj, 1e-300))

  # Classify genes
  df$category <- "NS"
  df$category[df$padj < 0.05 & df$log2FoldChange >  1] <- "Up"
  df$category[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"

  # Top 15 genes to label (guard against empty significant set)
  sig_df <- df[df$category != "NS", ]
  if (nrow(sig_df) > 0) {
    sig_df      <- sig_df[order(sig_df$padj), ]
    label_genes <- head(sig_df$gene_id, 15)
  } else {
    label_genes <- character(0)
  }
  df$label <- ifelse(df$gene_id %in% label_genes, df$gene_id, "")

  colours <- c("Up" = "#dc2626", "Down" = "#2563eb", "NS" = "#94a3b8")

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = category, label = label)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(size = 2.8, max.overlaps = 20, show.legend = FALSE,
                   segment.size = 0.3, color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", linewidth = 0.4) +
    scale_color_manual(values = colours,
                       labels = c("Up" = sprintf("Up (%d)", sum(df$category == "Up")),
                                  "Down" = sprintf("Down (%d)", sum(df$category == "Down")),
                                  "NS" = "NS")) +
    labs(
      title    = sprintf("Volcano Plot: %s vs %s", numerator, denominator),
      subtitle = "padj < 0.05 & |log2FC| > 1",
      x        = "log2 Fold Change",
      y        = "-log10(padj)",
      color    = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "top"
    )

  ggsave(outfile, plot = p, width = 7, height = 6, dpi = 150)
  invisible(p)
}
