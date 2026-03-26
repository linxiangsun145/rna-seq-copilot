# plot_ma.R — MA plot (log2FC vs mean expression)
# Note: accepts res_df (data.frame from write.csv output), not the raw DESeqResults
# object, so it works correctly after lfcShrink which drops the 'stat' column.

plot_ma <- function(res, outfile) {
  # Accept either a DESeqResults object or a plain data.frame
  df <- if (is.data.frame(res)) res else as.data.frame(res)
  df <- df[!is.na(df$padj), ]
  df$significant <- df$padj < 0.05

  p <- ggplot(df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.4) +
    scale_color_manual(values = c("TRUE" = "#dc2626", "FALSE" = "#94a3b8"),
                       labels = c("TRUE" = "padj < 0.05", "FALSE" = "NS")) +
    labs(
      title  = "MA Plot",
      x      = "log10(Mean normalised count + 1)",
      y      = "log2 Fold Change",
      color  = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "top"
    )

  ggsave(outfile, plot = p, width = 7, height = 6, dpi = 150)
  invisible(p)
}
