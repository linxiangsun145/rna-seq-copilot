# plot_pca.R — PCA plot using VST-transformed data

plot_pca <- function(vsd, meta, group_col, outfile) {
  pca_data <- plotPCA(vsd, intgroup = group_col, returnData = TRUE)
  pct_var  <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(x = PC1, y = PC2,
                             color = .data[[group_col]],
                             label = name)) +
    geom_point(size = 4, alpha = 0.85) +
    geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title  = "PCA of VST-normalised counts",
      x      = paste0("PC1: ", pct_var[1], "% variance"),
      y      = paste0("PC2: ", pct_var[2], "% variance"),
      color  = group_col
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title   = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(outfile, plot = p, width = 7, height = 6, dpi = 150)
  invisible(p)
}
