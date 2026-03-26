# plot_heatmap.R — Sample-to-sample distance heatmap using VST data

plot_heatmap <- function(vsd, meta, outfile) {
  sampleDists       <- dist(t(assay(vsd)))
  sampleDistMatrix  <- as.matrix(sampleDists)

  # Build column annotation: take the first column of meta as a data.frame
  # (pheatmap requires an explicit data.frame, not a matrix or tibble)
  ann_col <- NULL
  if (ncol(meta) > 0) {
    ann_col <- data.frame(
      row.names = rownames(meta),
      Group     = as.character(meta[, 1])
    )
  }

  colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

  png(outfile, width = 800, height = 700, res = 120)
  on.exit(dev.off(), add = TRUE)
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col                      = colours,
    annotation_col           = ann_col,
    main                     = "Sample-to-sample distances (VST)",
    fontsize                 = 11
  )
  invisible(NULL)
}
