options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

cran_pkgs <- c("optparse", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "jsonlite")
to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install)
} else {
  cat("All CRAN packages already installed.\n")
}

# BioConductor: DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_pkgs <- c("DESeq2")
bioc_missing <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(bioc_missing) > 0) {
  cat("Installing Bioconductor packages:", paste(bioc_missing, collapse = ", "), "\n")
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
} else {
  cat("All Bioconductor packages already installed.\n")
}

# Optional: ashr for lfcShrink(type="ashr")
if (!"ashr" %in% rownames(installed.packages())) {
  cat("Installing ashr...\n")
  install.packages("ashr")
}

cat("\n=== Verification ===\n")
for (pkg in c("optparse", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
              "jsonlite", "DESeq2", "ashr")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-18s %s\n", pkg, if (ok) "[OK]" else "[MISSING]"))
}
