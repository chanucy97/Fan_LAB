# Package inventory for the public hcc_DLC1 workflow.
# Install only what is needed for the analysis module you are running.

hcc_dlc1_cran_packages <- c(
  "data.table",
  "dplyr",
  "e1071",
  "forestplot",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "grid",
  "Matrix",
  "msigdbr",
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "readr",
  "remotes",
  "scales",
  "stringr",
  "survival",
  "survminer",
  "tibble",
  "tidyr"
)

hcc_dlc1_bioc_packages <- c(
  "AnnotationDbi",
  "Biobase",
  "BiocParallel",
  "biomaRt",
  "clusterProfiler",
  "ComplexHeatmap",
  "DOSE",
  "edgeR",
  "enrichplot",
  "fgsea",
  "GOSemSim",
  "GSVA",
  "limma",
  "monocle",
  "org.Hs.eg.db",
  "preprocessCore",
  "sva"
)

hcc_dlc1_github_packages <- c(
  "enblacar/SCpubr",
  "junjunlab/ClusterGVis",
  "omnideconv/immunedeconv@v2.0.3"
)

hcc_dlc1_optional_packages <- c(
  "CellChat",
  "EPIC",
  "Seurat",
  "xCell"
)

install_hcc_dlc1_cran_packages <- function(pkgs = hcc_dlc1_cran_packages) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(missing)
}

install_hcc_dlc1_bioc_packages <- function(pkgs = hcc_dlc1_bioc_packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
  invisible(missing)
}

install_hcc_dlc1_github_packages <- function(pkgs = hcc_dlc1_github_packages) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  for (pkg in pkgs) {
    remotes::install_github(pkg, dependencies = TRUE, upgrade = "never")
  }
  invisible(pkgs)
}
