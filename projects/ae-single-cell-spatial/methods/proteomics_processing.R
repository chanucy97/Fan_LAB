# Method extracts from the protein-intensity PCA/heatmap workflow.
# Input is a numeric protein-by-sample matrix supplied by the caller.
# No specimen labels, protein selection, observed values or saved results included.

protein_pca <- function(intensity) {
  x <- as.matrix(intensity)
  storage.mode(x) <- "double"
  if (ncol(x) < 2L) stop("At least two sample columns are required")
  keep <- complete.cases(x) & apply(x, 1, function(z) all(is.finite(z) & z > 0))
  if (sum(keep) < 2L) stop("At least two complete positive proteins are required")
  z <- log2(x[keep, , drop = FALSE])
  z <- sweep(z, 1, rowMeans(z), FUN = "-")
  pc <- prcomp(t(z), center = TRUE, scale. = FALSE)
  list(pca = pc, retained_rows = which(keep))
}

protein_heatmap_values <- function(intensity) {
  m <- as.matrix(intensity)
  storage.mode(m) <- "double"
  if (nrow(m) < 2L || ncol(m) < 2L) stop("At least two proteins and samples are required")
  m[!is.finite(m) | m <= 0] <- NA_real_
  m <- log2(m)
  row_mean <- rowMeans(m, na.rm = TRUE)
  row_sd <- apply(m, 1, sd, na.rm = TRUE)
  z <- sweep(sweep(m, 1, row_mean, FUN = "-"), 1, row_sd, FUN = "/")
  z[!is.finite(z)] <- NA_real_
  # Zero replacement is for clustering distance only; displayed missing values stay NA.
  cluster_matrix <- z
  cluster_matrix[is.na(cluster_matrix)] <- 0
  dendrogram <- as.dendrogram(hclust(dist(cluster_matrix), method = "complete"))
  list(z = z, row_dendrogram = dendrogram)
}
