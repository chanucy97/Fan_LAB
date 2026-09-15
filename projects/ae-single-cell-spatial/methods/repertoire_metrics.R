# Reusable clonotype summaries adapted from the study repertoire analysis.
# Inputs must already contain caller-defined clonotypes and one row per cell.
# Repertoire keys must separate libraries/assays; clones are counted within keys.
# Diversity and rarefaction formulas retain the source definitions.
# Adaptations: base-R grouping, explicit inputs, validation, and no data export.
# Pielou-derived clonality is NA for one observed clonotype because log(1) = 0.

.validate_clone_counts <- function(counts) {
  counts <- as.numeric(counts)
  if (!length(counts) || any(!is.finite(counts)) ||
      any(counts <= 0) || any(counts != floor(counts))) {
    stop("Provide positive integer counts for observed clonotypes.")
  }
  counts
}

repertoire_shannon <- function(counts) {
  counts <- .validate_clone_counts(counts)
  p <- counts / sum(counts)
  -sum(p * log(p))
}

repertoire_inverse_simpson <- function(counts) {
  counts <- .validate_clone_counts(counts)
  p <- counts / sum(counts)
  1 / sum(p^2)
}

repertoire_gini <- function(counts) {
  counts <- sort(.validate_clone_counts(counts))
  n <- length(counts)
  sum((2 * seq_len(n) - n - 1) * counts) / (n * sum(counts))
}

summarize_clonotypes <- function(cells, repertoire_column = "repertoire",
                                 clone_column = "clonotype",
                                 cell_column = "cell") {
  fields <- c(repertoire_column, clone_column, cell_column)
  if (!is.data.frame(cells) || !all(fields %in% names(cells)) ||
      anyDuplicated(fields)) {
    stop("Provide distinct repertoire, clonotype, and cell columns.")
  }
  key <- as.character(cells[[repertoire_column]])
  cell <- as.character(cells[[cell_column]])
  clone <- as.character(cells[[clone_column]])
  if (anyNA(key) || anyNA(cell) || any(!nzchar(key)) || any(!nzchar(cell)) ||
      anyDuplicated(data.frame(repertoire = key, cell = cell))) {
    stop("Each repertoire/cell pair must be present once and have nonempty keys.")
  }
  keep <- !is.na(clone) & nzchar(clone)
  if (!any(keep)) stop("No cells with assigned clonotypes were supplied.")
  tab <- stats::aggregate(
    rep.int(1L, sum(keep)),
    by = list(repertoire = key[keep], clonotype = clone[keep]), FUN = sum
  )
  names(tab)[3] <- "clone_cells"
  pieces <- split(tab, tab$repertoire, drop = TRUE)
  pieces <- lapply(pieces, function(x) {
    x$repertoire_cells <- sum(x$clone_cells)
    x$clone_frequency <- x$clone_cells / x$repertoire_cells
    x$clone_rank <- rank(-x$clone_cells, ties.method = "min")
    x
  })
  result <- do.call(rbind, unname(pieces))
  rownames(result) <- NULL
  result
}

repertoire_metrics <- function(counts, expansion_thresholds = c(2L, 4L),
                                hyperexpanded_above = 30L, top_n = 10L) {
  counts <- .validate_clone_counts(counts)
  thresholds <- as.numeric(expansion_thresholds)
  parameters <- c(thresholds, hyperexpanded_above, top_n)
  if (!length(thresholds) || any(!is.finite(parameters)) ||
      any(parameters < 1) || any(parameters != floor(parameters)) ||
      anyDuplicated(thresholds) || length(hyperexpanded_above) != 1L ||
      length(top_n) != 1L) {
    stop("Clone-size thresholds and top_n must be positive integers.")
  }
  n <- length(counts)
  total <- sum(counts)
  ordered <- sort(counts, decreasing = TRUE)
  entropy <- repertoire_shannon(counts)
  d50 <- which(cumsum(ordered) >= 0.5 * total)[1]
  result <- c(
    repertoire_cells = total,
    repertoire_richness = n,
    shannon = entropy,
    inverse_simpson_hill_q2 = repertoire_inverse_simpson(counts),
    gini = repertoire_gini(counts),
    clonality_1_minus_pielou = if (n > 1L) 1 - entropy / log(n) else NA_real_,
    d50_clone_count = d50,
    d50_clone_fraction = d50 / n,
    top1_fraction = max(counts) / total,
    top_n_fraction = sum(ordered[seq_len(min(top_n, n))]) / total,
    hyperexpanded_cell_fraction = sum(counts[counts > hyperexpanded_above]) / total,
    hyperexpanded_clone_fraction = mean(counts > hyperexpanded_above),
    sample_coverage_goods = 1 - sum(counts == 1L) / total
  )
  for (threshold in thresholds) {
    result[paste0("expanded_cell_fraction_n", threshold)] <-
      sum(counts[counts >= threshold]) / total
    result[paste0("expanded_clone_fraction_n", threshold)] <-
      mean(counts >= threshold)
  }
  result
}

expected_repertoire_richness <- function(counts, depth) {
  counts <- .validate_clone_counts(counts)
  total <- sum(counts)
  if (length(depth) != 1L || !is.finite(depth) || depth < 0 ||
      depth != floor(depth) || depth > total) {
    stop("Rarefaction depth must be an integer between zero and observed depth.")
  }
  if (depth == total) return(length(counts))
  # Hypergeometric expectation for sampling cells without replacement.
  sum(1 - exp(lchoose(total - counts, depth) - lchoose(total, depth)))
}

repertoire_rarefaction <- function(counts, depths) {
  counts <- .validate_clone_counts(counts)
  if (!length(depths)) stop("Provide at least one caller-selected depth.")
  data.frame(
    depth = depths,
    expected_richness = vapply(
      depths, function(depth) expected_repertoire_richness(counts, depth), numeric(1)
    ),
    observed_depth = sum(counts)
  )
}
