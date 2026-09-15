# Reusable composition methods extracted from the study analysis.
# Rows are independent biological samples; columns are caller-defined classes.
# Supply counts and verified group labels. No study data or group names are stored.
# These helpers implement the primary centroid-distance test, not every sensitivity.

make_count_matrix <- function(data, sample_col, label_col, sample_levels, label_levels) {
  valid_column <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  if (!is.data.frame(data) || !valid_column(sample_col) ||
      !valid_column(label_col) || anyDuplicated(names(data)) ||
      !all(c(sample_col, label_col) %in% names(data))) {
    stop("Provide a data frame with unambiguous sample and label columns")
  }
  sample_levels <- as.character(sample_levels)
  label_levels <- as.character(label_levels)
  valid_levels <- function(x) {
    length(x) > 0L && !anyNA(x) && !anyDuplicated(x) && all(nzchar(x))
  }
  if (!valid_levels(sample_levels) || !valid_levels(label_levels)) {
    stop("Supplied sample and label levels must be unique and nonmissing")
  }
  sample_values <- as.character(data[[sample_col]])
  label_values <- as.character(data[[label_col]])
  if (anyNA(sample_values) || anyNA(label_values) ||
      !all(sample_values %in% sample_levels) ||
      !all(label_values %in% label_levels)) {
    stop("Input sample and label values must be nonmissing and within supplied levels")
  }
  sample_factor <- factor(sample_values, levels = sample_levels)
  label_factor <- factor(label_values, levels = label_levels)
  tab <- table(sample_factor, label_factor)
  storage.mode(tab) <- "numeric"
  tab
}

# Add the same pseudocount to each class, normalize within sample, then center logs.
clr_from_counts <- function(counts, pseudocount = 0.5) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("counts must be a finite nonnegative numeric matrix")
  }
  if (nrow(counts) < 1L || ncol(counts) < 2L ||
      length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount <= 0) {
    stop("Provide samples, at least two classes and a positive pseudocount")
  }
  adjusted <- sweep(
    counts + pseudocount, 1,
    rowSums(counts) + pseudocount * ncol(counts), "/"
  )
  log_adjusted <- log(adjusted)
  log_adjusted - rowMeans(log_adjusted)
}

# Euclidean distances on CLR coordinates are Aitchison distances.
# Group centroids are recomputed for every supplied label assignment.
distance_to_centroid <- function(x, labels, bias_adjust = FALSE) {
  labels <- factor(labels)
  result <- numeric(nrow(x))
  for (g in levels(labels)) {
    idx <- which(labels == g)
    center <- colMeans(x[idx, , drop = FALSE])
    d <- sqrt(rowSums(
      (x[idx, , drop = FALSE] - rep(center, each = length(idx)))^2
    ))
    if (bias_adjust && length(idx) > 1L) {
      d <- d * sqrt(length(idx) / (length(idx) - 1))
    }
    result[idx] <- d
  }
  result
}

primary_values <- function(x, labels, reference, comparison) {
  labels <- as.character(labels)
  d <- distance_to_centroid(x, labels, bias_adjust = FALSE)
  reference_value <- mean(d[labels == reference])
  comparison_value <- mean(d[labels == comparison])
  c(
    reference_value = reference_value,
    comparison_value = comparison_value,
    difference = comparison_value - reference_value,
    ratio = comparison_value / reference_value
  )
}

# Enumerate all assignments while holding the observed group sizes fixed.
# This is intended for small cohorts; the number of assignments grows rapidly.
all_label_assignments <- function(n, n_comparison, reference, comparison) {
  comparison_sets <- utils::combn(seq_len(n), n_comparison, simplify = FALSE)
  lapply(comparison_sets, function(idx) {
    z <- rep(reference, n)
    z[idx] <- comparison
    z
  })
}

# x must contain sample-level CLR coordinates, for example from clr_from_counts().
# The statistic is the difference in mean distance to each group's own centroid.
# The full enumeration includes the observed assignment; use its exact tail fraction.
# Multiplicity correction belongs to a caller-defined family of tests.
exact_centroid_test <- function(x, labels, reference, comparison) {
  x <- as.matrix(x)
  labels <- as.character(labels)
  reference <- as.character(reference)
  comparison <- as.character(comparison)
  if (length(reference) != 1L || length(comparison) != 1L ||
      is.na(reference) || is.na(comparison) || identical(reference, comparison)) {
    stop("Specify two distinct group labels")
  }
  if (!is.numeric(x) || ncol(x) < 2L || any(!is.finite(x)) ||
      length(labels) != nrow(x) || anyNA(labels)) {
    stop("Provide finite sample coordinates and one group label per row")
  }
  if (!setequal(unique(labels), c(reference, comparison)) ||
      any(table(labels) < 2L)) {
    stop("Both specified groups must have at least two independent samples")
  }
  observed <- primary_values(x, labels, reference, comparison)
  assignments <- all_label_assignments(
    nrow(x), sum(labels == comparison), reference, comparison
  )
  null <- vapply(
    assignments,
    function(z) primary_values(x, z, reference, comparison),
    numeric(length(observed))
  )
  rownames(null) <- names(observed)
  data.frame(
    reference_value = unname(observed["reference_value"]),
    comparison_value = unname(observed["comparison_value"]),
    difference = unname(observed["difference"]),
    ratio = unname(observed["ratio"]),
    exact_permutation_p = mean(
      abs(null["difference", ]) >= abs(observed["difference"]) - 1e-12
    ),
    permutations = length(assignments)
  )
}
