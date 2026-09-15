# Reusable paired bulk-expression methods extracted from the study analysis.
# Supply count matrices, verified sample metadata and gene sets.
# No sample identities, exclusions, gene sets or study results are embedded.
# Dependencies: edgeR and fgsea. Installation and execution are caller-controlled.

# Generic robust quasi-likelihood core; design rows must match count columns.
run_edger <- function(count_matrix, metadata, design, coef_name) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR is required")
  if (ncol(count_matrix) != nrow(metadata) || nrow(design) != nrow(metadata)) {
    stop("Counts, metadata and design dimensions do not agree")
  }
  y <- edgeR::DGEList(counts = count_matrix)
  # Filter on biological condition groups, preserving the original paired method.
  # Design-based leverage in a subject-fixed model can admit additional sparse genes.
  keep <- edgeR::filterByExpr(y, group = metadata$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  coef_index <- match(coef_name, colnames(design))
  if (is.na(coef_index)) stop("Coefficient not found: ", coef_name)
  qlf <- edgeR::glmQLFTest(fit, coef = coef_index)
  tab <- as.data.frame(qlf$table)
  tab$gene_id <- rownames(tab)
  tab$FDR <- stats::p.adjust(tab$PValue, method = "BH")
  df_total <- qlf$df.total
  if (length(df_total) == 1L) df_total <- rep(df_total, nrow(tab))
  names(df_total) <- rownames(qlf$table)
  tab$df_total <- unname(df_total[tab$gene_id])
  tab$ql_se_approx <- ifelse(
    tab$F > 0, abs(tab$logFC) / sqrt(pmax(tab$F, 0)), NA_real_
  )
  tcrit <- stats::qt(0.975, pmax(tab$df_total, 1))
  tab$logFC_ci_low_approx <- tab$logFC - tcrit * tab$ql_se_approx
  tab$logFC_ci_high_approx <- tab$logFC + tcrit * tab$ql_se_approx
  tab$ci_method <- "approximate one-coefficient robust QL interval"
  tab <- tab[order(tab$PValue, -abs(tab$logFC)), ]
  list(
    y = y, fit = fit, qlf = qlf, table = tab,
    metadata = metadata, design = design
  )
}

# metadata has sample_id, subject_id and condition columns.
# condition_levels sets the reference first and comparison second.
# Each subject must contribute one sample in each condition.
# Cohort inclusion/exclusion decisions must be made by the caller before this call.
run_paired_bulk <- function(count_matrix, metadata, condition_levels) {
  required <- c("sample_id", "subject_id", "condition")
  if (!all(required %in% names(metadata))) {
    stop("metadata requires sample_id, subject_id and condition")
  }
  condition_levels <- as.character(condition_levels)
  if (length(condition_levels) != 2L || anyNA(condition_levels) ||
      anyDuplicated(condition_levels)) {
    stop("Provide two distinct condition levels, reference then comparison")
  }
  md <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (anyNA(md[, required]) || anyDuplicated(md$sample_id) ||
      is.null(colnames(count_matrix)) || anyDuplicated(colnames(count_matrix)) ||
      !all(md$sample_id %in% colnames(count_matrix))) {
    stop("Metadata must uniquely identify available samples without missing labels")
  }
  md$subject_id <- factor(md$subject_id)
  md$condition <- factor(as.character(md$condition), levels = condition_levels)
  if (anyNA(md$condition) || any(table(md$subject_id, md$condition) != 1L)) {
    stop("Each subject must have exactly one sample in each supplied condition")
  }
  if (nlevels(md$subject_id) < 2L) {
    stop("At least two complete subject pairs are required for residual degrees of freedom")
  }
  rownames(md) <- md$sample_id
  design <- stats::model.matrix(
    ~ subject_id + condition, data = md,
    contrasts.arg = list(condition = stats::contr.treatment(condition_levels, base = 1L))
  )
  if (qr(design)$rank != ncol(design)) stop("The paired design is not full rank")
  coef_name <- colnames(design)[ncol(design)]
  run_edger(count_matrix[, md$sample_id, drop = FALSE], md, design, coef_name)
}

# Preserve the full filtered ranking; do not preselect nominally significant genes.
rank_from_fit <- function(fit) {
  tab <- fit$table
  score <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
  names(score) <- tab$gene_id
  sort(score[is.finite(score)], decreasing = TRUE)
}

# pathways is one named list of gene vectors for one prespecified database family.
# Call separately for separate databases so adjusted values retain their families.
# Set the random seed explicitly before calling if a repeatable run is required.
run_ranked_gsea <- function(fit, pathways, min_size = 15L, max_size = 500L) {
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("fgsea is required")
  if (!is.list(pathways) || is.null(names(pathways)) ||
      anyNA(names(pathways)) || any(!nzchar(names(pathways))) ||
      anyDuplicated(names(pathways))) {
    stop("pathways must be a named list with unique nonempty pathway names")
  }
  ranks <- rank_from_fit(fit)
  if (is.null(names(ranks)) || anyNA(names(ranks)) ||
      any(!nzchar(names(ranks))) || anyDuplicated(names(ranks))) {
    stop("The full ranking requires unique, nonempty gene identifiers")
  }
  res <- as.data.frame(fgsea::fgseaMultilevel(
    pathways, ranks, minSize = min_size, maxSize = max_size, eps = 0, nproc = 1
  ))
  res$leadingEdge <- vapply(res$leadingEdge, paste, collapse = ";", character(1))
  res[order(res$padj, -abs(res$NES)), ]
}
