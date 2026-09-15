# Reusable fixed-QC and embedding methods adapted from the study analysis code.
# Inputs and batch labels are supplied by the caller; functions return objects.
# Numeric QC/normalization/dimension defaults follow the source workflow.
# Adaptations: input validation, explicit Harmony choice, and optional caller seed.
# Required packages: Seurat; harmony when use_harmony is TRUE.

build_qc_seurat <- function(count_matrices, batch_labels,
                            min_cells = 3L, min_features = 200L,
                            max_features = 7500L, max_mito_percent = 25,
                            mito_pattern = "^MT-", verbose = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
  if (!is.list(count_matrices) || !length(count_matrices)) {
    stop("Provide a named list of gene-by-cell count matrices.")
  }
  keys <- names(count_matrices)
  if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
    stop("Count-matrix names must be nonempty and unique.")
  }
  if (is.null(names(batch_labels)) || anyDuplicated(names(batch_labels)) ||
      !setequal(names(batch_labels), keys)) {
    stop("Provide one named batch label for each count matrix.")
  }
  batch_labels <- as.character(batch_labels[keys])
  if (anyNA(batch_labels) || any(!nzchar(batch_labels))) {
    stop("Batch labels must be nonmissing and nonempty.")
  }
  limits <- c(min_cells, min_features, max_features, max_mito_percent)
  if (any(!is.finite(limits)) || min_cells < 1 || min_features < 0 ||
      max_features < min_features || max_mito_percent < 0 ||
      max_mito_percent > 100) {
    stop("Invalid QC parameters.")
  }
  if (length(mito_pattern) != 1L || is.na(mito_pattern)) {
    stop("Provide one mitochondrial-gene regular expression.")
  }
  objects <- lapply(seq_along(count_matrices), function(i) {
    counts <- count_matrices[[i]]
    if (length(dim(counts)) != 2L || !nrow(counts) || !ncol(counts) ||
        is.null(rownames(counts)) || is.null(colnames(counts)) ||
        anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
      stop("Each count matrix needs unique gene and cell names.")
    }
    object <- Seurat::CreateSeuratObject(
      counts = counts, project = keys[i],
      min.cells = min_cells, min.features = min_features
    )
    object$library_id <- keys[i]
    object$batch <- batch_labels[i]
    object <- Seurat::RenameCells(object, add.cell.id = keys[i])
    object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
      object, pattern = mito_pattern
    )
    keep <- object$nFeature_RNA >= min_features &
      object$nFeature_RNA <= max_features &
      object$percent.mt <= max_mito_percent
    if (!any(keep)) stop("A supplied count matrix has no cells passing QC.")
    if (verbose) message("Applying fixed feature and mitochondrial filters.")
    subset(object, cells = colnames(object)[keep])
  })
  if (length(objects) == 1L) return(objects[[1]])
  Reduce(function(x, y) merge(x, y), objects)
}

cluster_single_cells <- function(object, use_harmony = TRUE, assay = "RNA",
                                  batch_column = "batch",
                                  variable_features = 2500L, n_pcs = 30L,
                                  resolution = 0.4, scale_factor = 10000,
                                  regress_variables = "percent.mt",
                                  seed = NULL, verbose = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
  if (!inherits(object, "Seurat")) stop("Provide a Seurat object.")
  if (length(assay) != 1L || !assay %in% Seurat::Assays(object)) {
    stop("Requested assay is absent from the Seurat object.")
  }
  Seurat::DefaultAssay(object) <- assay
  if (length(n_pcs) != 1L || !is.finite(n_pcs) ||
      n_pcs < 1 || n_pcs != as.integer(n_pcs) || ncol(object) <= n_pcs) {
    stop("n_pcs must be a positive integer smaller than the number of cells.")
  }
  if (length(variable_features) != 1L || !is.finite(variable_features) ||
      variable_features <= n_pcs || variable_features != as.integer(variable_features) ||
      length(resolution) != 1L || !is.finite(resolution) || resolution <= 0 ||
      length(scale_factor) != 1L || !is.finite(scale_factor) || scale_factor <= 0) {
    stop("Invalid variable-feature, resolution, or normalization parameters.")
  }
  if (!is.null(seed) && (length(seed) != 1L || !is.finite(seed) ||
                        seed < 0 || seed > .Machine$integer.max ||
                        seed != as.integer(seed))) {
    stop("seed must be NULL or a nonnegative integer.")
  }
  metadata <- object[[]]
  if (!all(regress_variables %in% names(metadata))) {
    stop("Requested regression variables are absent from cell metadata.")
  }
  if (use_harmony) {
    if (!requireNamespace("harmony", quietly = TRUE)) stop("harmony is required.")
    if (!batch_column %in% names(metadata) || anyNA(metadata[[batch_column]]) ||
        length(unique(metadata[[batch_column]])) < 2L) {
      stop("Harmony requires a metadata column with at least two batch levels.")
    }
  }
  if (!is.null(seed)) set.seed(seed)
  object <- Seurat::NormalizeData(
    object, normalization.method = "LogNormalize",
    scale.factor = scale_factor, verbose = verbose
  )
  object <- Seurat::FindVariableFeatures(
    object, selection.method = "vst", nfeatures = variable_features,
    verbose = verbose
  )
  features <- Seurat::VariableFeatures(object)
  if (length(features) <= n_pcs) stop("Too few variable genes for the requested PCA.")
  object <- Seurat::ScaleData(
    object, features = features, vars.to.regress = regress_variables,
    verbose = verbose
  )
  pca_args <- list(object = object, features = features, npcs = n_pcs, verbose = verbose)
  if (!is.null(seed)) pca_args$seed.use <- seed
  object <- do.call(Seurat::RunPCA, pca_args)
  dims <- seq_len(n_pcs)
  reduction <- "pca"
  if (use_harmony) {
    object <- harmony::RunHarmony(
      object = object, group.by.vars = batch_column,
      reduction.use = "pca", dims.use = dims, verbose = verbose
    )
    reduction <- "harmony"
  }
  object <- Seurat::FindNeighbors(
    object, reduction = reduction, dims = dims, verbose = verbose
  )
  cluster_args <- list(object = object, resolution = resolution, verbose = verbose)
  if (!is.null(seed)) cluster_args$random.seed <- seed
  object <- do.call(Seurat::FindClusters, cluster_args)
  umap_args <- list(object = object, reduction = reduction, dims = dims, verbose = verbose)
  if (!is.null(seed)) umap_args$seed.use <- seed
  do.call(Seurat::RunUMAP, umap_args)
}
