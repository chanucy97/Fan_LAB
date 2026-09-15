# Analysis methods

## Scope and inputs

This guide describes analysis procedures and their statistical units. It contains no study findings. Function extracts expose reusable parts of the original workflows; orchestration, final annotations, selected result panels and dataset-specific decisions are deliberately outside this release. Methods described without an accompanying helper remain descriptions, not newly implemented pipelines.

Required inputs are supplied by the user running the code: expression count matrices, sample metadata, independently verified pairing, approved annotations, receptor sequences/clone counts, spatial masks and calibrated coordinates, protein intensities and gene-set resources. The repository supplies no sample records or populated data templates.

## Single-cell RNA sequencing

Starting from filtered expression matrices, calculate detected-feature counts and mitochondrial fractions; apply specified QC cutoffs before normalization. Use library-size normalization and log transformation, variable-feature selection, scaling and PCA. Where specified, use Harmony for integration, followed by neighbour graph construction, clustering and UMAP. `single_cell_preprocessing.R` exposes the underlying operations and their caller-configurable parameters.

Interpret cluster identity using marker expression, lineage consistency and externally supplied reviewed annotations. Do not infer biological identity from UMAP geometry alone. Final cluster-to-cell-type maps and marker-confidence decisions are not supplied in this release. Optional doublet filtering, pseudobulk analysis, lineage trajectories, transcription-factor activity and cell-communication workflows require their own documented inputs, references and assumptions; the preprocessing helper does not implement those entire downstream workflows.

## Sample-level composition and comparisons

Aggregate cells into a sample-by-category count matrix using caller-provided labels. Add a declared pseudocount, compute proportions and use the centred log-ratio transformation for compositional comparisons. Group heterogeneity can be expressed as mean Euclidean distance to the current group centroid in CLR space. Permute sample labels while preserving group sizes to obtain an exact two-sided test of the group difference. Recompute centroids for every assignment.

`statistics.R` provides the reusable transformation and permutation calculations. Permutation assumes exchangeability under the null; the caller must account for pairing, batch structure and any non-independent samples. Do not permute individual cells as though they were biological replicates. Define the multiple-testing family in advance and preserve adjusted and unadjusted P values as separate fields.

## Paired bulk differential expression and enrichment

Use raw gene counts, match columns exactly to metadata, and filter low-expression genes under the declared analysis design. Apply TMM normalization and a paired `subject + condition` design in edgeR. Fit robust quasi-likelihood models and test the condition coefficient. Keep log fold change, raw P value and BH-adjusted values distinct.

For rank-based enrichment, the source workflow uses signed square-root quasi-likelihood F statistics as the ranking statistic and fgsea multilevel analysis with caller-provided pathway sets. `paired_bulk.R` exposes these steps. Cohort exclusions, actual pair identities and expected findings are not built into the functions. Database versions and identifier mapping are inputs that must be recorded by the caller.

## Immune receptor repertoire

Define clones consistently within each library from approved receptor calls, retain the unit used for clone counting, and summarize library-level clone sizes before comparing groups. Diversity metrics operate on positive clone counts rather than on patient-identifying sequences. `repertoire_metrics.R` extracts the source workflow's generic clone/diversity calculations; annotations, library identities and retained receptor records are not included.

Receptor-to-expression linkage, sequence-quality filtering, source-identity checks and external antigen-database matching are separate prerequisites. Shared sequence patterns or database matches do not by themselves establish antigen specificity.

## Spatial transcriptomics

Keep each independently acquired spatial specimen and coordinate system separate. Register pathology masks to expression coordinates and calibrate distances using the platform's physical bin scale. Exclude invalid or non-reportable areas using explicit masks. Reference-program projection measures relative support for supplied programs; it is not automatically a cell-proportion estimate.

`spatial_methods.py` contains two method extracts:

1. **Program scoring:** the caller selects a gene matrix; raw counts are normalized as log1p(count/library size × 10,000), scaled per gene by a positive-reference quantile, clipped to the unit interval and averaged. Detected-gene support is tracked separately.
2. **Directional extent:** threshold scores within the eligible mask, require gene support, identify eight-neighbour components in centroid-referenced sectors, retain components meeting size and boundary-seed criteria, and measure their farthest boundary distance. Nearest-region assignment separates competing exterior territories. Distance is bin-centre distance minus half a bin.

Connected-positive area and the geometric envelope area are different quantities. Field-edge censoring is indicated when retained signal touches the invalid-mask edge. The helper requires an invalid frame around the grid to make that boundary explicit. Directions describe the supplied image coordinates and must not be called anatomical directions without validated orientation. No region shapes, gene lists, measured extents or region-specific results are distributed.

## Proteomics

Preserve the distinction between vendor identification/quantification and downstream analysis. For protein-intensity PCA, retain complete positive rows, log2-transform intensities, centre each protein across samples and apply PCA without additional variance scaling. For heatmaps, compute row-wise z scores in log2 space. Keep missing values missing in display data; the helper uses a separate zero-filled matrix only for clustering distances.

`proteomics_processing.R` exposes these transformations. It does not supply differential proteins, observed intensities, host/parasite results, significance thresholds selected from outcomes, or vendor mass-spectrometry preprocessing software. Species-specific identifier mapping and pathway backgrounds remain caller-controlled inputs.

## Immune deconvolution and experimental-group analysis

Bulk immune deconvolution requires compatible gene identifiers, correctly normalized expression, versioned references and method-specific assumptions. Interpret abundance-like fractions separately from enrichment scores; agreement among methods is a sensitivity check, not independent experimental validation. The methods are described here without redistributing external reference matrices or licensed components.

For animal or other experimental-group analyses, keep the independent experimental unit in the metadata, distinguish descriptive cellular patterns from replicate-level inference, and analyse measured phenotypes using their original units and missingness rules. Human and animal layers are not interchangeable causal evidence. No group assignments, observed phenotypes, intervention outcomes or experimental-result figures are included.

## Validation boundary

The release is checked for syntax, self-contained function definitions, metadata-free content and basic algorithm behaviour on artificial fixtures. These fixtures are software checks, not simulated study results, and are not deposited as article data. The original datasets were not rerun and study findings were not reproduced for this methods-only publication.
