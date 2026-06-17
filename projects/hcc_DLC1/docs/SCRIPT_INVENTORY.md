# Script Inventory

## `scripts/DLC1_analysis_public.R`

Public-safe version of the user-provided DLC1 R workflow. The script preserves the analysis logic but replaces private absolute paths with `hcc_dlc1_root()` and `HCC_DLC1_ROOT`.

Major modules:

- Pan-cancer DLC1 expression and TCGA/GTEx comparison.
- TCGA survival, TMB, MSI, checkpoint, and immune deconvolution analyses.
- TCGA-LIHC optimal-cutoff survival analysis.
- TCGA-LIHC bulk DEG, GO/KEGG, and GSEA analyses.
- GSE149614 single-cell analysis and fibroblast-state annotation.
- Monocle and ClusterGVis trajectory/enrichment follow-up.
- CellChat C2/C3/C4-to-C1 ligand-receptor screening.
- Spatial transcriptomics DLC1/C1/C3 visualization and neighbor statistics.

The script is not a data package. It expects local controlled inputs under `HCC_DLC1_ROOT`.
