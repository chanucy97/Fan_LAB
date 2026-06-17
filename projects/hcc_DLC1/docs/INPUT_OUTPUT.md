# Input And Output Boundary

## Public Repository Inputs

This repository includes:

- Public-safe R analysis code.
- Package inventory.
- Runtime-root configuration example.
- Documentation and safety checklist.

## Controlled Inputs

Do not commit:

- TCGA, GTEx, GEO, or institutional expression matrices after download or preprocessing.
- Clinical/survival tables, mutation tables, MSI/TMB tables, sample annotations, or barcode-level merged tables.
- Seurat, CellChat, Monocle, ClusterGVis, spatial transcriptomics, or other serialized R objects.
- Source CSV/TSV/XLSX tables used for figures.
- Generated figures, tables, reports, or manuscript outputs.

## Generated Outputs

The script can create:

- CSV/TSV analysis tables.
- PDF/PNG figures.
- RDS/RData intermediate objects.
- Session information and module README files.

These outputs should remain under `HCC_DLC1_ROOT` or another private analysis directory. They are ignored by the project `.gitignore`.

## Release Rule

Commit code and templates. Do not commit data, identifiers, downloaded matrices, clinical tables, spatial objects, generated figures, generated tables, RDS/RData objects, or manuscript outputs.
