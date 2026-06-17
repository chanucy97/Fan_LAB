# hcc_DLC1

Public code skeleton for an HCC-focused DLC1 analysis project. The workflow collects bulk TCGA/GTEx analyses, TCGA-LIHC survival and enrichment analyses, single-cell fibroblast-state analyses, CellChat ligand-receptor screening, and spatial transcriptomics follow-up.

This directory contains code and documentation only. Raw matrices, RData/RDS objects, clinical tables, spatial objects, figures, tables, and generated outputs are not included.

## Included

- `scripts/DLC1_analysis_public.R`: public-safe version of the user-provided DLC1 R workflow.
- `requirements.R`: package inventory and optional installer helper.
- `configs/hcc_dlc1_root.example.R`: example runtime-root configuration.
- `docs/RUNNING.md`: execution notes.
- `docs/INPUT_OUTPUT.md`: controlled input and output boundary.
- `docs/DATA_CODE_AVAILABILITY.md`: release position for code and controlled data.
- `docs/SENSITIVE_CHECKLIST.md`: pre-commit safety checks.

## Runtime Root

Set `HCC_DLC1_ROOT` to a private directory that contains local inputs and receives generated outputs:

```r
Sys.setenv(HCC_DLC1_ROOT = "D:/controlled/hcc_DLC1")
```

If `HCC_DLC1_ROOT` is unset, the script uses the current working directory.

## Quick Start

```r
source("requirements.R")
install_hcc_dlc1_cran_packages()

Sys.setenv(HCC_DLC1_ROOT = "D:/controlled/hcc_DLC1")
source("scripts/DLC1_analysis_public.R")
```

The script is intentionally modular but still reflects the original working analysis notebook/script style. Review and run the relevant section for the data available in your private runtime root.

## Public-Safety Boundary

The original WeChat-delivered script was not committed verbatim. The public copy removes private absolute paths and routes local data through `HCC_DLC1_ROOT`. Do not commit raw data, RDS/RData objects, clinical tables, spatial transcriptomics objects, generated tables, figures, or manuscript outputs.
