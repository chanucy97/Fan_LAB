# Running hcc_DLC1

## 1. Prepare A Private Runtime Root

Keep controlled input files outside this public checkout:

```r
Sys.setenv(HCC_DLC1_ROOT = "D:/controlled/hcc_DLC1")
```

Suggested private layout:

```text
HCC_DLC1_ROOT/
  bulk_survival_optimal_cutoff/
  bulk_lihc/
    lihc.gdc_2022.rda
  single_cell/
    GSE149614/
    DLC1_sc_149614_restart_by_original_style/
  spatial_transcriptomics/
    GSE238264/
```

The script accepts this layout through `hcc_dlc1_root()` and `file.path(...)` calls. Adjust local filenames as needed before running a section.

## 2. Install Packages

```r
source("requirements.R")
install_hcc_dlc1_cran_packages()
install_hcc_dlc1_bioc_packages()
```

Install optional or GitHub packages only for modules that require them:

```r
install_hcc_dlc1_github_packages()
```

## 3. Run Relevant Modules

The public script is a consolidated analysis script. It contains:

- Pan-cancer DLC1 expression, survival, TMB/MSI, checkpoint, and optional immune deconvolution.
- TCGA-LIHC optimal-cutoff survival analysis.
- TCGA-LIHC bulk DEG and enrichment analysis.
- Single-cell fibroblast-state analysis.
- Monocle/ClusterGVis follow-up.
- CellChat ligand-receptor screening.
- Spatial transcriptomics HE/native coordinate and C1/C3 neighbor analyses.

Run only the section whose controlled input files are available:

```r
Sys.setenv(HCC_DLC1_ROOT = "D:/controlled/hcc_DLC1")
source("scripts/DLC1_analysis_public.R")
```

For safer interactive work, open the script and execute the needed section rather than sourcing all modules at once.
