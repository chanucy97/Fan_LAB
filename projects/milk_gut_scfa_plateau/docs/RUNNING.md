# Running The Code

## Environment

The scripts were written for a controlled local analysis workspace. Set `FANLAB_RESEARCH2_ROOT` to that workspace root before running any script.

```powershell
$env:FANLAB_RESEARCH2_ROOT = "D:\controlled\Research2"
```

The controlled root is expected to contain `final_manuscript_planning_20260605/` plus upstream derived table folders such as `direction1_publication_figures/`, `direction3_multiomics_blueprint/`, `p0_multiomics_extension/`, `p1_clinical_modifiers/`, `p2_bridge_analyses/`, and related analysis outputs.

## R scripts

Install the R packages used by the figure and sensitivity scripts:

```r
install.packages(c(
  "dplyr", "forcats", "ggplot2", "gridExtra", "patchwork", "png",
  "ragg", "readr", "scales", "stringr", "svglite", "tidyr"
))
```

Run from the repository root or from the controlled workspace root:

```powershell
Rscript projects\milk_gut_scfa_plateau\scripts\figures\fig1_cohort_design_overview\redraw_fig1_coverage_heatmap.R
Rscript projects\milk_gut_scfa_plateau\scripts\figures\fig5_milk_scfa_lagged_candidates\redraw_fig5_milk_scfa_lagged_candidates.R
Rscript projects\milk_gut_scfa_plateau\scripts\analysis\gut_microbes_lagged_sensitivity.R
```

## Python delivery scripts

Install Python dependencies:

```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r projects\milk_gut_scfa_plateau\requirements.txt
```

Then run, for example:

```powershell
python projects\milk_gut_scfa_plateau\scripts\manuscript_delivery\build_all_figures_pdf.py
python projects\milk_gut_scfa_plateau\scripts\manuscript_delivery\build_gut_microbes_manuscript_pdf.py
```

## Node workbook scripts

The workbook-generation scripts use the Codex workspace spreadsheet artifact runtime (`@oai/artifact-tool`). If that runtime is unavailable, use the scripts as provenance for the workbook build logic or port the sheet-writing calls to another Excel library.

```powershell
node projects\milk_gut_scfa_plateau\scripts\tables\build_manuscript_tables.mjs
node projects\milk_gut_scfa_plateau\scripts\tables\build_gut_microbes_supplementary_tables.mjs
```

## Public-data limitation

This repository intentionally omits real source data. External users can inspect the analysis logic and required input schema, but cannot reproduce private results without controlled data access.
