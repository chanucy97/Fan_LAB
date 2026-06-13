# Human milk bioactives, infant gut succession, and SCFA maturation

This project contains reproducible analysis and manuscript-delivery code for a Qinghai-Tibetan Plateau mother-infant cohort manuscript focused on human milk bioactive dynamics, infant gut Bifidobacterium succession, and fecal short-chain fatty acid (SCFA) maturation.

## Important data notice

**No real participant-level clinical data, sequencing files, raw assay exports, sample mapping tables, or manuscript circulation packages are included in this repository.** The public repository contains code and documentation only. The scripts expect a controlled local analysis workspace that is not uploaded to GitHub.

## What is included

- R scripts for Figure 1-7 redraws and lagged milk-to-SCFA sensitivity analysis.
- Node scripts for manuscript table and supplementary table workbook construction.
- Python scripts for manuscript figure PDF assembly and Word/PDF delivery builds.
- Running instructions, input/output descriptions, and data/code availability wording.
- A sensitive-information checklist for future commits.

## What is not included

This repository does not include raw metagenomic FASTQ files, sylph database files, stool or milk assay raw exports, clinical covariate tables, dyad/sample identifier mapping files, supplementary table workbooks, manuscript drafts, Word/PDF circulation files, ZIP packages, server paths, credentials, or any patient-level records.

## Repository layout

```text
milk_gut_scfa_plateau/
  README.md
  requirements.txt
  configs/
    research2_root.example.ps1
  docs/
    DATA_CODE_AVAILABILITY.md
    INPUT_OUTPUT.md
    RUNNING.md
    SENSITIVE_CHECKLIST.md
    SCRIPT_INVENTORY.md
  scripts/
    analysis/
    figures/
    manuscript_delivery/
    tables/
```

## Controlled workspace root

Most scripts use the environment variable `FANLAB_RESEARCH2_ROOT`. If it is not set, they use the current working directory. That root must contain the controlled source-data folders used during analysis, including `final_manuscript_planning_20260605/` and the upstream derived table directories described in `docs/INPUT_OUTPUT.md`.

PowerShell example:

```powershell
$env:FANLAB_RESEARCH2_ROOT = "D:\controlled\Research2"
Rscript projects\milk_gut_scfa_plateau\scripts\figures\fig5_milk_scfa_lagged_candidates\redraw_fig5_milk_scfa_lagged_candidates.R
```

## Code availability principle

The public code is maintained in the Fan_LAB GitHub repository under:

`projects/milk_gut_scfa_plateau`

Raw and participant-level data remain controlled because of privacy, ethics, and collaborating-institution data governance requirements.
