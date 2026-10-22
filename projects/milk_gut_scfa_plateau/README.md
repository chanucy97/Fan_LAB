# Qinghai-Tibetan Plateau mother-infant study: GlycanEvidence and analysis methods

This project contains reproducible analysis and manuscript-delivery code for a Qinghai-Tibetan Plateau mother-infant cohort manuscript focused on human milk bioactive dynamics, infant gut Bifidobacterium succession, and fecal short-chain fatty acid (SCFA) maturation.

## Current method source

[GlycanEvidence 0.1.0.dev10](glycan-evidence/) integrates separate marker/locus evidence with candidate-structure inference, joint fragment allocation, ordinary-assembly recovery and coding/base-level sequence checks. The public snapshot dated 22 October 2026 includes implementation, frozen rules, synthetic tests and detailed input/output contracts.

Start with the [methods code map](METHODS_CODE_MAP.md), [structure comparison](methods/structure-comparison/) and [sequence/population analysis](methods/sequence-and-population/) modules. These provide the reusable method components behind the current manuscript. Their provenance records distinguish original implementation from portable public adapters. The `scripts/` directory below retains the historical June 2026 cohort-analysis version and is not the source for every later figure or reanalysis.

## Important data notice

The repository contains method code, documentation, artificial test fixtures and selected aggregate regression records. Participant-level records, sequencing reads, assay exports and sample-identity tables follow the [study data-access route](docs/DATA_CODE_AVAILABILITY.md). Software users supply the inputs required by the documented interfaces.

## What is included

- GlycanEvidence dev10, method-comparison evaluators, sequence/population components and synthetic tests.
- Historical R scripts for Figure 1-7 redraws and lagged milk-to-SCFA sensitivity analysis.
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
