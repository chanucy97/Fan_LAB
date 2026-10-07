# Generic four-dimension survival modelling source

This package contains source code and an example configuration only. It includes no patient rows, patient IDs, outcome values, fitted model objects, site paths, or study-specific cohort counts. The package is designed for four information dimensions: clinical, conventional radiomics, habitat radiomics, and 2.5D deep features. Each imaging dimension has two source tables. Four fusion approaches are evaluated in parallel.

## Required inputs

Supply a private manifest.csv with the columns ID, group, cv_fold, rfs_event, and rfs_time_months. Choose any cohort labels; enter the development label in config.json as development_group. Development rows must have consecutive outer-fold values from zero to outer_folds-1. Validation rows use cv_fold=-1. The code checks duplicate IDs, missing fields, invalid endpoint values, and adequate events in each fold. It derives cohort sizes and event totals from the input; no patient ID or count is hard-coded.

Supply clinical.csv and six imaging feature tables. Each has one ID column and one row per manifest patient. Extra imaging feature rows are documented and ignored; missing manifest patients cause a failure. clinical_schema.json is user supplied, defining every clinical variable as continuous or categorical; categorical variables must list allowed numeric codes. The schema file in this package is only a toy example. Set pnm_models in the configuration only if conventional categorical benchmark models are needed. All such variable names must appear in the clinical schema.

Copy config.example.json to a private config.json, edit paths, group label and optional settings, and run:

    python run_rigorous_analysis.py --config config.json --input-root /path/to/private/inputs --outdir /path/to/new/output --preflight-only
    python run_rigorous_analysis.py --config config.json --input-root /path/to/private/inputs --outdir /path/to/new/output

Input specifications may include an optional sha256 field to lock private files. The public example omits hashes so it can be adapted to a new dataset. If included, a checksum mismatch stops the run. Keep real data, private configuration, generated manifest_used.csv, clinical_matrix_used.csv, predictions.csv, and model_artifacts/ out of public repositories.

The time horizon is 36 months and the endpoint column names reflect recurrence-free survival. This workflow assumes enough observed events and follow-up to estimate its nested Cox models and 36-month metrics. It is not a universal survival-analysis pipeline for every dataset shape or outcome definition. The software does not validate clinical endpoint adjudication or upstream image extraction provenance. Any study using it should document its own input definitions, feature-extraction and ICC procedures, model settings, and reporting decisions.
