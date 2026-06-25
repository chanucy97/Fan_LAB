# AE non-local metastasis pathomics code

This repository folder contains a de-identified reproduction code bundle for the manuscript on audited prediction of non-local metastasis in alveolar echinococcosis (AE) using H&E whole-slide images, clinical variables, and WSI-clinical fusion models.

The bundle is designed for manuscript review and method transparency. It does not include raw whole-slide images, patient-level clinical source tables, institutional file paths, model checkpoints, or protected health information.

## What is included

- `scripts/evaluate_mil_predictions.py`: summarize WSI/MIL prediction exports at slide and patient level.
- `scripts/run_clinical_extratrees.py`: reproduce clinical ExtraTrees-style tabular baselines under fixed patient folds.
- `scripts/run_fusion_extratrees.py`: combine clinical variables with WSI probability exports under the same patient folds.
- `scripts/make_roc_figure.R`: build publication-style ROC panels from prediction CSV files.
- `configs/cpathsoftware_clam_sb_conch_example.yaml`: path-free example of the CPathSoftWare CLAM-SB + CONCH configuration used as a template.
- `examples/input_schema/*.csv`: column schemas for required inputs.

## Data access

Raw WSIs and patient-level clinical data are controlled by the Affiliated Hospital of Qinghai University and are not publicly distributed because they may contain identifiable medical information. To run the scripts, users must provide de-identified tables following the schemas in `examples/input_schema/`.

## Minimal workflow

```bash
python scripts/evaluate_mil_predictions.py \
  --predictions-glob "results/mil/*/run_predictions.csv" \
  --output-dir outputs/wsi_metrics

python scripts/run_clinical_extratrees.py \
  --clinical-table data/clinical_features.csv \
  --fold-table data/folds.csv \
  --id-col patient_id \
  --label-col label \
  --numeric-cols Age,TB,DB,AST,ALT,Lesions \
  --categorical-cols Gender,Nationalities,Address,Therapy,HBsAg,Vascular_Invasion_binary \
  --output-dir outputs/clinical

python scripts/run_fusion_extratrees.py \
  --clinical-table data/clinical_features.csv \
  --wsi-predictions outputs/wsi_metrics/patient_predictions.csv \
  --fold-table data/folds.csv \
  --id-col patient_id \
  --label-col label \
  --numeric-cols Age,TB,DB,AST,ALT,Lesions \
  --categorical-cols Gender,Nationalities,Address,Therapy,HBsAg,Vascular_Invasion_binary \
  --output-dir outputs/fusion

Rscript scripts/make_roc_figure.R \
  --predictions outputs/fusion/patient_predictions.csv \
  --out outputs/figures/fusion_roc
```

## Reproducibility notes

The original WSI feature extraction and MIL training were performed with local scripts derived from the CPathSoftWare pipeline. This public bundle exposes the tabular, prediction-export, evaluation, and fusion layers and provides a sanitized configuration template for the WSI training layer. Exact raw WSI features and local paths are not included.

## No clinical-use statement

This code is for research reproducibility only. It is not validated for clinical diagnosis, treatment selection, or patient-level decision making.
