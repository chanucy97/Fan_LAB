# AE Metastasis Pathomics

Public code skeleton for a whole-slide-image pathomics study of non-local metastasis prediction in hepatic alveolar echinococcosis (AE).

This project is for audited AE metastasis modeling, not for a diagnostic differential-classification task. It keeps code and documentation public while the private manuscript text, author contact details, ethics and funding text, whole-slide images, feature files, source data tables, patient-level clinical data, final submission packages, QA screenshots, and large rendered figures remain controlled outside Git.

## What Is Included

- `scripts/build_ae_metastasis_submission_package.py`: public-safe manuscript/package builder that reads controlled local inputs at runtime.
- `scripts/render_review_docx_windows.py`: Windows DOCX-to-PDF QA renderer with configurable local paths.
- `docs/public_manuscript_template.md`: placeholder manuscript structure without unpublished text.
- `docs/RUNNING.md`: setup and execution notes.
- `docs/INPUT_OUTPUT.md`: input and output boundary.
- `docs/DATA_CODE_AVAILABILITY.md`: public release position for code and controlled data.
- `docs/SENSITIVE_CHECKLIST.md`: pre-release checklist for future updates.

## What Is Not Included

- Real WSI files or patch/feature tensors.
- Patient-level or slide-level clinical source tables.
- Source CSV/TSV/XLSX files used to make figures or supplementary tables.
- Unsubmitted manuscript text, author emails, ethics approvals, or grant details.
- Word/PDF submission packages, ZIP archives, rendered QA PNGs, and large figure assets.

## Runtime Root

Set `AE_METASTASIS_ROOT` to a controlled local directory that contains private inputs and receives generated outputs:

```powershell
$env:AE_METASTASIS_ROOT = "D:\controlled\ae_metastasis_pathomics"
```

If `AE_METASTASIS_ROOT` is not set, scripts use the current working directory. Generated outputs are ignored by Git.

## Quick Start

```powershell
cd projects\ae_metastasis_pathomics
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
.\.venv\Scripts\python scripts\build_ae_metastasis_submission_package.py --write-example-tables
```

For real manuscript assembly, pass a controlled manuscript markdown file that lives outside the public repository:

```powershell
.\.venv\Scripts\python scripts\build_ae_metastasis_submission_package.py `
  --manuscript "$env:AE_METASTASIS_ROOT\controlled_manuscript\ae_metastasis_submission_final.md" `
  --out "$env:AE_METASTASIS_ROOT\ae_metastasis_submission_package"
```

## Public-Safety Rule

Before committing future changes, run the sensitive-content checks in `docs/SENSITIVE_CHECKLIST.md`. The project is designed so code can be public while data, unpublished writing, and submission artifacts remain controlled.
