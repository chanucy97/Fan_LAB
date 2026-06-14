# Running The AE Metastasis Pathomics Skeleton

## 1. Create An Environment

```powershell
cd projects\ae_metastasis_pathomics
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
```

## 2. Set The Controlled Runtime Root

Keep real manuscript and data files outside this public checkout:

```powershell
$env:AE_METASTASIS_ROOT = "D:\controlled\ae_metastasis_pathomics"
```

Expected private layout:

```text
%AE_METASTASIS_ROOT%\
  controlled_manuscript\
    ae_metastasis_submission_final.md
  controlled_figures\
    figure1.png
    figure2.png
  ae_metastasis_submission_package\
```

The exact private layout can differ; pass paths with CLI flags when needed.

## 3. Build A Public Template Package

This uses the tracked placeholder manuscript and writes no controlled data:

```powershell
.\.venv\Scripts\python scripts\build_ae_metastasis_submission_package.py --write-example-tables
```

## 4. Build From Controlled Manuscript Text

```powershell
.\.venv\Scripts\python scripts\build_ae_metastasis_submission_package.py `
  --manuscript "$env:AE_METASTASIS_ROOT\controlled_manuscript\ae_metastasis_submission_final.md" `
  --out "$env:AE_METASTASIS_ROOT\ae_metastasis_submission_package" `
  --figures-dir "$env:AE_METASTASIS_ROOT\controlled_figures" `
  --embed-figures `
  --write-example-tables
```

Outputs remain under `AE_METASTASIS_ROOT` and are ignored by Git.

## 5. Render DOCX QA Pages On Windows

LibreOffice is required. The script uses `LIBREOFFICE_SOFFICE` when set, otherwise it tries the default Windows install path.

```powershell
$env:LIBREOFFICE_SOFFICE = "C:\Program Files\LibreOffice\program\soffice.com"

.\.venv\Scripts\python scripts\render_review_docx_windows.py `
  --docx "$env:AE_METASTASIS_ROOT\ae_metastasis_submission_package\AE_metastasis_submission.docx" `
  --out "$env:AE_METASTASIS_ROOT\ae_metastasis_submission_package\qa_pdf_pages"
```

If `pypdfium2` is not already installed in the environment, install dependencies from `requirements.txt` or pass `--install-pdfium` to create a local controlled cache under `AE_METASTASIS_ROOT`.
