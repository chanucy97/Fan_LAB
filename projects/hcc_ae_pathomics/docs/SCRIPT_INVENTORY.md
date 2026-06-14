# Script Inventory

## `scripts/build_hcc_ae_submission_package.py`

Public-safe replacement for the local final-submission builder. The local source script embedded the complete unpublished manuscript and submission metadata, so this repository version uses external markdown input and a placeholder template.

Key behavior:

- Reads `HCC_AE_ROOT` or uses the current working directory.
- Builds a DOCX from markdown headings, paragraphs, bullet lists, and simple pipe tables.
- Optionally embeds controlled local figures in the DOCX without copying them into Git.
- Optionally writes example supplementary table templates with headers only.
- Writes a manifest that records the release boundary.

## `scripts/render_review_docx_windows.py`

Public-safe Windows QA renderer for a controlled DOCX.

Key behavior:

- Reads `HCC_AE_ROOT` for controlled output locations.
- Allows `LIBREOFFICE_SOFFICE` to override the LibreOffice path.
- Converts DOCX to PDF, renders pages, and creates a contact sheet.
- Does not hard-code private project paths.
- Does not auto-install `pypdfium2` unless `--install-pdfium` is passed.
