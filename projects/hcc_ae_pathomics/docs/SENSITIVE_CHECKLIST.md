# Sensitive-Content Checklist

Run this checklist before every public commit that touches `projects/hcc_ae_pathomics/`.

## Files That Must Not Be Added

- WSI files: `*.svs`, `*.ndpi`, `*.tif`, `*.tiff`.
- Patch, feature, or embedding files: `*.h5`, `*.hdf5`, `*.pt`, `*.pth`, `*.npy`, `*.npz`.
- Source data tables: `*.csv`, `*.tsv`, `*.xlsx`, `*.xls`.
- Manuscript packages: `*.docx`, `*.doc`, `*.pdf`, `*.zip`, `*.tar`, `*.tar.gz`.
- QA pages or large figures: `*.png`, `*.jpg`, `*.jpeg`, `*.webp`.

## Text That Must Be Reviewed

- Local absolute paths such as `C:\Users\...`.
- Email addresses and personal contact details.
- Ethics approval numbers or committee names that have not been approved for release.
- Grant numbers, funder text, and author-contribution text.
- Real slide identifiers, patient identifiers, or center-specific source tables.
- Complete unpublished manuscript prose.
- API tokens, passwords, private keys, or service credentials.

## Suggested Local Scan

```powershell
rg -n -i "C:\\Users|Administrator|New project 5|[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}|ethic|irb|grant|fund|patient|TCGA-|\\.svs|\\.docx|\\.pdf|\\.zip|password|token|secret|api[_-]?key" projects\hcc_ae_pathomics
```

Some checklist and documentation hits are expected. Investigate any hit that contains real identifiers, private paths, unpublished text, or credentials.
