# Sensitive-Content Checklist

Run this checklist before every public commit that touches `projects/hcc_DLC1/`.

## Files That Must Not Be Added

- Serialized R objects: `*.rds`, `*.rda`, `*.RData`.
- Expression, clinical, mutation, or spatial tables: `*.csv`, `*.tsv`, `*.xlsx`, `*.xls`, `*.txt`.
- Large matrices and omics objects: `*.h5`, `*.h5ad`, `*.loom`, `*.mtx`.
- Generated figures and reports: `*.pdf`, `*.png`, `*.jpg`, `*.jpeg`, `*.tif`, `*.tiff`.
- Archives or compressed downloads: `*.zip`, `*.tar`, `*.tar.gz`, `*.gz`.

## Text That Must Be Reviewed

- Private absolute paths, including personal home directories and chat-cache paths.
- Patient identifiers, full barcodes beyond public-safe use, clinical row-level data, or sample-level metadata.
- Tokens, passwords, API keys, credentials, or private access URLs.
- Unreleased manuscript text or journal submission artifacts.

## Suggested Local Scan

```powershell
rg -n -i "/home/|C:\\Users|xwechat|wxid_|password|token|secret|api[_-]?key|BEGIN (RSA|OPENSSH|PRIVATE)|@[A-Z0-9._%+-]+\\.[A-Z]{2,}" projects\hcc_DLC1
git diff --cached --name-only | rg -n "\.(rds|rda|RData|h5|h5ad|loom|mtx|csv|tsv|xlsx|xls|txt|pdf|png|jpg|jpeg|tif|tiff|zip|tar|gz)$"
```

Documentation examples may intentionally mention generic patterns. Investigate any hit containing real private paths, identifiers, credentials, or data.
