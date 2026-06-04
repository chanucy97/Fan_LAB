# GitHub pre-commit sensitive-information checklist

Run this checklist before committing or pushing.

## 1. Check staged files

```bash
git status --short
git diff --cached --name-only
```

Only these file types should normally be staged:

- `.py`
- `.md`
- `.txt`
- `.yaml`
- `.yml`
- `.csv` under `data/example/` only
- `.gitkeep`

## 2. Never stage these materials

- patient raw data
- clinical tables
- sample number mapping tables
- hospital ID or visit ID files
- raw flow-cytometry exports
- sequencing raw data
- WSI images
- model weights
- compressed archives
- manuscript submission packages
- Word/PDF manuscript drafts
- ethics approval documents
- ICMJE disclosure forms
- account names, passwords, tokens, API keys, server paths

## 3. Search for common sensitive patterns

From the Fan_LAB repository root:

```bash
git grep -n -I -E "password|passwd|token|secret|api[_-]?key|Authorization|Bearer" -- projects/tbnk_lupus_nephritis
git grep -n -I -E "patient_id|hospital_id|visit_id|sample_no|sample_id|medical_record|phone|name|identifier" -- projects/tbnk_lupus_nephritis
git grep -n -I -E "H:|C:\\\\Users|/mnt/|/home/|ssh|scp|sftp" -- projects/tbnk_lupus_nephritis
```

Expected behavior:

- `patient` may appear only in documentation warnings or availability statements.
- Local Windows paths should not appear in scripts or configs.
- No token/password/API key should appear anywhere.

## 4. Check for large or disallowed files

```bash
git ls-files --stage projects/tbnk_lupus_nephritis
```

Do not commit:

- `*.zip`
- `*.docx`
- `*.pdf`
- `*.tiff`
- `*.svs`
- `*.h5`
- `*.bam`
- `*.fastq`
- `*.pt`
- `*.pth`
- `*.pkl`

## 5. Final rule

If a file was generated during manuscript submission or contains real patient or clinical information, do not upload it to GitHub.

中文提醒：只上传代码、说明文档、环境文件、配置模板和脱敏示例数据；不要上传患者原始数据、临床表、编号映射表、投稿文件、伦理文件、账号密码、token 或本地服务器路径。
