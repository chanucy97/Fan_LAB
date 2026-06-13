# Sensitive-Information Checklist

Before committing future changes to this project, confirm that the commit does not include:

- Real participant-level clinical records or identifiable metadata.
- Raw milk assay exports, stool assay exports, FASTQ files, sylph intermediate outputs, or GTDB database files.
- Sample, dyad, hospital, or instrument identifier mapping tables that can re-identify participants.
- Manuscript circulation packages, reviewer drafts, Word/PDF submissions, or ZIP archives.
- Credentials, tokens, private keys, server paths, or environment files.
- Large generated artifacts such as TIFF/PDF figure exports, workbooks, model files, and logs.

Use `git status -sb`, `git diff --cached --stat`, and a sensitive-word scan before pushing.
