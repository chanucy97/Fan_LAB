# Input And Output Boundary

## Public Repository Inputs

The public repository includes only:

- Python scripts.
- Placeholder markdown template.
- Documentation and configuration examples.
- Dependency list.

## Controlled Runtime Inputs

Keep these outside the public checkout:

- Real manuscript markdown or journal-specific text.
- Whole-slide images and patch outputs.
- Feature tensors or embeddings.
- Patient-level and slide-level clinical or provenance tables.
- Source data tables for figures.
- Reviewed Word/PDF/ZIP submission packages.
- Rendered QA pages and large figures.

## Generated Outputs

The builder can create:

- DOCX manuscript draft assembled from controlled markdown.
- Example supplementary table templates with column headers only.
- A build manifest describing input paths and release boundaries.

The QA renderer can create:

- PDF converted from a controlled DOCX.
- Page PNGs.
- A contact sheet.
- A render manifest.

These outputs should remain under `AE_METASTASIS_ROOT` or another private working directory. They are not intended for public Git commits.

## Release Boundary

Use this rule for future updates:

```text
Commit code and templates.
Do not commit real data, unpublished manuscript text, personal contact details,
ethics/funding text, final submission packages, rendered QA pages, or large assets.
```
