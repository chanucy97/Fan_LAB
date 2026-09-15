# AE single-cell and spatial analysis methods

This methods-only release shares methodological descriptions and reusable function extracts for single-cell, spatial and related multi-omics analysis. The functions were adapted from the study's working analysis code, with input data, group labels and annotations supplied by the caller.

**Release scope:** analysis methods. No study result tables, observed cell counts, final cluster annotations, confidence decisions, patient/sample mappings, source images, manuscript findings or expected-result assertions are included. This is not the complete original analysis/figure pipeline, and it should not be cited as an end-to-end reproduction package.

## Contents

| File | Method |
| --- | --- |
| [methods/single_cell_preprocessing.R](methods/single_cell_preprocessing.R) | Seurat preprocessing and optional Harmony integration |
| [methods/statistics.R](methods/statistics.R) | Compositional transformation and exact sample-label permutation |
| [methods/paired_bulk.R](methods/paired_bulk.R) | Paired edgeR quasi-likelihood analysis and ranked pathway enrichment |
| [methods/repertoire_metrics.R](methods/repertoire_metrics.R) | Repertoire diversity and clone summaries |
| [methods/spatial_methods.py](methods/spatial_methods.py) | Gene-program scoring and directional boundary-distance measurement |
| [methods/proteomics_processing.R](methods/proteomics_processing.R) | Protein-intensity PCA and heatmap transformation |
| [docs/METHODS.md](docs/METHODS.md) | Methodological overview, input requirements and interpretation rules |
| [docs/CODE_AVAILABILITY.md](docs/CODE_AVAILABILITY.md) | Manuscript wording for this limited deposit |
| [environment/README.md](environment/README.md) | Dependencies and validation scope |

R files define functions only. Source the chosen file and pass your own input objects; no analysis runs when they are sourced. The Python module likewise performs no file I/O or automatic analysis. Methodological numeric defaults describe algorithms, not observed study results.

The [Fan Lab website](https://chanucy97.github.io/Fan_LAB/ae-single-cell-spatial.html) provides a public entry. Use an exact repository commit when citing a particular version. The release date is 2026-09-15; no DOI or project-specific open-source licence is assigned here.

## 中文说明

此目录仅公开分析方法、通用函数及必要的输入说明。具体细胞数量、最终注释、研究结果、样本配对记录和结果校验值均未收入。函数需要运行者自行提供数据和参数；本目录不代表全文全部分析及图表的完整复现代码。
