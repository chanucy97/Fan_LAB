# 肝泡型包虫病术后复发多模态生存建模

**Multimodal survival modelling for postoperative recurrence in hepatic alveolar echinococcosis**

Version **1.0.0** · Release date **2026-10-07**

[Project page](https://chanucy97.github.io/Fan_LAB/hae-recurrence-multimodal.html) · [Repository directory](https://github.com/chanucy97/Fan_LAB/tree/main/projects/hae-recurrence-multimodal) · [Source manifest and SHA-256 hashes](source-manifest.json)

[Version 1.0.0 source snapshot](https://github.com/chanucy97/Fan_LAB/tree/hae-recurrence-v1.0.0/projects/hae-recurrence-multimodal) provides the named reference for this deposit.

## Purpose and scope / 用途与范围

This resource provides the Python core for survival modelling and evaluation using four tabular information dimensions: clinical variables, conventional radiomics, habitat radiomics, and 2.5D deep features. It implements training-partition preprocessing, nested out-of-fold stacking, four fusion approaches, Harrell C-index evaluation, and probability-based assessments at 36 months. The title above is the descriptive name of this software resource.

本项目提供基于已提取表格特征的生存建模与评价代码。输入包括临床变量、常规影像组学、生境影像组学和 2.5D 深度特征，输出包括模型预测、审计记录及评价汇总。公开范围为通用建模与评价核心、示例配置和依赖说明；患者级数据、拟合模型、上游图像处理与特征提取流程，以及其他文章特定分析不在本次发布范围内。

Version 1.0.0 preserves the ten supplied files in `code/` byte for byte. [source-manifest.json](source-manifest.json) records their SHA-256 hashes; [VERSION](VERSION) records the release version. Implementation behavior and usage limits are documented in [Implementation notes](docs/IMPLEMENTATION_NOTES.md). English and Chinese manuscript wording is provided in [Methods and code availability](docs/METHODS_AND_CODE_AVAILABILITY.md).

## Code map / 代码职责

| File | Responsibility |
| --- | --- |
| [`code/run_rigorous_analysis.py`](code/run_rigorous_analysis.py) | Command-line entry point; manifest and configuration checks, nested model fitting, output and provenance recording. |
| [`code/analysis_core.py`](code/analysis_core.py) | Training-partition preprocessing, feature screening, sparse Cox models, and inner out-of-fold base predictions. |
| [`code/attention_core.py`](code/attention_core.py) | Four-score attention Cox model, internal epoch selection, and Breslow conversion to 36-month event probability. |
| [`code/evaluate_survival.py`](code/evaluate_survival.py) | Harrell C-index, IPCW AUC and Brier score, descriptive calibration, decision curves, and paired bootstrap comparisons. |
| [`code/clinical_schema.py`](code/clinical_schema.py) | Explicit clinical-variable types, numeric category validation, and training-derived handling of missing or unseen categories. |

The remaining supplied files are the original `README.md`, two example JSON files, `requirements.txt`, and `requirements-tested-windows.txt`. The latter is the submitted environment record; it does not constitute a full-training test of this public release.

## Inputs / 输入

Prepare the following CSV files in a private input directory. Configure their filenames under `inputs` in your private configuration.

| Configuration key | Content |
| --- | --- |
| `manifest` | Patient identifiers, cohort labels, outer-fold assignments, event indicators, and follow-up times. |
| `clinical` | Clinical predictors matching the user-supplied schema exactly. |
| `ct_p_global` | Conventional radiomics features from the CT source. |
| `t2wi_global` | Conventional radiomics features from the T2-weighted MRI source. |
| `ct_p_habitat` | Habitat radiomics features from the CT source. |
| `t2wi_habitat` | Habitat radiomics features from the T2-weighted MRI source. |
| `ct_p_2p5d` | Pre-extracted 2.5D features from the CT source. |
| `t2wi_2p5d` | Pre-extracted 2.5D features from the T2-weighted MRI source. |

The manifest must contain exactly `ID`, `group`, `cv_fold`, `rfs_event`, and `rfs_time_months`. Identifiers must be unique and complete. `rfs_event` is 0 or 1, and follow-up times must be finite positive values in months. The endpoint definition and clinical adjudication are supplied by the user.

- Set `development_group` to the cohort used for training; at least one separate validation cohort is required.
- Provide development `cv_fold` assignments as consecutive integers from `0` to `outer_folds - 1`. Validation patients use `-1`. The software uses these supplied outer folds and generates its inner splits from the training partition.
- Every predictor table requires one `ID` column and a row for every manifest patient. Predictor values are numeric; missing predictor values are handled within training partitions. The manifest itself permits no missing values.
- Define every clinical predictor as `continuous` or `categorical` in a private clinical schema. For categorical predictors, declare numeric `allowed_codes`. The supplied schema is a toy example that must be replaced with the intended variable definitions.
- Set `clinical_schema` to a path resolved relative to the configuration file. Each `inputs` entry uses `relative_path` under either `--input-root` (`location: "root"`) or the source-code directory (`location: "package"`). Use `root` for private patient inputs.
- Optional input `sha256` entries lock file contents. Optional `code_sha256` entries lock implementation files. The public examples contain adaptable paths and settings, not a study-specific analysis lock.
- The prediction and evaluation horizon is fixed at **36 months** in the source. Changing a JSON field alone does not change this horizon.

Nested fitting requires adequate events and censoring within its smaller training and validation subsets. Preflight checks validate the input interface; the full run also applies model-level checks that can stop on an unsuitable partition. See the [implementation notes](docs/IMPLEMENTATION_NOTES.md).

## Run / 运行

Use a separate Python environment and install the declared dependencies. From this project directory:

```sh
cd code
python -m pip install -r requirements.txt
```

Copy `config.example.json` and `clinical_schema.example.json` to a private location. Edit the configuration, schema, input filenames, cohort label, and fold settings to match your data. The following paths are placeholders; replace them with absolute paths for your operating system.

```sh
python run_rigorous_analysis.py --config "/path/to/private/config.json" --input-root "/path/to/private/inputs" --outdir "/path/to/private/new-output" --preflight-only
python run_rigorous_analysis.py --config "/path/to/private/config.json" --input-root "/path/to/private/inputs" --outdir "/path/to/private/new-output"
```

The output directory must be new or empty. A full run writes predictions, model artifacts, fitted-preprocessing audits, provenance, and evaluation CSV files. It also copies the supplied manifest, clinical matrix, and configuration into the output directory. Keep these outputs in a controlled private location because they can contain patient identifiers, outcomes, and fitted models.

The example configuration uses `pnm_models: {}`. Optional categorical benchmark models require explicit variable lists that match the clinical schema. Epochs, fold counts, feature-selection caps, and bootstrap replicates are configurable; the recorded attention optimizer settings must match the constants implemented in the source.

## Reporting / 报告

Describe the actual input definitions, fold assignments, configuration, software environment, and source version used for an analysis. The repository provides [EN/CN code-availability wording](docs/METHODS_AND_CODE_AVAILABILITY.md) scoped to this release. Authors should write their data-availability statement separately using their actual access arrangements.
