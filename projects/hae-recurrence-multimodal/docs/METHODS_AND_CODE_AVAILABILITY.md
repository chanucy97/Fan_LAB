# Methods and code availability / 方法与代码可用性

The text below describes release **1.0.0**, dated **7 October 2026**, of the tabular survival-modelling and evaluation core. Use the actual analysis configuration and software version when describing a study run.

## Materials and Methods — English

**Concise version**

Python source code for the tabular survival-modelling and evaluation core, together with example configurations, is available in the Fan_LAB repository (version 1.0.0; 7 October 2026): https://github.com/chanucy97/Fan_LAB/tree/main/projects/hae-recurrence-multimodal.

**Two-sentence version**

The released Python implementation covers training-partition preprocessing, nested out-of-fold fusion of clinical and pre-extracted imaging features, Harrell C-index evaluation, and 36-month IPCW AUC, Brier score, grouped calibration and decision-curve analysis. Source code and example configurations are available in the Fan_LAB repository (version 1.0.0; 7 October 2026): https://github.com/chanucy97/Fan_LAB/tree/main/projects/hae-recurrence-multimodal.

## 材料与方法 — 中文

**简版**

基于表格特征的生存建模与评价核心 Python 源代码及示例配置已公开于 Fan_LAB 仓库（版本 1.0.0，2026 年 10 月 7 日）：https://github.com/chanucy97/Fan_LAB/tree/main/projects/hae-recurrence-multimodal。

**两句版**

公开的 Python 实现包括训练分区内预处理、临床变量与预提取影像特征的嵌套折外融合建模、Harrell C-index 评价，以及 36 个月 IPCW AUC、Brier 评分、分组校准和决策曲线分析。源代码和示例配置已公开于 Fan_LAB 仓库（版本 1.0.0，2026 年 10 月 7 日）：https://github.com/chanucy97/Fan_LAB/tree/main/projects/hae-recurrence-multimodal。

## Code availability — English

The tabular survival-modelling and evaluation source code and example configurations are publicly available in the Fan_LAB repository (version 1.0.0; released 7 October 2026) at https://github.com/chanucy97/Fan_LAB/tree/main/projects/hae-recurrence-multimodal. File-level SHA-256 hashes are provided in the source manifest. The release covers the modelling and evaluation core; patient-level data, fitted model artifacts, upstream image-processing and feature-extraction code, and additional article-specific analyses are not included in this public release.

## 代码可用性 — 中文

基于表格特征的生存建模与评价源代码及示例配置已公开于 Fan_LAB 仓库（版本 1.0.0，发布于 2026 年 10 月 7 日）：https://github.com/chanucy97/Fan_LAB/tree/main/projects/hae-recurrence-multimodal。源文件清单提供各文件的 SHA-256 校验值。本次发布覆盖建模与评价核心；患者级数据、拟合模型、上游图像处理与特征提取代码，以及其他文章特定分析未纳入本次公开版本。

## Author usage notes / 作者使用说明

- Select one Methods version and adapt it to the actual analysis. A public source release describes the available implementation; the study's Methods should separately state the configuration and workflow actually used.
- The repository title is a descriptive software-resource name. Use the manuscript's own title when identifying the article.
- Data availability requires a separate author statement reflecting actual arrangements. The text above makes no claim about data-access approval, request procedures, or permission to share patient data.
- A later change to the source or release scope requires a corresponding version and statement update.

[Project page](https://chanucy97.github.io/Fan_LAB/hae-recurrence-multimodal.html) · [Source manifest](../source-manifest.json) · [Implementation notes](IMPLEMENTATION_NOTES.md)

For a version-specific reference, use the [version 1.0.0 source snapshot](https://github.com/chanucy97/Fan_LAB/tree/hae-recurrence-v1.0.0/projects/hae-recurrence-multimodal) and record the associated Git commit.
