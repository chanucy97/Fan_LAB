# 实验室代码资产初步盘点

盘点时间：2026-05-13

本文件来自本机只读扫描，未移动、未复制、未上传任何既有项目文件。

## 总览

| 路径 | 文件数 | 代码/文档/配置数 | 估计体积 | Git 状态 |
| --- | ---: | ---: | ---: | --- |
| `C:\Users\Administrator\Documents\New project 2` | 233 | 112 | 60.69 MB | 已初始化，本地未提交 |
| `C:\Users\Administrator\Documents\New project 3` | 186 | 34 | 78.69 MB | 已初始化，本地未提交 |
| `C:\Users\Administrator\Documents\New project 4` | 5 | 3 | 0.14 MB | 已初始化，本地未提交 |
| `C:\Users\Administrator\Documents\New project 5` | 0 | 0 | 0 MB | 已初始化，空项目 |
| `C:\Users\Administrator\Desktop\binglizuxue` | 153 | 55 | 2.53 MB | 未见 `.git` |
| `C:\Users\Administrator\Desktop\he` | 858 | 13 | 103.82 MB | 未见 `.git` |

## 识别到的主要资产

### AE 病理组学/WSI-临床融合

相关路径：

- `C:\Users\Administrator\Desktop\binglizuxue`
- `C:\Users\Administrator\Documents\New project 2`

主要内容：

- Python、Shell、YAML 配置和 Markdown 文档。
- AE clean patient-level split/config、MIL 训练配置、HPC 提交脚本。
- 病理组学、临床关联、WSI-临床融合相关输出和报告。

建议：

- 建议作为独立 private 仓库：`ae-pathomics-fusion`。
- GitHub 中只保留代码、配置模板、说明文档、小型示例表。
- 不上传真实临床原始表、患者标识、WSI 大文件、服务器训练输出、模型权重。

### Research2 母婴肠道菌群/SCFA

相关路径：

- `C:\Users\Administrator\Desktop\he`
- `C:\Users\Administrator\Documents\New project 3`

主要内容：

- R 脚本、Markdown 稿件、Word/PPT/补充表交付包。
- 宏基因组、SCFA、投稿材料和图表更新输出。

建议：

- 建议作为独立 private 仓库：`research2-infant-microbiome-scfa`。
- GitHub 中保存分析脚本、图表生成脚本、稿件 Markdown、投稿清单。
- 原始测序数据、临床信息、未脱敏表格和大体积结果留在本地/服务器/网盘。

### 中文基金/申请书文档

相关路径：

- `C:\Users\Administrator\Documents\New project 4`

主要内容：

- 中文项目申请书 Markdown 草稿。
- 省自然基金/科技厅项目设计相关文档。

建议：

- 如果只是个人写作版本管理，可以先放在本仓库 `projects/grants/`。
- 如果后续多人协作和敏感内容较多，可单独建 private 文档仓库。

## 当前风险点

- 多个项目已经初始化 Git，但还没有第一次提交；需要先补 `.gitignore` 再提交。
- 现有目录里存在 zip、docx、pptx、xlsx、pdf、输出目录和可能包含隐私的数据表，不能直接全量上传。
- `Desktop\he` 和 `Desktop\binglizuxue` 是关键数据/结果目录，不能简单 `git add .`。

## 推荐纳入顺序

1. 先提交当前 `lab-code-hub` 索引仓库。
2. 为 `New project 2`、`New project 3`、`New project 4` 分别加 `.gitignore` 并做本地首次提交。
3. 将 `Desktop\binglizuxue` 中的代码和配置复制到一个干净仓库目录，不直接把整个桌面目录变成仓库。
4. 将 `Desktop\he` 中的分析脚本和文档产物整理到干净仓库目录，数据目录只保留索引说明。

