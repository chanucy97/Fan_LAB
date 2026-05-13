# Fan_LAB

这是樊海宁课题组代码与项目交接文档的总入口仓库。它的目标不是一次性塞进所有数据和结果，而是先把代码资产、项目边界、数据规则和协作流程管理起来。

## 推荐定位

- 管理实验室可复用代码、项目说明、运行环境、服务器命令和交接文档。
- 不直接保存患者原始数据、测序原始数据、WSI 图像、模型权重、压缩包、投稿版敏感材料。
- 对已有独立项目仓库，优先保留原仓库历史；这里作为总索引和规范入口。

## 当前已识别项目

| 项目 | 本机路径 | 建议管理方式 |
| --- | --- | --- |
| AE 病理组学/WSI-临床融合 | `C:\Users\Administrator\Desktop\binglizuxue`、`C:\Users\Administrator\Documents\New project 2` | 独立私有仓库或子模块；严格排除临床表、WSI、预测大文件 |
| Research2 母婴肠道菌群/SCFA | `C:\Users\Administrator\Desktop\he`、`C:\Users\Administrator\Documents\New project 3` | 独立私有仓库；保存分析脚本和投稿文档模板，排除原始数据 |
| 省自然基金/项目申请书 | `C:\Users\Administrator\Documents\New project 4` | 可放文档仓库或本仓库 `projects/grants/`，按版本管理 |
| 空白新项目 | `C:\Users\Administrator\Documents\New project 5` | 暂不纳入 |

详细盘点见 [docs/project_inventory.md](docs/project_inventory.md)。

## 建议目录

```text
Fan_LAB/
  README.md
  docs/
    project_inventory.md
    data_policy.md
    github_setup.md
  projects/
    README.md
  shared/
    README.md
    scripts/
  templates/
    project_README_template.md
```

## 下一步

1. 在 GitHub 创建一个 private 仓库，例如 `Fan_LAB`。
2. 先推送本仓库的索引、规范和模板。
3. 逐个项目清理 `.gitignore`，确认没有敏感数据后，再决定是迁入代码还是保留独立仓库。
4. 为每个项目补一个简短 `README.md`：研究目标、输入数据位置、运行命令、输出位置、负责人、敏感数据规则。
