# 数据与隐私边界

本仓库默认按医学/生物信息项目的敏感数据标准处理。原则是：代码可托管，数据不默认托管。

## 可以进入 GitHub 的内容

- 分析脚本：`.py`、`.R`、`.sh`、`.ps1`。
- 配置模板：`.yaml`、`.json`、`.toml`，前提是不含服务器密码、token、患者标识。
- 文档：`README.md`、运行说明、项目交接文档、投稿清单。
- 小型示例数据：必须脱敏，且仅用于演示代码格式。
- 环境文件：`requirements.txt`、`environment.yml`、`renv.lock`。

## 默认不进入 GitHub 的内容

- 患者姓名、住院号、电话、身份证号、病理号、样本条码映射表。
- 原始测序数据、WSI 全切片图像、临床原始表。
- 模型权重、训练日志、服务器大体积输出。
- `.env`、API key、服务器账号、SSH 密钥。
- 投稿前敏感材料、未公开审稿意见、合作方未授权数据。
- zip/rar/7z 压缩包、docx/pptx/xlsx/pdf 等大体积交付包，除非明确需要版本管理。

## 处理大文件

如果确实需要托管大文件，优先顺序：

1. 本地/服务器/网盘保存数据，GitHub 只保存数据路径和校验说明。
2. 使用 Git LFS 管理少量必须版本化的大文件。
3. 对复杂数据流水线使用 DVC 或类似工具，GitHub 保存元数据。

## 提交前检查

每次上传前至少检查：

```powershell
git status --short
git diff --cached --name-only
git diff --cached
```

如涉及表格或配置，再额外搜索：

```powershell
rg -n "password|token|secret|身份证|住院号|姓名|电话|病理号|patient" .
```

