# GitHub 托管落地方案

## 最简单方案

1. GitHub 上创建一个 private 仓库：`Fan_LAB`。
2. 当前本地目录作为总入口仓库：

   `C:\Users\Administrator\Documents\New project 6`

3. 提交本仓库的 README、规范和模板。
4. 以后每个课题先在本仓库登记，再决定是否单独建仓库。

## 推荐组织方式

如果实验室后续会有多人参与，建议创建 GitHub Organization：

- `Fan_LAB`：樊海宁课题组总索引、规范、模板。
- `ae-pathomics-fusion`：AE 病理组学/WSI-临床融合代码。
- `research2-infant-microbiome-scfa`：Research2 母婴肠道菌群/SCFA 代码。
- `grant-writing`：项目申请书与文档草稿，可选。

## 我可以托管的部分

可以交给 Codex 处理：

- 本机项目盘点。
- `.gitignore`、`README.md`、项目模板。
- 检查是否有明显敏感文件准备进入提交。
- 初始化仓库、整理目录、创建提交。
- 在本机 `gh` 已登录时推送到 GitHub。
- 后续维护分支、标签、发布说明和项目交接文档。

需要你确认或授权：

- GitHub 仓库名和组织/个人账号归属。
- 仓库必须 private 还是可以 public。
- 哪些成员可以访问。
- 是否允许推送到远程。
- 是否允许把某些文档、表格或压缩包纳入版本管理。

## 本机推送命令模板

首次提交：

```powershell
git add README.md docs projects shared templates .gitignore
git commit -m "Initialize lab code hub"
```

如果使用 GitHub CLI 创建私有仓库：

```powershell
gh repo create Fan_LAB --private --source . --remote origin --push
```

如果 GitHub 上已经创建远程仓库：

```powershell
git remote add origin https://github.com/<OWNER>/Fan_LAB.git
git branch -M main
git push -u origin main
```
