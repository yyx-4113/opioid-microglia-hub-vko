# GitHub 复现包发布操作手册（中文）

本手册说明如何把本项目的复现包发布到 GitHub，使稿件 §6 Data Availability 链接在投稿时可解析。

## 一、前置条件
- 已安装 Git 且 `git` 在 PATH。
- 已登录 GitHub CLI：`gh auth login`（账号 `yyx-4113`，scope 含 `repo`）。
- 确认本机 `02_scripts/`、`03_results/`、`04_figures/`、`06_docking/` 为最新权威产物（以 `manuscript_plosone_v6.docx` 引用为准）。

## 二、本地仓库准备（仅首次）
1. 在本地复现包目录初始化：`git init -b main`
2. 编写标准配置：`README.md`、`CITATION.cff`、`LICENSE`(MIT)、`.github/workflows/release.yml`、`GITHUB_DEPOSIT_SOP.md`、`author_verification_statement.md`。
3. **不要提交** `01_data/`（原始公共输入，体积大、超 100 MB，且已在稿件中给出 GEO 链接）。
4. 加入产物目录并提交：`git add 00_pipeline 02_scripts 03_results 04_figures 06_docking 07_submission README.md CITATION.cff LICENSE .github GITHUB_DEPOSIT_SOP.md author_verification_statement.md && git commit -m "v6.0.0 reproduction package"`

## 三、创建远程仓库并推送
- 创建公开仓库：`gh repo create yyx-4113/opioid-microglia-hub-vko --public --description "..." --source . --remote origin --push`
- 或手动：`gh repo create yyx-4113/opioid-microglia-hub-vko --public`，再 `git remote add origin https://github.com/yyx-4113/opioid-microglia-hub-vko.git && git push -u origin main`

## 四、打标签（与稿件 §6 一致）
- `git tag v6.0.0 && git push origin v6.0.0`
- 标签触发 `.github/workflows/release.yml`，自动打包 `results-bundle.zip` 并生成 GitHub Release。

## 五、核对清单（投稿前）
- [ ] 仓库为 **Public**。
- [ ] 存在 tag **v6.0.0**，且 `02_scripts/`、`03_results/`、`06_docking/`(含 `step8c_decoy_summary.json`)、`04_figures/` 均在仓库根目录可见。
- [ ] 稿件 §6 链接 `https://github.com/yyx-4113/opioid-microglia-hub-vko`（tag v6.0.0）可打开。
- [ ] 接受后再补 Zenodo 归档快照（带 MANIFEST 校验和）。

## 六、版本号约定
- 每次稿件大版本（v4/v5/v6…）对应一个 GitHub tag（v4.0.0/v5.0.0/v6.0.0…）。
- tag 必须与稿件 §6 中写明的一致，否则 Data Availability 链接无法解析。
