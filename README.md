# 格物 · 生信分析包（gewu-pack-bio）

格物桌面（Gewu Desktop）的**生信垂直领域包**。格物核心是通用 agent 工作台——生信分析作为可安装的领域包存在，装则生信能力齐备，卸则回到纯通用态，核心零残留。

## 为什么是独立仓库

BLOCK-83 裁决（2026-09-12）：领域能力不进核心，核心只提供领域中性的扩展点。本包是**租户 #1 / 参考实现**——它证明平台够用，也示范后续领域包（化学、影像、金融等）的接入方式。

## 内容

| 目录 | 内容 | 迁移自（来源仓） |
| --- | --- | --- |
| `plugins/` | 15 个分析插件（seurat-standard / deseq2 / scanpy-advanced / cellchat / scvi-tools / gsea / go-kegg / volcano-plot / wgcna / pca-explorer / celltypist / expression-qc / deg-standardizer / gene-correlation / sample-metadata-validator）+ 插件清单 | `resources/built-in-plugins/` |
| `scripts/` | 69 个 R/Python 分析脚本 | `core/biof3-server/tool-scripts/` |
| `demo-data/` | 演示数据集 | `core/biof3-server/tool-demo-data/` |
| `pack.json` | 包契约：声明提供物（插件/脚本/语言/确认钩子/产物类型）与所需核心能力 | 本包新建 |

## 安装

1. 把本包的 `plugins/` 目录登记为插件的 installedDir：
   - 环境变量：`BIOF3_PLUGINS_DIR=<本包路径>/plugins`
   - 或应用设置页的插件目录
2. 重启应用；分析入口与插件商店即出现生信能力。

## 平台契约（依赖的核心扩展点）

| 核心扩展点 | 用途 |
| --- | --- |
| `plugin-source:installed-dir` | 核心从 installedDir 发现插件（不内置领域内容） |
| `script-host-languages` | 语言以数据描述（扩展名/标签/命令名），宿主按描述符执行 |
| `task-confirmation-hooks` | 领域确认门（本包注册 `scientific-task-confirmation`） |
| `tool-catalog-projection` | 工具目录由已装插件投影（`GET /api/tools/catalog`） |
| `script-host:executor` | 领域中性的脚本执行宿主 |

## 迁移来源

本包内容自 `ShengXinF3/biof3-desktop` 迁出（BLOCK-83 P4），保留 `git log` 中的历史可追溯性。核心仓不再包含生信内容，`npm run` 门禁 `test/domain-neutrality.test.mjs` 保证这一点不被回退。
