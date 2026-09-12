# gene-correlation

> **科学子问题**: 高变异基因之间的共表达关系如何？是否存在功能模块？

## 职责

从表达矩阵中选取方差最高的 top N 基因，计算基因间相关性矩阵（Pearson 或 Spearman），绘制相关性热图并输出相关系数矩阵。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| expr_matrix | file | ✓ | — | 表达矩阵 (CSV: gene × sample) |
| top_n | integer | — | `30` | 按方差选取的高变异基因数 |
| method | text | — | `pearson` | 相关性方法 (pearson / spearman) |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `heatmap.png` | image | 基因相关性热图 |
| `heatmap.pdf` | image | 热图 PDF 版本 |
| `corr_matrix.csv` | table | 相关系数矩阵 |
| `manifest.json` | manifest | 工具 manifest |

## R 包依赖

- pheatmap (>=1.0.12)
- jsonlite (>=1.8.0)

## 平台

- ✅ macos-arm64 / macos-x64 / windows-x64

## 注意

此插件使用 legacy `manifest.json` 格式（非 `artifact-manifest.json`）。
升级到完整 Artifact lineage 格式已跟踪为 HRT-1 任务。
