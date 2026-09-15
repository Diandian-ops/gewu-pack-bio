# volcano-plot

> **科学子问题**: 差异分析结果中哪些基因在统计显著性和效应量上最值得关注？

## 职责

从差异分析结果（log2FoldChange + pvalue/padj）绘制 publication-ready 火山图，支持阈值标注与 top 基因高亮。

## DAG 位置

```
deg-standardizer → volcano-plot
```

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| deg_results | file | ✓ | — | 差异分析结果表 (CSV/TSV) |
| gene_col | text | — | `gene` | 基因 ID 列名 |
| lfc_col | text | — | `log2FoldChange` | log2FC 列名 |
| pval_col | text | — | `padj` | 用于判定和纵轴的显著性列；显式改为 `pvalue` 属科学变化 |
| pval_threshold | number | — | `1.3` | `-log10(所选显著性列)` 阈值 |
| lfc_threshold | number | — | `1` | |log2FC| 阈值 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `volcano_plot.png` | image | 火山图 |
| `report.html` | report | HTML 报告 |
| `manifest.json` | manifest | 工具 manifest |

## R 包依赖

- ggplot2
- ggrepel

## 平台

- ✅ macos-arm64 / macos-x64 / windows-x64

## 注意

此插件使用 legacy `manifest.json` 格式（非 `artifact-manifest.json`）。
升级到完整 Artifact lineage 格式已跟踪为 HRT-1 任务。
