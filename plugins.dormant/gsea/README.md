# gsea

> **科学子问题**: 差异表达基因的 ranked list 是否得到 GSEA 基因集富集证据支持？

## 职责

对差异分析产生的 ranked gene list（按统计显著性排序）执行 GSEA（Gene Set Enrichment Analysis），检测基因集在排序列表两端是否富集。与 GO/KEGG 超几何检验互补：GSEA 不需要预先定义"显著"阈值，使用全量 ranked list。

## DAG 位置

```
deg-standardizer → gsea
```

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| ranked_list | file | ✓ | — | 排序基因列表 (CSV: gene, rank 或 gene, metric) |
| organism | text | — | `hsa` | 物种 (hsa=human, mmu=mouse) |
| gene_id_type | text | — | `symbol` | ID 类型 (symbol/ensembl/entrez) |
| padj_cutoff | number | — | `0.05` | 富集 padj 阈值 |
| min_set_size | integer | — | `10` | 基因集最小大小 |
| max_set_size | integer | — | `500` | 基因集最大大小 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `gsea_results.csv` | table | GSEA 富集结果 |
| `gsea_plot.png` | image | GSEA 富集图 |
| `report.html` | report | HTML 报告 |
| `manifest.json` | manifest | 工具 manifest |

## R 包依赖

- clusterProfiler (>=4.18.0)
- enrichplot (>=1.30.0)
- DOSE (>=4.4.0)
- org.Hs.eg.db / org.Mm.eg.db
- msigdbr (>=26.1.0)
- ggplot2

## Capability Bundle

- `r-enrichment` (>=1.0.0) — 包含物种注释数据库和 MSigDB 基因集

## 平台

- ✅ macos-arm64
- ⏳ macos-x64 / windows-x64 / linux-x64 (planned)

## 注意

此插件使用 legacy `manifest.json` 格式（非 `artifact-manifest.json`）。
升级到完整 Artifact lineage 格式已跟踪为 HRT-1 任务。
