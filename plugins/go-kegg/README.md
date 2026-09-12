# go-kegg

> **科学子问题**: 差异表达基因是否得到 GO、KEGG 功能通路证据支持？

## 职责

对差异分析产生的显著基因列表执行 GO（BP/CC/MF）和 KEGG 通路富集分析，输出富集结果表和可视化报告。

## DAG 位置

```
deg-standardizer → go-kegg
```

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| genes | file | ✓ | — | 基因列表 (CSV: 一列基因 ID) |
| organism | text | — | `hsa` | 物种 (hsa=human, mmu=mouse) |
| gene_id_type | text | — | `symbol` | ID 类型 (symbol/ensembl/entrez) |
| padj_cutoff | number | — | `0.05` | 富集 padj 阈值 |
| ont | text | — | `ALL` | GO 本体 (BP/CC/MF/ALL) |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `go_results.csv` | table | GO 富集结果 |
| `kegg_results.csv` | table | KEGG 通路富集结果 |
| `report.html` | report | HTML 报告（含图） |
| `manifest.json` | manifest | 工具 manifest |

## R 包依赖

- clusterProfiler (>=4.18.0)
- enrichplot (>=1.30.0)
- DOSE (>=4.4.0)
- org.Hs.eg.db / org.Mm.eg.db
- KEGGREST
- ggplot2

## Capability Bundle

- `r-enrichment` (>=1.0.0) — 包含物种注释数据库和 KEGG REST 缓存

## 平台

- ✅ macos-arm64
- ⏳ macos-x64 / windows-x64 / linux-x64 (planned)

## 注意

此插件使用 legacy `manifest.json` 格式（非 `artifact-manifest.json`）。
升级到完整 Artifact lineage 格式已跟踪为 HRT-1 任务。
