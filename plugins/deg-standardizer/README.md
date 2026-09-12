# deg-standardizer

> **科学子问题**: 差异分析结果是否能够转换成稳定、可供下游插件消费的标准 DEG Artifact？

## 职责

将 DESeq2、limma 或兼容差异结果表转换为标准 DEG Artifact，供下游 volcano-plot、pca-explorer、go-kegg、gsea 消费。

## 关键原则

- **不覆盖上游原始结果** — 只生成标准化副本
- **不伪造缺失 p 值或 padj** — NA padj 标记为非显著（padj=1），不编造
- **缺少 padj 列时阻断** — 不能没有 FDR 控制就标准化
- **ranked_genes 保留全部基因** — 不只保留显著基因，供 GSEA 使用
- **阈值属于参数和科学决策** — padj_threshold 和 effect_threshold 必须记录

## DAG 位置

```
deseq2 → deg-standardizer ├── volcano-plot
                           ├── pca-explorer
                           ├── go-kegg
                           └── gsea
```

## 标准输出列

| 列 | 说明 |
|---|---|
| `gene` | 基因 ID |
| `gene_id_type` | ID 类型 (symbol/ensembl/entrez) |
| `log2FoldChange` | 效应量 |
| `pvalue` | 原始 p 值 |
| `padj` | FDR 校正后 p 值 |
| `direction` | up / down / neutral |
| `significant` | TRUE / FALSE (基于 padj + effect 阈值) |
| `source_method` | 来源方法 (DESeq2/limma/edgeR) |

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| deg_results | file | ✓ | — | 差异分析结果表 (CSV/TSV) |
| deg_artifact | text | — | — | 上游 DEG Artifact 引用 |
| gene_col | text | — | `gene` | gene ID 列名 |
| gene_id_type | text | — | `symbol` | ID 类型 |
| effect_col | text | — | `log2FoldChange` | effect size 列名 |
| pval_col | text | — | `pvalue` | p value 列名 |
| padj_col | text | — | `padj` | adjusted p value 列名 |
| padj_threshold | number | — | `0.05` | 显著性阈值 |
| effect_threshold | number | — | `1` | |log2FC| 阈值 |
| source_method | text | — | `DESeq2` | 来源方法 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `deg_standardized.csv` | table | 完整标准化 DEG 表 |
| `deg_up.csv` | table | 显著上调基因 |
| `deg_down.csv` | table | 显著下调基因 |
| `ranked_genes.csv` | table | 全量排序基因列表 (含 rank 列) |
| `filter_summary.json` | json | 过滤摘要 |
| `report.html` | report | HTML 标准化报告 |
| `artifact-manifest.json` | manifest | Artifact 清单 |

## 状态

- `valid` — 标准化成功
- `warning` — 存在非致命问题（NA 值、重复基因、无显著基因）
- `blocked` — 缺少必需列或表为空

## 空结果处理

无显著基因是**有效科学结果**，不是执行失败：
- 仍生成完整的 `deg_standardized.csv`、`ranked_genes.csv`
- `deg_up.csv` 和 `deg_down.csv` 为空表（仅有表头）
- 状态为 `warning`，明确标注 "valid empty result"

## R 包依赖

- jsonlite (必需)
- digest (必需)

## 平台

- ✅ macos-arm64
- ⏳ macos-x64 / windows-x64 / linux-x64 (planned)
