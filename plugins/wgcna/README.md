# wgcna

> **科学子问题**: 基因表达数据中是否存在共表达模块？这些模块与外部性状（如处理条件、临床指标）是否显著关联？模块内的 Hub 基因是谁？

## 职责

执行 WGCNA 标准流程（Tutorial I 7 步）：

1. **QC** — `goodSamplesGenes()` + 样本聚类树
2. **软阈值选择** — `pickSoftThreshold()` 自动选 power（R² ≥ 0.85）
3. **网络构建** — `blockwiseModules(networkType="signed")`
4. **模块聚类** — `plotDendroAndColors()` 动态合并
5. **模块-性状关联** — `moduleEigengenes()` + `cor()` + p 值
6. **Hub 基因** — `intramodularConnectivity()` top kIM
7. **GS vs MM 散点** — 最显著模块内基因显著性 vs 模块成员度

## DAG 位置

```
expression-qc → wgcna → (模块基因子集) → go-kegg / gsea (二跳绑定)
```

WGCNA 产出的模块基因列表可作为下游 GO/KEGG 富集的输入（通过 bind_artifact 二跳绑定）。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| expression | file | ✓ | — | 表达矩阵 (CSV: gene × sample, log2 标准化) |
| traits | file | ✓ | — | 性状矩阵 (CSV: sample × numeric traits) |
| min_module_size | number | — | `30` | 最小模块大小 |
| merge_cut_height | number | — | `0.25` | 合并相似模块的高度阈值 |
| network_type | select | — | `signed` | signed / unsigned |
| soft_power | number | — | 自动 | 软阈值 power（留空自动选） |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `sample_dendrogram.png` | image | 样本聚类树 |
| `soft_threshold.png` | image | 软阈值选择图 |
| `module_dendrogram.png` | image | 模块聚类树 + 颜色标注 |
| `module_trait_heatmap.png` | image | 模块-性状相关性热图 |
| `gs_mm_scatter.png` | image | GS vs MM 散点图 |
| `module_genes.csv` | table | 模块基因列表 |
| `module_trait_cor.csv` | table | 模块-性状相关系数 + p 值 |
| `hub_genes.csv` | table | Hub 基因 (top kIM) |
| `network.rds` | rds | 完整网络 RDS 对象 |
| `report.html` | report | HTML 报告 |
| `manifest.json` | manifest | 工具 manifest |

## R 包依赖

- WGCNA (>=1.74)
- flashClust (>=1.1.4)
- impute (>=1.84.0)
- ggplot2 (>=4.0.0)
- dplyr (>=1.1.0)
- jsonlite (>=1.8.0)
- base64enc (>=0.1-3)

## Capability Bundle

- `r-wgcna` (>=1.0.0) — 包含 WGCNA 及其依赖包

## 平台

- ✅ macos-arm64
- ⏳ macos-x64 / windows-x64 / linux-x64 (planned)

## 规模限制

- 上限 ~10000 基因 × ~500 样本
- 建议先用 expression-qc 筛选高变基因再输入

## 注意

此插件使用 legacy `manifest.json` 格式（非 `artifact-manifest.json`）。
升级到完整 Artifact lineage 格式已跟踪为 HRT-1 任务。
