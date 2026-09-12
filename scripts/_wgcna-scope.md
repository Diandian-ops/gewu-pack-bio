# WGCNA 工具范围卡（步骤 0 落地）

> 按 `.kiro/steering/tool-dev.md` "步骤 0：工具范围与官方文档" 流程产出。
> Seurat-standard 之后的第二个落地实例。

## Tool ID

`wgcna`

## Version targeted

- BioF3 R runtime: `R 4.5.2`
- 关键包版本（已 `Rscript packageVersion()` 验证）:
  - `WGCNA 1.74`
  - `flashClust 1.1-4`
  - `impute 1.84.0`
  - `preprocessCore 1.72.0`
  - `GO.db 3.22.0`

## Canonical workflow source

**官方主 vignette**: WGCNA Tutorial I — Network construction and module detection
(Horvath lab, https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)

**主 reference**: Langfelder & Horvath, 2008. *WGCNA: an R package for weighted correlation network analysis*. BMC Bioinformatics 9:559.

## Steps (in-scope) — Tutorial I 的 7 步

1. **Load + QC** — 读表达矩阵 (rows=samples × cols=genes after transpose) + traits 矩阵; `goodSamplesGenes()` 检测离群基因 / 样本; sample dendrogram 树检测异常样本
2. **Soft threshold** — `pickSoftThreshold(powerVector = c(1:10, seq(12, 30, 2)))` 自动选 power; 输出 scale-free fit (R²) + mean connectivity 双图; 推荐 R² ≥ 0.85 的最小 power
3. **Network construction (one-step)** — `blockwiseModules(power=auto, TOMType="signed", minModuleSize=30, mergeCutHeight=0.25, networkType="signed")`
4. **Module dendrogram + colors** — `plotDendroAndColors()` 展示动态合并前后的模块
5. **Module-trait correlation** — `moduleEigengenes(MEs)` + `cor(MEs, traits)` + `corPvalueStudent()` ; 输出 ME-trait 热图（cells: r 值 + p 值文本）
6. **Hub genes per module** — `intramodularConnectivity()` 计算 kIM (intramodular connectivity); 每模块 top 10 hub 基因
7. **Gene-trait significance scatter** — 对最显著的模块画 GS vs MM 散点图（gene significance vs module membership），验证模块基因与性状关联

## Out-of-scope (不做)

- ❌ **Tutorial II Trait Association 进阶**（GO enrichment 内置 — 用户可以把 module 基因导出到 BioF3 GO/KEGG 工具）
- ❌ **Tutorial III Consensus Networks**（多数据集共识，独立工具）
- ❌ **Block-wise 大数据集模式**（>20k 基因；本工具限 5000-10000 基因）
- ❌ **Module preservation analysis**（需要 validation 数据集）
- ❌ **Signed hybrid network type**（多用 signed）
- ❌ **Bicor robust correlation**（默认 pearson；用户可在 RunModal 改）
- ❌ **CytoScape 网络导出**（exportNetworkToCytoscape 工具）

> `_data.ts` description: "WGCNA 标准流程（Tutorial I）：QC → 软阈值 → 网络构建 → 模块识别 → 模块-性状关联 → Hub 基因。不含: 多数据集共识、模块保守性分析、Cytoscape 导出。最大 ~10000 基因 × ~500 样本。"

## Default params (官方默认 / Tutorial I 推荐)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `corType` | `pearson` | 官方默认（bicor 更鲁棒但慢，留给高级用户） |
| `networkType` | `signed` | 官方推荐（unsigned 会让正负相关混在一起）|
| `TOMType` | `signed` | 配合 networkType=signed |
| `minModuleSize` | `30` | 官方默认（小数据集可降到 20）|
| `mergeCutHeight` | `0.25` | 官方默认（合并相似模块的高度阈值）|
| `power` | **auto** | 跑 `pickSoftThreshold(c(1:10, seq(12,30,2)))` 找 R² ≥ 0.85 的最小 power; 找不到则取 R² 最大的 |
| `deepSplit` | `2` | 官方默认（0-4，越大越细分）|
| `maxBlockSize` | `5000` | BioF3 限定（防止超大网络 OOM）|
| `randomSeed` | `54321` | 官方默认 |

**偏离官方默认的注解**:
- `maxBlockSize=5000` 偏离官方默认 5000（实际就是默认值，写出来是为了防止后续误改）。物理意义：> 5000 基因会自动分块，影响一致性

## Expected outputs (8 张图 + 6 份中间数据)

### 图

| # | 文件名 | 步骤 | 函数 | 中间数据 |
|---|--------|------|------|---------|
| 1 | `sample_dendrogram.png` | QC | `hclust(dist(datExpr))` | qc_summary.csv |
| 2 | `soft_threshold.png` | Step 2 | `pickSoftThreshold` (双图: SFT R² + mean k) | sft_results.csv |
| 3 | `gene_dendrogram.png` | Step 3-4 | `plotDendroAndColors(geneTree, mergedColors)` | module_assignment.csv |
| 4 | `module_eigengenes.png` | Step 5 | ME barplot 各样本 | module_eigengenes.csv |
| 5 | `module_trait_heatmap.png` | Step 5 | r 值 heatmap (cells: r + p) | module_trait_corr.csv |
| 6 | `module_size_barplot.png` | Step 4 | 每模块基因数条形图 | module_assignment.csv |
| 7 | `gs_mm_scatter.png` | Step 7 | 最显著模块 GS vs MM 散点 | gs_mm_data.csv |
| 8 | `hub_gene_table.png` | Step 6 | top hub 基因 dot/bar 图（per module）| hub_genes.csv |

### 中间数据

| 文件 | 内容 | 用途 |
|------|------|------|
| `module_assignment.csv` | 每个基因 → 模块颜色 | 重画 dendrogram / 模块大小 / 用户筛选 |
| `module_eigengenes.csv` | 每个模块的 ME 在每个样本的值 | 重画 ME barplot / 用户做下游 |
| `module_trait_corr.csv` | 模块 × 性状 r + p 矩阵 | 重画热图 |
| `hub_genes.csv` | 每模块 top 10 hub + kIM | 用户筛重点 |
| `sft_results.csv` | pickSoftThreshold 全部 power 拟合结果 | 重画 sft 图 |
| `gs_mm_data.csv` | 每基因的 GS (vs trait) + MM (vs ME) | 重画 GS-MM scatter / 用户挑核心基因 |
| `network.rds` | 完整 blockwiseModules 输出 | 高级用户 / 后续 Cytoscape 导出 |

## Demo 数据准备

- 上游: TCGA-BRCA (UCSC Xena)
  - 表达: HiSeqV2 log2(RSEM+1)
  - 临床: BRCA_clinicalMatrix
- 预处理: 
  - 218 样本（PAM50 平衡，每亚型 ≤ 50）
  - 5000 高变基因（mean > 1, top variance）
  - 已编码 traits (ER/PR/HER2 → 0/1, Stage → 1-4, Age → 数值, PAM50 → 5 binary)
- 文件:
  - `wgcna_expr.csv` 7.0 MB (rows=genes, cols=samples, 含 gene 列)
  - `wgcna_traits.csv` 7.9 KB (rows=samples, 含 sample 列 + 10 traits)
- 预计跑完时间: 1-2 分钟（5000 基因 218 样本 + signed network）

## 工具间数据联动

- ❌ 上游不导入（WGCNA 是入口工具，需要 raw expression 矩阵）
- ✅ 下游可被导入: 用户可以把任意模块的基因列表（从 `module_assignment.csv` 筛颜色）导出到 BioF3 的 GO/KEGG 工具或 GSEA 工具做功能富集

## 已读官方文档（T1 检查清单）

- [x] WGCNA Tutorial I (network construction + module detection) — Horvath lab UCLA
- [x] BioF3 自己的教程 `docs/integration/module04` 对照口径（用 power=8 / unsigned，是 BioF3 自己的选择，**工具默认跟官方更新的 signed**，文档说明可调）
- [x] `blockwiseModules` rdocumentation 默认参数（power=6 是占位，用户必须用 pickSoftThreshold 自动选）
- [x] BioF3 runtime 版本验证: WGCNA 1.74 / R 4.5.2（2026-07-19 demo 真流程）

## 已知风险

1. **样本聚类有离群点**：用户上传的真实数据可能有"主板坏一块"导致某个样本 PCA 飞出，现行流程要决策"自动剔除还是仅警告"。我选**仅警告 + 在 dendrogram 显示阈值线**，不自动剔除（避免静默改用户数据）
2. **scale-free 拟合 R² 上不到 0.85**：小数据集 / 稀疏网络的常见问题。`pickSoftThreshold` 找不到 R² ≥ 0.85 时取**最大 R² 对应的 power**作为 fallback（不报错），用户在报告里能看到"请注意 R² 仅 X.XX"提示
3. **模块全归到 grey**：power 太高 / mergeCutHeight 太低 → 没识别出模块。脚本里加 sanity check：如果非 grey 模块 < 2，报错让用户调参（不勉强出图）
4. **OOM**：5000 基因 × 218 样本 signed TOM ≈ 200MB matrix，4GB 容器可承受。但用户传 10000 基因 × 500 样本会撞内存——脚本加 dim check，超出限制直接拒绝

## 下一步

按 14 步流程：
- 步骤 1: ID 已定 `wgcna`
- 步骤 2: `_data.ts` 已存在但参数和输出过简，按本范围卡更新
- 步骤 3: `routes/tools.js` 已注册 `wgcna.R`
- 步骤 4: `index.tsx` `TOOL_COLORS` + `REPORT_URLS` 之前已注册（保留）
- 步骤 5: 编写 `wgcna.R`（顶部注释从本文档抄）
- 步骤 6-14: 按流程走
