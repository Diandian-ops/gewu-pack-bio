# Seurat-standard 工具范围卡（步骤 0 落地）

> 本文档是 Seurat-standard 工具的 **ground truth**, 写代码前由此推导所有工程决策。
> 按 `.kiro/steering/tool-dev.md` "步骤 0：工具范围与官方文档" 流程产出。
> 完成后将转录到 `seurat.R` 顶部注释。

## Tool ID

`seurat-standard`

## Version targeted

- BioF3 R 容器: `R 4.5.1` (2025-06-13)
- 关键包版本（已 `Rscript packageVersion()` 验证）:
  - `Seurat 5.4.0`
  - `SeuratObject 5.3.0`
  - `ggplot2 3.5.2` (这是当前锁定的版本，配合 enrichplot 的兼容性 fix)
  - `dplyr 1.1.4`
  - `patchwork 1.3.2`
  - `Matrix 1.7.4`

## Canonical workflow source

**官方主 vignette**: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html ("Guided Clustering Tutorial", Satija Lab)

**标准数据集**: PBMC 3k (10X Genomics 公开发布，2,700 cells × 13,714 genes，可在 https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz 下载)

> 注：BioF3 自己的 `docs/single-cell/module03.md` 教程也以 PBMC 3k 为标准数据集，工具与教程口径一致。

## Steps (in-scope) — Canonical Workflow 9 步

按官方 PBMC 3k tutorial 顺序：

1. **Load** — `Read10X()` 读 10X 三件套 (matrix.mtx + barcodes.tsv + features.tsv) → `CreateSeuratObject(min.cells=3, min.features=200)`
2. **QC + filter** — `PercentageFeatureSet(pattern="^MT-")` 计算线粒体比例 → `VlnPlot()` 三联图 → `subset(nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)` 过滤
3. **Normalize** — `NormalizeData(normalization.method="LogNormalize", scale.factor=10000)`
4. **HVGs** — `FindVariableFeatures(selection.method="vst", nfeatures=2000)` → `VariableFeaturePlot()` 标记 top 10
5. **Scale** — `ScaleData(features = all.genes)` 全基因 scale (官方默认)
6. **PCA** — `RunPCA(features = VariableFeatures(...))` → `ElbowPlot()` 查肘点 → 决定 dims
7. **Cluster** — `FindNeighbors(dims=1:10)` → `FindClusters(resolution=0.5)`
8. **UMAP** — `RunUMAP(dims=1:10)` → `DimPlot(reduction="umap", label=TRUE)`
9. **Markers** — `FindAllMarkers(only.pos=TRUE)` → top10 → `DoHeatmap()` + `VlnPlot()` 验证 top markers

## Out-of-scope (不做)

`_data.ts` description 字段必须显式标注:

- ❌ **SCTransform 归一化** — 进阶，留给后续工具或用户在 RunModal 里改
- ❌ **Integration / Harmony / RPCA** — 多样本整合属于不同工具范围
- ❌ **WNN 多组学** — CITE-seq / 10X Multiome 不在范围
- ❌ **Spatial 分析** — `docs/single-cell/module09` 涉及，单独工具
- ❌ **Pseudotime / RNA velocity** — 后续可能开 monocle3 / scvelo 工具
- ❌ **Cell-type annotation** — 仅给 cluster ID，不做自动注释 (SingleR 等)
- ❌ **Doublet 检测** — 工具有相关教程但本工具不集成 (避免依赖额外包)
- ❌ **手动 marker → 细胞类型映射** — 留给用户在 RunModal 里改

> 在 `_data.ts` description 写: "Seurat 标准流程: QC → 归一化 → 高变基因 → PCA → 聚类 → UMAP → Marker 识别。不含: 整合、SCTransform、空间组学、轨迹分析、自动细胞类型注释。"

## Default params (从官方 PBMC 3k tutorial)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `min.cells` | 3 | 官方 `CreateSeuratObject` |
| `min.features` | 200 | 官方 `CreateSeuratObject` |
| `mt_pattern` | `^MT-` | 官方 `PercentageFeatureSet`（人类）。⚠️ 鼠用 `^mt-` |
| `nFeature_max` | 2500 | 官方 PBMC3k 教程 subset 上限 (BioF3 暴露给用户调) |
| `nFeature_min` | 200 | 同上 |
| `percent_mt_max` | 5 | 同上 |
| `normalization.method` | `LogNormalize` | 官方默认 |
| `scale.factor` | 10000 | 官方默认 |
| `selection.method` | `vst` | 官方默认 (FindVariableFeatures) |
| `nfeatures` | 2000 | 官方默认 (FindVariableFeatures) |
| `dims` | `1:10` | 官方默认。注：BioF3 docs/single-cell/module03 用 `1:20`（dataset-specific 选择）。工具默认跟官方，文档说明可调 |
| `resolution` | 0.5 | 官方默认 (FindClusters) |
| `only.pos` | TRUE | 官方默认 (FindAllMarkers) |
| `top_n_markers` | 10 | 用于 DoHeatmap 展示 top markers (官方教程用 top10) |
| `species` | `human` | 默认 human, 切换 `mouse` 时改 mt_pattern 为 `^mt-` |

**所有偏离官方默认的注解:**
- 无（本工具完全跟随官方 PBMC 3k tutorial 默认值）

## Expected outputs (10 张图 + 6 份中间数据)

### 图 (按 canonical 9 步顺序对应)

| # | 文件名 | 步骤 | 函数 | 中间数据 |
|---|--------|------|------|---------|
| 1 | `qc_violin.png` | QC | VlnPlot(c("nFeature_RNA","nCount_RNA","percent.mt")) | qc_metrics.csv |
| 2 | `qc_scatter.png` | QC | FeatureScatter × 2 | qc_metrics.csv |
| 3 | `hvg_plot.png` | HVGs | VariableFeaturePlot + LabelPoints(top10) | hvg_table.csv |
| 4 | `pca_dimplot.png` | PCA | DimPlot(reduction="pca") | seurat_obj.rds |
| 5 | `elbow_plot.png` | PCA | ElbowPlot() | pca_stdev.csv |
| 6 | `umap_clusters.png` | UMAP | DimPlot(reduction="umap", label=TRUE) | umap_coords.csv |
| 7 | `markers_heatmap.png` | Markers | DoHeatmap(top10 per cluster) | markers_table.csv |
| 8 | `markers_violin.png` | Markers | VlnPlot(top 6 markers) | markers_table.csv |
| 9 | `feature_plot.png` | Markers | FeaturePlot(top 6 markers, reduction="umap") | umap_coords.csv + counts |
| 10 | `pca_dim_loadings.png` | PCA | VizDimLoadings(dims=1:2) | pca_loadings.csv |

### 中间数据

| 文件 | 内容 | 用途 |
|------|------|------|
| `qc_metrics.csv` | 每个 cell 的 nFeature_RNA / nCount_RNA / percent.mt | 重画 QC 三联 + scatter |
| `hvg_table.csv` | 高变基因列表 + 平均表达 + 标准化方差 | 重画 HVG plot |
| `pca_stdev.csv` | PCA 各 PC 的 stdev | 重画 elbow plot |
| `pca_loadings.csv` | PC1-2 top loading 基因和值 | 重画 VizDimLoadings |
| `umap_coords.csv` | 每 cell 的 UMAP_1 / UMAP_2 / cluster | 重画 UMAP（任意配色 / 任意标 cluster）|
| `markers_table.csv` | FindAllMarkers 完整输出 | 重画 heatmap / violin / 用户筛选基因 |
| `seurat_obj.rds` | 完整 Seurat 对象 | 复杂二次绘图（如 FeaturePlot 任意基因）|

## Demo 数据准备

- 上游数据集: PBMC 3k (官方 link)
- 预处理: 不预处理，直接用三件套
- 用户上传形式: 三个文件 (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`) 或单个 .h5
- 预期跑完时间: 约 1-2 分钟（PBMC 3k 是小数据集，scaling 和 PCA 是瓶颈）

## 工具间数据联动

- ❌ 上游不导入 (Seurat 是单细胞流程入口)
- ✅ 下游可被导入: `seurat_obj.rds` 可作为后续 CellChat / monocle3 等工具的输入（届时实现）

## 已读官方文档（T1 检查清单）

- [x] 官方主 vignette: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
- [x] BioF3 容器版本验证: Seurat 5.4.0 / R 4.5.1
- [x] BioF3 自己的教程 `docs/single-cell/module01-03` 对照口径
- [ ] (T2 用时查) Seurat v5 cheat sheet: https://satijalab.org/seurat/articles/seurat5_essential_commands.html

## 已知风险

1. **大数据集容器 OOM**: 当前 R 容器 4G memory 上限，PBMC 3k (2700 cells) 没问题，但用户传 50k+ cells 会撞上限。需要 R 脚本里加 cell 数 sanity check + 给出明确报错，而不是让 OOM 杀容器
2. **dims 选择**: tutorial 说"7-12 都合理"。工具默认 `dims=1:10` 跟官方，但如果 demo 数据 elbow plot 看起来不像 PBMC3k 经典形状（用户上传的非标准数据），UMAP 可能不理想。文档说明这是可调的
3. **Mitochondrial pattern**: 物种切换时（mouse vs human）需要换 `^MT-` vs `^mt-`，UI 必须暴露 species 选项

## 下一步

按 13 步流程：
- 步骤 1: ID 已定 `seurat-standard`
- 步骤 2: `_data.ts` 添加工具定义（in-scope / out-of-scope 写入 description）
- 步骤 3: `routes/tools.js` 注册（ID + 显示名 "Seurat 标准流程"）
- 步骤 4: `index.tsx` `TOOL_COLORS` + `REPORT_URLS`
- 步骤 5: 编写 `seurat.R` (顶部注释从本文档 0.6 段抄过去)
- ... (按 13 步走完)
