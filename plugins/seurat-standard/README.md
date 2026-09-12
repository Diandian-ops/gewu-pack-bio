# seurat-standard

> **科学子问题**: 在给定 10X 三件套 / CSV / TSV / MTX / H5 输入下，单细胞 RNA-seq 标准 Seurat 流程能否在 biof3-r-runtime (dev Electron 边界内) 跑通 QC → Normalize → HVG → Scale → PCA → Cluster → UMAP → Markers？

## 职责

使用 Seurat 5.x 对单细胞 RNA-seq 输入执行官方 PBMC 3k tutorial 标准分析管线：

1. **Load** — 自动探测 10X 三件套 (matrix.mtx + barcodes.tsv + features.tsv) 或单个 CSV/TSV/MTX/H5/RDS
2. **QC + filter** — `PercentageFeatureSet(pattern="^MT-")` + `subset(nFeature_RNA > min & < max & percent.mt < max_mt)`
3. **Normalize** — `NormalizeData(LogNormalize, scale.factor=10000)`
4. **HVGs** — `FindVariableFeatures(selection.method="vst", nfeatures=N)`
5. **Scale + PCA** — `ScaleData(features=all.genes)` + `RunPCA` + `ElbowPlot`
6. **Cluster** — `FindNeighbors(dims=1:N)` + `FindClusters(resolution=R)`
7. **UMAP** — `RunUMAP(dims=1:N)` + `DimPlot(reduction="umap", label=TRUE)`
8. **Markers** — `FindAllMarkers(only.pos=TRUE)` → `DoHeatmap` + `VlnPlot` + `FeaturePlot` 验证 top markers
9. **落盘** — `seurat_obj.rds` + figures + tables + report.html + manifest.json

## DAG 位置

```
sample-metadata-validator → seurat-standard → (future: cellchat / singleR / monocle3)
```

`seurat-standard` 是 single-cell 流程的 R 侧入口；`scripts/seurat-standard.R` 接受 10X 三件套或单个矩阵文件（RDS / CSV / TSV / MTX / H5）。`seurat_obj.rds` 可作为后续 CellChat / monocle3 工具的输入。

## 复用既有 ground truth

本插件的 R 脚本 `scripts/seurat-standard.R` 与既有 `core/biof3-server/tool-scripts/seurat-standard.R` 同源，遵循 `core/biof3-server/tool-scripts/_seurat-standard-scope.md` 步骤 0 的全部 in-scope 步骤（9 步）+ out-of-scope 边界声明。**原 `core/biof3-server/tool-scripts/seurat-standard.R` 保留**作为向后兼容；本插件是独立拷贝，作为 builtin-canonical plugin 暴露给 BioF3 调度器。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| `matrix_file` | file (`.csv/.tsv/.mtx/.h5/.rds`) | ✓ | — | 10X 三件套或单文件；自动探测 |
| `species` | select | — | `human` | `human` → `^MT-`；`mouse` → `^mt-` |
| `min_features` | number | — | `200` | CreateSeuratObject min.features |
| `max_features` | number | — | `2500` | subset nFeature_RNA 上限 |
| `max_mt` | number | — | `5` | subset percent.mt 上限（%）|
| `nfeatures` | number | — | `2000` | FindVariableFeatures nfeatures |
| `dims` | number | — | `10` | FindNeighbors/RunUMAP dims=1:N |
| `resolution` | number | — | `0.5` | FindClusters resolution |
| `top_n_markers` | number | — | `10` | 每个 cluster top N marker (heatmap/vln/feature) |

## 输出（8 件套）

| 文件 | 类型 | 说明 |
|------|------|------|
| `seurat_obj.rds` | data | 完整 Seurat 对象（PCA/UMAP/leiden/FindAllMarkers）|
| `markers_table.csv` | table | FindAllMarkers 完整输出 |
| `qc_metrics.csv` | table | 每 cell 的 nFeature / nCount / percent.mt |
| `umap_clusters.png` | plot | UMAP 按 seurat_cluster 着色（双格式 PNG + PDF）|
| `umap_coords.csv` | table | 每 cell 的 UMAP_1 / UMAP_2 / cluster |
| `hvg_plot.png` | plot | VariableFeaturePlot + LabelPoints(top10) |
| `report.html` | report | 完整分析报告（SCI 主题）|
| `manifest.json` | manifest | Job manifest（含 stats）|

辅助中间产物（10 张 PNG/PDF 图 + 4 份 CSV）：

| 文件 | 说明 |
|------|------|
| `qc_violin.png` + `.pdf` | QC 三联小提琴图 |
| `qc_scatter.png` + `.pdf` | QC scatter（nCount vs nFeature / percent.mt）|
| `pca_dimplot.png` + `.pdf` | PCA DimPlot |
| `pca_dim_loadings.png` + `.pdf` | VizDimLoadings (dims=1:2) |
| `elbow_plot.png` + `.pdf` | PCA ElbowPlot |
| `markers_heatmap.png` + `.pdf` | DoHeatmap(top N per cluster) |
| `markers_violin.png` + `.pdf` | VlnPlot(top 6 markers) |
| `feature_plot.png` + `.pdf` | FeaturePlot(top 6 markers) |
| `hvg_table.csv` | 高变基因列表 + mean / variance / standardized |
| `pca_stdev.csv` | PCA 各 PC 的 stdev |
| `pca_loadings.csv` | PC1-2 top loadings |
| `summary.txt` | 分析摘要 |

## 状态

**dev Electron 边界内 verified**（G0 阶段）：

- biof3-r-runtime v1.0.0 已装 `Seurat>=5.0.0` + `SeuratObject>=5.0.0` + `ggplot2` + `patchwork` + `dplyr` + `Matrix` + `jsonlite` + `base64enc` + `png`（见 `resources/built-in-plugins/r-runtime-packages.txt`）
- 不依赖额外 R 包安装（PHASE 1 SUB-1.2 任务边界内 runtime **不扩展**）
- 离线可执行；产出 progress.json 驱动 UI 进度条

> **G0 阶段纪律（按 AGENTS.md "当前阶段纪律"）**：本插件仅在 dev Electron 边界内 verified；packaged / Provider / GA / dist staleness 属 G1+ 阶段事项，按本阶段纪律**不外推**。`qualificationEligible: false` 是当前 G0 开发期常态，不是缺陷、不被视为缺口。

## R 包依赖

- Seurat ≥ 5.0.0（biof3-r-runtime 声明）
- SeuratObject ≥ 5.0.0
- ggplot2 ≥ 3.4.0
- patchwork ≥ 1.3.0（VlnPlot / FeaturePlot 拼图）
- dplyr ≥ 1.1.0（top markers 取 top N per cluster）
- Matrix ≥ 1.7.0（CsparseMatrix 转换）
- jsonlite ≥ 1.8.0（params.json 解析 + manifest 写盘）
- base64enc ≥ 0.1-3（report.html 内嵌 PNG base64）
- png ≥ 0.1-8（PNG hash 重写稳定 sha256）

## 平台

- ✅ macos-arm64 / windows-x64（dev Electron runtime 就位即可跑）
- ⏳ packaged 边界验证属 G1+ 阶段，本批次不外推

## 关联 Skill

- **PHASE 1 SUB-1.6**：`sc-analysis-standard` skill 在 SUB-1.6 子批才创建，**不动**（单工作包串行纪律）
- 本插件 standalone 可被任意 `biof3-r-runtime` R 调用链直接消费

## 注意

`scripts/seurat-standard.R` 末尾对所有输出 PNG 调 `canonicalize_png()`，用 `png::writePNG(png::readPNG(f), f)` 重写以稳定 sha256（与 Python 侧 `Pillow Image.save` 同义），保证下游 artifact-lineage / Result Studio hash 不漂移。