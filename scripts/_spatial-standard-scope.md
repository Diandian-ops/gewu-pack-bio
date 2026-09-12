# Spatial Standard 空间转录组标准分析工具范围卡

## Tool ID

`spatial-standard`（py-server task type: `spatial_standard`）

## 版本锁定

- 本地 Python runtime: Python 3.11+
- 关键包: `scanpy` / `squidpy` / `anndata` / `matplotlib`

## Canonical workflow source

**官方文档**: Squidpy tutorials (Palla et al., 2022 Nat Commun) + Scanpy best practices

## Steps (in-scope) — 空间转录组标准 4 件套

1. **数据加载 + QC** — h5ad 输入，检测 platform / has_image
2. **标准化 + HVG + PCA** — normalize_total → log1p → top HVG → scale → PCA
3. **邻居图 + Leiden 聚类 + UMAP** — sq.gr.spatial_neighbors / sc.tl.leiden / sc.tl.umap
4. **空间可视化** — 4 张图：空间聚类 / UMAP / 邻域富集 z-score / Top SVGs 空间分布
5. **空间高变基因 (Moran's I)** — sq.gr.spatial_autocorr top-N
6. **结果保存** — result_spatial.h5ad + metadata.csv + features.csv + summary.json + manifest.json + report.html

## Out-of-scope (不做)

- ❌ **多方法 SVG 鉴定**（Moran + Geary + Sepal 交叉验证，见 spatial-svg 工具）
- ❌ **细胞类型解卷积**（见 spatial-deconv 工具）
- ❌ **配体-受体分析**（见 spatial-lr 工具）
- ❌ **子聚类**（当前工具集不包含）
- ❌ **整合批次**（当前工具集不包含）
- ❌ **HE 图像分割 / 对齐**（用户预处理）

## Default params (脚本默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `nPCs` | `30` | Scanpy 常用 |
| `resolution` | `0.5` | Leiden 常用 |
| `coordType` | `generic` | Squidpy 默认（兼容非网格）|
| `nNeighbors` | `6` | Visium 六邻域 |
| `topSVGs` | `50` | 可视化常用 |
| `nHVGs` | `2000` | Scanpy 默认 |

## Expected outputs

图（PNG）：plot_001_spatial_clusters / plot_002_umap / plot_003_nhood_enrichment / plot_004_top_svgs
数据表：top_svgs.csv / metadata.csv / features.csv
文件：result_spatial.h5ad / summary.json / manifest.json / report.html

## 工具间数据联动

- ✅ 上游：用户提供带 obsm['spatial'] 的 h5ad
- ✅ 下游：result_spatial.h5ad 可作为本地文件输入 spatial-svg / spatial-deconv / spatial-lr
