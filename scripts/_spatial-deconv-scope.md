# Spatial Deconv 空间转录组解卷积工具范围卡

## Tool ID

`spatial-deconv`（py-server task type: `spatial_deconv`）

## 版本锁定

- 本地 Python runtime: Python 3.11+
- 关键包: `scanpy` / `squidpy` / `tangram` / `torch` / `pandas`

## Canonical workflow source

**官方文档**: Tangram tutorial (Biancalani et al., 2021 Nat Commun) + Squidpy deconvolution

## Steps (in-scope) — Tangram 细胞类型解卷积

1. **加载空间数据** — 从本地 tool job 目录读取上传的 h5ad
2. **自动补全标准化** — 检测 log1p 缺失则 normalize_total + log1p
3. **加载单细胞参考** — sq.datasets.sc_mouse_cortex（~21k cells，自动下采样防 OOM）
4. **挑选 training genes** — sc 与 st 共有基因，按 sc mean expression top-N
5. **Tangram 映射** — tg.pp_adatas + tg.map_cells_to_space（cells 模式，CPU）
6. **投影细胞注释** — tg.project_cell_annotations → 每 spot 细胞类型丰度矩阵
7. **可视化** — 4 张图：每 spot 主要细胞类型 / Top-4 细胞类型空间丰度 / 整体组成 / 共定位 heatmap
8. **结果保存** — cell_type_abundance.csv + result_spatial.h5ad + manifest + report.html

## Out-of-scope (不做)

- ❌ **cell2location / RCTD 解卷积**（其他方法，独立工具）
- ❌ **spot-level 聚类重新注释**（见 spatial-standard）
- ❌ **配体-受体**（见 spatial-lr，解卷积结果可做下游但本工具不含）
- ❌ **自定义单细胞参考**（当前固定用 sc_mouse_cortex；用户参考需改代码）

## Default params (脚本默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `refLabelKey` | `cell_subclass` | sc_mouse_cortex 细分类型；fallback cell_class |
| `maxRefCells` | `5000` | 防 OOM 下采样上限 |
| `nTrainingGenes` | `1000` | Tangram 常用 |
| `num_epochs` | `300` | Tangram cells 模式常用 |

## Expected outputs

图（PNG）：plot_001_dominant_celltype / plot_002_top_celltypes_grid / plot_003_celltype_proportion / plot_004_colocalization
数据表：cell_type_abundance.csv / metadata.csv / features.csv
文件：result_spatial.h5ad / summary.json / manifest.json / report.html

## 工具间数据联动

- ✅ 上游：spatial-standard 的 result_spatial.h5ad
- ✅ 下游：细胞类型丰度可进 spatial-lr 做细胞类型间通讯
