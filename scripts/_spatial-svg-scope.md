# Spatial SVG 空间高变基因多方法鉴定工具范围卡

## Tool ID

`spatial-svg`（py-server task type: `spatial_svg`）

## 版本锁定

- 本地 Python runtime: Python 3.11+
- 关键包: `scanpy` / `squidpy` / `matplotlib`

## Canonical workflow source

**官方文档**: Squidpy spatial autocorr vignette (Moran's I + Geary's C) + Sepal

## Steps (in-scope) — 三方法交叉 SVG 鉴定

1. **加载输入** — 从本地 tool job 目录读取上传的 h5ad
2. **自动补全预处理** — 检测 log1p / HVG / spatial_connectivities 缺失则补
3. **Method 1: Moran's I** — sq.gr.spatial_autocorr(mode='moran')，全局空间自相关
4. **Method 2: Geary's C** — sq.gr.spatial_autocorr(mode='geary')，局部空间差异
5. **Method 3: Sepal** — sq.gr.sepal，扩散时间法，偏好模式化分布
6. **跨方法交集** — robust SVG = ≥2 方法显著
7. **可视化** — 4 张图：父聚类上下文 / Top-4 robust SVGs 空间分布 / 三方法 Venn / 分数相关性散点
8. **结果保存** — top_svgs.csv（per-gene: moran_I/geary_C/sepal_score/n_methods/robust）+ result_spatial.h5ad + manifest

## Out-of-scope (不做)

- ❌ **单方法 Moran's I**（见 spatial-standard 的 top SVGs 段，本工具是多方法深化）
- ❌ **细胞类型解卷积**（见 spatial-deconv）
- ❌ **配体-受体**（见 spatial-lr）
- ❌ **Spatial domain detection (STAGATE / SpaGCN)**（独立工具）
- ❌ **SVG 功能富集**（用户自行做 GO/KEGG）

## Default params (脚本默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `topN` | `50` | per-method top-N 用于交集 |
| `nPerms` | `100` | Moran/Geary 排列次数 |
| `minMethods` | `2` | robust 阈值（≥2 方法显著）|
| `nGenesTest` | `1000` | 限制测试基因数控时长 |

## Expected outputs

图（PNG）：plot_001_spatial_clusters / plot_002_top_robust_svgs / plot_003_method_venn / plot_004_score_corr
数据表：top_svgs.csv / metadata.csv / features.csv
文件：result_spatial.h5ad / summary.json / manifest.json

## 工具间数据联动

- ✅ 上游：spatial-standard 的 result_spatial.h5ad（推荐通过 /spatial 平台链）
- ✅ 下游：robust SVG 列表可回导做功能富集
