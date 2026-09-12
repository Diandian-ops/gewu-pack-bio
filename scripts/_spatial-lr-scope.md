# Spatial LR 空间配体-受体分析工具范围卡

## Tool ID

`spatial-lr`（py-server task type: `spatial_lr`）

## 版本锁定

- 本地 Python runtime: Python 3.11+
- 关键包: `scanpy` / `squidpy` / `liana` / `pandas`

## Canonical workflow source

**官方文档**: LIANA bivariate tutorial (Dimitrov et al., 2024 Nat Commun) + Squidpy spatial

## Steps (in-scope) — 双变量 LR 分析

1. **加载空间数据** — 从本地 tool job 目录读取上传的 h5ad
2. **自动补全** — 检测 log1p / spatial_connectivities 缺失则补
3. **LIANA bivariate** — li.mt.bivariate(resource + local metric)，直接在 spot 物理邻居上算 LR 共现
4. **per-LR-pair 统计** — mean_score / frac_active（>0.5 阈值）
5. **可视化** — 4 张图：Top-15 LR 显著性 / Top-1 LR 空间分布 / Top-N LR 分布 boxplot / Top-4 LR 2×2 空间散点
6. **结果保存** — top_lr_pairs.csv + lr_score_per_spot.csv + result_spatial.h5ad + manifest + report.html

## Out-of-scope (不做)

- ❌ **cluster-level 细胞通讯**（见 cellchat R 工具，需要预注释细胞类型）
- ❌ **CellPhoneDB / NicheNet**（其他方法，独立工具）
- ❌ **需要解卷积的 LR**（本工具是 spot-level bivariate，不依赖细胞类型注释）
- ❌ **时间序列 LR 动态**（独立工具）

## Default params (脚本默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `resourceName` | `mouseconsensus` | 小鼠 LR 数据库（人类用 consensus）|
| `localName` | `cosine` | bivariate 相似度度量 |
| `nPerms` | `100` | 显著性排列次数 |
| `topN` | `20` | 提取/绘图 top-N LR pairs |
| `nzProp` | `0.05` | L/R 最小非零比例 |

## Expected outputs

图（PNG）：plot_001_top_lr_significance / plot_002_top_lr_spatial / plot_003_lr_score_distribution / plot_004_top_lr_pairs_grid
数据表：top_lr_pairs.csv / lr_score_per_spot.csv / metadata.csv / features.csv
文件：result_spatial.h5ad / summary.json / manifest.json / report.html

## 工具间数据联动

- ✅ 上游：spatial-standard 的 result_spatial.h5ad
- ✅ 下游：LR pairs 可结合 spatial-deconv 的细胞类型丰度做细胞类型间通讯推断
