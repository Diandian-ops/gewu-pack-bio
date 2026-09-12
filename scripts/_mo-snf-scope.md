# SNF 相似性网络融合工具范围卡

## Tool ID

`mo-snf`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `SNFtool` / `mclust`

## Canonical workflow source

**官方文档**: SNFtool vignette (Wang et al., 2014 Sci Rep)

## Steps (in-scope) - 多组学网络融合

1. **读 2-3 层组学** - CSV, rows=samples, cols=features
2. **各层相似性矩阵** - SNFtool::affinityMatrix
3. **网络融合** - SNF(K, alpha, T_iter)
4. **亚型发现** - spectral clustering on fused network
5. **可视化** - fused network heatmap / eigen-gap / NMI+ARI 验证

## Out-of-scope (不做)

- ❌ **MOFA 因子分解**（见 mo-mofa 工具，SNF 是聚类导向）
- ❌ **Multi-view NMF**（独立工具）
- ❌ **iCluster**（独立工具）
- ❌ **Single-sample network**（用户自行）
- ❌ **Spatial multi-omics**（见 spatial 工具）

## Default params (官方默认)

| 参数 | 默认值 | 来源 |
|------|--------|------|
| `K` | `20` | SNF 默认（邻居数）|
| `alpha` | `0.5` | SNF 默认（超参）|
| `T_iter` | `20` | SNF 默认（迭代数）|
| `max_k` | `8` | 聚类上限 |

## Expected outputs

图（双格式 PNG+PDF）：fused network heatmap / eigen-gap / subtype PCA / validation metrics
中间数据：fused_network.csv / subtype_assignment.csv / metrics_per_k.csv

## 工具间数据联动

- ✅ 上游：用户提供 2-3 层组学矩阵（样本对齐）
- ✅ 下游：subtype_assignment.csv 可进 ml-classifier / km-survival
