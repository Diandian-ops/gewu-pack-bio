# Consensus Clustering 工具范围卡

## Tool ID

`consensus-cluster`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `ConsensusClusterPlus` / `NMF` / `mclust` / `cluster`

## Canonical workflow source

**官方文档**: ConsensusClusterPlus vignette (Monti et al., 2003)

## Steps (in-scope) - 无监督亚型发现

1. **ConsensusClusterPlus** - k = 2 to max_k，resampling + voting consensus
2. **Auto k selection** - delta area elbow + CDF plateau
3. **可视化** - CDF / delta area / consensus matrix panel
4. **Silhouette score** - per k
5. **PCA projection** - by subtype
6. **ARI validation** - 如提供 known_labels

## Out-of-scope (不做)

- ❌ **有监督分类**（见 ml-classifier）
- ❌ **单细胞聚类**（见 seurat-standard）
- ❌ **NMF 分解**（独立工具，本工具仅借用其距离）
- ❌ **多组学联合聚类**（见 mo-snf / mo-mofa）
- ❌ **Spatial domain detection**（见 spatial 工具）

## Default params (官方默认)

| 参数 | 默认值 | 来源 |
|------|--------|------|
| `max_k` | `8` | CCP 常用上限 |
| `reps` | `50` | CCP 默认 1000，BioF3 降到 50 平衡速度 |
| `pItem` | `0.8` | CCP 默认 |
| `clusterAlg` | `hc` | CCP 默认 |
| `distance` | `pearson` | CCP 默认 |

## Expected outputs

图（双格式 PNG+PDF）：CDF / delta area / consensus matrix / silhouette / PCA
中间数据：consensus_matrix.csv / subtype_assignment.csv / metrics_per_k.csv

## 工具间数据联动

- ✅ 上游：用户提供 expression 矩阵
- ✅ 下游：subtype_assignment.csv 可进 ml-classifier 当 label；或 km-survival 比较亚型预后
