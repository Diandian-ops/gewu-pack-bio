# MOFA2 多组学因子分析工具范围卡

## Tool ID

`mo-mofa`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `MOFA2`（Bioconductor）

## Canonical workflow source

**官方文档**: MOFA2 vignette (Argelaguet et al., 2020 Mol Syst Biol)

## Steps (in-scope) - 多组学因子分解

1. **读 2-4 层组学** - CSV, rows=features, cols=samples
2. **训练 MOFA** - num_factors / convergence_mode
3. **方差解释** - heatmap per factor × view
4. **Top weights** - 每 factor top features
5. **Factor scatter** - 样本在 factor 空间投影
6. **Factor correlation** - factor 间关联
7. **Data overview** - 各层缺失/分布

## Out-of-scope (不做)

- ❌ **MOFA+ 训练后下游（imputation / timecourse）**（高级）
- ❌ **SPIRE / DIABLO**（mixOmics 多组学，独立工具）
- ❌ **Single-cell MOFA**（独立工具）
- ❌ **Cluster of cells on factors**（用户自行）
- ❌ **Annotation of factors**（用户自行做富集）

## Default params (官方默认)

| 参数 | 默认值 | 来源 |
|------|--------|------|
| `num_factors` | `15` | MOFA 默认 |
| `convergence_mode` | `slow` | MOFA 默认（最稳）|
| `scale_views` | `FALSE` | MOFA 默认（各层量纲不同）|
| `ard_weights` | `TRUE` | MOFA 默认 |

## Expected outputs

图（双格式 PNG+PDF）：variance explained heatmap / top weights / factor scatter / factor correlation / data overview
中间数据：factor_scores.csv / weight_tables.csv / mofa_model.rds

## 工具间数据联动

- ✅ 上游：用户提供 2-4 层组学矩阵（样本对齐）
- ✅ 下游：factor_scores.csv 可进 consensus-cluster / ml-classifier
