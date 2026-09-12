# SHAP 可解释性工具范围卡

## Tool ID

`shap-explain`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `iml` (kernel SHAP) / 各模型包

## Canonical workflow source

**官方文档**: iml vignette + Lundberg & Lee 2017 (SHAP)

## Steps (in-scope)

1. **读上游 ml-classifier model** - upstream_jobid
2. **自动模型类型检测** - LR / SVM / RF / XGBoost
3. **Global SHAP importance** - 全局特征重要性
4. **Beeswarm** - 全样本 SHAP 分布
5. **Force plot** - 单样本解释
6. **Dependence plots** - top N 特征

## Out-of-scope (不做)

- ❌ **特征选择**（见 ml-feature-select，SHAP 是解释非选择）
- ❌ **TreeSHAP 原生**（用 iml kernel SHAP，通用但较慢）
- ❌ **LIME**（独立工具）
- ❌ **Counterfactual examples**（高级）
- ❌ **SHAP for survival models**（rsf 暂不支持 iml）

## Default params (官方默认)

| 参数 | 默认值 | 来源 |
|------|--------|------|
| `top_n` | `20` | 可视化常用 |
| `n_explain_samples` | `30` | kernel SHAP 采样数 |
| `sample_index` | `1` | force plot 默认样本 |

## Expected outputs

图（双格式 PNG+PDF）：global importance barplot / beeswarm / force plot / dependence plots
中间数据：shap_values.csv / shap_summary.csv

## 工具间数据联动

- ✅ 上游：必须接 ml-classifier 的 trained_model.rds + features/labels
- ❌ 下游：解释性结果，无直接下游
