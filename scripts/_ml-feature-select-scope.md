# ML 特征选择工具范围卡

## Tool ID

`ml-feature-select`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `caret` / `Boruta` / `mRMRe` / `glmnet` / `e1071`

## Canonical workflow source

**官方文档**: caret vignette + Boruta (Kursa & Rudnicki 2010) + mRMR (Peng 2005)

## Steps (in-scope) - 5 种方法对比

1. **Variance filter** - top variance 粗筛
2. **Univariate Wilcoxon / t-test** - 单变量统计
3. **mRMR** - max relevance / min redundancy
4. **Boruta** - wrapper, RF-based
5. **RFE** - recursive feature elimination, LR-based
6. **共识** - Venn consensus + method correlation heatmap + volcano (univariate)

## Out-of-scope (不做)

- ❌ **LASSO 特征选择**（见 lasso-cox 工具）
- ❌ **Deep learning feature learning**（本工具不包含）
- ❌ **SHAP 特征重要性**（见 shap-explain，是后验解释而非选择）
- ❌ **Stability selection**（高级，留给用户）
- ❌ **特征工程 / 交互项构造**（用户自行预处理）

## Default params (官方默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `top_variance_n` | `2000` | 粗筛常用 |
| `univariate_p` | `0.05` | 常规 |
| `boruta_p` | `0.01` | Boruta 默认 |
| `rfe_sizes` | `c(10,20,50,100)` | caret 常用 |
| `consensus_top_n` | `50` | 下游可用 |

## Expected outputs

图（双格式 PNG+PDF）：Venn 共识图 / method correlation heatmap / volcano plot
中间数据：5 份 ranking table + consensus_top_genes.csv + summary

## 工具间数据联动

- ✅ 上游：用户提供 expression + label
- ✅ 下游：consensus_top_genes.csv 可导入 ml-classifier 做建模；或导入 lasso-cox
