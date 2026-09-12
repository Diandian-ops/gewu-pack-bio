# LASSO Cox 预后 Signature 构建工具范围卡

## Tool ID

`lasso-cox`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `glmnet 4.1+` / `survival 3.x` / `survminer 0.4+`

## Canonical workflow source

**官方文档**: glmnet vignette + Simon et al., 2011 JSS

## Steps (in-scope) - 预后 Signature 构建主流程

1. **读 expression + clinical** - 候选基因池（可选）
2. **train/test 分割** - train_ratio（默认 0.7），固定 seed
3. **LASSO CV** - cv.glmnet 选 lambda（alpha=1）
4. **Signature 构建** - 非零系数基因 -> risk score
5. **三联图验证** - KM (train+test) / ROC / Risk score 分布
6. **HTML 报告 + manifest**

## Out-of-scope (不做)

- ❌ **Elastic Net (alpha<1) / Ridge**（默认 LASSO，用户可手改）
- ❌ **Stepwise Cox**（独立工具）
- ❌ **Nomogram / Calibration**（见 rsf-survival 或独立工具）
- ❌ **外部验证队列**（用户自行跑 test 集）
- ❌ **time-dependent AUC 一致性指数**（见 rsf-survival）

## Default params (官方默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `train_ratio` | `0.7` | 临床常规 |
| `seed` | `42` | 可复现 |
| `alpha` | `1` (LASSO) | 官方默认 LASSO |
| `nfolds` | `10` | cv.glmnet 默认 |

## Expected outputs

图（双格式 PNG+PDF）：LASSO CV 曲线 / 系数路径 / KM (train+test) / ROC / Risk score 三联图
中间数据：signature_genes.csv / risk_scores.csv / lasso_model.rds

## 工具间数据联动

- ✅ 上游：可导入 deseq2 deg_results.csv 作为候选池；km-survival 单基因结果
- ✅ 下游：signature 基因可进 rsf-survival 做非线性建模；risk score 可进 shap-explain
