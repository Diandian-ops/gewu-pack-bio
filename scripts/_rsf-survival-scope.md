# Random Survival Forest (RSF) 工具范围卡

## Tool ID

`rsf-survival`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `randomForestSRC` / `timeROC` / `pec`

## Canonical workflow source

**官方文档**: randomForestSRC vignette + Ishwaran & Kogalur 2008

## Steps (in-scope) - 非线性预后建模

1. **C-index** - train / OOB / test
2. **Time-dependent ROC** - 用户指定时间点
3. **VIMP** - Variable Importance ranking
4. **Partial dependence** - top N features
5. **Risk score 三联图** - KM + scatter + heatmap
6. **Trained model** - .rds

## Out-of-scope (不做)

- ❌ **线性 Cox / LASSO Cox**（见 lasso-cox 工具）
- ❌ **DeepSurv / Deep learning survival**（独立工具）
- ❌ **Competing risks RSF**（独立工具）
- ❌ **Nomogram / Calibration**（线性工具更合适）
- ❌ **time-varying coefficients**（高级）

## Default params (官方默认)

| 参数 | 默认值 | 来源 |
|------|--------|------|
| `ntree` | `500` | RF 标准 |
| `mtry` | `sqrt(p)` | RSF 默认 |
| `nodesize` | `3` | RSF 默认（survival 比 class 小）|
| `splitrule` | `logrank` | RSF 默认 |
| `importance` | `TRUE` | VIMP 开启 |

## Expected outputs

图（双格式 PNG+PDF）：C-index 表 / timeROC / VIMP barplot / partial dependence / risk 三联图
中间数据：vimp_ranking.csv / risk_scores.csv / rsf_model.rds

## 工具间数据联动

- ✅ 上游：可导入 lasso-cox 的 signature 基因 / ml-feature-select 结果
- ✅ 下游：rsf_model.rds 可进 SHAP 或外部验证
