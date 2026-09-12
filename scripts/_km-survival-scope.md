# KM 生存分析工具范围卡

## Tool ID

`km-survival`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `survival 3.x` / `survminer 0.4+` / `timeROC 0.4+`

## Canonical workflow source

**官方文档**: survminer vignette + timeROC vignette

## Steps (in-scope) - 单基因预后探索 4 件套

1. **读 expression + clinical**
2. **cutoff 分组** - median / optimal (maxstat) / tertile
3. **KM survival curve + log-rank**
4. **univariate Cox + forest plot**
5. **time-dependent ROC** - 1/3/5 year
6. **expression distribution**

## Out-of-scope (不做)

- ❌ **多基因 signature**（见 lasso-cox 工具）
- ❌ **多变量 Cox / nomogram**（见 rsf-survival 或独立工具）
- ❌ **竞争风险模型**（finegray，独立工具）
- ❌ **landmark analysis**（时变协变量，高级）
- ❌ **batch 多基因批量筛查**（性能/多重检验，独立工具）

## Default params (官方默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `cutoff_method` | `median` | 临床最常用（optimal 易过拟合）|
| `time_col` | `OS.time` | TCGA 标准 |
| `event_col` | `OS` | TCGA 标准（1=死）|
| `roc_times` | `1,3,5` year | 临床常规 |

## Expected outputs

图（双格式 PNG+PDF）：KM curve / forest plot / timeROC / expression distribution
中间数据：km_results.csv / cox_results.csv / roc_results.csv

## 工具间数据联动

- ✅ 上游：用户提供 expression + clinical
- ✅ 下游：单基因结果可整合进 lasso-cox 候选池
