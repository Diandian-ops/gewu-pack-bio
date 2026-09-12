# ML 分类器工具范围卡

## Tool ID

`ml-classifier`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `caret` / `glmnet` / `e1071` / `randomForest` / `xgboost` / `pROC`

## Canonical workflow source

**官方文档**: caret vignette + 各算法主文档

## Steps (in-scope) - 4 算法对比

1. **4 算法** - LR / SVM / RF / XGBoost
2. **5-fold CV** - ROC / PR curves per algorithm
3. **Calibration plot** - 校准度
4. **Confusion matrices** - 混淆矩阵
5. **Feature importance comparison** - 跨算法
6. **Learning curve** - 学习曲线
7. **Trained model** - .rds 保存供下游 SHAP

## Out-of-scope (不做)

- ❌ **特征选择**（见 ml-feature-select）
- ❌ **SHAP 解释**（见 shap-explain）
- ❌ **深度学习分类器**（本工具不包含）
- ❌ **Bayesian optimization 超参**（留高级用户）
- ❌ **Ensemble stacking**（独立工具）
- ❌ **生存预测**（见 rsf-survival）

## Default params (官方默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `cv_folds` | `5` | caret 常用 |
| `metric` | `ROC` | 二分类标准 |
| `preProcess` | `center-scale` | caret 标准（中心化+方差+Z）|

## Expected outputs

图（双格式 PNG+PDF）：ROC / PR / calibration / confusion matrix / feature importance / learning curve
中间数据：metrics_summary.csv / trained_model.rds / predictions.csv

## 工具间数据联动

- ✅ 上游：可导入 ml-feature-select 的 consensus_top_genes.csv
- ✅ 下游：trained_model.rds 可被 shap-explain 读取
