# pca-explorer

> **科学子问题**: 样本在降维空间中的分布是否与实验设计分组一致？是否存在批次效应或离群样本？

## 职责

从表达矩阵执行 PCA 降维，绘制 PC1/PC2 散点图（按分组着色），输出方差贡献表和 PCA 坐标表。

## DAG 位置

```
deg-standardizer → pca-explorer
```

也可独立使用 expression-qc 的 `validated_counts.csv` 作为输入。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| expr | file | ✓ | — | 表达矩阵 (CSV: gene × sample) |
| groups | file | — | — | 样本分组表 (CSV: sample, condition) |
| sample_col | text | — | `sample` | 样本列名 |
| group_col | text | — | `condition` | 分组列名 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `pca_plot.png` | image | PCA 散点图 |
| `variance_table.csv` | table | 各 PC 方差贡献 |
| `report.html` | report | HTML 报告 |
| `manifest.json` | manifest | 工具 manifest |

## 运行时

- Python >= 3.10
- pandas, scikit-learn, matplotlib

## 平台

- ✅ macos-arm64 / macos-x64 / windows-x64

## 注意

此插件使用 Python（非 R），使用 legacy `manifest.json` 格式。
升级到完整 Artifact lineage 格式已跟踪为 HRT-1 任务。
