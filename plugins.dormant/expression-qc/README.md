# expression-qc

> **科学子问题**: 表达/counts 矩阵是否具有足够的数据质量，可以进入差异表达建模？

## 职责

对 counts 矩阵执行质量门检查，**不执行差异表达**：

- gene × sample 方向检查
- 非数值列检测
- NA / Inf 检测
- 负 counts 检测
- 小数 counts 检测（可能是 TPM/FPKM）
- 重复 gene ID 检测
- 全零基因检测
- 低表达基因过滤（规则写入参数和报告）
- Library size 统计与均匀性检查
- 样本相关性矩阵 + 热图
- PCA 坐标 + 图
- 样本离群检测（warning，**不自动删除样本**）
- metadata 对齐状态检查

## DAG 位置

```
sample-metadata-validator → expression-qc → deseq2
```

## 关键原则

- **不根据预期结论删除样本**
- **不自动删除用户样本** — 离群样本只输出 warning，要求确认
- **过滤规则必须写入参数和报告**
- 只做 QC 和质量门，不执行差异表达

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| counts | file | ✓ | — | 表达/counts 矩阵 (CSV/TSV: gene × sample) |
| metadata | file | — | — | 样本信息表（可选） |
| validated_metadata | text | — | — | 上游 validated_metadata Artifact 引用 |
| min_count | number | — | `10` | 最低 count 阈值 |
| min_samples | integer | — | `2` | 最少表达样本数 |
| transformation | text | — | `none` | 转换模式 (none / log2) |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `qc_summary.json` | json | QC 摘要（含所有检查项） |
| `gene_filter_statistics.csv` | table | 基因过滤统计 |
| `sample_library_size.csv` | table | 每样本 library size |
| `sample_correlation.csv` | table | 样本相关性矩阵 |
| `sample_correlation.png` | image | 相关性热图 |
| `pca_coordinates.csv` | table | PCA 坐标 |
| `pca.png` | image | PCA 散点图 |
| `validated_counts.csv` | table | 过滤后 counts |
| `report.html` | report | HTML QC 报告 |
| `artifact-manifest.json` | manifest | Artifact 清单 |

## 状态

- `valid` — 质量门通过
- `warning` — 存在潜在问题（离群样本、library size 差异等），需用户确认
- `blocked` — 数据质量不足以继续

## R 包依赖

- jsonlite (必需)
- digest (必需)
- pheatmap (可选，用于热图；缺失时 fallback 到 base R heatmap)

## 平台

- ✅ macos-arm64
- ⏳ macos-x64 / windows-x64 / linux-x64 (planned)
