# deseq2

> **科学子问题**: 在给定实验设计和样本质量前提下，处理条件是否导致可靠的基因表达变化？

## 职责

使用 DESeq2 对 counts 矩阵执行差异表达基因检测，包含：
- DESeqDataSet 构建（design formula）
- DESeq() 标准流程（estimation → dispersion → Wald test）
- lfcShrink（apeglm）收缩效应量
- 输出 DEG 结果表、normalized counts、VST 矩阵
- MA plot、PCA plot、dispersion plot

## DAG 位置

```
expression-qc → deseq2 → deg-standardizer
```

上游 `expression-qc` 的 `validated_counts.csv` 和 `sample-metadata-validator` 的 `validated_metadata.csv` 作为输入。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| counts | file | ✓ | — | counts 矩阵 (CSV: gene × sample) |
| coldata | file | ✓ | — | 样本信息表 (CSV: sample, condition, ...) |
| design | text | — | `condition` | 设计公式变量名 |
| contrast_ref | text | — | `Control` | 对照水平 |
| contrast_treat | text | — | `Treatment` | 处理水平 |
| padj_cutoff | number | — | `0.05` | padj 阈值 |
| lfc_cutoff | number | — | `1` | |log2FC| 阈值 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `deg_results.csv` | table | DEG 结果表 (gene, baseMean, log2FoldChange, lfcSE, pvalue, padj) |
| `normalized_counts.csv` | table | DESeq2 normalized counts |
| `vst_matrix.csv` | table | VST transformed matrix |
| `report.html` | report | HTML 报告（含 MA plot、PCA、dispersion） |
| `manifest.json` | manifest | 工具 manifest |

## R 包依赖

- DESeq2 (>=1.40.0)
- ggplot2 (>=3.4.0)
- ggrepel (>=0.9.0)
- pheatmap (>=1.0.12)
- apeglm (>=1.20.0)
- jsonlite (>=1.8.0)
- base64enc (>=0.1-3)

## 平台

- ✅ macos-arm64
- ⏳ macos-x64 / windows-x64 / linux-x64 (planned)

## 注意

此插件使用 legacy `manifest.json` 格式（非 `artifact-manifest.json`）。
升级到完整 Artifact lineage 格式已跟踪为 HRT-1 任务。
