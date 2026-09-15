# sample-metadata-validator

> **科学子问题**: counts 矩阵和 metadata 是否能够支持指定的实验设计与 contrast？

## 职责

在进入表达 QC 和差异分析之前，验证输入数据的结构和实验设计是否可靠：

- counts 矩阵是否有样本列
- metadata 是否有样本标识列
- counts 样本名与 metadata 样本名是否一致
- 是否有重复样本
- 是否有缺失分组
- reference/treatment 是否真实存在
- 每组样本数是否足够
- 是否只有单个样本
- 是否出现完全混淆的设计
- 对齐时是否丢失样本
- 是否可以继续建模

## DAG 位置

```
输入文件探测 → sample-metadata-validator → expression-qc → deseq2 → ...
```

状态为 `valid` 或 `warning` 的 Artifact 可被下游 `expression-qc` 和 `deseq2` 消费。
状态为 `blocked` 的 Artifact **不能**进入下游分析。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| counts | file | ✓ | — | 基因表达 counts 矩阵 (CSV/TSV: gene × sample) |
| metadata | file | ✓ | — | 样本信息表 (CSV/TSV) |
| sample_column | text | — | `sample` | metadata 中标识样本的列名 |
| design_column | text | — | `condition` | metadata 中用于对比的分组列名 |
| reference_group | text | — | `Control` | 对照组 |
| treatment_group | text | — | `Treatment` | 实验组 |
| min_replicates | integer | — | `2` | 每组最少样本数 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `validated_metadata.csv` | table | 对齐后的 metadata |
| `sample_alignment.csv` | table | 样本对齐详情 |
| `design_summary.json` | json | 设计摘要 |
| `validation_diagnostics.json` | json | 校验诊断 |
| `report.html` | report | HTML 校验报告 |
| `artifact-manifest.json` | manifest | Artifact 清单 |

## 状态

- `valid` — 所有检查通过，可以继续建模
- `warning` — 存在非致命问题（如部分样本不匹配、低重复数），需要用户确认
- `blocked` — 存在致命问题，不能进入下游分析

## 失败模式

| code | 说明 |
|------|------|
| `missing_required_file_inputs` | counts 或 metadata 文件未提供 |
| `sample_id_column_missing` | metadata 中找不到 sample 列 |
| `sample_mismatch` | 样本名严重不一致 |
| `duplicated_samples` | 存在重复样本 ID |
| `missing_design_column` | metadata 中找不到 design 列 |
| `contrast_level_missing` | reference/treatment 分组不存在 |
| `insufficient_replicates` | 重复数不足 |
| `invalid_design` | 设计完全混淆 |
| `runtime_not_ready` | R 运行时或包未就绪 |
| `platform_not_supported` | 平台不支持 |

## 运行

```bash
Rscript scripts/sample-metadata-validator.R <job_dir>
```

`<job_dir>` 需包含：
- `params.json` — 参数
- `counts.csv` (或 `counts.tsv`) — counts 矩阵
- `metadata.csv` (或 `metadata.tsv` / `coldata.csv`) — metadata

## R 包依赖

- jsonlite
- digest

## 平台

- ✅ macos-arm64
- ⏳ macos-x64 (planned)
- ⏳ windows-x64 (planned)
- ⏳ linux-x64 (planned)

## Demo 数据

`demo-data/` 目录包含：
- `demo-counts.csv` + `demo-metadata.csv` — valid fixture
- `demo-counts-single-rep.csv` + `demo-metadata-single-rep.csv` — insufficient replicates fixture
- `demo-counts-mismatch.csv` + `demo-metadata-mismatch.csv` — sample mismatch fixture
- `demo-metadata-no-sample-col.csv` — sample column missing fixture
