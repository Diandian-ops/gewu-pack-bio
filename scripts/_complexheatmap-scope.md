# ComplexHeatmap 工具范围卡（步骤 0 落地）

> tool-dev "步骤 0：工具范围与官方文档" 第三个落地实例
> 作为 BioF3 SCI 视觉风格规范的首个工具

## Tool ID

`complexheatmap`

## Version targeted

- BioF3 R 容器: R 4.5.1
- 关键包：
  - `ComplexHeatmap 2.24.1`
  - `circlize 0.4.17`
  - `RColorBrewer 1.1.3`
  - `viridis 0.6.5`
  - `dendextend 1.19.1`

## Canonical workflow source

**官方 vignette**: ComplexHeatmap Complete Reference (https://jokergoo.github.io/ComplexHeatmap-reference/book/)

**Reference**: Gu Z. (2022). Complex heatmap visualization. *iMeta* 1(3):e43.

## Steps (in-scope) — 标准热图工作流

1. **Load + 检查输入** — CSV 数值矩阵 + 可选 sample annotation
2. **Z-score 标准化** — `t(scale(t(mat)))`，可选项（用户可关闭）
3. **聚类决策** — 行聚类 / 列聚类的开关 + 距离方法 (euclidean / pearson / spearman) + 链接方法 (ward.D2 / complete / average)
4. **基础热图** — 单层 Heatmap + 颜色渐变（蓝-白-红 或 viridis）+ 图例
5. **顶部 annotation** — 用 sample annotation CSV 画 column annotation bar（categorical → discrete colors / numeric → gradient）
6. **行/列分组** — 用 `row_split` / `column_split` + `cluster_within_group` 让相似样本聚到一起
7. **多组热图叠加（可选）** — 第二个 annotation matrix（如 mutation status）作为 `+ Heatmap()` 拼接
8. **导出** — 主热图 + 一个简化版 + 一个用户友好版（带行名 / 标记 top genes）

## Out-of-scope (不做)

- ❌ **OncoPrint** — 突变分析专用，留给 maftools-mutation 工具
- ❌ **UpSet plot** — `UpSet()` 函数，独立工具更合适
- ❌ **Density / scatter heatmap** — `densityHeatmap()`, 边角特性
- ❌ **3D / interactive heatmap** — 离 sci publication 标准远
- ❌ **跨大数据集的 modulePreservation 风格热图** — 太特化
- ❌ **自定义色卡复杂注释** — 用户在 RunModal 改

> `_data.ts` description: "ComplexHeatmap 高级热图：z-score 标准化 + 行/列聚类 + sample annotation + 行分组。SCI 风格双格式输出（PNG + PDF）。不含: OncoPrint / UpSet / 突变可视化（用其他工具）。"

## Default params (按官方推荐 + BioF3 SCI 默认)

| 参数 | 默认 | 来源 / 说明 |
|------|------|-------|
| `z_score` | TRUE | 行 z-score 标准化（最常见的 use case） |
| `cluster_rows` | TRUE | |
| `cluster_columns` | TRUE | |
| `clustering_distance` | `pearson` | 表达数据用 pearson 比 euclidean 更合理 |
| `clustering_method` | `ward.D2` | 官方推荐（也是 pheatmap / heatmap.2 默认） |
| `show_row_names` | "auto" | 行数 ≤ 50 时 show，否则隐藏（多基因时不可读） |
| `show_column_names` | TRUE | 样本通常 ≤ 200 |
| `top_n_label` | 20 | 行数 > 50 时仅在右侧标记 top 20 高方差基因 |
| `color_scheme` | `div` | div = 蓝-白-红（z-score 用），seq = viridis（原始值用） |
| `color_range_quantile` | 0.99 | 用 [-quantile(|x|, 0.99), +quantile(|x|, 0.99)] 截断颜色范围，避免极端值压扁主梯度 |

## Expected outputs (5 张图 + 4 中间数据)

### 图

| # | 文件名 | 说明 |
|---|--------|------|
| 1 | `heatmap_main.{png,pdf}` | 主热图：z-score + 聚类 + annotation + 标记 top genes |
| 2 | `heatmap_unscaled.{png,pdf}` | 未标准化热图（用户自行比较） |
| 3 | `dendrogram_rows.{png,pdf}` | 行（基因）聚类树独立图（debug 用） |
| 4 | `dendrogram_columns.{png,pdf}` | 列（样本）聚类树独立图 |
| 5 | `correlation_heatmap.{png,pdf}` | 样本间 correlation matrix（QC 用，验证聚类合理性） |

### 中间数据

| 文件 | 内容 |
|------|------|
| `matrix_zscore.csv` | z-score 后的矩阵（用户可重画 / 自定义热图）|
| `matrix_clustered.csv` | 按聚类顺序 reorder 后的矩阵 |
| `column_annotation.csv` | 处理后的 sample annotation（colors 已 binding） |
| `clustering_summary.txt` | 聚类参数 + 行/列树高度统计 |

## Demo 数据

**用 TCGA-BRCA 子集**（与 WGCNA 共用上游）：
- 矩阵: 100 高变基因 × 80 样本（更小子集让热图清晰可读）
- Annotation: PAM50 子型 + ER / HER2 状态

复用现有的 `wgcna_expr.csv`，但子选 100 基因 + 80 样本 + 加 annotation：

```r
# 在 build-complexheatmap-demo.R 里:
expr  <- read.csv(".../wgcna/expression.csv")
trts  <- read.csv(".../wgcna/traits.csv")

# 选 80 样本 (按 PAM50 平衡, 每类 16 个)
selected_samples <- ...
# 选 100 高变基因
top100_var <- ...

mat_demo <- expr[top100_var, selected_samples]
anno_demo <- data.frame(
  Sample = selected_samples,
  PAM50 = ifelse(trts$LumA == 1, "LumA",
          ifelse(trts$LumB == 1, "LumB",
          ifelse(trts$Her2 == 1, "Her2",
          ifelse(trts$Basal == 1, "Basal", "Normal")))),
  ER = ifelse(trts$ER == 1, "Positive", "Negative"),
  HER2 = ifelse(trts$HER2 == 1, "Positive", "Negative")
)
```

## 工具间数据联动

- ✅ **从 DESeq2 导入**：可读取 DESeq2 的 `vst_matrix.csv`（top N DEGs 由用户在前端选）
- ✅ **从 WGCNA 导入**：可读取 WGCNA 的 module 基因子集（特定 module 颜色的基因画热图）
- ❌ 上游不需要其他工具

## 已读官方文档（T1）

- [x] ComplexHeatmap Complete Reference 第 1-3 章（基础 Heatmap + Annotation + 颜色控制）
- [x] BioF3 自己的教程 `docs/integration/module03.md`（多组学整合用 ComplexHeatmap）口径
- [x] 容器版本验证：ComplexHeatmap 2.24.1 / circlize 0.4.17

## 已知风险

1. **大矩阵 OOM**：> 2000 行 × 200 列时 hierarchical clustering 慢（O(n²)）。脚本里加 dim check：> 5000 × 500 时拒绝
2. **Pearson 距离对零方差行报错**：z-score 后零方差行会变 NaN → `ward.D2` 报错。脚本里先过滤零方差行，warn 日志
3. **Annotation 列太多视觉爆炸**：> 6 个 annotation 列时图变窄。脚本里 hard-cap 到 6（取前 6），warn 日志
4. **PDF 输出对超大矩阵慢**：1000 行的 PDF 矢量化要几秒。可接受，不优化

## 下一步

按 14 步流程：
- 步骤 0: 本范围卡 ✅
- 步骤 1: ID = `complexheatmap` ✅
- 步骤 2: 加 `_data.ts` 工具定义
- 步骤 3: 加 `routes/tools.js` 注册
- 步骤 4: 加 `index.tsx` TOOL_COLORS / REPORT_URLS
- 步骤 5: 编写 `complexheatmap.R`（用 `_biof3-theme.R`，做第一个 sci 风格落地）
- 步骤 6-14: 按流程走
