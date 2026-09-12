# maftools 工具范围卡（步骤 0 落地）

> tool-dev "步骤 0：工具范围与官方文档" 第五个落地实例
> 第三个使用 SCI 视觉风格规范的工具（base R 系，全部用 `save_grid_biof3()` 双格式）

## Tool ID

`maftools-mutation`

## Version targeted

- `maftools 2.24.0`（最新 Bioc release）
- `BSgenome.Hsapiens.UCSC.hg19 1.4.3` / `BSgenome.Hsapiens.UCSC.hg38 1.4.5`
- `NMF 0.28`（mutational signature 用）
- `R.utils 2.13.0`

## Canonical workflow source

**官方 vignette**: https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html

**Reference**: Mayakonda A et al. (2018). Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Research* 28:1747-1756.

## Steps (in-scope) — Tutorial 主流程 7 步

1. **Load** — `read.maf(maf, clinicalData)` 读 MAF 文件 + 临床表（可选）
2. **MAF 概览** — `plotmafSummary()` 输出 6 panel 概览（样本变异分布 + 突变类型 + SNV class + Top 10 基因）
3. **Oncoplot** — `oncoplot(top=20)` 突变景观图（细胞学注释 + Pathway 注释 + 临床 annotation track）
4. **TiTv 比例** — `titv()` + `plotTiTv()` 转换/颠换比例
5. **Lollipop 图** — 对 top 5-10 基因画 protein domain + mutation 位置（lollipop）
6. **Somatic Interactions** — `somaticInteractions()` 互斥/共现矩阵 + p value
7. **Mutational Signatures** — `trinucleotideMatrix()` + `extractSignatures()` 提取 mutation signatures（COSMIC v3 比较）

## Out-of-scope (不做)

- ❌ **TCGA cohort comparison** — `tcgaCompare()` 需要联网下载，离线不稳
- ❌ **OncodriveCLUST** — 高级 driver gene 识别，依赖 NMF 多次拟合，慢
- ❌ **Drug-Gene Interactions** — `drugInteractions()` 依赖 DGIdb 联网
- ❌ **PathwayAnalyzer** — `OncogenicPathways()` 输出有限，留给 ORA 工具
- ❌ **MAF 转 ICGC / 自定义格式转换** — 数据格式问题，留给用户预处理
- ❌ **survival analysis on mutation** — `mafSurvival()` 留给 KM 工具

> `_data.ts` description: "maftools v2 突变景观分析: MAF 概览 → Oncoplot → TiTv → Lollipop → 互斥共现 → Mutational Signature。SCI 风格双格式输出（PNG + PDF）。不含: TCGA cohort 对比、Drug-Gene 互作、生存分析（用 KM 工具）。"

## Default params (按官方 + BioF3 SCI 默认)

| 参数 | 默认 | 说明 |
|------|------|------|
| `top_n_oncoplot` | 20 | Oncoplot 显示 top N 基因 |
| `lollipop_top_n` | 5 | Lollipop 图为 top N 突变基因画 |
| `genome` | `hg19` | 与 LAML demo 一致；用户可切换 hg38 |
| `n_signatures` | 3 | mutational signature 数量（小数据集 3 个稳定，大数据集 5-6） |
| `min_mut_for_signature` | 5 | extractSignatures 需要每样本最少 5 个突变 |

## Expected outputs (8 张图 × PNG + PDF + 7 中间数据)

### 图

| # | 文件名 | 函数 | 说明 |
|---|--------|------|------|
| 1 | `maf_summary` | `plotmafSummary()` | 6 panel 总览 |
| 2 | `oncoplot` | `oncoplot()` | 突变景观（Top N 基因）|
| 3 | `titv` | `plotTiTv()` | 转换/颠换比例 |
| 4 | `lollipop_top1` | `lollipopPlot()` | top 1 基因 lollipop |
| 5 | `lollipop_top2` | `lollipopPlot()` | top 2 基因 lollipop |
| 6 | `somatic_interactions` | `somaticInteractions()` | 互斥/共现矩阵 |
| 7 | `mutation_signature` | `plotSignatures()` | 提取的突变 signature 谱 |
| 8 | `signature_cosmic_compare` | `compareSignatures` | 与 COSMIC v3 signatures 对比 |

### 中间数据

| 文件 | 内容 |
|------|------|
| `gene_summary.csv` | 每个基因的突变数 / 样本数 / 突变类型分布 |
| `sample_summary.csv` | 每个样本的总变异数 + 突变类型分布 |
| `clinical_data.csv` | 临床数据透传（如有）|
| `titv_data.csv` | TiTv 矩阵 |
| `interaction_matrix.csv` | somatic interactions p value 矩阵 |
| `signature_contributions.csv` | 每个样本各 signature 的贡献度 |
| `maf.rds` | 完整 maftools MAF 对象 |

## Demo 数据

**用本地 R runtime 中 maftools 自带的 TCGA-LAML**：
- `system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')` 64KB
- `system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')` 4.5KB

可复制到本地插件的 `demo-data/` 目录。文件极小（< 70KB），是真实 TCGA 数据，193 AML 样本 + FAB 分类 + 生存数据。

## 工具间数据联动

- ❌ 上游不需要（MAF 是入口）
- ⏳ 下游可被 KM 工具 / GO/KEGG 工具消费（用户筛 top 突变基因导出）

## 已读官方文档（T1）

- [x] maftools Bioconductor vignette 主流程
- [x] customizing oncoplots vignette
- [x] BioF3 自己的教程 `docs/genomics/module04` 突变分析
- [x] 容器版本验证：maftools 2.24.0

## 已知风险

1. **maftools 全是 base R 绘图**：所有图必须用 `save_grid_biof3()` 双格式包装，**不能套 theme_biof3()**（base R 不接 ggplot2 主题）。这是 `_biof3-theme.R` helper 的另一个验证场景
2. **mutational signature 对小数据集不稳**：< 100 样本时 `extractSignatures()` 可能拒绝 / 输出空。脚本里加 `tryCatch`，失败时跳过 signature 部分
3. **lollipop 基因要有蛋白域注释**：`lollipopPlot()` 内部查 protein domains，对一些非编码 / 罕见基因会失败。脚本里 `tryCatch` 兜住
4. **OOM**：MAF 数据集本身很小，`oncoplot` / `titv` 都不耗内存；mutational signature 是 NMF，但 LAML 193 样本 < 1 分钟可跑完

## 下一步

按 14 步流程：
- 步骤 0: 范围卡 ✅
- 步骤 1: ID = `maftools-mutation` ✅
- 步骤 2: 更新 `_data.ts`
- 步骤 3: `routes/tools.js` 已注册（`maftools.R`）
- 步骤 4: TOOL_COLORS / REPORT_URLS 已注册
- 步骤 5: 编写 `maftools.R`
- 步骤 6-14: 按流程
