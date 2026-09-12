# DESeq2 工具范围卡

## Tool ID

`deseq2`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `DESeq2 1.46+`

## Canonical workflow source

**官方主 vignette**: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## Steps (in-scope) - Bioconductor 主 vignette 主流程

1. **Load + filter** - 读 counts + coldata；low-count filter
2. **DESeq() 拟合** - 标准三步（DESeq / results / lfcShrink(apeglm)）
3. **结果可视化** - plotPCA / plotDispEsts / plotMA / 火山图 / 热图 / sample distance / library size
4. **中间数据保存** - vst_matrix / coldata_info / dispersion_data
5. **HTML 报告 + manifest**

## Out-of-scope (不做)

- ❌ **时间序列分析**（LRT 多阶设计，独立工具）
- ❌ **batch effect 校正**（limma/removeBatchEffect，用户自行预处理）
- ❌ **ashr shrinkage**（默认 apeglm，官方推荐）
- ❌ **TXimport 上游导入**（要求用户直接提供 gene-level counts）
- ❌ **多因子交互设计**（~A*B，留给高级用户手改设计矩阵）

## Default params (官方默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `lfcShrink type` | `apeglm` | 官方推荐（Zhu 2019）|
| `alpha` | `0.05` | 官方默认 |
| `cooks cutoff` | `TRUE` | 官方默认（自动剔除离群）|
| `betaPrior` | `FALSE` | 官方默认（DESeq2 v1.28+）|

## Expected outputs

图（双格式 PNG+PDF）：plotPCA / plotDispEsts / plotMA / 火山图 / 热图 / sample distance / library size
中间数据：vst_matrix.csv / coldata_info.csv / dispersion_data.csv / deg_results.csv

## 工具间数据联动

- ✅ 上游：用户提供 raw counts（gene × sample）+ coldata
- ✅ 下游：deg_results.csv 可被 go-kegg / gsea / lasso-cox 工具导入
