# GSEA 基因集富集分析工具范围卡

## Tool ID

`gsea`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `clusterProfiler 4.18+` / `enrichplot 1.30+` / `msigdbr 26+` / `ggplot2 4.0+`

## Canonical workflow source

**官方文档**: clusterProfiler GSEA vignette + Subramanian et al., 2005 PNAS

## Steps (in-scope)

1. **读 ranked gene list** - gene + log2FC，或 DESeq2 result
2. **ID 转换** - SYMBOL -> ENTREZID
3. **gseGO BP** - GO 生物过程
4. **gseKEGG** - 通路
5. **Hallmark** - msigdbr H 集合
6. **可视化** - dotplot / waterfall / running score / ridge / emapplot / cnetplot / treeplot
7. **中间数据** - gsea_res.rds + HTML 报告

## Out-of-scope (不做)

- ❌ **ORA 过表达分析**（见 go-kegg 工具）
- ❌ **ssGSEA / GSVA 单样本富集**（见 bindea-immune 工具的 GSVA 路径）
- ❌ **iGSEA / GSEA-P 软件**（外部）
- ❌ **leading-edge 基因导出到 Cytoscape**（手动）
- ❌ **自定义 GMT**（留给高级用户）

## Default params (官方默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `minSize` | `15` | 官方默认（msigdbr 常用 15-500）|
| `maxSize` | `500` | 官方默认 |
| `pvalueCutoff` | `0.25` | GSEA 官方默认（比 ORA 宽松）|
| `pAdjustMethod` | `BH` | 官方默认 |
| `seed` | `TRUE` | 官方默认（结果可复现）|

## Expected outputs

图（双格式 PNG+PDF）：dotplot / waterfall / running_score / ridgeplot / emapplot / cnetplot / treeplot
中间数据：gsea_go_results.csv / gsea_kegg_results.csv / hallmark_results.csv / gsea_res.rds

## 工具间数据联动

- ✅ 上游：可导入 deseq2 的 deg_results.csv（按 log2FC 排序）
- ✅ 下游：leading-edge 基因可下载
