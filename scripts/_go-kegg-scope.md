# GO/KEGG 富集分析工具范围卡

## Tool ID

`go-kegg`

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `clusterProfiler 4.18+` / `enrichplot 1.30+` / `ggplot2 4.0+`

## Canonical workflow source

**官方文档**: clusterProfiler 主 vignette (Yu et al., 2012 OMICS)

## Steps (in-scope)

1. **读输入** - gene list 或 DESeq2 result
2. **ID 转换** - SYMBOL -> ENTREZID
3. **enrichGO** - BP / MF / CC 三大类
4. **enrichKEGG** - 通路富集
5. **compareCluster** - Up vs Down 对比（如输入是 DESeq2 result）
6. **可视化** - dotplot / barplot / cnetplot / treeplot / heatplot / emapplot
7. **中间数据** - ego_*.rds + *_data.csv + HTML 报告

## Out-of-scope (不做)

- ❌ **GSEA**（有序全基因富集，见 gsea 工具）
- ❌ **ReactomePA / WikiPathways**（用户可自行扩展）
- ❌ **自定义 GMT 富集**（留给高级用户）
- ❌ **Enrichment map 导出 Cytoscape**（手动）
- ❌ **DAVID / Metascape 对接**（外部平台）

## Default params (官方默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `pvalueCutoff` | `0.05` | 官方默认 |
| `qvalueCutoff` | `0.2` | 官方默认 |
| `pAdjustMethod` | `BH` | 官方默认 |
| `organism` | `hsa` | BioF3 默认人类 |
| `ont` | `ALL (BP/MF/CC)` | 一次跑全 |

## Expected outputs

图（双格式 PNG+PDF）：dotplot / barplot / cnetplot / treeplot / heatplot / emapplot / compareCluster
中间数据：go_bp_results.csv / go_mf_results.csv / go_cc_results.csv / kegg_results.csv / ego_*.rds

## 工具间数据联动

- ✅ 上游：可导入 deseq2 的 deg_results.csv
- ✅ 下游：富集结果表可下载，gene 列表可回导
