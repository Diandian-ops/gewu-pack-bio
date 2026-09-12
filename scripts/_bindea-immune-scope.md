# 免疫浸润评估工具范围卡

## Tool ID

`bindea-immune`（脚本文件 `bindea-immune.R`）

## Version targeted

- BioF3 R 容器: `R 4.5.x`
- 关键包: `GSVA 1.54+`

## Canonical workflow source

**官方文档**: GSVA vignette + Bindea et al., 2013 Immunity（28 immune cell types）

## Steps (in-scope)

1. **读 TPM 矩阵** - 行=基因 Symbol，列=样本
2. **ssGSEA / Bindea 28 immune cells** - method 参数切换
3. **免疫细胞比例矩阵**
4. **可视化** - 堆叠条形图 / 相关性热图 / 组间箱线图
5. **HTML 报告 + manifest**

## Out-of-scope (不做)

- ❌ **CIBERSORT**（需外部许可/代码，独立工具）
- ❌ **xCell / TIMER / EPIC**（其他反卷积方法，独立工具）
- ❌ **Bulk 细胞类型 deconvolution 全谱**（见 spatial-deconv 工具）
- ❌ **免疫微环境打分**（ESTIMATE，独立工具）
- ❌ **单细胞参考映射**（见 scvi-annotate）

## Default params (官方默认)

| 参数 | 默认值 | 来源 / 说明 |
|------|--------|---------|
| `method` | `ssgsea` | GSVA 官方 ssgsea 模式 |
| `kcdf` | `Gaussian` | 官方默认（连续 TPM 数据）|
| `tau` | `0.4` | ssgsza 默认 |
| `min_size` | `1` | Bindea 集合已精选 |

## Expected outputs

图（双格式 PNG+PDF）：堆叠条形图 / 相关性热图 / 组间箱线图
中间数据：immune_scores.csv / immune_correlation.csv

## 工具间数据联动

- ✅ 上游：用户提供 TPM 矩阵（deseq2 vst 或上游归一化结果）
- ✅ 下游：免疫评分可整合进多组学 / ML 工具
