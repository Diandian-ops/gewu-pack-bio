# rsf-survival demo data

TCGA-LIHC 子集(BioF3 ML W2-D 派生):
- features.csv: top 200 高方差基因 × 365 肿瘤样本
- clinical.csv: 365 样本, OS_event=1: 130 (35.6%)

## 字段
- sample_id: TCGA barcode
- OS_time: 总生存时间(单位:days)
- OS_event: 0=alive, 1=dead

## 来源
UCSC Xena Hub TCGA-LIHC HiSeqV2 + clinicalMatrix

## 预期输出
C-index + 时间 ROC + VIMP + KM + Risk 三联图(< 90s 跑完)

生成时间: 2026-05-27
