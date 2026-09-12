# consensus-cluster demo data

TCGA-LIHC top-variance features 子集(派生自 ml-feature-select demo):

- `matrix.csv`:200 features × 150 samples
- `known_labels.csv`:Tumor / Normal 已知标签(用于 ARI 验证)

## 用法

consensus-cluster 工具应该 **不知道**已知标签的情况下找出最优 k,
理想情况下:
1. Δ(area) 法选 k* = 2(因为这是 Tumor vs Normal 二分)
2. mean silhouette > 0.3(数据有清晰二分结构)
3. ARI vs known labels > 0.5(发现的 2 个亚型与 Tumor/Normal 高度对应)

## 数据来源

复用 ml-feature-select demo 的 vst_hvg5000 子集。原始数据:TCGA-LIHC(`~/biof3-data/tcga-lihc/`)

本 demo 用于 BioF3 ML W3 module07 教程演示 + consensus-cluster 工具公开 demo。
