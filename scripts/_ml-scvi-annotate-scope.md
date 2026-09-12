# scVI 细胞类型注释工具范围卡

## Tool ID

`ml-scvi-annotate`（task type: `scvi_annotation`）

## Version targeted

- BioF3 Python 容器: `Python 3.11+`
- 关键包: `scvi-tools` / `scanpy` / `scikit-learn` / `anndata`

## Canonical workflow source

**官方文档**: scvi-tools tutorials (Lopez et al., 2018 Nat Methods)

## Steps (in-scope) — scVI 参考映射注释

1. **读参考 + 查询 h5ad** — 参考数据已预训练或现场训练
2. **scVI 训练** — 参考集 latent 空间学习
3. **查询集映射** — SCANVI / transfer learning
4. **kNN 注释** — latent 空间 kNN 投票
5. **可视化** — UMAP (参考+查询) / 注释置信度 / 混淆矩阵
6. **中间数据** — latent_embedding.csv / annotation_results.csv

## Out-of-scope (不做)

- ❌ **Seurat label transfer**（见 seurat-standard 工具）
- ❌ **手动 marker 注释**（用户自行）
- ❌ **CellTypist / scArches**（其他注释工具，独立）
- ❌ **参考集构建/训练**（用户自行准备或用预训练模型）
- ❌ **空间细胞类型注释**（见 spatial-deconv）

## Default params (官方默认)

| 参数 | 默认值 | 来源 |
|------|--------|------|
| `n_latent` | `30` | scVI 默认 |
| `n_layers` | `2` | scVI 默认 |
| `max_epochs` | `100` | scVI 常用 |
| `knn_k` | `15` | kNN 常用 |

## Expected outputs

图（PNG）：UMAP (参考+查询) / 注释置信度 / 混淆矩阵
中间数据：latent_embedding.csv / annotation_results.csv / summary.json

## 工具间数据联动

- ✅ 上游：用户提供参考 h5ad + 查询 h5ad
- ✅ 下游：annotation_results.csv 可回导 seurat-standard 做下游分析
