# scanpy-advanced

> **科学子问题**: 在给定 AnnData 输入下，单细胞 RNA-seq 标准 高级分析管线能否在 biof3-py-runtime (dev Electron 边界内) 跑通 QC → normalize → HVG → PCA → UMAP → Leiden → markers？

## 职责

使用 Scanpy 1.9+ 对 AnnData (`.h5ad`) 执行单细胞 RNA-seq 高级分析管线：

1. **QC**：sc.pp.filter_cells / sc.pp.filter_genes (基于 min_genes / min_cells)
2. **Normalize**：sc.pp.normalize_total + sc.pp.log1p
3. **HVG**：sc.pp.highly_variable_genes (seurat flavor)
4. **Scale + PCA**：sc.pp.scale + sc.tl.pca
5. **Neighbors + UMAP**：sc.pp.neighbors + sc.tl.umap
6. **Leiden 聚类**：sc.tl.leiden (igraph flavor)
7. **Marker genes**：sc.tl.rank_genes_groups (默认 wilcoxon)
8. **落盘 figures + tables + report.html + 持久化 AnnData**

## DAG 位置

```
sample-metadata-validator → scanpy-advanced → (future: celltypist / cellchat / scvi-tools)
```

`sample-metadata-validator` 的 `validated_counts.csv` 经本地 h5ad 包装或上游 `scanpy-cluster` 链式产物进入本插件。`scripts/scanpy-advanced.py` 不直接读 CSV — 只接受 `.h5ad`。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| `anndata_file` | file (`.h5ad`) | ✓ | — | AnnData 输入 |
| `min_genes` | number | — | `200` | 每细胞最少检出基因数 |
| `min_cells` | number | — | `3` | 每基因最少在多少细胞中检出 |
| `n_top_genes` | number | — | `2000` | HVG 数量 |
| `n_pcs` | number | — | `30` | PCA 主成分数 |
| `n_neighbors` | number | — | `15` | k-NN 邻居数 |
| `resolution` | number | — | `0.5` | Leiden resolution |
| `method` | select | — | `wilcoxon` | `rank_genes_groups` 方法 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `preprocessed.h5ad` | data | 预处理后 AnnData (含 PCA / UMAP / leiden / rank_genes_groups) |
| `marker_genes.csv` | table | 每 cluster 的 top markers (gene, score, pvalue, padj) |
| `qc_summary.csv` | table | QC 摘要 (n_obs_in/post, clusters, markers...) |
| `umap.png` | plot | UMAP 可视化（未着色） |
| `umap_clustered.png` | plot | UMAP 按 leiden 着色 |
| `hvg_plot.png` | plot | HVG 分布 |
| `report.html` | report | 简洁分析报告 |
| `manifest.json` | manifest | 工具 manifest (含 stats) |

## 状态

**dev Electron 边界内 verified**（G0 阶段）：

- biof3-py-runtime v1.0.0 已装 `scanpy>=1.9` + `anndata>=0.9` + `numpy/pandas/scipy/matplotlib/seaborn/Pillow`（见 `resources/built-in-plugins/py-runtime-requirements.txt`）
- 不依赖额外 pip 安装（PHASE 1 SUB-1.1 任务边界内 runtime **不扩展**）
- 离线可执行；产出 progress.json 驱动 UI 进度条

> **G0 阶段纪律（按 AGENTS.md "当前阶段纪律"）**：本插件仅在 dev Electron 边界内 verified；packaged / Provider / GA / dist staleness 属 G1+ 阶段事项，按本阶段纪律**不外推**。`qualificationEligible: false` 是当前 G0 开发期常态，不是缺陷、不被视为缺口。

## Python 包依赖

- scanpy ≥ 1.9（biof3-py-runtime 声明）
- anndata ≥ 0.9
- numpy ≥ 1.24
- pandas ≥ 2.0
- scipy ≥ 1.10
- matplotlib ≥ 3.7
- seaborn ≥ 0.12（间接）
- Pillow ≥ 9.0（PNG hash 重写）

## 平台

- ✅ macos-arm64 / windows-x64（dev Electron runtime 就位即可跑）
- ⏳ packaged 边界验证属 G1+ 阶段，本批次不外推

## 关联 Skill

- **PHASE 1 SUB-1.6**：`sc-analysis-standard` skill 在 SUB-1.6 子批才创建，**不动**（单工作包串行纪律）
- 本插件 standalone 可被任意 `biof3-py-runtime` Python 调用链直接消费

## 注意

`scripts/scanpy-advanced.py` 末尾对所有输出 PNG 调 `_canonicalize_png()`，用 Pillow read+save 一次以稳定 sha256（与 R 侧 `png::writePNG(png::readPNG(f), f)` 同义），保证下游 artifact-lineage / Result Studio hash 不漂移。
