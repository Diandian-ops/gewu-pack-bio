# celltypist

> **科学子问题**: 在给定 AnnData 输入下，Celltypist 1.6+ 自动 + 手动 cell type annotation 能否在 biof3-py-runtime (dev Electron 边界内) 跑通 Model.predict + majority voting？

## 职责

使用 Celltypist 1.6+ 对 AnnData (`.h5ad`) 执行 single-cell RNA-seq cell type 自动注释：

1. **QC**：sc.pp.filter_cells / sc.pp.filter_genes (基于 min_genes / min_cells)
2. **normalize + log1p + HVG**：上游或本脚本自跑
3. **PCA + UMAP + Leiden**：上游或本脚本自跑（celltypist 注释通常基于已有 leiden clusters）
4. **Celltypist Model.predict**：
   - `params.model` (本地 .pkl 路径) → 直接 `Model.load(path)` 离线友好
   - `params.model_name` (e.g. Immune_All_Low.pkl) → 需联网下载
   - 默认走 `Immune_All_Low.pkl`（PBMC 68 类型）
5. **majority voting**：在 over_clustering 指定的 cluster 列上做多数投票
6. **fallback**：若 celltypist 不可用 / model 加载失败，使用 scanpy-based 启发式 fallback（cluster-level 多数投票 + cluster 占比作为 conf_score；`fallback_used=1`）
7. **落盘 `predictions.csv` / `probabilities.csv` / `confusion_matrix.csv` / `annotated.h5ad` / `umap_annotation.png` / `umap_confidence.png` / `report.html` / `manifest.json`**

## DAG 位置

```
scanpy-advanced → celltypist → (后续: cellchat / scvi-tools)
```

Celltypist 接收 `scanpy-advanced.outputs.preprocessed.h5ad`（或任意上游 QC 后的 .h5ad 输入），输出 `.obs.predicted_labels` / `conf_score` / `majority_voting`，供下游 cell-cell communication / latent representation 等分析。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| `anndata_file` | file (`.h5ad`) | ✓ | — | AnnData 输入 |
| `model` | string | — | `null` | 本地 celltypist model .pkl 路径（离线友好，import 优先） |
| `model_name` | string | — | `Immune_All_Low.pkl` | celltypist 内置 model 名（需联网下载） |
| `majority_voting` | boolean | — | `true` | 是否启用 over_clustering 多数投票 |
| `over_clustering` | string | — | `"leiden"` | 多数投票 cluster 列（obs 列） |
| `p_thres` | number | — | `0.5` | multi-label 概率阈值（仅 mode='prob match' 生效） |
| `mode` | select | — | `best match` | celltypist annotate mode: `best match` / `prob match` / `binary` |
| `min_prop` | number | — | `0.0` | majority voting 最小 cluster 占比 |
| `min_genes` | number | — | `200` | 每细胞最少检出基因数 |
| `min_cells` | number | — | `3` | 每基因最少在多少细胞中检出 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `predictions.csv` | table | per-cell predicted_labels / conf_score / majority_voting |
| `probabilities.csv` | table | top-N predicted labels per cell |
| `confusion_matrix.csv` | table | real vs predicted (仅当 obs.cell_type 存在) |
| `annotated.h5ad` | data | 含 obs.predicted_labels / conf_score 的 AnnData |
| `umap_annotation.png` | plot | UMAP 按 predicted_labels 着色 |
| `umap_confidence.png` | plot | UMAP 按 conf_score 着色 |
| `report.html` | report | 简洁报告（含 stats + 模型来源 + fallback 标记） |
| `manifest.json` | manifest | Job manifest (含 stats + fallback_used ∈ {0,1}) |

## 状态

**dev Electron 边界内 verified**（G0 阶段）：

- biof3-py-runtime v1.0.0 已装 `scanpy>=1.9` + `anndata>=0.9` + `numpy/pandas/scipy/matplotlib/seaborn/Pillow`（见 `resources/built-in-plugins/py-runtime-requirements.txt`）
- celltypist 本体（PyPI ~10MB）**不在** biof3-py-runtime；install.json 触发 `pip install celltypist>=1.6` 局部补装（**不**动 biof3-py-runtime 主体）
- celltypist 不可用 / model 缺失时自动走 scanpy-based heuristic fallback（`fallback_used=1`），保证 valid fixture 始终跑通
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
- **celltypist ≥ 1.6**（PyPI bootstrap，install.json 触发）

## 平台

- ✅ macos-arm64 / windows-x64 (dev Electron runtime + celltypist bootstrap 即跑)
- ⏳ packaged 边界验证属 G1+ 阶段，本批次不外推

## 关联 Skill

- **PHASE 1 SUB-1.6**：`sc-analysis-standard` skill 在 SUB-1.6 子批才创建，**不动**（单工作包串行纪律）
- 本插件 standalone 可被任意 `biof3-py-runtime` Python 调用链直接消费

## 模型选择

| 来源 | 参数 | 是否需要联网 | 备注 |
|------|------|-------------|------|
| 本地 .pkl | `model: /path/to/model.pkl` | ❌ | 离线友好，PHASE 1 推荐 |
| celltypist 内置名 | `model_name: Immune_All_Low.pkl` | ✅ | 需联网下载（约 50-100MB） |
| 默认 | 都不提供 | ✅ | 默认 `Immune_All_Low.pkl`（PBMC 68 类型） |

G0 离线优先：本插件 PHASE 1 默认假定用户提供本地 `.pkl`（`model` 路径）或接受 fallback（`fallback_used=1` + 启发式注释）。**不**自动触发下载。

## 注意

`scripts/celltypist.py` 末尾对所有输出 PNG 调 `_canonicalize_png()`，用 Pillow read+save 一次以稳定 sha256（与 R 侧 `png::writePNG(png::readPNG(f), f)` 同义），保证下游 artifact-lineage / Result Studio hash 不漂移。
