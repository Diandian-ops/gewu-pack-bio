# cellchat

> **科学子问题**: 在给定 Seurat 对象（带细胞类型注释）下，CellChat 能否在 biof3-r-runtime (dev Electron 边界内) 跑通讯概率推断 → 通路级通讯 → 通讯结果表 + 热图 + 弦图 + 气泡图？

## 职责

使用 CellChat 1.6+ 对 Seurat 对象执行单细胞细胞-细胞通讯推断：

1. **DB 准备**：`subsetDB(CellChatDB = CellChatDB.human/mouse/zebrafish, search = db_category)`
2. **createCellChat**：从 Seurat 对象构建 CellChat 实例，`group.by = cell_type_col`
3. **computeCommunProb**：`type = "truncatedMean"`, `trim = 0.1`
4. **filterCommunication**：`min.cells = min_cells`
5. **computeCommunProbPathway**：通路级（CellChatDB pathway 注释）
6. **aggregateNet**：聚合 cell group × cell group 通讯数量与强度
7. **可视化**：`netVisual_heatmap` + `netVisual_chord_cell` + `netVisual_bubble`
8. **可选层级推断**：`selectK` + `identifyCommunicationPatterns`（`compute_hierarchy = "true"` 时）
9. **落盘**：communications_table.csv / pathway_enrichment.csv / interaction_heatmap.png / chord_diagram.png / bubble_plot.png / cellchat_object.rds / report.html / manifest.json

## DAG 位置

```
scanpy-advanced / seurat → cellchat → (future: 下游 netAnalysis_* 复用 cellchat_object.rds)
```

上游 Seurat 对象由 `scanpy-advanced`（Python → Seurat conversion via sceasy / anndata2ri）或 `seurat` 标准 R 管线产出。本插件只接受 `.rds` Seurat 对象。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| `seurat_rds` | file (`.rds`) | ✓ | — | Seurat 对象 RDS |
| `cell_type_col` | text | — | `cell_type` | meta.data 中的细胞类型列名 |
| `species` | select | — | `human` | `human` / `mouse` / `zebrafish` |
| `db_category` | select | — | `Secreted Signaling` | `Secreted Signaling` / `Cell-Cell Contact` / `ECM-Receptor` / `Non-protein Signaling` / `all` |
| `min_cells` | number | — | `10` | 每组最少细胞数（computeCommunProb + filterCommunication） |
| `pvalue_threshold` | number | — | `0.05` | 置换检验 p-value 阈值 |
| `n_patterns` | number | — | `5` | NMF latent pattern 数量（仅 `compute_hierarchy=true` 时使用） |
| `compute_hierarchy` | select | — | `false` | 是否执行 selectK + identifyCommunicationPatterns |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `communications_table.csv` | table | 长表：每行一个 LR pair（source / target / ligand / receptor / prob / pval / pathway_name） |
| `pathway_enrichment.csv` | table | 按 pathway 聚合（pathway / n_lr_pairs / mean_prob / pvalue） |
| `interaction_heatmap.png` | plot | cell group × cell group 通讯数量热图（netVisual_heatmap） |
| `chord_diagram.png` | plot | 弦图（netVisual_chord_cell，top-1 pathway） |
| `bubble_plot.png` | plot | 通路气泡图（netVisual_bubble，sources × targets） |
| `cellchat_object.rds` | data | CellChat 对象 RDS（供下游 `netAnalysis_*` 复用） |
| `report.html` | report | 简洁分析报告 |
| `manifest.json` | manifest | 工具 manifest（含 stats） |

## 运行时

- **R runtime**: biof3-r-runtime v1.0.0（R 4.5.2 + Seurat 5.5.1 + Bioconductor 3.22）
- **额外 R 包**: CellChat 1.6+（Bioconductor 包）+ igraph / ComplexHeatmap / circlize / NMF / ggalluvial / BiocManager
- **安装**: `Rscript -e 'BiocManager::install("CellChat", ask=FALSE, update=FALSE)'`（一次性，~120MB）
- **数据库**: `CellChatDB.human` / `CellChatDB.mouse` / `CellChatDB.zebrafish` 已随 CellChat 包内置，运行时无需联网
- **状态**: dev Electron 边界内 verified（需 CellChat install）
- **不外推**: packaged / Provider / GA / dist staleness（G1+ 阶段事项）

## 已知边界

- **CellChat 本体不在 biof3-r-runtime 显式清单**（r-runtime-packages.txt 26 包）：通过 `install.json` 中的 BiocManager::install 一次性补齐；如安装失败（如离线 CRAN 镜像阻塞）acceptance.json 中 runtime-not-ready fixture 保持 `skipOnSupportedPlatform:true`，插件仍可 8 件套合规但运行时 blocked
- **compute_hierarchy = true 时执行 selectK（NMF rank selection）+ identifyCommunicationPatterns**：CPU-heavy，2-5 min on PBMC 3k 默认 fixture
- **chord_diagram.png 仅画 top-1 pathway**：避免长 list 让 chord 不可读；如需多 pathway 可后续扩展 chord_grid
- **bubble_plot 限定 top-20 pathways**：避免 sources × targets × pathways 笛卡尔爆
- **CellChat 不输出通路级 p-value**：pathway_enrichment.csv 的 pvalue 字段恒为 NA，下游如需真 p-value 需 aggregate 后自定 hypergeometric

## 相关文档

- `docs/block-11-phase-1-batch-plan.md` §3.4（cellchat 工作包规划）
- `docs/block-11-bioSkills-runtime-inventory.md` §5 / §7（biof3-r-runtime 已装包清单 + PHASE 1-4 runtime 扩展需求）
- `docs/block-11-bioSkills-implementation-patterns.md` §3-6（plugin / environment / tool-definition / acceptance 模板）