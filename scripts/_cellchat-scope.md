# CellChat 工具范围卡（步骤 0 落地）

> tool-dev "步骤 0：工具范围与官方文档" 第四个落地实例
> 第二个使用 SCI 视觉风格规范 (theme_biof3 + ggsave_biof3) 的工具

## Tool ID

`cellchat`

## Version targeted

- BioF3 R 容器: R 4.5.1
- `CellChat 2.2.0` / `NMF 0.28` / `ggalluvial 0.12.5` / `future 1.69.0`
- `circlize 0.4.17` / `ComplexHeatmap 2.24.1`（CellChat 内部依赖）

## Canonical workflow source

**官方 vignette**: CellChat single-dataset analysis tutorial
(https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html)

**Reference**: Jin S et al. (2025). CellChat for systematic analysis of cell-cell communication from single-cell transcriptomics. *Nature Protocols*.

## Steps (in-scope) — Official CellChat Single-Dataset Workflow

1. **Load** — 接收 normalized expression matrix (CSV) + cell metadata (CSV with cell_type column)
2. **CellChat 对象创建** — `createCellChat()` + 选择 LR database（CellChatDB.human / mouse）
3. **预处理** — `subsetData()` → `identifyOverExpressedGenes()` → `identifyOverExpressedInteractions()`
4. **通讯概率推断** — `computeCommunProb(raw.use=TRUE)` + `filterCommunication(min.cells=10)`
5. **通路水平汇总** — `computeCommunProbPathway()` + `aggregateNet()`
6. **可视化输出**:
   - 通讯网络图（圆形 + 热图）— `netVisual_circle` + `netVisual_heatmap`
   - 显著 LR 气泡图 — `netVisual_bubble`
   - 通路 chord/violin — `netVisual_aggregate`(top pathways)
   - 信号角色（incoming / outgoing） — `netAnalysis_signalingRole_scatter`

## Out-of-scope (不做)

- ❌ **多数据集对比** — `mergeCellChat` + `compareCellChat`，独立工具更合理
- ❌ **Mouse 之外的物种** — 默认 human，mouse 切换；不内建其他物种
- ❌ **Spatial CellChat** — 空间组学专用，独立工具
- ❌ **NMF pattern analysis** — `identifyCommunicationPatterns()` 进阶分析，留给 RunModal
- ❌ **配体-受体水平的复杂热图** — `netVisual_heatmap(cellchat, signaling=...)` 单 pathway 切片图，留给 RunModal

> `_data.ts` description: "CellChat 单数据集细胞-细胞通讯分析（基于 LR 数据库）：通讯网络图 + 显著 LR 气泡 + 信号角色分布。需要：归一化表达矩阵 + 细胞类型标注。不含: 多数据集对比、空间 CellChat、NMF pattern。"

## Default params

| 参数 | 默认 | 说明 |
|------|------|------|
| `species` | `human` | 决定用 CellChatDB.human 还是 CellChatDB.mouse |
| `signaling_db` | `Secreted Signaling` | 三个子库之一，最常用 |
| `min_cells_per_group` | 10 | filterCommunication 的最小细胞数阈值 |
| `top_pathways_to_visualize` | 5 | 出 chord/violin 的通路数（pathway 多时）|
| `compute_method` | `truncatedMean` | computeCommunProb 默认 (`triMean` 用于细胞数 ≥30/group) |

## Expected outputs (8 张图 + 4 中间数据)

### 图（每张 PNG + PDF 双格式 = 16 文件）

| # | 文件名 | 说明 |
|---|--------|------|
| 1 | `network_circle_count` | 圆形网络图（边=互作数） |
| 2 | `network_circle_weight` | 圆形网络图（边=互作强度）|
| 3 | `network_heatmap` | 通讯热图（行/列 = cell type）|
| 4 | `bubble_top_LR` | top 显著 LR 对气泡图 |
| 5 | `signaling_role_scatter` | 信号角色（outgoing vs incoming） |
| 6 | `pathway_aggregate_top1` | 最显著 pathway 的 chord 图 |
| 7 | `pathway_violin_top` | top pathway 在各 cell type 的 violin |
| 8 | `pathway_count_barplot` | 各 cell type 的 outgoing/incoming pathway 数 |

### 中间数据

| 文件 | 内容 |
|------|------|
| `lr_significant.csv` | 显著的配体-受体对（含 cell_type pair + prob + pval） |
| `pathway_significance.csv` | 每个通路的总体强度排序 |
| `signaling_role_data.csv` | 每个 cell type 的 outgoing / incoming 强度 |
| `cellchat.rds` | 完整 CellChat 对象 |

## Demo 数据

- 来源: PBMC 3k 全基因 normalized matrix（基于 BioF3 seurat-standard demo 派生）
- 维度: 9821 基因 × 605 cells
- Cell types: 9 个 PBMC canonical（Naive_CD4_T / Memory_CD4 / CD14_Mono / B / CD8_T / FCGR3A_Mono / NK / DC / Platelet）
- 文件:
  - `expression.csv` 19.2 MB（normalized log expression）
  - `meta.csv` 14.6 KB（cell_id, cell_type）
- 验证跑通时间: 2 分钟（120 秒），16 个显著通路（MHC-I / CD99 / MHC-II / MIF / GALECTIN 等 PBMC 经典）

## 工具间数据联动

- ✅ **后续可从 Seurat-standard 导入**（待实现）：用户跑完 Seurat 后导出 normalized matrix + cluster 标注 → 自动填到 CellChat 输入
- ❌ 上游不需要其他工具

## 已读官方文档（T1）

- [x] CellChat single-dataset analysis vignette
- [x] CellChat 2.2.0 NEWS / API doc
- [x] BioF3 自己的教程 `docs/single-cell/module07.md`

## 已知风险

1. **Demo 19.2MB → base64 26MB JSON**：仍在 200MB body limit 内，但前端 `textToBase64` 处理 19MB 要 1-2 秒（已 chunked，OK）
2. **OOM**：CellChat 对 8000+ cells 矩阵会撞 4GB 容器内存。脚本加 dim check：> 5000 cells → 拒绝
3. **物种支持**：仅 human / mouse（CellChatDB 限制）
4. **Cell type 太少（< 3）会无意义**：脚本检查 `length(unique(cell_type)) >= 3`，否则报错

## 下一步

按 14 步流程：
- 步骤 0: 范围卡 ✅
- 步骤 1: ID = `cellchat` ✅
- 步骤 2: 更新 `_data.ts` 工具定义
- 步骤 3: `routes/tools.js` 已注册 ✅
- 步骤 4: `index.tsx` `TOOL_COLORS` + `REPORT_URLS` 已注册 ✅
- 步骤 5: 编写 `cellchat.R`（按 SCI 风格规范用 _biof3-theme.R）
- 步骤 6-14: 按流程走
