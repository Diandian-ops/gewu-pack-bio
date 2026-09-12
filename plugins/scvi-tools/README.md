# scvi-tools

> **科学子问题**：在给定 AnnData 输入下，scvi-tools 的 scVI 变分自编码器能否在 biof3-py-runtime + 扩展依赖（scvi-tools>=1.0 + torch>=2.0）边界内训练 + 提取 latent embedding + 标准化表达 + UMAP 用于批次整合？

## 职责

使用 scvi-tools 1.0+ 对 AnnData (`.h5ad`) 执行 scVI 深度生成模型管线：

1. **preflight**：校验 h5ad 输入 + counts layer
2. **注册模型**：`scvi.model.SCVI.setup_anndata` (含 batch_key 可选)
3. **训练**：`scvi.model.SCVI.train(max_epochs, batch_size, lr, accelerator='cpu')`
4. **Latent embedding**：`model.get_latent_representation()` → 落盘 CSV
5. **Normalized expression**：`model.get_normalized_expression(library_size=1e4)` → 落盘 CSV
6. **UMAP**：基于 latent 走 scanpy 邻居图 + UMAP
7. **落盘**：8 outputs + latent_umap_batch.png + persisted AnnData

## DAG 位置

```
sample-metadata-validator → scanpy-advanced → scvi-tools（深度生成批次整合）
```

`scanpy-advanced` 的 `preprocessed.h5ad`（含 PCA/UMAP/leiden）可作为本插件上游入参；`scvi-tools` 在隐空间做批次整合 + 深度生成式标准化。

## 输入

| 参数 | 类型 | 必需 | 默认 | 说明 |
|------|------|------|------|------|
| `anndata_file` | file (`.h5ad`) | ✓ | — | AnnData 输入 |
| `batch_key` | string | — | `batch` | `obs` 列名（如果在 obs 中存在则注入 batch_key） |
| `layer` | string | — | `counts` | counts layer 名（scvi-tools 需要 raw counts） |
| `n_latent` | number | — | `10` | scVI latent dim |
| `n_hidden` | number | — | `128` | scVI hidden dim |
| `n_layers` | number | — | `2` | scVI hidden layer 数 |
| `max_epochs` | number | — | `5` | 训练 epoch 数 |
| `batch_size` | number | — | `64` | 训练 batch_size |
| `device` | select | — | `cpu` | `cpu` / `auto`（auto 在 G0 阶段降级为 cpu） |
| `seed` | number | — | `42` | 随机种子 |

## 输出

| 文件 | 类型 | 说明 |
|------|------|------|
| `latent_embedding.csv` | table | scVI 隐空间嵌入 (cells × n_latent) |
| `normalized_expression.csv` | table | 标准化表达矩阵 (cells × genes) |
| `integrated.h5ad` | data | 整合后 AnnData（含 `obsm['X_scVI']`） |
| `latent_umap.png` | plot | UMAP 隐空间（未着色） |
| `latent_umap_batch.png` | plot | UMAP 隐空间（按 batch 着色） |
| `training_loss.csv` | table | epoch vs train_loss |
| `qc_summary.csv` | table | QC 摘要 |
| `report.html` | report | 简洁分析报告 |
| `manifest.json` | manifest | 工具 manifest（含 stats） |

## 状态

**dev Electron 边界内 verified**（G0 阶段）：

- biof3-py-runtime v1.0.0 基础栈 + `pip install scvi-tools torch` 扩展（`install.json` 步骤）
- 实测版本：scvi-tools 1.4.2 + torch 2.13.0（macOS arm64, Python 3.11.9, CPU）
- 270 cells × 100 vars demo fixture 训练 5 epochs → `(270, 10)` latent + `(270, 100)` normalized expression
- 设备策略：**CPU**（避免引入 GPU 依赖链；GPU/MPS 验证属 G1+ 阶段）
- 不修改 biof3-py-runtime 制品本身（扩展走 install.json 步骤在 biof3-py-runtime venv 内部执行）

> **G0 阶段纪律（按 AGENTS.md "当前阶段纪律"）**：本插件仅在 dev Electron 边界内 verified；packaged / Provider / GA / dist staleness 属 G1+ 阶段事项，按本阶段纪律**不外推**。`qualificationEligible: false` 是当前 G0 开发期常态，不是缺陷、不被视为缺口。

## Python 包依赖

| Package | 版本 | 来源 |
|---|---|---|
| scanpy | ≥ 1.9 | biof3-py-runtime 基础栈 |
| anndata | ≥ 0.9 | biof3-py-runtime 基础栈 |
| numpy | ≥ 1.24 | biof3-py-runtime 基础栈 |
| pandas | ≥ 2.0 | biof3-py-runtime 基础栈 |
| scipy | ≥ 1.10 | biof3-py-runtime 基础栈 |
| matplotlib | ≥ 3.7 | biof3-py-runtime 基础栈 |
| seaborn | ≥ 0.12 | biof3-py-runtime 基础栈 |
| Pillow | ≥ 9.0 | biof3-py-runtime 基础栈 |
| **scvi-tools** | **≥ 1.0** | **额外 pip install（install.json）** |
| **torch** | **≥ 2.0** | **额外 pip install（install.json）** |

> scvi-tools + torch 不在 biof3-py-runtime 基础栈中（PHASE 0 runtime inventory §6.1 + §3 已确认）。在 venv 内部执行 `pip install scvi-tools torch` 即可（实测 macOS arm64 ~5min，scvi-tools 1.4.2 + torch 2.13.0）。

## Runtime 扩展决策

按 `agent-tools/block-11-bioSkills-integration-plan-v0.2.md` §6.1，本插件在 SUB-1.5 启动时已**显式裁决**：

- **不修改 biof3-py-runtime 制品本身**（避免 147MB → ~600MB 体积膨胀）
- **走 install.json 步骤**：biof3-py-runtime venv 内部 `pip install scvi-tools torch`
- **CPU 模式**：避免触发 GPU 依赖链（GPU/MPS 验证属 G1+ 阶段）
- **未来扩展决策保留**：G1+ 阶段可考虑（a）增量加 PyTorch 到 biof3-py-runtime；（b）建独立 `biof3-scvi-tools-runtime` 子集

## 平台

- ✅ macos-arm64 / windows-x64（dev Electron runtime + install.json 步骤就位即可跑）
- ⏳ packaged 边界验证属 G1+ 阶段，本批次不外推

## 关联 Skill

- **PHASE 1 SUB-1.6**：`sc-analysis-standard` skill 在 SUB-1.6 子批才创建，**不动**（单工作包串行纪律）
- 本插件在 `sc-analysis-standard` skill 的 `plugins.optional` 中（与 celltypist / cellchat 并列）

## 注意

- `scripts/scvi-tools.py` 末尾对所有输出 PNG 调 `_canonicalize_png()`，用 Pillow read+save 一次以稳定 sha256（与 R 侧 `png::writePNG(png::readPNG(f), f)` 同义），保证下游 artifact-lineage / Result Studio hash 不漂移
- `runtime-not-ready` 失败模式在 `acceptance.json` 标 `skipOnSupportedPlatform: true`（dev Electron 扩展后视为支持，跳过）
- 不在 closeout 报告列举"未做 packaged / 未做 Provider"等 G1+ 阶段事项（按 AGENTS.md G0 阶段纪律）
