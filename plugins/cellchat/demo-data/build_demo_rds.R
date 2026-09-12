#!/usr/bin/env Rscript
# ============================================================
# build_demo_rds.R — 生成 cellchat demo fixture Seurat RDS
#
# 270 cells × 100 PBMC marker + LR-family genes (real symbol names)
# 3 cell types: CD4_T / CD8_T / B_cell
# Deterministic via set.seed(20260818)
#
# 设计：基因列表优先选 CellChatDB.human 中能形成 LR pair 的 symbols（单基
# 因受体 ligand），加上经典 PBMC markers，确保 demo fixture 跑通
# computeCommunProb 时至少有若干 LR pair（valid fixture 不落入 empty-result
# 分支）。完整 1939 个 LR pair 中 100 个 gene 仍只覆盖约 30 个 ligand +
# 30 个 receptor = 约 50-100 LR pair 命中。
#
# Run:
#   Rscript resources/built-in-plugins/cellchat/demo-data/build_demo_rds.R
#
# Idempotent: 重新运行会覆盖相同文件。
# ============================================================

set.seed(20260818)
suppressMessages(library(Seurat))

# PBMC markers + CellChatDB LR-family symbols (single-gene receptor preferred)
pbmc_genes <- c(
  # CD4_T markers (15)
  "CD3D","CD3E","CD3G","CD4","IL7R","CCR7","TCF7","LEF1","MAL","CD2",
  "CD5","CD28","ICOS","CTLA4","CD40LG",
  # CD8_T markers (15)
  "CD8A","CD8B","CCL5","GZMK","GZMA","PRF1","NKG7","GNLY","KLRD1","FASLG",
  "FAS","IFNG","TBX21","CCL3","CCL4",
  # B_cell markers (10)
  "MS4A1","CD19","CD22","CD74","HLA-DRA","HLA-DRB1","CD79A","CD79B","PAX5","BLK",
  # Chemokines (10)
  "CCL19","CCL21","CXCL12","CXCL8","CXCL1","CCL2","CCL20","NAMPT","MDK","GAS6",
  # Notch / Wnt (8)
  "NOTCH1","NOTCH2","JAG1","JAG2","DLL1","WNT5A","FZD4","FZD5",
  # TNF superfamily (10)
  "TNFSF4","TNFRSF4","TNFSF9","TNFRSF9","CD70","CD27","TNFSF11","TNFRSF11A",
  "LGALS9","HAVCR2",
  # Immune checkpoints (8)
  "CD80","CD86","CD274","PDCD1","CD40","PROS1","AXL","MERTK",
  # Cytokine receptors / growth factors (14)
  "TGFB1","TGFB2","TGFBR1","TGFBR2","VEGFA","FLT1","KDR","BMP4","BMP7",
  "ACVR1","BMPR2","CCR7","CXCR4","TLR4",
  # Housekeeping (10)
  "ACTB","GAPDH","B2M","HPRT1","RPL13A","RPS18","TUBB","UBC","PPIA","SDHA"
)

# 去重 + 截断到 100
pbmc_genes <- unique(pbmc_genes)
if (length(pbmc_genes) > 100) pbmc_genes <- pbmc_genes[1:100]
if (length(pbmc_genes) < 100) {
  pbmc_genes <- c(pbmc_genes, sprintf("g%03d", seq_len(100 - length(pbmc_genes))))
}

n_obs <- 270
n_vars <- length(pbmc_genes)
n_types <- 3
cells_per_type <- n_obs %/% n_types

# 模拟 counts: 3 clusters with shifted mean on a subset of genes
counts <- matrix(
  stats::rnbinom(n_obs * n_vars, mu = 5, size = 0.3),
  nrow = n_vars, ncol = n_obs
)
# cluster 1 (CD4_T) bumps marker genes + TGF/VEGFA pathway (genes 1-40)
counts[1:40, 1:cells_per_type] <- counts[1:40, 1:cells_per_type] + 8L
# cluster 2 (CD8_T) bumps CD8 markers + FASLG/FAS (genes 21-60)
counts[21:60, (cells_per_type + 1):(2 * cells_per_type)] <-
  counts[21:60, (cells_per_type + 1):(2 * cells_per_type)] + 8L
# cluster 3 (B_cell) bumps B_cell markers (genes 51-70)
counts[51:70, (2 * cells_per_type + 1):n_obs] <-
  counts[51:70, (2 * cells_per_type + 1):n_obs] + 8L

rownames(counts) <- pbmc_genes
colnames(counts) <- sprintf("cell%03d", seq_len(n_obs))

seurat <- CreateSeuratObject(
  counts = counts,
  project = "pbmc_3k_mini",
  min.cells = 0,
  min.features = 0
)
seurat <- NormalizeData(
  seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
seurat <- FindVariableFeatures(
  seurat,
  selection.method = "vst",
  nfeatures = n_vars,
  verbose = FALSE
)

cell_types <- c("CD4_T", "CD8_T", "B_cell")[((seq_len(n_obs) - 1L) %/% cells_per_type) + 1L]
seurat$cell_type <- factor(cell_types)

DefaultAssay(seurat) <- "RNA"

# Resolve out_path from --file= commandArgs
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}
out_path <- file.path(script_dir, "pbmc_3k_mini.rds")

saveRDS(seurat, out_path)
cat(sprintf("[demo-data build] Wrote %s (%d cells x %d features, 3 cell types)\n",
            out_path, ncol(seurat), nrow(seurat)))