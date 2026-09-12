#!/usr/bin/env Rscript
# ============================================================
# build_demo_rds.R — 生成 seurat-standard demo Seurat RDS fixture
#
# 输出: resources/built-in-plugins/seurat-standard/demo-data/pbmc_3k_mini.rds
#
# 设计原则：
#   - 270 cells × 100 genes × 3 simulated clusters (与 scanpy-advanced demo 形态对齐)
#   - 确定性: set.seed(42) 锁死
#   - 体积小: < 100 KB（避免污染 repo size）
#   - 满足：min.cells=3 + min.features=1（小 fixture 用 1 让全部 270 cells 保留）
#
# 运行:
#   /path/to/biof3-r-runtime/bin/Rscript build_demo_rds.R
# ============================================================

suppressMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(42)

n_obs <- 270
n_vars <- 100
n_types <- 3
cells_per_type <- n_obs %/% n_types

# IMPORTANT: Seurat expects matrix where rows = genes, cols = cells.
# But CsparseMatrix creates a (nrow, ncol) matrix. Need to build as (n_vars × n_obs)
# then transpose. Below: X is genes × cells = n_vars × n_obs.
X <- matrix(
  rnbinom(n_vars * n_obs, mu = 5, size = 0.3),
  nrow = n_vars,
  ncol = n_obs
)
storage.mode(X) <- "double"

# Cluster mean shifts — apply per cell (column)
bump_A <- matrix(rnbinom(cells_per_type * 30, mu = 8, size = 0.3), nrow = 30, ncol = cells_per_type)
bump_B <- matrix(rnbinom(cells_per_type * 30, mu = 8, size = 0.3), nrow = 30, ncol = cells_per_type)
X[1:30, 1:cells_per_type] <- X[1:30, 1:cells_per_type] + bump_A
X[31:60, (cells_per_type + 1):(2 * cells_per_type)] <-
  X[31:60, (cells_per_type + 1):(2 * cells_per_type)] + bump_B

# Gene names & cell barcodes
gene_names <- sprintf("gene_%03d", seq_len(n_vars))
cell_barcodes <- sprintf("cell_%04d", seq_len(n_obs))
cluster_labels <- c(rep("A", cells_per_type),
                    rep("B", cells_per_type),
                    rep("C", n_obs - 2 * cells_per_type))

X <- as(X, "CsparseMatrix")
rownames(X) <- gene_names
colnames(X) <- cell_barcodes

# Build Seurat object — use min.features = 1 so all 270 cells stay
obj <- CreateSeuratObject(counts = X, project = "biof3-pbmc3k-mini",
                          min.cells = 3, min.features = 1)
obj$simulated_cluster <- cluster_labels

# Save RDS — script_dir-relative
script_dir <- tryCatch({
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) ".")
out_path <- file.path(script_dir, "pbmc_3k_mini.rds")

saveRDS(obj, out_path)
cat(sprintf("[demo-data] wrote %s (%d bytes, %d genes × %d cells)\n",
            out_path,
            file.info(out_path)$size,
            nrow(obj),
            ncol(obj)))