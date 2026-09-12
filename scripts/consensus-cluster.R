#!/usr/bin/env Rscript
# ============================================================
# Consensus Clustering — Unsupervised subtype discovery
#
# Based on ConsensusClusterPlus (Bioconductor)
#   k = 2 to max_k, resampling + voting consensus
#   Auto k selection: delta area elbow + CDF plateau
#
# Outputs:
#   - CDF / delta area / consensus matrix panel
#   - Silhouette score per k
#   - PCA projection by subtype
#   - ARI validation (if known_labels provided)
#   - Subtype assignment table
#
# BioF3 SCI visual style: theme_biof3(), biof3_palette(), ggsave_biof3()
# ============================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(cluster)
})

need_pkgs <- c("ConsensusClusterPlus", "NMF", "mclust")
# BLOCK-58 — 依赖执行路径零安装：缺包毫秒级 fail-fast（结构化错误，绝不联网）。
# 依赖单一事实源：_dependencies.json；显式补装：scripts/install-r-tool-dependencies.mjs
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("MISSING_R_PACKAGES: %s\n", paste(missing_pkgs, collapse = ", ")))
  cat("Install first via: node scripts/install-r-tool-dependencies.mjs --tool consensus-cluster\n")
  quit(status = 1)
}

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(job_dir, "_biof3-theme.R"))

cat("[consensus-cluster] Starting consensus clustering (SCI style)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
max_k <- as.integer(params$max_k %||% 8)
reps <- as.integer(params$reps %||% 1000)
distance <- params$distance %||% "pearson"

# Read data
read_matrix <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f)) read.delim(f, row.names = 1, check.names = FALSE)
  else read.csv(f, row.names = 1, check.names = FALSE)
}

matrix_file <- list.files(job_dir, pattern = "^matrix", full.names = TRUE)[1]
known_labels_file <- list.files(job_dir, pattern = "^known_labels", full.names = TRUE)[1]

mat <- as.matrix(read_matrix(matrix_file))

# Remove constant / near-constant features
var_row <- apply(mat, 1, var, na.rm = TRUE)
mat <- mat[var_row > 0.001, , drop = FALSE]

sample_ids <- colnames(mat)
cat(sprintf("[consensus-cluster] Samples: %d, Features: %d\n", ncol(mat), nrow(mat)))

# ============================================================
# ConsensusClusterPlus
# ============================================================
cat(sprintf("[consensus-cluster] Running ConsensusClusterPlus (k=2:%d, reps=%d)...\n", max_k, reps))
suppressMessages(library(ConsensusClusterPlus))

# Distance mapping
dist_map <- list(pearson = "pearson", spearman = "spearman", euclidean = "euclidean")
ccp_dist <- dist_map[[distance]] %||% "pearson"

ccp_res <- ConsensusClusterPlus(
  d = mat,
  maxK = max_k,
  reps = reps,
  pItem = 0.8,
  pFeature = 0.8,
  clusterAlg = "hc",
  distance = ccp_dist,
  title = file.path(output_dir, "ccp"),
  plot = NULL,  # We generate our own plots
  seed = 42
)

# ============================================================
# Extract metrics for k selection
# ============================================================
cat("[consensus-cluster] Extracting consensus metrics...\n")
k_vals <- 2:max_k
cdf_vals <- sapply(k_vals, function(k) ccp_res[[k]]$cdf)
delta_area <- diff(cdf_vals)

metrics_df <- data.frame(
  k = k_vals,
  cdf = cdf_vals,
  delta_area = c(NA, delta_area)
)
metrics_df <- metrics_df[order(metrics_df$k), ]
write.csv(metrics_df, file.path(output_dir, "metrics.csv"), row.names = FALSE)

# ============================================================
# Delta area plot (optimal k)
# ============================================================
cat("[consensus-cluster] Plotting delta area...\n")
# Elbow: find largest drop in delta_area
if (length(delta_area) >= 2) {
  elbow_idx <- which.max(-delta_area) + 1  # +1 because delta is k[i] - k[i-1]
  optimal_k <- k_vals[elbow_idx]
} else {
  optimal_k <- 2
}

p_delta <- ggplot(metrics_df, aes(x = k, y = delta_area)) +
  geom_line(color = "#0f766e", linewidth = 0.8) +
  geom_point(size = 2, color = "#0f766e") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "#ef4444") +
  annotate("text", x = optimal_k, y = max(delta_area, na.rm = TRUE) * 0.9,
           label = paste("Optimal k =", optimal_k), color = "#ef4444", hjust = -0.1) +
  labs(title = "Delta Area (Consensus CDF)", x = "k (clusters)", y = "Δ Area") +
  theme_biof3()
ggsave_biof3(p_delta, file.path(output_dir, "delta_area"), width = 6, height = 4.5)

# CDF plot
p_cdf <- ggplot(metrics_df, aes(x = k, y = cdf)) +
  geom_line(color = "#0f766e", linewidth = 0.8) +
  geom_point(size = 2, color = "#0f766e") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "#ef4444") +
  labs(title = "CDF (Cumulative Distribution Function)", x = "k (clusters)", y = "CDF") +
  theme_biof3()
ggsave_biof3(p_cdf, file.path(output_dir, "cdf_plot"), width = 6, height = 4.5)

# ============================================================
# Consensus matrix heatmaps (for k = 2:max_k, panel)
# ============================================================
cat("[consensus-cluster] Plotting consensus matrices...\n")

# Plot each k as a separate panel, then combine
library(gridExtra)
consensus_plots <- list()
for (k in k_vals) {
  cm <- ccp_res[[k]]$consensusMatrix
  # Reorder by cluster assignment
  cls <- ccp_res[[k]]$clrs$consensusClass
  ord <- order(cls)
  cm_ord <- cm[ord, ord]

  # Convert to long format for ggplot
  cm_df <- expand.grid(x = seq_len(ncol(cm_ord)), y = seq_len(nrow(cm_ord)))
  cm_df$value <- as.vector(cm_ord)

  p_cm <- ggplot(cm_df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#0f766e", limits = c(0, 1)) +
    labs(title = paste("k =", k), x = NULL, y = NULL, fill = "Consensus") +
    coord_equal() +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5))
  consensus_plots[[as.character(k)]] <- p_cm
}

save_grid_biof3(file.path(output_dir, "consensus_matrices"), width = 12, height = 10, expr = {
  grid.arrange(grobs = lapply(consensus_plots, ggplotGrob), ncol = 3,
               top = textGrob("Consensus Matrices", gp = gpar(fontsize = 14, fontface = "bold")))
})

# ============================================================
# Silhouette score per k
# ============================================================
cat("[consensus-cluster] Computing silhouette scores...\n")
sil_scores <- sapply(k_vals, function(k) {
  cls <- ccp_res[[k]]$clrs$consensusClass
  if (length(unique(cls)) < 2) return(0)
  dm <- as.dist(1 - ccp_res[[k]]$consensusMatrix)
  sil <- silhouette(cls, dm)
  mean(sil[, 3], na.rm = TRUE)
})

sil_df <- data.frame(k = k_vals, silhouette = sil_scores)
write.csv(sil_df, file.path(output_dir, "silhouette.csv"), row.names = FALSE)

p_sil <- ggplot(sil_df, aes(x = k, y = silhouette)) +
  geom_line(color = "#0f766e", linewidth = 0.8) +
  geom_point(size = 2, color = "#0f766e") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "#ef4444") +
  labs(title = "Silhouette Score", x = "k (clusters)", y = "Mean Silhouette") +
  theme_biof3()
ggsave_biof3(p_sil, file.path(output_dir, "silhouette"), width = 6, height = 4.5)

# ============================================================
# Subtype assignment (for optimal k)
# ============================================================
cat("[consensus-cluster] Assigning subtypes for optimal k =", optimal_k, "...\n")
opt_cls <- ccp_res[[optimal_k]]$clrs$consensusClass
subtype_df <- data.frame(
  sample = sample_ids,
  cluster = opt_cls,
  stringsAsFactors = FALSE
)
write.csv(subtype_df, file.path(output_dir, "subtype_assignment.csv"), row.names = FALSE)

# ============================================================
# PCA projection by subtype
# ============================================================
cat("[consensus-cluster] PCA projection...\n")
pca_res <- prcomp(t(mat), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  cluster = factor(opt_cls),
  sample = sample_ids,
  stringsAsFactors = FALSE
)

pca_pal <- biof3_palette(optimal_k)
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, label = sample)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = pca_pal) +
  labs(title = paste("PCA by Subtype (k =", optimal_k, ")"),
       x = sprintf("PC1 (%.1f%%)", summary(pca_res)$importance[2, 1] * 100),
       y = sprintf("PC2 (%.1f%%)", summary(pca_res)$importance[2, 2] * 100),
       color = "Subtype") +
  theme_biof3()
ggsave_biof3(p_pca, file.path(output_dir, "pca_projection"), width = 6.5, height = 5)

# ============================================================
# ARI validation (if known_labels provided)
# ============================================================
ari_text <- "No known labels provided, skipping ARI validation.\n"
if (file.exists(known_labels_file)) {
  cat("[consensus-cluster] Computing ARI with known labels...\n")
  suppressMessages(library(mclust))
  known_df <- read_matrix(known_labels_file)
  if (ncol(known_df) >= 2) {
    known_ids <- as.character(known_df[[1]])
    known_labels <- as.character(known_df[[2]])
  } else {
    known_ids <- rownames(known_df)
    known_labels <- as.character(known_df[[1]])
  }
  names(known_labels) <- known_ids

  common_samp <- intersect(sample_ids, known_ids)
  if (length(common_samp) > 0) {
    found_cls <- opt_cls[match(common_samp, sample_ids)]
    known_cls <- known_labels[common_samp]
    ari_val <- adjustedRandIndex(found_cls, known_cls)

    ari_text <- sprintf(
      "ARI Validation (k = %d)\nSamples with known labels: %d\nAdjusted Rand Index: %.4f\nInterpretation: %s\n",
      optimal_k, length(common_samp), ari_val,
      ifelse(ari_val > 0.7, "Excellent agreement",
             ifelse(ari_val > 0.4, "Moderate agreement", "Weak agreement"))
    )
    writeLines(ari_text, file.path(output_dir, "ari_validation.txt"))
  }
}

# ============================================================
# Summary
# ============================================================
summary_text <- paste0(
  "Consensus Clustering Subtype Discovery\n\n",
  sprintf("Samples: %d, Features: %d\n", ncol(mat), nrow(mat)),
  sprintf("k range: 2 to %d, reps: %d, distance: %s\n", max_k, reps, distance),
  sprintf("Optimal k (delta area elbow): %d\n", optimal_k),
  sprintf("Mean silhouette (optimal k): %.3f\n", sil_scores[optimal_k - 1]),
  "\nConsensus distribution (optimal k):\n",
  paste(sprintf("  Subtype %d: %d", seq_len(optimal_k), table(opt_cls)), collapse = "\n"),
  "\n\n", ari_text
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = list(
    list(name = "metrics.csv", type = "table", label = "Metrics (CDF + delta area)"),
    list(name = "silhouette.csv", type = "table", label = "Silhouette scores"),
    list(name = "subtype_assignment.csv", type = "table", label = "亚型注释表"),
    list(name = "delta_area.png", type = "plot", label = "Δarea 最优 k"),
    list(name = "delta_area.pdf", type = "file", label = "Δarea 最优 k (PDF)"),
    list(name = "cdf_plot.png", type = "plot", label = "CDF 累积分布"),
    list(name = "cdf_plot.pdf", type = "file", label = "CDF 累积分布 (PDF)"),
    list(name = "consensus_matrices.png", type = "plot", label = "一致性矩阵热图"),
    list(name = "consensus_matrices.pdf", type = "file", label = "一致性矩阵热图 (PDF)"),
    list(name = "silhouette.png", type = "plot", label = "Silhouette score"),
    list(name = "silhouette.pdf", type = "file", label = "Silhouette score (PDF)"),
    list(name = "pca_projection.png", type = "plot", label = "PCA 亚型投影"),
    list(name = "pca_projection.pdf", type = "file", label = "PCA 亚型投影 (PDF)"),
    list(name = "ari_validation.txt", type = "text", label = "ARI 验证"),
    list(name = "summary.txt", type = "text", label = "分析摘要")
  ),
  summary = list(
    samples = ncol(mat),
    features = nrow(mat),
    optimal_k = optimal_k,
    distance = distance
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

cat("[consensus-cluster] Pipeline complete.\n")
