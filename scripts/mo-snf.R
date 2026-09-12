#!/usr/bin/env Rscript
# ============================================================
# SNF (Similarity Network Fusion) — Multi-omics network fusion
#
# Input: 2-3 omics layers (CSV, rows=samples, cols=features)
# Output: fused network heatmap / subtype assignment / eigen-gap /
#         NMI+ARI validation (if known_labels provided)
#
# BioF3 SCI visual style: theme_biof3(), biof3_palette(), ggsave_biof3()
# ============================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(pheatmap)
})

need_pkgs <- c("SNFtool", "mclust")
# BLOCK-58 — 依赖执行路径零安装：缺包毫秒级 fail-fast（结构化错误，绝不联网）。
# 依赖单一事实源：_dependencies.json；显式补装：scripts/install-r-tool-dependencies.mjs
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("MISSING_R_PACKAGES: %s\n", paste(missing_pkgs, collapse = ", ")))
  cat("Install first via: node scripts/install-r-tool-dependencies.mjs --tool mo-snf\n")
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

cat("[mo-snf] Starting SNF multi-omics fusion (SCI style)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
K <- as.integer(params$K %||% 20)
alpha <- as.numeric(params$alpha %||% 0.5)
T_iter <- as.integer(params$T_iter %||% 20)
max_k <- as.integer(params$max_k %||% 8)

# Read omics layers
read_matrix <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f)) read.delim(f, row.names = 1, check.names = FALSE)
  else read.csv(f, row.names = 1, check.names = FALSE)
}

layer_files <- list.files(job_dir, pattern = "^omics_layer", full.names = TRUE)
if (length(layer_files) < 2) {
  stop("Need at least 2 omics layers. Found: ", length(layer_files))
}

omics_list <- list()
for (i in seq_along(layer_files)) {
  m <- as.matrix(read_matrix(layer_files[i]))
  # Remove constant / near-constant features
  vr <- apply(m, 2, var, na.rm = TRUE)
  m <- m[, vr > 0.001, drop = FALSE]
  omics_list[[i]] <- m
}

# Intersect samples (rownames must match across layers)
sample_ids <- Reduce(intersect, lapply(omics_list, rownames))
if (length(sample_ids) == 0) {
  stop("No common samples across layers. Check row names (sample IDs).")
}
cat(sprintf("[mo-snf] Common samples: %d, Layers: %d\n", length(sample_ids), length(omics_list)))

for (i in seq_along(omics_list)) {
  omics_list[[i]] <- omics_list[[i]][sample_ids, , drop = FALSE]
}

# Read optional known labels
known_labels_file <- list.files(job_dir, pattern = "^known_labels", full.names = TRUE)[1]
known_labels <- NULL
if (file.exists(known_labels_file)) {
  kdf <- read_matrix(known_labels_file)
  if (ncol(kdf) >= 2) {
    sids <- as.character(kdf[[1]])
    labs <- as.character(kdf[[2]])
  } else {
    sids <- rownames(kdf)
    labs <- as.character(kdf[[1]])
  }
  known_labels <- setNames(labs, sids)
}

# ============================================================
# SNF pipeline
# ============================================================
cat("[mo-snf] Building similarity networks...\n")
suppressMessages(library(SNFtool))

# Standardize each layer
omics_std <- lapply(omics_list, standardNormalization)

# Affinity matrices (KNN similarity networks)
W_list <- lapply(omics_std, function(d) {
  affinityMatrix(d, K = K, sigma = alpha)
})

# Fuse networks
if (length(W_list) >= 2) {
  W_fused <- SNF(W_list, K = K, t = T_iter)
} else {
  W_fused <- W_list[[1]]
}

# ============================================================
# Eigen-gap for optimal k
# ============================================================
cat("[mo-snf] Estimating optimal k via eigen-gap...\n")
optimal_k <- estimateNumberOfClustersGivenGraph(W_fused, 2:max_k)$
  optimalK
if (is.null(optimal_k) || optimal_k < 2 || optimal_k > max_k) {
  optimal_k <- 2
}
cat(sprintf("[mo-snf] Optimal k (eigen-gap): %d\n", optimal_k))

# Spectral clustering
clusters <- spectralClustering(W_fused, optimal_k)
names(clusters) <- sample_ids

# ============================================================
# Subtype assignment table
# ============================================================
subtype_df <- data.frame(
  sample = sample_ids,
  cluster = clusters,
  stringsAsFactors = FALSE
)
write.csv(subtype_df, file.path(output_dir, "subtype_assignment.csv"), row.names = FALSE)

# ============================================================
# Fused network heatmap
# ============================================================
cat("[mo-snf] Plotting fused network heatmap...\n")
# Reorder by cluster
ord <- order(clusters)
W_ord <- W_fused[ord, ord]

hm_pal <- biof3_palette_seq(100, option = "ylorbr")
save_grid_biof3(file.path(output_dir, "fused_network_heatmap"), width = 7, height = 6, expr = {
  pheatmap(W_ord,
           color = hm_pal,
           border_color = NA,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           show_colnames = FALSE,
           main = "Fused Similarity Network (SNF)")
})

# ============================================================
# Eigen-gap plot
# ============================================================
cat("[mo-snf] Plotting eigen-gap...\n")
egap <- sapply(2:max_k, function(k) {
  eig <- eigen(W_fused, symmetric = TRUE, only.values = TRUE)$values
  # Eigen-gap = difference between k-th and (k+1)-th eigenvalues
  if (length(eig) >= k + 1) eig[k] - eig[k + 1] else 0
})
egap_df <- data.frame(k = 2:max_k, gap = egap)

p_egap <- ggplot(egap_df, aes(x = k, y = gap)) +
  geom_line(color = "#0f766e", linewidth = 0.8) +
  geom_point(size = 2, color = "#0f766e") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "#ef4444") +
  annotate("text", x = optimal_k, y = max(egap, na.rm = TRUE) * 0.9,
           label = paste("Optimal k =", optimal_k), color = "#ef4444", hjust = -0.1) +
  labs(title = "Eigen-gap (Spectral Clustering)", x = "k (clusters)", y = "Eigen-gap") +
  theme_biof3()
ggsave_biof3(p_egap, file.path(output_dir, "eigen_gap"), width = 6, height = 4.5)

# ============================================================
# NMI / ARI validation (if known_labels provided)
# ============================================================
nmi_ari_text <- "No known labels provided, skipping NMI/ARI validation.\n"
if (!is.null(known_labels)) {
  cat("[mo-snf] Computing NMI and ARI...\n")
  suppressMessages(library(mclust))

  common_samp <- intersect(sample_ids, names(known_labels))
  if (length(common_samp) > 0) {
    found_cls <- clusters[match(common_samp, sample_ids)]
    known_cls <- known_labels[common_samp]

    nmi_val <- NA
    tryCatch({
      nmi_val <- NMI(found_cls, known_cls)
    }, error = function(e) {
      cat("[mo-snf] NMI computation failed:", conditionMessage(e), "\n")
    })

    ari_val <- adjustedRandIndex(found_cls, known_cls)

    nmi_ari_text <- sprintf(
      "NMI / ARI Validation (k = %d)\nSamples with known labels: %d\nNMI: %.4f\nARI: %.4f\n",
      optimal_k, length(common_samp),
      ifelse(is.na(nmi_val), 0, nmi_val), ari_val
    )
    writeLines(nmi_ari_text, file.path(output_dir, "nmi_ari_validation.txt"))
  }
}

# ============================================================
# Summary
# ============================================================
summary_text <- paste0(
  "SNF Multi-Omics Network Fusion\n\n",
  sprintf("Layers: %d, Common samples: %d\n", length(omics_list), length(sample_ids)),
  sprintf("SNF parameters: K=%d, alpha=%.2f, T=%d\n", K, alpha, T_iter),
  sprintf("Optimal k (eigen-gap): %d\n", optimal_k),
  "\nSubtype distribution:\n",
  paste(sprintf("  Subtype %d: %d", seq_len(optimal_k), table(clusters)), collapse = "\n"),
  "\n\n", nmi_ari_text
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = list(
    list(name = "subtype_assignment.csv", type = "table", label = "亚型注释表"),
    list(name = "fused_network_heatmap.png", type = "plot", label = "融合网络热图"),
    list(name = "fused_network_heatmap.pdf", type = "file", label = "融合网络热图 (PDF)"),
    list(name = "eigen_gap.png", type = "plot", label = "Eigen-gap 最优 k"),
    list(name = "eigen_gap.pdf", type = "file", label = "Eigen-gap (PDF)"),
    list(name = "nmi_ari_validation.txt", type = "text", label = "NMI/ARI 验证"),
    list(name = "summary.txt", type = "text", label = "分析摘要")
  ),
  summary = list(
    layers = length(omics_list),
    samples = length(sample_ids),
    optimal_k = optimal_k,
    K = K,
    T_iter = T_iter
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

cat("[mo-snf] Pipeline complete.\n")
