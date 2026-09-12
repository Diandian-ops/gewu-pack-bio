#!/usr/bin/env Rscript
# ============================================================
# MOFA2 Multi-Omics Factor Analysis
#
# Input: 2-4 omics layers (CSV, rows=features, cols=samples)
# Output: variance explained heatmap / top weights / factor scatter /
#         factor correlation / data overview / factor scores table
#
# BioF3 SCI visual style: theme_biof3(), biof3_palette(), ggsave_biof3()
# ============================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(pheatmap)
})

need_pkgs <- c("MOFA2")
# BLOCK-58 — 依赖执行路径零安装：缺包毫秒级 fail-fast（结构化错误，绝不联网）。
# 依赖单一事实源：_dependencies.json；显式补装：scripts/install-r-tool-dependencies.mjs
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("MISSING_R_PACKAGES: %s\n", paste(missing_pkgs, collapse = ", ")))
  cat("Install first via: node scripts/install-r-tool-dependencies.mjs --tool mo-mofa\n")
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

cat("[mo-mofa] Starting MOFA2 multi-omics analysis (SCI style)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
num_factors <- as.integer(params$num_factors %||% 15)
conv_mode <- params$convergence_mode %||% "slow"

# Convergence settings
conv_iters <- switch(conv_mode, slow = 1000, medium = 500, fast = 100, 1000)

# Read omics layers
read_matrix <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f)) read.delim(f, row.names = 1, check.names = FALSE)
  else read.csv(f, row.names = 1, check.names = FALSE)
}

layer_files <- list.files(job_dir, pattern = "^omics_layer", full.names = TRUE)
layer_names <- paste0("layer", seq_along(layer_files))

if (length(layer_files) < 2) {
  stop("Need at least 2 omics layers. Found: ", length(layer_files))
}

omics_list <- list()
for (i in seq_along(layer_files)) {
  m <- as.matrix(read_matrix(layer_files[i]))
  # Remove constant features
  vr <- apply(m, 1, var, na.rm = TRUE)
  m <- m[vr > 0.001, , drop = FALSE]
  omics_list[[layer_names[i]]] <- m
}

# Intersect samples across layers
common_samples <- Reduce(intersect, lapply(omics_list, colnames))
if (length(common_samples) == 0) {
  stop("No common samples across layers. Check column names.")
}
cat(sprintf("[mo-mofa] Common samples: %d\n", length(common_samples)))

for (nm in names(omics_list)) {
  omics_list[[nm]] <- omics_list[[nm]][, common_samples, drop = FALSE]
}

# Read optional group file
group_file <- list.files(job_dir, pattern = "^group_file", full.names = TRUE)[1]
groups <- NULL
if (file.exists(group_file)) {
  gdf <- read_matrix(group_file)
  if (ncol(gdf) >= 2) {
    sids <- as.character(gdf[[1]])
    glabs <- as.character(gdf[[2]])
  } else {
    sids <- rownames(gdf)
    glabs <- as.character(gdf[[1]])
  }
  groups <- setNames(glabs[sids %in% common_samples], sids[sids %in% common_samples])
  # Ensure all samples have a group
  missing <- setdiff(common_samples, names(groups))
  if (length(missing) > 0) {
    groups[missing] <- "Unknown"
  }
  groups <- groups[common_samples]
}

# ============================================================
# Build MOFA2 model
# ============================================================
cat("[mo-mofa] Building MOFA2 model (", length(omics_list), "layers,", num_factors, "factors)...\n")
suppressMessages(library(MOFA2))

mofa_data <- create_mofa_object(data = omics_list)
mofa_data <- set_data_options(mofa_data, scale_views = TRUE)

# Set model options
model_opts <- get_default_model_options(mofa_data)
model_opts$num_factors <- num_factors

# Training options
train_opts <- get_default_training_options(mofa_data)
train_opts$maxiter <- conv_iters
train_opts$verbose <- FALSE

tryCatch({
  mofa_model <- run_mofa(mofa_data, model_options = model_opts,
                          training_options = train_opts, outfile = file.path(output_dir, "mofa_model.hdf5"))
}, error = function(e) {
  cat("[mo-mofa] MOFA2 training failed:", conditionMessage(e), "\n")
  # Fallback: write empty outputs and exit gracefully
  writeLines("MOFA2 training failed. Data may be too sparse or dimensions mismatch.\n", file.path(output_dir, "summary.txt"))
  writeLines('{"files":[],"summary":{"status":"failed"}}', file.path(output_dir, "manifest.json"))
  quit(status = 0)
})

cat("[mo-mofa] Model trained.\n")

# Extract factors
factors_df <- get_factors(mofa_model, as.data.frame = TRUE)
factors_wide <- reshape(factors_df, direction = "wide", idvar = "sample", timevar = "factor")
write.csv(factors_wide, file.path(output_dir, "factor_scores.csv"), row.names = FALSE)

# Extract variance explained
var_df <- get_variance_explained(mofa_model, as.data.frame = TRUE)
write.csv(var_df, file.path(output_dir, "variance_table.csv"), row.names = FALSE)

# ============================================================
# Plot 1: Variance explained heatmap (factor x view)
# ============================================================
cat("[mo-mofa] Plotting variance explained heatmap...\n")
var_mat <- get_variance_explained(mofa_model, per_factor = TRUE)$R2PerFactor
if (!is.null(var_mat) && ncol(var_mat) > 0) {
  hm_pal <- biof3_palette_seq(100, option = "ylgn")
  save_grid_biof3(file.path(output_dir, "variance_heatmap"), width = 7, height = 5, expr = {
    pheatmap(var_mat,
             color = hm_pal,
             border_color = NA,
             display_numbers = TRUE,
             number_format = "%.1f",
             fontsize = 9,
             main = "Variance Explained (%) per Factor & View")
  })
}

# ============================================================
# Plot 2: Total variance per view
# ============================================================
cat("[mo-mofa] Plotting total variance per view...\n")
total_r2 <- get_variance_explained(mofa_model)$R2Total
total_df <- data.frame(
  view = names(total_r2),
  r2 = as.numeric(total_r2),
  stringsAsFactors = FALSE
)

p_total <- ggplot(total_df, aes(x = reorder(view, r2), y = r2)) +
  geom_col(fill = "#0f766e", width = 0.7) +
  labs(title = "Total Variance Explained per View", x = "View", y = "R² (%)") +
  coord_flip() +
  theme_biof3()
ggsave_biof3(p_total, file.path(output_dir, "variance_total"), width = 5.5, height = 4)

# ============================================================
# Plot 3: Top weights per factor
# ============================================================
cat("[mo-mofa] Plotting top weights...\n")
weights_list <- list()
for (v in names(omics_list)) {
  w <- get_weights(mofa_model, view = v, as.data.frame = TRUE)
  if (nrow(w) > 0) {
    # Top 5 absolute weights per factor
    w$abs_w <- abs(w$value)
    w_top <- do.call(rbind, by(w, w$factor, function(df) {
      head(df[order(-df$abs_w), ], 5)
    }))
    weights_list[[v]] <- w_top
  }
}

if (length(weights_list) > 0) {
  w_all <- do.call(rbind, weights_list)
  # Pick top 3 factors by total variance
  factor_order <- names(sort(rowSums(var_mat, na.rm = TRUE), decreasing = TRUE))[1:min(3, nrow(var_mat))]
  w_sub <- w_all[w_all$factor %in% factor_order, ]

  p_weights <- ggplot(w_sub, aes(x = reorder(feature, abs_w), y = value, fill = factor)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_grid(view ~ factor, scales = "free_y") +
    labs(title = "Top Weights per Factor", x = "Feature", y = "Weight", fill = "Factor") +
    coord_flip() +
    theme_biof3() +
    theme(strip.text = element_text(size = 9))
  ggsave_biof3(p_weights, file.path(output_dir, "top_weights"), width = 10, height = 8)
}

# ============================================================
# Plot 4: Factor scatter (F1 vs F2)
# ============================================================
cat("[mo-mofa] Plotting factor scatter...\n")
f1 <- factors_wide$factor.1
f2 <- factors_wide$factor.2
scatter_df <- data.frame(
  sample = factors_wide$sample,
  F1 = f1,
  F2 = f2,
  stringsAsFactors = FALSE
)
if (!is.null(groups)) {
  scatter_df$group <- groups[as.character(scatter_df$sample)]
}

p_scatter <- ggplot(scatter_df, aes(x = F1, y = F2)) +
  geom_point(size = 2, alpha = 0.7, color = "#0f766e") +
  geom_text_repel(aes(label = sample), size = 2.5, max.overlaps = 10, color = "black") +
  labs(title = "Factor Scatter (F1 vs F2)", x = "Factor 1", y = "Factor 2") +
  theme_biof3()

if (!is.null(groups)) {
  p_scatter <- p_scatter + aes(color = group) +
    scale_color_manual(values = biof3_palette(length(unique(groups)))) +
    labs(color = "Group")
}
ggsave_biof3(p_scatter, file.path(output_dir, "factor_scatter"), width = 6.5, height = 5)

# ============================================================
# Plot 5: Factor correlation
# ============================================================
cat("[mo-mofa] Plotting factor correlation...\n")
fac_mat <- as.matrix(factors_wide[, -1])
cor_fac <- cor(fac_mat, use = "pairwise.complete.obs")
hm_pal2 <- biof3_palette_div(100)
save_grid_biof3(file.path(output_dir, "factor_cor"), width = 6, height = 5, expr = {
  pheatmap(cor_fac,
           color = hm_pal2,
           border_color = NA,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize = 9,
           main = "Factor Correlation")
})

# ============================================================
# Plot 6: Data overview (missing pattern per view)
# ============================================================
cat("[mo-mofa] Plotting data overview...\n")
overview_list <- list()
for (v in names(omics_list)) {
  m <- omics_list[[v]]
  overview_list[[v]] <- data.frame(
    view = v,
    samples = ncol(m),
    features = nrow(m),
    missing_pct = round(sum(is.na(m)) / prod(dim(m)) * 100, 2)
  )
}
overview_df <- do.call(rbind, overview_list)
write.csv(overview_df, file.path(output_dir, "data_overview.csv"), row.names = FALSE)

p_overview <- ggplot(overview_df, aes(x = view, y = features, fill = view)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(samples, " samples\n", missing_pct, "% NA")), vjust = -0.2, size = 3) +
  labs(title = "Data Overview", x = "View", y = "Features", fill = "View") +
  theme_biof3() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave_biof3(p_overview, file.path(output_dir, "data_overview"), width = 5.5, height = 4.5)

# ============================================================
# Summary
# ============================================================
summary_text <- sprintf(
  "MOFA2 Multi-Omics Factor Analysis\n\nLayers: %s\nCommon samples: %d\nFactors: %d\nConvergence: %s (%d iters)\n\nTotal variance explained:\n%s\n\nTop 3 factors by total variance:\n%s\n",
  paste(names(omics_list), collapse = ", "),
  length(common_samples),
  num_factors,
  conv_mode, conv_iters,
  paste(sprintf("  %s: %.2f%%", total_df$view, total_df$r2), collapse = "\n"),
  paste(sprintf("  Factor %s: %.2f%%", factor_order, rowSums(var_mat, na.rm = TRUE)[factor_order]), collapse = "\n")
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = list(
    list(name = "variance_table.csv", type = "table", label = "方差分解表"),
    list(name = "factor_scores.csv", type = "table", label = "因子得分表"),
    list(name = "data_overview.csv", type = "table", label = "数据总览"),
    list(name = "variance_heatmap.png", type = "plot", label = "方差解释热图"),
    list(name = "variance_heatmap.pdf", type = "file", label = "方差解释热图 (PDF)"),
    list(name = "variance_total.png", type = "plot", label = "各层总方差"),
    list(name = "variance_total.pdf", type = "file", label = "各层总方差 (PDF)"),
    list(name = "top_weights.png", type = "plot", label = "因子权重 Top 特征"),
    list(name = "top_weights.pdf", type = "file", label = "因子权重 Top 特征 (PDF)"),
    list(name = "factor_scatter.png", type = "plot", label = "因子散点(F1 vs F2)"),
    list(name = "factor_scatter.pdf", type = "file", label = "因子散点 (PDF)"),
    list(name = "factor_cor.png", type = "plot", label = "因子相关性"),
    list(name = "factor_cor.pdf", type = "file", label = "因子相关性 (PDF)"),
    list(name = "data_overview.png", type = "plot", label = "数据总览"),
    list(name = "data_overview.pdf", type = "file", label = "数据总览 (PDF)"),
    list(name = "summary.txt", type = "text", label = "分析摘要")
  ),
  summary = list(
    layers = length(omics_list),
    samples = length(common_samples),
    factors = num_factors,
    convergence = conv_mode
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

cat("[mo-mofa] Pipeline complete.\n")
