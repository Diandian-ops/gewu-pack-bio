#!/usr/bin/env Rscript
# ============================================================
# SHAP Explainability — Global + single-sample + dependence
#
# Input: upstream ml-classifier job ID (reads model.rds + features/labels)
# Output: global SHAP importance / beeswarm / force plot / dependence plots
#
# Uses iml::Shapley (kernel SHAP) with automatic model-type detection.
# BioF3 SCI visual style: theme_biof3(), biof3_palette(), ggsave_biof3()
# ============================================================

suppressMessages({
  library(ggplot2)
  library(jsonlite)
})

need_pkgs <- c("iml", "e1071", "randomForest", "xgboost", "glmnet")
# BLOCK-58 — 依赖执行路径零安装：缺包毫秒级 fail-fast（结构化错误，绝不联网）。
# 依赖单一事实源：_dependencies.json；显式补装：scripts/install-r-tool-dependencies.mjs
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("MISSING_R_PACKAGES: %s\n", paste(missing_pkgs, collapse = ", ")))
  cat("Install first via: node scripts/install-r-tool-dependencies.mjs --tool shap-explain\n")
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

cat("[shap-explain] Starting SHAP explainability (SCI style)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
upstream_jobid <- params$upstream_jobid
algorithms <- params$algorithms %||% "all"
top_n <- as.integer(params$top_n %||% 20)
n_explain <- as.integer(params$n_explain_samples %||% 30)
sample_index <- as.integer(params$sample_index %||% 1)

# Upstream job directory (jobs are siblings under the same parent)
up_dir <- file.path(dirname(job_dir), upstream_jobid)
if (!dir.exists(up_dir)) {
  stop("Upstream job directory not found: ", up_dir)
}

# Read upstream data
read_matrix <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f)) read.delim(f, row.names = 1, check.names = FALSE)
  else read.csv(f, row.names = 1, check.names = FALSE)
}

features_file <- list.files(up_dir, pattern = "^features", full.names = TRUE)[1]
labels_file   <- list.files(up_dir, pattern = "^labels",   full.names = TRUE)[1]
features <- as.matrix(read_matrix(features_file))
labels_df <- read_matrix(labels_file)

if (ncol(labels_df) >= 2) {
  sample_ids <- rownames(labels_df)
  labels_vec <- labels_df[[2]]
} else {
  sample_ids <- rownames(labels_df)
  labels_vec <- labels_df[[1]]
}
names(labels_vec) <- sample_ids

common <- intersect(colnames(features), sample_ids)
features <- features[, common, drop = FALSE]
labels_vec <- labels_vec[common]
labels_vec <- as.factor(labels_vec)

X <- t(features)
var_x <- apply(X, 2, var, na.rm = TRUE)
X <- X[, var_x > 0, drop = FALSE]
X_df <- as.data.frame(X)

cat(sprintf("[shap-explain] Samples: %d, Features: %d\n", nrow(X), ncol(X)))

# Load upstream model
model_path <- file.path(up_dir, "output", "model.rds")
if (!file.exists(model_path)) {
  stop("Upstream model not found: ", model_path)
}
model <- readRDS(model_path)
cat("[shap-explain] Model loaded: ", class(model)[1], "\n")

# Detect model type
model_type <- NULL
if (inherits(model, "xgb.Booster")) {
  model_type <- "xgb"
  suppressMessages(library(xgboost))
} else if (inherits(model, "randomForest")) {
  model_type <- "rf"
  suppressMessages(library(randomForest))
} else if (inherits(model, "cv.glmnet") || inherits(model, "glmnet")) {
  model_type <- "lr"
  suppressMessages(library(glmnet))
} else if (inherits(model, "svm")) {
  model_type <- "svm"
  suppressMessages(library(e1071))
} else {
  stop("Unsupported model class: ", paste(class(model), collapse = ", "))
}

# ============================================================
# Select top_n features for SHAP (to control runtime)
# ============================================================
cat("[shap-explain] Selecting top", top_n, "features...\n")

imp_scores <- switch(model_type,
  xgb = {
    imp <- xgb.importance(model = model)
    setNames(imp$Gain, imp$Feature)
  },
  rf = {
    imp <- importance(model)[, "MeanDecreaseGini"]
    setNames(imp, names(imp))
  },
  lr = {
    coef_mat <- as.matrix(coef(model, s = "lambda.1se"))[-1, , drop = TRUE]
    abs(coef_mat)
  },
  svm = {
    # SVM has no built-in importance — use linear model approximation (Ridge)
    suppressMessages(library(glmnet))
    y_num <- as.numeric(labels_vec) - 1
    cv_ridge <- cv.glmnet(X, y_num, family = "binomial", alpha = 0, nfolds = 5)
    coef_ridge <- as.matrix(coef(cv_ridge, s = "lambda.1se"))[-1, , drop = TRUE]
    abs(coef_ridge)
  }
)
imp_scores[is.na(imp_scores)] <- 0
imp_scores <- imp_scores[names(imp_scores) %in% colnames(X)]

top_features <- names(sort(imp_scores, decreasing = TRUE))[seq_len(min(top_n, length(imp_scores)))]
X_sub <- X[, top_features, drop = FALSE]
X_df_sub <- as.data.frame(X_sub)

cat(sprintf("[shap-explain] Explaining with %d features\n", ncol(X_sub)))

# ============================================================
# Build iml Predictor
# ============================================================
cat("[shap-explain] Building iml predictor...\n")
suppressMessages(library(iml))

# Predict function that returns probability of positive class
predict_proba <- function(model, newdata) {
  nd <- as.matrix(newdata)
  probs <- switch(model_type,
    xgb = predict(model, xgb.DMatrix(data = nd)),
    rf = predict(model, as.data.frame(nd), type = "prob")[, levels(labels_vec)[2]],
    lr = {
      p <- predict(model, newx = nd, s = "lambda.1se", type = "response")
      if (is.matrix(p) && ncol(p) == 1) as.vector(p) else p[, 1]
    },
    svm = {
      df <- as.data.frame(nd)
      pr <- predict(model, df, probability = TRUE)
      attr(pr, "probabilities")[, levels(labels_vec)[2]]
    }
  )
  data.frame(prob = as.numeric(probs))
}

predictor <- Predictor$new(
  model = model,
  data = X_df_sub,
  y = labels_vec,
  predict.fun = function(m, nd) predict_proba(m, nd)$prob
)

# ============================================================
# Compute SHAP for N samples
# ============================================================
n_explain <- min(n_explain, nrow(X_sub))
sample_idx <- sample(seq_len(nrow(X_sub)), n_explain)

cat(sprintf("[shap-explain] Computing SHAP for %d samples...\n", n_explain))
all_shap <- list()
for (i in seq_along(sample_idx)) {
  sid <- sample_idx[i]
  shap_obj <- Shapley$new(predictor,
                            x.interest = X_df_sub[sid, , drop = FALSE],
                            sample.size = 50)
  res <- shap_obj$results
  res$sample <- sid
  res$feature.value <- as.numeric(X_df_sub[sid, as.character(res$feature)])
  all_shap[[i]] <- res
  if (i %% 10 == 0) cat(sprintf("  %d / %d done\n", i, n_explain))
}
shap_df <- do.call(rbind, all_shap)

# ============================================================
# Global SHAP importance (mean |phi|)
# ============================================================
cat("[shap-explain] Plotting global importance...\n")
shap_global <- aggregate(phi ~ feature, data = shap_df, FUN = function(x) mean(abs(x)))
names(shap_global)[2] <- "mean_abs_shap"
shap_global <- shap_global[order(-shap_global$mean_abs_shap), ]
write.csv(shap_global, file.path(output_dir, "shap_summary.csv"), row.names = FALSE)

p_global <- ggplot(shap_global, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(fill = "#0f766e", width = 0.7) +
  labs(title = "SHAP Global Importance", x = "Feature", y = "Mean |SHAP|") +
  coord_flip() +
  theme_biof3()
ggsave_biof3(p_global, file.path(output_dir, "global_importance"), width = 6, height = 7)

# ============================================================
# Beeswarm plot
# ============================================================
cat("[shap-explain] Plotting beeswarm...\n")
p_bee <- ggplot(shap_df, aes(x = phi, y = reorder(feature, phi, FUN = function(x) mean(abs(x))), color = feature.value)) +
  geom_point(size = 1.2, alpha = 0.6, position = position_jitter(height = 0.2)) +
  scale_color_gradient2(low = "#3b82f6", mid = "white", high = "#ef4444", midpoint = 0) +
  labs(title = "SHAP Beeswarm", x = "SHAP value", y = "Feature", color = "Feature value") +
  theme_biof3()
ggsave_biof3(p_bee, file.path(output_dir, "beeswarm"), width = 7, height = 6)

# ============================================================
# Force plot (single sample)
# ============================================================
cat("[shap-explain] Plotting force plot...\n")
sid <- min(sample_index, n_explain)
sample_shap <- shap_df[shap_df$sample == sample_idx[sid], ]
sample_shap <- sample_shap[order(abs(sample_shap$phi), decreasing = TRUE), ]
sample_shap$direction <- ifelse(sample_shap$phi > 0, "Push up", "Push down")

p_force <- ggplot(sample_shap, aes(x = reorder(feature, abs(phi)), y = phi, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Push up" = "#ef4444", "Push down" = "#3b82f6")) +
  labs(title = sprintf("SHAP Force Plot (sample %d)", sid), x = "Feature", y = "SHAP value", fill = NULL) +
  coord_flip() +
  theme_biof3()
ggsave_biof3(p_force, file.path(output_dir, "force_plot"), width = 6, height = 7)

# ============================================================
# Dependence plots (top 4 features)
# ============================================================
cat("[shap-explain] Plotting dependence plots...\n")
top4 <- head(shap_global$feature, 4)
pd_plots <- list()
for (feat in top4) {
  feat_shap <- shap_df[shap_df$feature == feat, ]
  p_dep <- ggplot(feat_shap, aes(x = feature.value, y = phi)) +
    geom_point(alpha = 0.5, color = "#0f766e") +
    geom_smooth(method = "loess", se = FALSE, color = "#ef4444", linewidth = 0.8) +
    labs(title = feat, x = "Feature value", y = "SHAP value") +
    theme_biof3()
  pd_plots[[feat]] <- p_dep
}

library(gridExtra)
save_grid_biof3(file.path(output_dir, "dependence"), width = 10, height = 8, expr = {
  grid.arrange(grobs = lapply(pd_plots, ggplotGrob), ncol = 2,
               top = textGrob("SHAP Dependence (top 4 features)", gp = gpar(fontsize = 14, fontface = "bold")))
})

# ============================================================
# Save per-sample SHAP matrix
# ============================================================
shap_wide <- reshape(shap_df[, c("sample", "feature", "phi")],
                     direction = "wide", idvar = "sample", timevar = "feature")
write.csv(shap_wide, file.path(output_dir, "shap_global_csv.csv"), row.names = FALSE)

# ============================================================
# Summary
# ============================================================
summary_text <- sprintf(
  "SHAP Explainability Report\n\nUpstream job: %s\nModel type: %s\nSamples explained: %d\nFeatures: %d\nTop feature by mean |SHAP|: %s (%.4f)\n",
  upstream_jobid, model_type, n_explain, ncol(X_sub),
  shap_global$feature[1], shap_global$mean_abs_shap[1]
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = list(
    list(name = "shap_summary.csv", type = "table", label = "SHAP 全局排序"),
    list(name = "shap_global_csv.csv", type = "table", label = "每样本 SHAP 矩阵"),
    list(name = "global_importance.png", type = "plot", label = "全局 SHAP 重要性"),
    list(name = "global_importance.pdf", type = "file", label = "全局 SHAP 重要性 (PDF)"),
    list(name = "beeswarm.png", type = "plot", label = "SHAP beeswarm"),
    list(name = "beeswarm.pdf", type = "file", label = "SHAP beeswarm (PDF)"),
    list(name = "force_plot.png", type = "plot", label = "单样本 force plot"),
    list(name = "force_plot.pdf", type = "file", label = "单样本 force plot (PDF)"),
    list(name = "dependence.png", type = "plot", label = "Top 4 依赖图"),
    list(name = "dependence.pdf", type = "file", label = "Top 4 依赖图 (PDF)"),
    list(name = "summary.txt", type = "text", label = "分析摘要")
  ),
  summary = list(
    upstream_job = upstream_jobid,
    model_type = model_type,
    n_samples = n_explain,
    n_features = ncol(X_sub),
    top_feature = shap_global$feature[1]
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

cat("[shap-explain] Pipeline complete.\n")
