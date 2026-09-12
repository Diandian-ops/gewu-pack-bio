#!/usr/bin/env Rscript
# ============================================================
# Random Survival Forest (RSF) — Non-linear prognostic modeling
#
# Outputs:
#   - C-index (train / OOB / test)
#   - Time-dependent ROC at user-specified time points
#   - VIMP (Variable Importance) ranking
#   - Partial dependence plots (top N features)
#   - Risk score triple plot (KM + scatter + heatmap)
#   - Trained model (.rds)
#
# BioF3 SCI visual style: theme_biof3(), biof3_palette(), ggsave_biof3()
# ============================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(survival)
})

need_pkgs <- c("randomForestSRC", "timeROC", "pec")
# BLOCK-58 — 依赖执行路径零安装：缺包毫秒级 fail-fast（结构化错误，绝不联网）。
# 依赖单一事实源：_dependencies.json；显式补装：scripts/install-r-tool-dependencies.mjs
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("MISSING_R_PACKAGES: %s\n", paste(missing_pkgs, collapse = ", ")))
  cat("Install first via: node scripts/install-r-tool-dependencies.mjs --tool rsf-survival\n")
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

cat("[rsf-survival] Starting RSF pipeline (SCI style)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
ntree <- as.integer(params$ntree %||% 500)
mtry_in <- params$mtry %||% "auto"
time_points_str <- params$time_points %||% "1,3,5"
top_n_features <- as.integer(params$top_n_features %||% 20)

# Parse time points (years)
time_points <- as.numeric(strsplit(time_points_str, ",")[[1]])

# Read data
read_matrix <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f)) read.delim(f, row.names = 1, check.names = FALSE)
  else read.csv(f, row.names = 1, check.names = FALSE)
}

features_file <- list.files(job_dir, pattern = "^features", full.names = TRUE)[1]
clinical_file <- list.files(job_dir, pattern = "^clinical", full.names = TRUE)[1]

features <- as.matrix(read_matrix(features_file))
clinical <- read_matrix(clinical_file)

# Standardize clinical column names (case-insensitive)
clinical_names <- tolower(colnames(clinical))
cols <- list(
  sample = grep("sample", clinical_names, value = TRUE)[1],
  time = grep("time|os_time|time_to_event|survival_time", clinical_names, value = TRUE)[1],
  event = grep("event|status|os_event|death|vital_status", clinical_names, value = TRUE)[1]
)

# If not found by pattern, use position
if (is.na(cols$sample)) cols$sample <- colnames(clinical)[1]
if (is.na(cols$time)) cols$time <- colnames(clinical)[2]
if (is.na(cols$event)) cols$event <- colnames(clinical)[3]

sample_ids <- as.character(clinical[[cols$sample]])
os_time <- as.numeric(clinical[[cols$time]])
os_event <- as.numeric(clinical[[cols$event]])

# Intersect samples
common <- intersect(colnames(features), sample_ids)
features <- features[, common, drop = FALSE]
idx <- match(common, sample_ids)
os_time <- os_time[idx]
os_event <- os_event[idx]

# Transpose: rows = samples, cols = features
X <- t(features)
var_x <- apply(X, 2, var, na.rm = TRUE)
X <- X[, var_x > 0, drop = FALSE]

# Auto-detect time unit: median > 30 → days → convert to years
if (median(os_time, na.rm = TRUE) > 30) {
  os_time <- os_time / 365.25
  cat("[rsf-survival] Time unit auto-detected: days → converted to years\n")
} else {
  cat("[rsf-survival] Time unit auto-detected: years\n")
}

# Remove NA
valid <- !is.na(os_time) & !is.na(os_event) & os_time > 0
X <- X[valid, , drop = FALSE]
os_time <- os_time[valid]
os_event <- os_event[valid]

cat(sprintf("[rsf-survival] Samples: %d, Features: %d, Events: %d\n", nrow(X), ncol(X), sum(os_event)))

# ============================================================
# Build survival data frame
# ============================================================
surv_df <- as.data.frame(X)
surv_df$time <- os_time
surv_df$event <- os_event

# mtry setting
p <- ncol(X)
if (mtry_in == "auto") {
  mtry <- max(1, floor(sqrt(p)))
} else {
  mtry <- as.integer(mtry_in)
}

# ============================================================
# Train RSF
# ============================================================
cat(sprintf("[rsf-survival] Training RSF (ntree=%d, mtry=%d)...\n", ntree, mtry))
suppressMessages(library(randomForestSRC))

set.seed(42)
rsf_fit <- rfsrc(Surv(time, event) ~ ., data = surv_df, ntree = ntree, mtry = mtry,
                 importance = "permute", splitrule = "logrank", 
                 cause = 1, verbose = FALSE)

cat(sprintf("[rsf-survival] OOB Error: %.4f\n", rsf_fit$err.rate[ntree]))

# ============================================================
# C-index (OOB + train)
# ============================================================
cat("[rsf-survival] Computing C-index...\n")

# OOB C-index
pred_oob <- rsf_fit$predicted.oob
cindex_oob <- survConcordance(Surv(time, event) ~ pred_oob, data = surv_df)$concordance

# Train C-index (in-bag)
pred_train <- predict(rsf_fit, newdata = surv_df, importance = "none")$predicted
cindex_train <- survConcordance(Surv(time, event) ~ pred_train, data = surv_df)$concordance

# Test C-index (if we do a quick 80/20 split)
set.seed(123)
train_idx <- sample(nrow(surv_df), floor(0.8 * nrow(surv_df)))
test_df <- surv_df[-train_idx, , drop = FALSE]
rsf_test <- predict(rsf_fit, newdata = test_df, importance = "none")
pred_test <- rsf_test$predicted
cindex_test <- tryCatch(
  survConcordance(Surv(time, event) ~ pred_test, data = test_df)$concordance,
  error = function(e) NA
)

cindex_df <- data.frame(
  Set = c("Train (in-bag)", "OOB", "Test (20% holdout)"),
  C_index = c(round(cindex_train, 4), round(cindex_oob, 4), round(ifelse(is.na(cindex_test), NA, cindex_test), 4)),
  N = c(length(train_idx), length(train_idx), nrow(test_df))
)
write.csv(cindex_df, file.path(output_dir, "cindex_table.csv"), row.names = FALSE)

# ============================================================
# VIMP (Variable Importance)
# ============================================================
cat("[rsf-survival] Extracting VIMP...\n")
vimp_raw <- rsf_fit$importance
vimp_df <- data.frame(
  feature = names(vimp_raw),
  vimp = as.numeric(vimp_raw),
  stringsAsFactors = FALSE
)
vimp_df <- vimp_df[order(-vimp_df$vimp), ]
vimp_df$rank <- seq_len(nrow(vimp_df))
write.csv(vimp_df, file.path(output_dir, "vimp_table.csv"), row.names = FALSE)

# VIMP plot (top N)
vimp_top <- head(vimp_df, top_n_features)
if (nrow(vimp_top) > 0) {
  p_vimp <- ggplot(vimp_top, aes(x = reorder(feature, vimp), y = vimp)) +
    geom_col(fill = "#0f766e", width = 0.7) +
    labs(title = sprintf("VIMP (top %d)", top_n_features), x = "Feature", y = "Permutation Importance") +
    coord_flip() +
    theme_biof3()
  ggsave_biof3(p_vimp, file.path(output_dir, "vimp"), width = 6, height = 7)
}

# ============================================================
# Time-dependent ROC
# ============================================================
cat("[rsf-survival] Computing time-dependent ROC...\n")
suppressMessages(library(timeROC))

roc_results <- list()
for (tp in time_points) {
  roc_obj <- timeROC(T = surv_df$time, delta = surv_df$event,
                     marker = pred_train, cause = 1,
                     times = tp, iid = FALSE, ROC = TRUE)
  auc_val <- roc_obj$AUC[1]
  roc_results[[as.character(tp)]] <- list(
    time_point = tp,
    auc = auc_val,
    tp = roc_obj$TP,
    fp = roc_obj$FP
  )
}

# Plot time-dependent ROC
roc_data <- do.call(rbind, lapply(names(roc_results), function(nm) {
  r <- roc_results[[nm]]
  data.frame(FPR = r$fp, TPR = r$tp, time = paste0(r$time_point, "yr"), stringsAsFactors = FALSE)
}))

roc_pal <- biof3_palette(length(time_points))
p_troc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = time)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = setNames(roc_pal, paste0(time_points, "yr"))) +
  labs(title = "Time-dependent ROC", x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = "Time") +
  coord_equal() +
  theme_biof3()
ggsave_biof3(p_troc, file.path(output_dir, "time_roc"), width = 6, height = 5)

# ROC summary table
roc_summary <- data.frame(
  Time_point_years = sapply(roc_results, function(r) r$time_point),
  AUC = sapply(roc_results, function(r) round(r$auc, 4)),
  stringsAsFactors = FALSE
)
write.csv(roc_summary, file.path(output_dir, "time_roc_table.csv"), row.names = FALSE)

# ============================================================
# Partial Dependence Plots (top 4 features)
# ============================================================
cat("[rsf-survival] Generating partial dependence plots...\n")
top4 <- head(vimp_df$feature, 4)

pd_plots <- list()
for (feat in top4) {
  pd <- partial(rsf_fit, partial.xvar = feat, partial.type = "mort",
                partial.values = quantile(surv_df[[feat]], probs = seq(0.05, 0.95, 0.05)),
                seed = 42)
  pd_df <- data.frame(
    value = pd$partial.values,
    pred = pd$predicted
  )
  p_pd <- ggplot(pd_df, aes(x = value, y = pred)) +
    geom_line(color = "#0f766e", linewidth = 0.8) +
    geom_rug(sides = "b", alpha = 0.3) +
    labs(title = feat, x = "Feature value", y = "Predicted mortality") +
    theme_biof3()
  pd_plots[[feat]] <- p_pd
}

# Combine 4 PD plots into one figure
library(gridExtra)
save_grid_biof3(file.path(output_dir, "partial_dependence"), width = 10, height = 8, expr = {
  grid.arrange(grobs = lapply(pd_plots, ggplotGrob), ncol = 2,
               top = textGrob("Partial Dependence (top 4 VIMP features)", gp = gpar(fontsize = 14, fontface = "bold")))
})

# ============================================================
# Risk Score Triple Plot
# ============================================================
cat("[rsf-survival] Risk score triple plot...\n")

# Risk score = predicted mortality
risk_score <- pred_train
risk_df <- data.frame(
  sample = seq_along(risk_score),
  risk = risk_score,
  time = surv_df$time,
  event = surv_df$event,
  stringsAsFactors = FALSE
)
risk_df$group <- ifelse(risk_df$risk > median(risk_df$risk), "High", "Low")

# KM curves by risk group
suppressMessages(library(survminer))
fit <- survfit(Surv(time, event) ~ group, data = risk_df)
p_km <- ggsurvplot(fit, data = risk_df, pval = TRUE, risk.table = TRUE,
                   palette = c("#ef4444", "#10b981"),
                   title = "KM by Risk Group",
                   xlab = "Time (years)", ylab = "Survival probability")

# Risk scatter
p_scatter <- ggplot(risk_df, aes(x = reorder(sample, risk), y = risk, color = group)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(High = "#ef4444", Low = "#10b981")) +
  labs(title = "Risk Score Distribution", x = "Sample (ordered)", y = "Risk score") +
  theme_biof3() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Risk heatmap (sample x top 10 features, colored by expression)
top10 <- head(vimp_df$feature, 10)
if (length(top10) > 0) {
  heat_mat <- scale(X[, top10, drop = FALSE])
  heat_df <- as.data.frame(heat_mat)
  heat_df$sample <- seq_len(nrow(heat_mat))
  heat_df$risk <- risk_score
  heat_df <- heat_df[order(heat_df$risk), ]
  heat_long <- reshape(as.data.frame(heat_df[, c("sample", top10)]),
                       direction = "long", varying = top10, v.names = "expr", times = top10, idvar = "sample")
  names(heat_long)[3] <- "feature"

  p_heat <- ggplot(heat_long, aes(x = sample, y = feature, fill = expr)) +
    geom_tile() +
    scale_fill_gradient2(low = "#3b82f6", mid = "white", high = "#ef4444", midpoint = 0) +
    labs(title = "Top 10 Feature Expression (by risk)", x = "Sample (ordered by risk)", y = NULL, fill = "Z-score") +
    theme_biof3() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# Save combined triple plot
save_grid_biof3(file.path(output_dir, "risk_triple"), width = 12, height = 10, expr = {
  if (exists("p_km") && exists("p_scatter") && exists("p_heat")) {
    km_plot <- p_km$plot + theme(legend.position = "bottom")
    grid.arrange(
      arrangeGrob(km_plot, p_scatter, ncol = 2),
      p_heat, ncol = 1,
      heights = c(1, 1),
      top = textGrob("Risk Score Triple Plot", gp = gpar(fontsize = 14, fontface = "bold"))
    )
  }
})

# ============================================================
# Save model
# ============================================================
saveRDS(rsf_fit, file.path(output_dir, "model.rds"))

# ============================================================
# Summary
# ============================================================
summary_text <- sprintf(
  "Random Survival Forest (RSF) Prognostic Model\n\nSamples: %d, Features: %d, Events: %d\nTime unit: years\nntree: %d, mtry: %d\n\nC-index:\n  Train (in-bag): %.3f\n  OOB: %.3f\n  Test (20%% holdout): %.3f\n\nTime-dependent AUC:\n%s\n\nTop 5 VIMP features:\n%s\n",
  nrow(X), ncol(X), sum(os_event), ntree, mtry,
  cindex_train, cindex_oob, ifelse(is.na(cindex_test), NA, cindex_test),
  paste(sapply(roc_results, function(r) sprintf("  %.1f yr: %.3f", r$time_point, r$auc)), collapse = "\n"),
  paste(head(sprintf("  %s: %.4f", vimp_df$feature, vimp_df$vimp), 5), collapse = "\n")
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = list(
    list(name = "cindex_table.csv", type = "table", label = "C-index 表"),
    list(name = "time_roc_table.csv", type = "table", label = "时间 ROC 表"),
    list(name = "vimp_table.csv", type = "table", label = "VIMP 排序表"),
    list(name = "time_roc.png", type = "plot", label = "时间依赖 ROC"),
    list(name = "time_roc.pdf", type = "file", label = "时间依赖 ROC (PDF)"),
    list(name = "vimp.png", type = "plot", label = "VIMP 排序"),
    list(name = "vimp.pdf", type = "file", label = "VIMP 排序 (PDF)"),
    list(name = "partial_dependence.png", type = "plot", label = "部分依赖图"),
    list(name = "partial_dependence.pdf", type = "file", label = "部分依赖图 (PDF)"),
    list(name = "risk_triple.png", type = "plot", label = "Risk 三联图"),
    list(name = "risk_triple.pdf", type = "file", label = "Risk 三联图 (PDF)"),
    list(name = "model.rds", type = "file", label = "RSF 模型 (.rds)"),
    list(name = "summary.txt", type = "text", label = "分析摘要")
  ),
  summary = list(
    samples = nrow(X),
    features = ncol(X),
    events = sum(os_event),
    cindex_oob = round(cindex_oob, 4)
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

cat("[rsf-survival] Pipeline complete.\n")
