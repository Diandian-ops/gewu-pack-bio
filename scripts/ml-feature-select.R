#!/usr/bin/env Rscript
# ============================================================
# ML Feature Selection — 5 methods comparison
#   1. Variance filter (top variance)
#   2. Univariate Wilcoxon / t-test
#   3. mRMR (max relevance / min redundancy)
#   4. Boruta (wrapper, RF-based)
#   5. RFE (recursive feature elimination, LR-based)
#
# Output: 5 ranking tables + Venn consensus + method correlation heatmap
#         + volcano plot (univariate) + consensus top-N table + summary
#
# BioF3 SCI visual style: theme_biof3(), biof3_palette(), ggsave_biof3()
# ============================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(VennDiagram)
})

# Optional packages — install if missing
need_pkgs <- c("caret", "Boruta", "mRMRe", "glmnet", "e1071")
# BLOCK-58 — 依赖执行路径零安装：缺包毫秒级 fail-fast（结构化错误，绝不联网）。
# 依赖单一事实源：_dependencies.json；显式补装：scripts/install-r-tool-dependencies.mjs
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("MISSING_R_PACKAGES: %s\n", paste(missing_pkgs, collapse = ", ")))
  cat("Install first via: node scripts/install-r-tool-dependencies.mjs --tool ml-feature-select\n")
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

cat("[ml-feature-select] Starting 5-method comparison (SCI style)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
top_n <- as.integer(params$top_n %||% 50)

# Read data
read_matrix <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f)) read.delim(f, row.names = 1, check.names = FALSE)
  else read.csv(f, row.names = 1, check.names = FALSE)
}

features_file <- list.files(job_dir, pattern = "^features", full.names = TRUE)[1]
labels_file   <- list.files(job_dir, pattern = "^labels",   full.names = TRUE)[1]

features <- as.matrix(read_matrix(features_file))
labels_df <- read_matrix(labels_file)

# Ensure labels is a data.frame with sample + label columns
if (ncol(labels_df) >= 2) {
  sample_ids <- rownames(labels_df)
  labels_vec <- labels_df[[2]]
} else {
  sample_ids <- rownames(labels_df)
  labels_vec <- labels_df[[1]]
}
names(labels_vec) <- sample_ids

# Intersect samples
common_samples <- intersect(colnames(features), sample_ids)
features <- features[, common_samples, drop = FALSE]
labels_vec <- labels_vec[common_samples]

# Ensure binary label (0/1 or factor with 2 levels)
labels_vec <- as.factor(labels_vec)
if (length(levels(labels_vec)) != 2) {
  stop("Labels must be binary (2 classes). Found: ", length(levels(labels_vec)), " classes.")
}

# Transpose: rows = samples, cols = features
X <- t(features)
# Remove constant features
var_x <- apply(X, 2, var, na.rm = TRUE)
X <- X[, var_x > 0, drop = FALSE]

cat(sprintf("[ml-feature-select] Samples: %d, Features: %d\n", nrow(X), ncol(X)))

# ============================================================
# Method 1: Variance Filter
# ============================================================
cat("[ml-feature-select] Method 1: Variance filter...\n")
var_scores <- apply(X, 2, var, na.rm = TRUE)
var_rank <- data.frame(
  feature = names(var_scores),
  score = var_scores,
  rank = rank(-var_scores, ties.method = "min")
)
var_rank <- var_rank[order(var_rank$rank), ]
write.csv(var_rank, file.path(output_dir, "variance_table.csv"), row.names = FALSE)

# ============================================================
# Method 2: Univariate Wilcoxon / t-test
# ============================================================
cat("[ml-feature-select] Method 2: Univariate test...\n")
g1 <- labels_vec == levels(labels_vec)[1]
g2 <- labels_vec == levels(labels_vec)[2]

uni_p <- sapply(seq_len(ncol(X)), function(j) {
  x1 <- X[g1, j]
  x2 <- X[g2, j]
  # Use t-test if both groups look normal-ish, else Wilcoxon
  if (length(x1) > 3 && length(x2) > 3) {
    tryCatch(t.test(x1, x2)$p.value, error = function(e) wilcox.test(x1, x2)$p.value)
  } else {
    wilcox.test(x1, x2)$p.value
  }
})
# Fold change (log2 mean ratio)
uni_fc <- sapply(seq_len(ncol(X)), function(j) {
  m1 <- mean(X[g1, j], na.rm = TRUE)
  m2 <- mean(X[g2, j], na.rm = TRUE)
  ifelse(m1 > 0 && m2 > 0, log2(m2 / m1), NA)
})

uni_df <- data.frame(
  feature = colnames(X),
  p_value = uni_p,
  log2FC = uni_fc,
  neg_log10_p = -log10(uni_p)
)
uni_df$padj <- p.adjust(uni_df$p_value, method = "BH")
uni_df$rank <- rank(uni_df$p_value, ties.method = "min")
uni_df <- uni_df[order(uni_df$rank), ]
write.csv(uni_df, file.path(output_dir, "univariate_table.csv"), row.names = FALSE)

# Volcano plot (univariate)
sig_thresh <- 0.05
fc_thresh <- 1
uni_df$sig <- ifelse(uni_df$padj < sig_thresh && abs(uni_df$log2FC) > fc_thresh, "Significant",
              ifelse(uni_df$padj < sig_thresh, "padj<0.05", "NS"))
top_volc <- head(uni_df[uni_df$sig == "Significant", ], 15)

volc_pal <- c(NS = "gray70", "padj<0.05" = "#fbbf24", Significant = "#ef4444")
p_volc <- ggplot(uni_df, aes(log2FC, neg_log10_p, color = sig)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_text_repel(data = top_volc, aes(label = feature), size = 3, max.overlaps = 15, color = "black") +
  scale_color_manual(values = volc_pal) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(sig_thresh), linetype = "dashed", color = "gray40") +
  labs(title = "Univariate Test Volcano", x = "log2 Fold Change", y = "-log10(p)", color = NULL) +
  theme_biof3()
ggsave_biof3(p_volc, file.path(output_dir, "volcano"), width = 6, height = 5)

# ============================================================
# Method 3: mRMR (max relevance min redundancy)
# ============================================================
cat("[ml-feature-select] Method 3: mRMR...\n")
mrmr_rank <- NULL
tryCatch({
  suppressMessages(library(mRMRe))
  # mRMRe needs a data.frame with target as last column
  mrmr_df <- as.data.frame(X)
  mrmr_df$label <- as.numeric(labels_vec) - 1  # 0/1

  # Use mRMR.ensemble for robustness; if too few samples, fall back to classic
  if (nrow(mrmr_df) >= 10) {
    obj <- mRMR.ensemble(data = mrmr_df, target_indices = ncol(mrmr_df),
                         solution_count = 1, feature_count = min(50, ncol(X)))
    sel <- as.vector(solutions(obj)[[1]])
    mrmr_rank <- data.frame(
      feature = colnames(X)[sel],
      rank = seq_along(sel)
    )
  }
}, error = function(e) {
  cat("[ml-feature-select] mRMR failed (", conditionMessage(e), "), using Spearman correlation fallback.\n")
})

if (is.null(mrmr_rank)) {
  # Fallback: rank by absolute Spearman correlation with label
  cors <- apply(X, 2, function(x) abs(cor(x, as.numeric(labels_vec), method = "spearman", use = "complete.obs")))
  cors[is.na(cors)] <- 0
  mrmr_rank <- data.frame(
    feature = names(cors),
    rank = rank(-cors, ties.method = "min")
  )
  mrmr_rank <- mrmr_rank[order(mrmr_rank$rank), ]
}
write.csv(mrmr_rank, file.path(output_dir, "mrmr_table.csv"), row.names = FALSE)

# ============================================================
# Method 4: Boruta (wrapper, RF-based)
# ============================================================
cat("[ml-feature-select] Method 4: Boruta...\n")
boruta_rank <- NULL
tryCatch({
  suppressMessages(library(Boruta))
  bor_df <- as.data.frame(X)
  bor_df$label <- labels_vec

  # Limit features for speed if too many
  if (ncol(X) > 5000) {
    pre_var <- apply(X, 2, var, na.rm = TRUE)
    keep_pre <- names(sort(pre_var, decreasing = TRUE))[1:5000]
    bor_df <- bor_df[, c(keep_pre, "label")]
  }

  bor_res <- Boruta(label ~ ., data = bor_df, doTrace = 0, maxRuns = 100)
  bor_imp <- attStats(bor_res)
  bor_imp$feature <- rownames(bor_imp)
  # Confirmed > Tentative > Rejected; within same decision, sort by meanImp desc
  bor_imp$decision_score <- ifelse(bor_imp$decision == "Confirmed", 3,
                            ifelse(bor_imp$decision == "Tentative", 2, 1))
  bor_imp$rank <- rank(-bor_imp$decision_score * 1e6 - bor_imp$meanImp, ties.method = "min")
  boruta_rank <- bor_imp[, c("feature", "decision", "meanImp", "rank")]
  boruta_rank <- boruta_rank[order(boruta_rank$rank), ]
}, error = function(e) {
  cat("[ml-feature-select] Boruta failed (", conditionMessage(e), "), using RF importance fallback.\n")
})

if (is.null(boruta_rank)) {
  # Fallback: Random Forest importance via ranger (lightweight)
  # BLOCK-58 — ranger 缺包时由上方 MISSING_R_PACKAGES fail-fast 统一拦截（不在此联网安装）
  suppressMessages(library(ranger))
  rf_df <- as.data.frame(X)
  rf_df$label <- labels_vec
  rf_fit <- ranger(label ~ ., data = rf_df, importance = "impurity", num.trees = 500, verbose = FALSE)
  imp <- rf_fit$variable.importance
  boruta_rank <- data.frame(
    feature = names(imp),
    decision = "Confirmed",
    meanImp = imp,
    rank = rank(-imp, ties.method = "min")
  )
  boruta_rank <- boruta_rank[order(boruta_rank$rank), ]
}
write.csv(boruta_rank, file.path(output_dir, "boruta_table.csv"), row.names = FALSE)

# ============================================================
# Method 5: RFE (recursive feature elimination, LR-based)
# ============================================================
cat("[ml-feature-select] Method 5: RFE...\n")
rfe_rank <- NULL
tryCatch({
  suppressMessages(library(caret))
  # Use glmnet (LASSO) as the base model for RFE — fast and handles high-dim well
  ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 3, verbose = FALSE)

  # For speed, subsample if > 2000 features
  if (ncol(X) > 2000) {
    pre_var <- apply(X, 2, var, na.rm = TRUE)
    keep_pre <- names(sort(pre_var, decreasing = TRUE))[1:2000]
    X_rfe <- X[, keep_pre, drop = FALSE]
  } else {
    X_rfe <- X
  }

  profile <- rfe(X_rfe, labels_vec, sizes = c(10, 20, 50, 100),
                 rfeControl = ctrl,
                 method = "glmnet",
                 trControl = trainControl(method = "none"),
                 metric = "Accuracy")

  # Extract variable importance from final model
  imp_rfe <- varImp(profile, scale = FALSE)
  rfe_rank <- data.frame(
    feature = rownames(imp_rfe),
    overall = imp_rfe$Overall,
    rank = rank(-imp_rfe$Overall, ties.method = "min")
  )
  rfe_rank <- rfe_rank[order(rfe_rank$rank), ]
}, error = function(e) {
  cat("[ml-feature-select] RFE failed (", conditionMessage(e), "), using LASSO coefficient fallback.\n")
})

if (is.null(rfe_rank)) {
  # Fallback: LASSO coefficients as importance
  suppressMessages(library(glmnet))
  cv_fit <- cv.glmnet(X, as.numeric(labels_vec) - 1, family = "binomial", nfolds = 5)
  coefs <- as.matrix(coef(cv_fit, s = "lambda.1se"))[-1, , drop = TRUE]  # drop intercept
  coefs <- abs(coefs)
  coefs[is.na(coefs)] <- 0
  rfe_rank <- data.frame(
    feature = names(coefs),
    overall = coefs,
    rank = rank(-coefs, ties.method = "min")
  )
  rfe_rank <- rfe_rank[order(rfe_rank$rank), ]
}
write.csv(rfe_rank, file.path(output_dir, "rfe_table.csv"), row.names = FALSE)

# ============================================================
# Consensus Top-N (5-method voting)
# ============================================================
cat("[ml-feature-select] Computing consensus top-", top_n, "...\n")
get_top <- function(df, n) {
  head(df$feature[order(df$rank)], n)
}

sets <- list(
  Variance = get_top(var_rank, top_n),
  Univariate = get_top(uni_df, top_n),
  mRMR = get_top(mrmr_rank, top_n),
  Boruta = get_top(boruta_rank, top_n),
  RFE = get_top(rfe_rank, top_n)
)

# Count votes
all_genes <- unique(unlist(sets))
votes <- sapply(all_genes, function(g) sum(sapply(sets, function(s) g %in% s)))
consensus_df <- data.frame(
  feature = all_genes,
  votes = as.integer(votes),
  stringsAsFactors = FALSE
)
consensus_df <- consensus_df[order(-consensus_df$votes, consensus_df$feature), ]
consensus_df$rank <- seq_len(nrow(consensus_df))
consensus_top <- head(consensus_df, top_n)
write.csv(consensus_df, file.path(output_dir, "consensus_table.csv"), row.names = FALSE)
write.csv(consensus_top, file.path(output_dir, "consensus_top_n.csv"), row.names = FALSE)

# ============================================================
# Venn Diagram (5-method overlap)
# ============================================================
cat("[ml-feature-select] Plotting Venn diagram...\n")
venn_pal <- biof3_palette(5)
venn_obj <- venn.diagram(
  x = sets,
  filename = NULL,
  fill = venn_pal,
  alpha = 0.6,
  cex = 0.8,
  cat.cex = 0.7,
  cat.dist = c(0.15, 0.15, 0.12, 0.12, 0.1),
  margin = 0.08,
  main = sprintf("Feature Selection Consensus (top %d)", top_n),
  main.cex = 1.2,
  fontfamily = "sans",
  cat.fontfamily = "sans"
)

save_grid_biof3(file.path(output_dir, "venn"), width = 8, height = 7, expr = {
  grid.draw(venn_obj)
})

# ============================================================
# Method correlation heatmap (how similar are the rankings?)
# ============================================================
cat("[ml-feature-select] Method correlation heatmap...\n")
# Create rank vectors for all features
methods <- c("Variance", "Univariate", "mRMR", "Boruta", "RFE")
all_feat <- colnames(X)

rank_mat <- matrix(NA, nrow = length(all_feat), ncol = 5,
                   dimnames = list(all_feat, methods))

fill_ranks <- function(df, col_name = "rank") {
  r <- setNames(df[[col_name]], df$feature)
  # unranked features get max rank + 1
  default <- max(r, na.rm = TRUE) + 1
  sapply(all_feat, function(f) ifelse(is.na(r[f]), default, r[f]))
}

rank_mat[, "Variance"]    <- fill_ranks(var_rank)
rank_mat[, "Univariate"]  <- fill_ranks(uni_df)
rank_mat[, "mRMR"]        <- fill_ranks(mrmr_rank)
rank_mat[, "Boruta"]      <- fill_ranks(boruta_rank)
rank_mat[, "RFE"]         <- fill_ranks(rfe_rank)

# Spearman correlation of rankings
cor_mat <- cor(rank_mat, method = "spearman", use = "pairwise.complete.obs")
cor_df <- as.data.frame(cor_mat)
cor_df$method <- rownames(cor_df)
write.csv(cor_df, file.path(output_dir, "method_correlation.csv"), row.names = FALSE)

# Heatmap
hm_pal <- biof3_palette_seq(100, option = "ylgn")
save_grid_biof3(file.path(output_dir, "method_correlation"), width = 5.5, height = 4.5, expr = {
  pheatmap(cor_mat,
           color = hm_pal,
           border_color = NA,
           display_numbers = TRUE,
           number_color = "black",
           fontsize = 10,
           main = "Method Rank Correlation (Spearman)")
})

# ============================================================
# Summary
# ============================================================
summary_text <- sprintf(
  "ML Feature Selection — 5-Method Comparison\n\nSamples: %d\nFeatures: %d\nTop-N consensus: %d\n\nMethod coverage (top %d):\n  Variance:    %d\n  Univariate:  %d\n  mRMR:        %d\n  Boruta:      %d\n  RFE:         %d\n\nConsensus distribution:\n  5/5 methods: %d\n  4/5 methods: %d\n  3/5 methods: %d\n  2/5 methods: %d\n  1/5 methods: %d\n",
  nrow(X), ncol(X), top_n, top_n,
  length(sets$Variance), length(sets$Univariate), length(sets$mRMR),
  length(sets$Boruta), length(sets$RFE),
  sum(votes == 5), sum(votes == 4), sum(votes == 3),
  sum(votes == 2), sum(votes == 1)
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = list(
    list(name = "variance_table.csv",     type = "table", label = "方差排序表"),
    list(name = "univariate_table.csv",  type = "table", label = "单因素 Wilcoxon 排序表"),
    list(name = "mrmr_table.csv",        type = "table", label = "mRMR 排序表"),
    list(name = "boruta_table.csv",      type = "table", label = "Boruta 排序表"),
    list(name = "rfe_table.csv",         type = "table", label = "RFE 排序表"),
    list(name = "consensus_table.csv",   type = "table", label = "共识投票表"),
    list(name = "consensus_top_n.csv",   type = "table", label = "共识 Top N 基因"),
    list(name = "method_correlation.csv",type = "table", label = "方法间相关性"),
    list(name = "venn.png",              type = "plot",  label = "5法韦恩图"),
    list(name = "venn.pdf",              type = "file",  label = "5法韦恩图 (PDF)"),
    list(name = "volcano.png",           type = "plot",  label = "单因素火山图"),
    list(name = "volcano.pdf",           type = "file",  label = "单因素火山图 (PDF)"),
    list(name = "method_correlation.png",type = "plot",  label = "方法间相关性热图"),
    list(name = "method_correlation.pdf",type = "file",  label = "方法间相关性热图 (PDF)"),
    list(name = "summary.txt",           type = "text",  label = "分析摘要")
  ),
  summary = list(
    samples = nrow(X),
    features = ncol(X),
    top_n = top_n,
    consensus_5 = sum(votes == 5),
    consensus_4 = sum(votes == 4)
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

cat("[ml-feature-select] Pipeline complete.\n")
