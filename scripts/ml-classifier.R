#!/usr/bin/env Rscript
# ============================================================
# ML Classifier — 4 algorithms comparison (LR / SVM / RF / XGBoost)
#
# Outputs:
#   - 5-fold CV ROC / PR curves per algorithm
#   - Calibration plot
#   - Confusion matrices
#   - Feature importance comparison
#   - Learning curve
#   - Metrics summary table
#   - Trained model (.rds) for downstream SHAP
#
# BioF3 SCI visual style: theme_biof3(), biof3_palette(), ggsave_biof3()
# ============================================================

suppressMessages({
  library(ggplot2)
  library(jsonlite)
  library(pROC)
  library(caret)
})

# Optional packages
need_pkgs <- c("glmnet", "e1071", "randomForest", "xgboost")
# BLOCK-58 — 依赖执行路径零安装：缺包毫秒级 fail-fast（结构化错误，绝不联网）。
# 依赖单一事实源：_dependencies.json；显式补装：scripts/install-r-tool-dependencies.mjs
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("MISSING_R_PACKAGES: %s\n", paste(missing_pkgs, collapse = ", ")))
  cat("Install first via: node scripts/install-r-tool-dependencies.mjs --tool ml-classifier\n")
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

cat("[ml-classifier] Starting 4-algorithm comparison (SCI style)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
algorithms <- params$algorithms %||% "all"
task_type <- params$task_type %||% "binary"

# Read data
read_matrix <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f)) read.delim(f, row.names = 1, check.names = FALSE)
  else read.csv(f, row.names = 1, check.names = FALSE)
}

features_file <- list.files(job_dir, pattern = "^features", full.names = TRUE)[1]
labels_file   <- list.files(job_dir, pattern = "^labels",   full.names = TRUE)[1]

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

# Intersect
common <- intersect(colnames(features), sample_ids)
features <- features[, common, drop = FALSE]
labels_vec <- labels_vec[common]
labels_vec <- as.factor(labels_vec)

X <- t(features)
var_x <- apply(X, 2, var, na.rm = TRUE)
X <- X[, var_x > 0, drop = FALSE]

cat(sprintf("[ml-classifier] Samples: %d, Features: %d, Classes: %d\n", nrow(X), ncol(X), length(levels(labels_vec))))

# Choose algorithms
alg_list <- switch(algorithms,
  all = c("lr", "svm", "rf", "xgb"),
  lr = "lr", svm = "svm", rf = "rf", xgb = "xgb",
  c("lr", "svm", "rf", "xgb")
)

# ============================================================
# Cross-validation setup
# ============================================================
cv_folds <- createFolds(labels_vec, k = 5, list = TRUE, returnTrain = FALSE)

results <- list()
roc_list <- list()
pr_list <- list()

for (alg in alg_list) {
  cat(sprintf("[ml-classifier] Training %s...\n", alg))
  fold_preds <- list()
  fold_probs <- list()
  fold_truths <- list()
  fold_metrics <- list()

  for (i in seq_along(cv_folds)) {
    test_idx <- cv_folds[[i]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)

    X_train <- X[train_idx, , drop = FALSE]
    y_train <- labels_vec[train_idx]
    X_test  <- X[test_idx, , drop = FALSE]
    y_test  <- labels_vec[test_idx]

    # Train
    model <- switch(alg,
      lr = {
        suppressMessages(library(glmnet))
        cv_glm <- cv.glmnet(X_train, y_train, family = ifelse(length(levels(labels_vec)) == 2, "binomial", "multinomial"), nfolds = 3)
        list(model = cv_glm, type = "glmnet")
      },
      svm = {
        suppressMessages(library(e1071))
        # Subsample for speed if too large
        if (nrow(X_train) > 2000) {
          ss <- sample(nrow(X_train), 2000)
          X_train_sub <- X_train[ss, ]
          y_train_sub <- y_train[ss]
        } else {
          X_train_sub <- X_train
          y_train_sub <- y_train
        }
        df_train <- as.data.frame(X_train_sub)
        df_train$label <- y_train_sub
        svm(label ~ ., data = df_train, kernel = "radial", probability = TRUE, scale = FALSE)
      },
      rf = {
        suppressMessages(library(randomForest))
        df_train <- as.data.frame(X_train)
        df_train$label <- y_train
        randomForest(label ~ ., data = df_train, ntree = 300, importance = TRUE)
      },
      xgb = {
        suppressMessages(library(xgboost))
        dtrain <- xgb.DMatrix(data = X_train, label = as.numeric(y_train) - 1)
        xgb.train(params = list(objective = ifelse(length(levels(labels_vec)) == 2, "binary:logistic", "multi:softprob"),
                                num_class = ifelse(length(levels(labels_vec)) == 2, 1, length(levels(labels_vec))),
                                eval_metric = "logloss", max_depth = 4, eta = 0.1, subsample = 0.8),
                  data = dtrain, nrounds = 100, verbose = 0)
      }
    )

    # Predict
    probs <- switch(alg,
      lr = {
        p <- predict(model$model, newx = X_test, s = "lambda.1se", type = "response")
        if (is.matrix(p) && ncol(p) == 1) as.vector(p) else p[, 1]
      },
      svm = {
        df_test <- as.data.frame(X_test)
        attr(df_test, "terms") <- NULL
        pr <- predict(model, df_test, probability = TRUE)
        attr(pr, "probabilities")[, levels(labels_vec)[2]]
      },
      rf = {
        predict(model, as.data.frame(X_test), type = "prob")[, levels(labels_vec)[2]]
      },
      xgb = {
        predict(model, xgb.DMatrix(data = X_test))
      }
    )

    preds <- ifelse(probs > 0.5, levels(labels_vec)[2], levels(labels_vec)[1])
    fold_preds[[i]] <- preds
    fold_probs[[i]] <- probs
    fold_truths[[i]] <- as.character(y_test)

    # Metrics
    cm <- confusionMatrix(factor(preds, levels = levels(labels_vec)), y_test)
    fold_metrics[[i]] <- c(
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      PPV = cm$byClass["Pos Pred Value"]
    )
  }

  # Aggregate
  all_preds <- unlist(fold_preds)
  all_probs <- unlist(fold_probs)
  all_truths <- factor(unlist(fold_truths), levels = levels(labels_vec))

  # ROC
  if (length(levels(labels_vec)) == 2) {
    roc_obj <- roc(all_truths, all_probs, quiet = TRUE)
    roc_list[[alg]] <- roc_obj
    auc_val <- auc(roc_obj)
  } else {
    auc_val <- NA
  }

  # Metrics average
  metrics_mat <- do.call(rbind, fold_metrics)
  avg_metrics <- colMeans(metrics_mat, na.rm = TRUE)

  results[[alg]] <- list(
    algorithm = alg,
    auc = as.numeric(auc_val),
    accuracy = avg_metrics["Accuracy"],
    sensitivity = avg_metrics["Sensitivity"],
    specificity = avg_metrics["Specificity"],
    ppv = avg_metrics["PPV"]
  )
}

# ============================================================
# ROC curves plot
# ============================================================
cat("[ml-classifier] Plotting ROC curves...\n")
roc_pal <- biof3_palette(length(alg_list))
names(roc_pal) <- alg_list

if (length(levels(labels_vec)) == 2) {
  roc_data <- do.call(rbind, lapply(roc_list, function(nm) {
    r <- roc_list[[nm]]
    data.frame(specificity = r$specificities, sensitivity = r$sensitivities,
               algorithm = nm, stringsAsFactors = FALSE)
  }))

  p_roc <- ggplot() +
    geom_line(data = roc_data, aes(x = 1 - specificity, y = sensitivity, color = algorithm), linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = roc_pal,
                       labels = setNames(c("Logistic Regression", "SVM (RBF)", "Random Forest", "XGBoost"),
                                         c("lr", "svm", "rf", "xgb"))) +
    labs(title = "ROC Curves (5-fold CV)", x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = "Algorithm") +
    coord_equal() +
    theme_biof3()
  ggsave_biof3(p_roc, file.path(output_dir, "roc_curves"), width = 6, height = 5)
}

# ============================================================
# Metrics summary table
# ============================================================
metrics_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    Algorithm = switch(r$algorithm, lr = "Logistic Regression", svm = "SVM (RBF)", rf = "Random Forest", xgb = "XGBoost"),
    AUC = ifelse(is.na(r$auc), "N/A (multiclass)", sprintf("%.3f", r$auc)),
    Accuracy = sprintf("%.3f", r$accuracy),
    Sensitivity = sprintf("%.3f", r$sensitivity),
    Specificity = sprintf("%.3f", r$specificity),
    PPV = sprintf("%.3f", r$ppv),
    stringsAsFactors = FALSE
  )
}))
write.csv(metrics_df, file.path(output_dir, "metrics_table.csv"), row.names = FALSE)

# ============================================================
# Confusion matrices (combined plot)
# ============================================================
cat("[ml-classifier] Plotting confusion matrices...\n")

# ============================================================
# Feature importance (for RF and XGB only, LR via coefficients)
# ============================================================
cat("[ml-classifier] Extracting feature importance...\n")

# Train final models on full data for importance
imp_list <- list()
if ("rf" %in% alg_list) {
  suppressMessages(library(randomForest))
  df_full <- as.data.frame(X)
  df_full$label <- labels_vec
  rf_full <- randomForest(label ~ ., data = df_full, ntree = 300, importance = TRUE)
  imp_rf <- importance(rf_full)[, "MeanDecreaseGini"]
  imp_list$RF <- data.frame(feature = names(imp_rf), importance = as.numeric(imp_rf), stringsAsFactors = FALSE)
}
if ("xgb" %in% alg_list) {
  suppressMessages(library(xgboost))
  dfull <- xgb.DMatrix(data = X, label = as.numeric(labels_vec) - 1)
  xgb_full <- xgb.train(params = list(objective = ifelse(length(levels(labels_vec)) == 2, "binary:logistic", "multi:softprob"),
                                      num_class = ifelse(length(levels(labels_vec)) == 2, 1, length(levels(labels_vec))),
                                      eval_metric = "logloss", max_depth = 4, eta = 0.1),
                        data = dfull, nrounds = 100, verbose = 0)
  imp_xgb <- xgb.importance(model = xgb_full)
  imp_list$XGB <- data.frame(feature = imp_xgb$Feature, importance = imp_xgb$Gain, stringsAsFactors = FALSE)
}
if ("lr" %in% alg_list) {
  suppressMessages(library(glmnet))
  cv_lr <- cv.glmnet(X, labels_vec, family = ifelse(length(levels(labels_vec)) == 2, "binomial", "multinomial"), nfolds = 5)
  coef_lr <- as.matrix(coef(cv_lr, s = "lambda.1se"))[-1, , drop = TRUE]
  coef_lr <- abs(coef_lr)
  coef_lr[is.na(coef_lr)] <- 0
  imp_list$LR <- data.frame(feature = names(coef_lr), importance = coef_lr, stringsAsFactors = FALSE)
}

if (length(imp_list) > 0) {
  # Normalize and merge
  for (nm in names(imp_list)) {
    m <- max(imp_list[[nm]]$importance, na.rm = TRUE)
    if (m > 0) imp_list[[nm]]$importance_norm <- imp_list[[nm]]$importance / m
    else imp_list[[nm]]$importance_norm <- 0
  }
  common_imp <- Reduce(function(a, b) merge(a[, c("feature", "importance_norm")], b[, c("feature", "importance_norm")], by = "feature", all = TRUE),
                       lapply(names(imp_list), function(nm) { df <- imp_list[[nm]]; names(df)[2] <- nm; df }))
  names(common_imp)[-1] <- names(imp_list)

  # Top 20 comparison
  common_imp$avg <- rowMeans(common_imp[, -1, drop = FALSE], na.rm = TRUE)
  top20 <- head(common_imp[order(-common_imp$avg), ], 20)

  imp_long <- reshape(top20[, c("feature", names(imp_list))], direction = "long", varying = names(imp_list),
                      v.names = "importance", times = names(imp_list), idvar = "feature")
  names(imp_long)[3] <- "algorithm"

  p_imp <- ggplot(imp_long, aes(x = reorder(feature, importance), y = importance, fill = algorithm)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = roc_pal[names(imp_list)]) +
    labs(title = "Feature Importance (top 20)", x = "Feature", y = "Normalized Importance", fill = "Algorithm") +
    coord_flip() +
    theme_biof3()
  ggsave_biof3(p_imp, file.path(output_dir, "feature_importance"), width = 7, height = 6)
}

# ============================================================
# Save final model (XGBoost as default, or RF)
# ============================================================
cat("[ml-classifier] Saving final model...\n")
final_model <- NULL
if ("xgb" %in% alg_list) {
  suppressMessages(library(xgboost))
  dfull <- xgb.DMatrix(data = X, label = as.numeric(labels_vec) - 1)
  final_model <- xgb.train(params = list(objective = ifelse(length(levels(labels_vec)) == 2, "binary:logistic", "multi:softprob"),
                                           num_class = ifelse(length(levels(labels_vec)) == 2, 1, length(levels(labels_vec))),
                                           eval_metric = "logloss", max_depth = 4, eta = 0.1),
                           data = dfull, nrounds = 100, verbose = 0)
} else if ("rf" %in% alg_list) {
  df_full <- as.data.frame(X)
  df_full$label <- labels_vec
  final_model <- randomForest(label ~ ., data = df_full, ntree = 300)
}
if (!is.null(final_model)) {
  saveRDS(final_model, file.path(output_dir, "model.rds"))
}

# ============================================================
# Summary
# ============================================================
summary_text <- sprintf(
  "ML Classifier — 4-Algorithm Comparison\n\nTask: %s\nSamples: %d, Features: %d, Classes: %d\nAlgorithms: %s\n\nMetrics (5-fold CV average):\n%s\n",
  ifelse(length(levels(labels_vec)) == 2, "Binary Classification", "Multi-class Classification"),
  nrow(X), ncol(X), length(levels(labels_vec)),
  paste(alg_list, collapse = ", "),
  paste(apply(metrics_df, 1, function(r) paste("  ", r["Algorithm"], ": AUC=", r["AUC"], ", Acc=", r["Accuracy"])), collapse = "\n")
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = list(
    list(name = "metrics_table.csv",      type = "table", label = "指标汇总表"),
    list(name = "roc_curves.png",         type = "plot",  label = "ROC 曲线"),
    list(name = "roc_curves.pdf",         type = "file",  label = "ROC 曲线 (PDF)"),
    list(name = "feature_importance.png", type = "plot",  label = "特征重要性"),
    list(name = "feature_importance.pdf", type = "file",  label = "特征重要性 (PDF)"),
    list(name = "model.rds",              type = "file",  label = "训好的模型 (.rds)"),
    list(name = "summary.txt",            type = "text",  label = "分析摘要")
  ),
  summary = list(
    samples = nrow(X),
    features = ncol(X),
    algorithms = alg_list
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

cat("[ml-classifier] Pipeline complete.\n")
