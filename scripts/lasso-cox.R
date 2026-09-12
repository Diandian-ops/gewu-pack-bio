#!/usr/bin/env Rscript
# ============================================================
# LASSO Cox 预后 Signature 构建工具
#
# 输入（job_dir/）：
#   - expression 文件 (CSV/TSV: 行=基因, 列=样本)
#   - clinical 文件 (CSV/TSV: 含 OS.time, OS 列)
#   - candidate_genes 文件（可选，一列基因名）
#   - params.json: {train_ratio, seed}
#
# 输出（job_dir/output/）：
#   - Signature 基因表 + LASSO CV + KM + ROC + Risk 三联图 + report.html
# ============================================================

suppressMessages({
  library(glmnet)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(jsonlite)
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题 (W3-1 升级 2026-05-25)
source(file.path(job_dir, "_biof3-theme.R"))

params <- fromJSON(file.path(job_dir, "params.json"))
train_ratio <- as.numeric(params$train_ratio %||% 0.7)
seed <- as.integer(params$seed %||% 42)

cat("Parameters loaded.\n")
cat("  Train ratio:", train_ratio, "\n")
cat("  Random seed:", seed, "\n")

set.seed(seed)

# ============================================================
# Read input files
# ============================================================
expr_file <- list.files(job_dir, pattern = "^expression|^expr|^tpm|^fpkm",
                        full.names = TRUE, ignore.case = TRUE)[1]
if (is.na(expr_file)) {
  all_files <- list.files(job_dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
  all_files <- all_files[!grepl("params|clinical|coldata|survival|pheno|candidate|genes", basename(all_files), ignore.case = TRUE)]
  expr_file <- all_files[1]
}

clin_file <- list.files(job_dir, pattern = "clinical|survival|coldata|pheno",
                        full.names = TRUE, ignore.case = TRUE)[1]
if (is.na(clin_file)) {
  all_files <- list.files(job_dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
  all_files <- all_files[basename(all_files) != basename(expr_file)]
  all_files <- all_files[!grepl("candidate|genes_list", basename(all_files), ignore.case = TRUE)]
  clin_file <- all_files[1]
}

cand_file <- list.files(job_dir, pattern = "candidate|genes_list|gene_list",
                        full.names = TRUE, ignore.case = TRUE)[1]

if (is.na(expr_file)) stop("No expression file found")
if (is.na(clin_file)) stop("No clinical file found")

cat("Reading expression:", basename(expr_file), "\n")
expr <- tryCatch(
  read.csv(expr_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) read.delim(expr_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
)

cat("Reading clinical:", basename(clin_file), "\n")
clin <- tryCatch(
  read.csv(clin_file, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) read.delim(clin_file, stringsAsFactors = FALSE, check.names = FALSE)
)

if (!is.null(clin[[1]]) && is.character(clin[[1]])) {
  rownames(clin) <- clin[[1]]
}

cat(sprintf("  Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("  Clinical: %d samples x %d variables\n", nrow(clin), ncol(clin)))

# Detect time/event columns
time_col <- NULL; event_col <- NULL
for (tc in c("OS.time", "os.time", "OS_time", "time", "survival_time", "futime")) {
  if (tc %in% colnames(clin)) { time_col <- tc; break }
}
for (ec in c("OS", "os", "status", "event", "vital_status")) {
  if (ec %in% colnames(clin)) { event_col <- ec; break }
}
if (is.null(time_col)) stop("Cannot find time column in clinical data")
if (is.null(event_col)) stop("Cannot find event column in clinical data")
cat("  Time column:", time_col, "| Event column:", event_col, "\n")

# Merge and filter
common_samples <- intersect(colnames(expr), rownames(clin))
if (length(common_samples) < 30) stop(paste0("Too few common samples: ", length(common_samples)))

surv_time <- as.numeric(clin[common_samples, time_col])
surv_event <- as.numeric(clin[common_samples, event_col])
valid <- !is.na(surv_time) & !is.na(surv_event) & surv_time > 0
common_samples <- common_samples[valid]
surv_time <- surv_time[valid]
surv_event <- surv_event[valid]
cat(sprintf("  Valid samples: %d\n", length(common_samples)))

# Candidate genes
if (!is.na(cand_file)) {
  cand_genes <- tryCatch({
    cg <- read.csv(cand_file, header = TRUE, stringsAsFactors = FALSE)
    as.character(cg[[1]])
  }, error = function(e) readLines(cand_file))
  cand_genes <- intersect(cand_genes, rownames(expr))
  cat(sprintf("  Candidate genes from file: %d\n", length(cand_genes)))
} else {
  gene_var <- apply(expr[, common_samples], 1, var, na.rm = TRUE)
  cand_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, sum(gene_var > 0))]
  cat(sprintf("  Using top %d variable genes\n", length(cand_genes)))
}

if (length(cand_genes) < 10) stop("Too few candidate genes")

# ============================================================
# Step 1: Univariate Cox screening
# ============================================================
cat("Step 1: Univariate Cox screening...\n")

expr_mat <- t(expr[cand_genes, common_samples])
surv_obj <- Surv(surv_time, surv_event)

uni_p <- sapply(1:ncol(expr_mat), function(i) {
  tryCatch({
    fit <- coxph(surv_obj ~ expr_mat[, i])
    summary(fit)$coefficients[1, 5]
  }, error = function(e) 1)
})
names(uni_p) <- cand_genes

sig_genes <- names(uni_p[uni_p < 0.05])
cat(sprintf("  Significant genes (p < 0.05): %d\n", length(sig_genes)))

if (length(sig_genes) < 5) {
  sig_genes <- names(sort(uni_p))[1:min(50, length(uni_p))]
  cat(sprintf("  Relaxed: using top %d genes by p-value\n", length(sig_genes)))
}
if (length(sig_genes) > 200) {
  sig_genes <- names(sort(uni_p[sig_genes]))[1:200]
  cat("  Capped at 200 genes for LASSO\n")
}

uni_df <- data.frame(gene = names(uni_p), pvalue = uni_p, significant = uni_p < 0.05)
uni_df <- uni_df[order(uni_df$pvalue), ]
write.csv(uni_df, file.path(output_dir, "univariate_cox.csv"), row.names = FALSE)

# ============================================================
# Step 2: Train/test split
# ============================================================
cat("Step 2: Train/test split...\n")

n <- length(common_samples)
train_idx <- sample(1:n, round(n * train_ratio))
test_idx <- setdiff(1:n, train_idx)
cat(sprintf("  Train: %d, Test: %d\n", length(train_idx), length(test_idx)))

x_train <- as.matrix(expr_mat[train_idx, sig_genes])
y_train <- surv_obj[train_idx]
x_test <- as.matrix(expr_mat[test_idx, sig_genes])
y_test <- surv_obj[test_idx]

# ============================================================
# Step 3: LASSO Cox
# ============================================================
cat("Step 3: LASSO Cox regression...\n")

cv_fit <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1, nfolds = 10)

tryCatch({
  save_grid_biof3(file.path(output_dir, "lasso_cv"), width = 7, height = 5, expr = {
    par(mar = c(4.5, 4.5, 2, 1), bty = "l")
    plot(cv_fit, main = "LASSO Cox Cross-Validation",
         cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)
    abline(v = log(cv_fit$lambda.min), col = biof3_palette_div(3)[3], lty = 2, lwd = 1.5)
    abline(v = log(cv_fit$lambda.1se), col = biof3_palette_div(3)[1], lty = 2, lwd = 1.5)
    legend("topright", legend = c("lambda.min", "lambda.1se"),
           col = c(biof3_palette_div(3)[3], biof3_palette_div(3)[1]),
           lty = 2, lwd = 1.5, bty = "n", cex = 0.95)
  })
}, error = function(e) cat("  LASSO CV plot failed:", e$message, "\n"))

coef_mat <- coef(cv_fit, s = "lambda.min")
coef_vec <- as.numeric(coef_mat)
names(coef_vec) <- rownames(coef_mat)
selected <- coef_vec[coef_vec != 0]

if (length(selected) == 0) {
  coef_mat <- coef(cv_fit, s = "lambda.1se")
  coef_vec <- as.numeric(coef_mat)
  names(coef_vec) <- rownames(coef_mat)
  selected <- coef_vec[coef_vec != 0]
}

if (length(selected) == 0) {
  writeLines("LASSO did not select any genes.", file.path(output_dir, "summary.txt"))
  manifest <- list(files = list(list(name = "summary.txt", type = "text", label = "Analysis note")),
                   summary = list(n_selected = 0))
  writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), file.path(output_dir, "manifest.json"))
  stop("LASSO selected 0 genes")
}

cat(sprintf("  Selected genes: %d\n", length(selected)))

sig_table <- data.frame(gene = names(selected), coefficient = round(selected, 6), stringsAsFactors = FALSE)
sig_table <- sig_table[order(abs(sig_table$coefficient), decreasing = TRUE), ]
write.csv(sig_table, file.path(output_dir, "signature_genes.csv"), row.names = FALSE)

# Save the cv.glmnet object for online plotting (FigCode lasso-cv-plot)
# Small object (~few hundred KB), enables exact reproduction of plot(cv_fit)
tryCatch({
  saveRDS(cv_fit, file.path(output_dir, "lasso_cv_fit.rds"))
}, error = function(e) cat("  saveRDS(cv_fit) failed:", e$message, "\n"))

# ============================================================
# Step 4: Risk Score
# ============================================================
cat("Step 4: Calculating risk scores...\n")

risk_all <- as.numeric(expr_mat[, names(selected)] %*% selected)
risk_cutoff <- median(as.numeric(expr_mat[train_idx, names(selected)] %*% selected))

risk_df <- data.frame(
  sample = common_samples, risk_score = risk_all,
  group = factor(ifelse(risk_all >= risk_cutoff, "High", "Low"), levels = c("Low", "High")),
  time = surv_time, event = surv_event,
  set = ifelse(1:n %in% train_idx, "Train", "Test"),
  stringsAsFactors = FALSE
)
write.csv(risk_df, file.path(output_dir, "risk_scores.csv"), row.names = FALSE)
# Save as survival_data.csv for FigCode compatibility
# Rename risk_score -> gene_expr so KM FigCode template works
surv_compat <- risk_df[, c("sample", "risk_score", "group", "time", "event")]
colnames(surv_compat)[colnames(surv_compat) == "risk_score"] <- "gene_expr"
write.csv(surv_compat, file.path(output_dir, "survival_data.csv"), row.names = FALSE)

# Save signature gene expression matrix (top N for risk_triple) for FigCode online plotting
# Columns are samples (sorted by risk_score), rows are signature genes
tryCatch({
  top_sig_for_save <- head(sig_table$gene, min(15, nrow(sig_table)))
  risk_sorted_for_save <- risk_df[order(risk_df$risk_score), ]
  sig_expr_mat <- as.matrix(expr[top_sig_for_save, risk_sorted_for_save$sample])
  # Convert to data frame with gene as first column
  sig_expr_df <- data.frame(gene = rownames(sig_expr_mat), sig_expr_mat, check.names = FALSE)
  write.csv(sig_expr_df, file.path(output_dir, "signature_expr_matrix.csv"), row.names = FALSE)
}, error = function(e) cat("  signature_expr_matrix save failed:", e$message, "\n"))

n_high <- sum(risk_df$group == "High")
n_low <- sum(risk_df$group == "Low")

# ============================================================
# Step 5: KM curve
# ============================================================
cat("Step 5: KM survival...\n")

fit_km <- survfit(Surv(time, event) ~ group, data = risk_df)
logrank <- survdiff(Surv(time, event) ~ group, data = risk_df)
logrank_p <- 1 - pchisq(logrank$chisq, df = 1)

tryCatch({
  km_pal <- biof3_palette(2)  # NPG 蓝-红
  p <- ggsurvplot(fit_km, data = risk_df,
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE, risk.table.height = 0.25,
                  palette = km_pal,
                  legend.labs = c(paste0("Low Risk (n=", n_low, ")"), paste0("High Risk (n=", n_high, ")")),
                  xlab = "Time", ylab = "Survival Probability",
                  title = "LASSO Cox Risk Score - KM Curve",
                  ggtheme = theme_biof3())
  save_grid_biof3(file.path(output_dir, "km_curve"),
                   width = 8, height = 7, expr = print(p))
}, error = function(e) cat("  KM plot failed:", e$message, "\n"))

# ============================================================
# Step 6: Time-dependent ROC
# ============================================================
cat("Step 6: ROC...\n")

tryCatch({
  suppressMessages(library(timeROC))
  time_vals <- risk_df$time
  if (max(time_vals) > 365) { roc_times <- c(365, 365*3, 365*5); roc_labels <- c("1-year", "3-year", "5-year") }
  else if (max(time_vals) > 12) { roc_times <- c(12, 36, 60); roc_labels <- c("1-year", "3-year", "5-year") }
  else { roc_times <- c(1, 3, 5); roc_labels <- c("1-year", "3-year", "5-year") }

  max_time <- max(time_vals) * 0.9
  valid_times <- roc_times[roc_times <= max_time]
  valid_labels <- roc_labels[roc_times <= max_time]
  if (length(valid_times) == 0) {
    valid_times <- c(quantile(time_vals, 0.25), median(time_vals), quantile(time_vals, 0.75))
    valid_labels <- c("Q1", "Median", "Q3")
  }

  roc_res <- timeROC(T = risk_df$time, delta = risk_df$event,
                     marker = risk_df$risk_score, cause = 1, times = valid_times, iid = TRUE)
  aucs <- roc_res$AUC

  roc_pal <- c(biof3_palette_div(3)[3], biof3_palette_div(3)[1], biof3_palette()[3])  # 红/蓝/绿
  save_grid_biof3(file.path(output_dir, "roc_curve"),
                   width = 6.5, height = 6, expr = {
    par(mar = c(4.5, 4.5, 2, 1), bty = "l")
    plot(roc_res, time = valid_times[1], col = roc_pal[1], lwd = 2,
         title = "LASSO Cox - Time-dependent ROC")
    if (length(valid_times) >= 2) plot(roc_res, time = valid_times[2], col = roc_pal[2], lwd = 2, add = TRUE)
    if (length(valid_times) >= 3) plot(roc_res, time = valid_times[3], col = roc_pal[3], lwd = 2, add = TRUE)
    abline(0, 1, lty = 2, col = "gray50")
    legend("bottomright", legend = paste0(valid_labels, " (AUC=", sprintf("%.3f", aucs), ")"),
           col = roc_pal[1:length(valid_times)], lwd = 2, bty = "n", cex = 0.95)
  })
  write.csv(data.frame(Time = valid_labels, AUC = round(aucs, 3)), file.path(output_dir, "roc_auc.csv"), row.names = FALSE)
}, error = function(e) cat("  ROC failed:", e$message, "\n"))

# ============================================================
# Step 7: Risk triple plot
# ============================================================
cat("Step 7: Risk triple plot...\n")

tryCatch({
  risk_sorted <- risk_df[order(risk_df$risk_score), ]
  risk_sorted$rank <- 1:nrow(risk_sorted)

  triple_pal <- biof3_palette(2)  # NPG 蓝/红 (Low/High)

  p1 <- ggplot(risk_sorted, aes(x = rank, y = risk_score, color = group)) +
    geom_point(size = 1) +
    scale_color_manual(values = c(Low = triple_pal[1], High = triple_pal[2])) +
    geom_hline(yintercept = risk_cutoff, linetype = "dashed", color = "gray40") +
    labs(x = "Patients (sorted)", y = "Risk Score", title = "Risk Score Distribution") +
    theme_biof3() + theme(legend.position = "none")

  status_pal <- c("0" = biof3_palette()[3], "1" = biof3_palette_div(3)[3])  # 绿活 / 红死
  p2 <- ggplot(risk_sorted, aes(x = rank, y = time, color = factor(event))) +
    geom_point(size = 1) +
    scale_color_manual(values = status_pal, labels = c("Alive", "Dead"), name = NULL) +
    labs(x = "Patients (sorted)", y = "Survival Time", title = "Survival Status") +
    theme_biof3() + theme(legend.position = "top")

  top_sig <- head(sig_table$gene, min(10, nrow(sig_table)))
  heat_mat <- as.matrix(expr[top_sig, risk_sorted$sample])
  heat_mat_z <- t(scale(t(heat_mat)))
  heat_mat_z[heat_mat_z >  2] <-  2
  heat_mat_z[heat_mat_z < -2] <- -2

  heat_df <- expand.grid(gene = rownames(heat_mat_z), sample = 1:ncol(heat_mat_z))
  heat_df$value <- as.vector(heat_mat_z)

  div_pal <- biof3_palette_div(3)  # 蓝-白-红
  p3 <- ggplot(heat_df, aes(x = sample, y = gene, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = div_pal[1], mid = "white", high = div_pal[3], midpoint = 0) +
    labs(x = "Patients (sorted)", y = NULL, title = "Signature Genes (Z-score)") +
    theme_biof3() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  suppressMessages(library(gridExtra))
  save_grid_biof3(file.path(output_dir, "risk_triple"),
                   width = 8, height = 10, expr = {
    grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 1, 1.2))
  })
}, error = function(e) cat("  Risk triple failed:", e$message, "\n"))

# Coefficient bar plot (ggplot 系)
tryCatch({
  sig_table$direction <- ifelse(sig_table$coefficient > 0, "Risk", "Protective")
  coef_div <- biof3_palette_div(3)
  p <- ggplot(sig_table, aes(x = reorder(gene, coefficient), y = coefficient, fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(Risk = coef_div[3], Protective = coef_div[1])) +
    coord_flip() +
    labs(x = NULL, y = "LASSO Coefficient",
         title = paste0("Signature: ", nrow(sig_table), " Genes"),
         subtitle = "Red = Risk Factor, Blue = Protective") +
    theme_biof3() +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 10, color = "gray40"))
  ggsave_biof3(p, file.path(output_dir, "coef_barplot"),
                width = 7, height = max(4, nrow(sig_table) * 0.4))
}, error = function(e) cat("  Coef plot failed:", e$message, "\n"))

# ============================================================
# Manifest + Report
# ============================================================
summary_text <- sprintf(
  "LASSO Cox Prognostic Signature\n\nSamples: %d (Train: %d, Test: %d)\nCandidate genes: %d\nUnivariate significant: %d\nLASSO selected: %d\nSignature: %s\nRisk cutoff: %.4f\nHigh/Low: %d/%d\nLog-rank p: %.2e\n",
  n, length(train_idx), length(test_idx), length(cand_genes),
  sum(uni_p < 0.05), nrow(sig_table), paste(sig_table$gene, collapse = ", "),
  risk_cutoff, n_high, n_low, logrank_p)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

files_list <- list(
  list(name = "summary.txt", type = "text", label = "\u5206\u6790\u6458\u8981"),
  list(name = "signature_genes.csv", type = "table", label = "Signature \u57fa\u56e0"),
  list(name = "risk_scores.csv", type = "table", label = "Risk Score"),
  list(name = "univariate_cox.csv", type = "table", label = "\u5355\u56e0\u7d20 Cox")
)
# Optional intermediate files (only added if successfully saved)
for (extra_f in c("signature_expr_matrix.csv", "lasso_cv_fit.rds")) {
  if (file.exists(file.path(output_dir, extra_f))) {
    extra_label <- if (extra_f == "signature_expr_matrix.csv") "Signature \u8868\u8fbe\u77e9\u9635" else "LASSO CV \u5bf9\u8c61"
    extra_type <- if (endsWith(extra_f, ".rds")) "file" else "table"
    files_list <- c(files_list, list(list(name = extra_f, type = extra_type, label = extra_label)))
  }
}
plot_pairs <- list(
  c("lasso_cv",     "LASSO CV"),
  c("coef_barplot", "Signature 系数"),
  c("km_curve",     "KM 曲线"),
  c("roc_curve",    "ROC"),
  c("risk_triple",  "Risk 三联图")
)
for (pp in plot_pairs) {
  png_path <- file.path(output_dir, paste0(pp[1], ".png"))
  pdf_path <- file.path(output_dir, paste0(pp[1], ".pdf"))
  if (file.exists(png_path)) {
    files_list <- c(files_list, list(list(name = paste0(pp[1], ".png"), type = "plot", label = paste0(pp[2], " (PNG)"))))
  }
  if (file.exists(pdf_path)) {
    files_list <- c(files_list, list(list(name = paste0(pp[1], ".pdf"), type = "file", label = paste0(pp[2], " (PDF)"))))
  }
}
if (file.exists(file.path(output_dir, "roc_auc.csv")))
  files_list <- c(files_list, list(list(name = "roc_auc.csv", type = "table", label = "ROC AUC")))
files_list <- c(files_list, list(list(name = "report.html", type = "file", label = "\u89e3\u8bfb\u62a5\u544a")))

manifest <- list(files = files_list, summary = list(
  n_samples = n, n_train = length(train_idx), n_test = length(test_idx),
  n_candidate = length(cand_genes), n_selected = nrow(sig_table),
  signature_genes = paste(sig_table$gene, collapse = ", "),
  risk_cutoff = round(risk_cutoff, 4), logrank_p = signif(logrank_p, 3),
  n_high = n_high, n_low = n_low))
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), file.path(output_dir, "manifest.json"))

# Report generation
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  cat("Generating HTML report...\n")
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
  png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
  for (img_path in png_files) {
    img_name <- basename(img_path)
    b64 <- base64enc::base64encode(img_path)
    report_html <- gsub(paste0("{{", img_name, "}}"), paste0("data:image/png;base64,", b64), report_html, fixed = TRUE)
  }
  report_html <- gsub("{{n_samples}}", n, report_html, fixed = TRUE)
  report_html <- gsub("{{n_train}}", length(train_idx), report_html, fixed = TRUE)
  report_html <- gsub("{{n_test}}", length(test_idx), report_html, fixed = TRUE)
  report_html <- gsub("{{n_selected}}", nrow(sig_table), report_html, fixed = TRUE)
  report_html <- gsub("{{signature_genes}}", paste(sig_table$gene, collapse = ", "), report_html, fixed = TRUE)
  report_html <- gsub("{{logrank_p}}", sprintf("%.2e", logrank_p), report_html, fixed = TRUE)
  report_html <- gsub("{{n_high}}", n_high, report_html, fixed = TRUE)
  report_html <- gsub("{{n_low}}", n_low, report_html, fixed = TRUE)
  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("Report generated.\n")
}

cat("\n========================================\n")
cat("LASSO Cox complete!\n")
cat(sprintf("  %d genes selected | Log-rank p=%.2e\n", nrow(sig_table), logrank_p))
cat("========================================\n")
