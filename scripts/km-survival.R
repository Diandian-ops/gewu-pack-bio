#!/usr/bin/env Rscript
# ============================================================
# KM 生存分析工具脚本
# Version targeted: survival 3.x / survminer 0.4+ / timeROC 0.4+
#
# Steps (in-scope) — 单基因预后探索 4 件套:
#   1. read expression + clinical
#   2. cutoff 分组 (median / optimal / tertile)
#   3. KM survival curve + log-rank
#   4. univariate Cox + forest plot
#   5. time-dependent ROC (1/3/5 year)
#   6. expression distribution
#
# 2026-05-25: 升级到 BioF3 SCI 视觉风格规范
#   - ggsurvplot (survminer 复合): save_grid_biof3() (含 risk table 子图)
#   - forest_plot / expression_dist (ggplot 系): theme_biof3 + ggsave_biof3
#   - timeROC (base R plot): save_grid_biof3
#   - 配色: biof3_palette(2) Low/High; biof3_palette_div() ROC 多线条
# ============================================================

suppressMessages({
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

# ★ Load BioF3 SCI 主题
source(file.path(job_dir, "_biof3-theme.R"))

cat("[km] Starting pipeline (SCI style v2026-05-25)\n")

params <- fromJSON(file.path(job_dir, "params.json"))
gene <- params$gene
cutoff_method <- params$cutoff_method %||% "median"
time_col <- params$time_col %||% "OS.time"
event_col <- params$event_col %||% "OS"

cat(sprintf("  gene=%s, cutoff=%s, time_col=%s, event_col=%s\n",
            gene, cutoff_method, time_col, event_col))

# ============================================================
# Read input files
# ============================================================
expr_file <- list.files(job_dir, pattern = "^expression|^expr|^tpm|^fpkm|^counts",
                        full.names = TRUE, ignore.case = TRUE)[1]
if (is.na(expr_file)) {
  all_files <- list.files(job_dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
  all_files <- all_files[!grepl("params|clinical|coldata|survival|pheno", basename(all_files), ignore.case = TRUE)]
  expr_file <- all_files[1]
}

clin_file <- list.files(job_dir, pattern = "clinical|survival|coldata|pheno",
                        full.names = TRUE, ignore.case = TRUE)[1]
if (is.na(clin_file)) {
  all_files <- list.files(job_dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
  all_files <- all_files[basename(all_files) != basename(expr_file)]
  clin_file <- all_files[1]
}

if (is.na(expr_file)) stop("No expression file found")
if (is.na(clin_file)) stop("No clinical file found")

cat(sprintf("[km] Reading expression: %s\n", basename(expr_file)))
expr <- tryCatch(
  read.csv(expr_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) read.delim(expr_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
)

cat(sprintf("[km] Reading clinical: %s\n", basename(clin_file)))
clin <- tryCatch(
  read.csv(clin_file, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) read.delim(clin_file, stringsAsFactors = FALSE, check.names = FALSE)
)

if (!is.null(clin[[1]]) && is.character(clin[[1]])) {
  rownames(clin) <- clin[[1]]
}

cat(sprintf("  Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("  Clinical: %d samples x %d variables\n", nrow(clin), ncol(clin)))

# ============================================================
# Validate inputs
# ============================================================
if (!gene %in% rownames(expr)) {
  match_idx <- which(toupper(rownames(expr)) == toupper(gene))
  if (length(match_idx) > 0) {
    gene <- rownames(expr)[match_idx[1]]
    cat(sprintf("  Gene matched (case-insensitive): %s\n", gene))
  } else {
    stop(paste0("Gene '", gene, "' not found. First 10 available: ",
                paste(head(rownames(expr), 10), collapse = ", ")))
  }
}

if (!time_col %in% colnames(clin)) {
  alt_time <- c("OS.time", "os.time", "OS_time", "time", "survival_time", "futime", "days_to_death")
  found <- intersect(alt_time, colnames(clin))
  if (length(found) > 0) {
    time_col <- found[1]
    cat(sprintf("  Auto-detected time column: %s\n", time_col))
  } else {
    stop(paste0("Time column '", time_col, "' not found. Available: ", paste(colnames(clin), collapse = ", ")))
  }
}

if (!event_col %in% colnames(clin)) {
  alt_event <- c("OS", "os", "status", "event", "vital_status", "dead")
  found <- intersect(alt_event, colnames(clin))
  if (length(found) > 0) {
    event_col <- found[1]
    cat(sprintf("  Auto-detected event column: %s\n", event_col))
  } else {
    stop(paste0("Event column '", event_col, "' not found. Available: ", paste(colnames(clin), collapse = ", ")))
  }
}

# ============================================================
# Merge expression + clinical
# ============================================================
common_samples <- intersect(colnames(expr), rownames(clin))
if (length(common_samples) < 10) {
  stop(paste0("Too few common samples: ", length(common_samples)))
}

surv_data <- data.frame(
  sample = common_samples,
  gene_expr = as.numeric(expr[gene, common_samples]),
  time = as.numeric(clin[common_samples, time_col]),
  event = as.numeric(clin[common_samples, event_col]),
  stringsAsFactors = FALSE
)

surv_data <- surv_data[!is.na(surv_data$time) & !is.na(surv_data$event) &
                       !is.na(surv_data$gene_expr) & surv_data$time > 0, ]
cat(sprintf("  Valid samples: %d\n", nrow(surv_data)))

if (nrow(surv_data) < 20) stop("Too few valid samples (need >= 20)")

# ============================================================
# Grouping by cutoff
# ============================================================
cat(sprintf("[km] Grouping by %s...\n", cutoff_method))

if (cutoff_method == "optimal") {
  tryCatch({
    cut_res <- surv_cutpoint(surv_data, time = "time", event = "event", variables = "gene_expr")
    cutoff_val <- cut_res$cutpoint$cutpoint
    surv_data$group <- ifelse(surv_data$gene_expr >= cutoff_val, "High", "Low")
    cat(sprintf("  Optimal cutoff: %.4f\n", cutoff_val))
  }, error = function(e) {
    cat("  surv_cutpoint failed, falling back to median.\n")
    cutoff_val <<- median(surv_data$gene_expr)
    surv_data$group <<- ifelse(surv_data$gene_expr >= cutoff_val, "High", "Low")
  })
} else if (cutoff_method == "tertile") {
  q33 <- quantile(surv_data$gene_expr, 1/3)
  q67 <- quantile(surv_data$gene_expr, 2/3)
  surv_data <- surv_data[surv_data$gene_expr <= q33 | surv_data$gene_expr >= q67, ]
  surv_data$group <- ifelse(surv_data$gene_expr >= q67, "High", "Low")
  cutoff_val <- q67
  cat(sprintf("  Tertile cutoff: %.4f / %.4f, n=%d\n", q33, q67, nrow(surv_data)))
} else {
  cutoff_val <- median(surv_data$gene_expr)
  surv_data$group <- ifelse(surv_data$gene_expr >= cutoff_val, "High", "Low")
  cat(sprintf("  Median cutoff: %.4f\n", cutoff_val))
}

surv_data$group <- factor(surv_data$group, levels = c("Low", "High"))
n_high <- sum(surv_data$group == "High")
n_low <- sum(surv_data$group == "Low")
cat(sprintf("  High: %d, Low: %d\n", n_high, n_low))

# ============================================================
# KM Survival Curve (survminer 复合, 用 save_grid_biof3)
# ============================================================
cat("[km] Fitting KM model...\n")

fit <- survfit(Surv(time, event) ~ group, data = surv_data)
logrank <- survdiff(Surv(time, event) ~ group, data = surv_data)
logrank_p <- 1 - pchisq(logrank$chisq, df = 1)

cat(sprintf("  Log-rank p: %.2e\n", logrank_p))

km_pal <- biof3_palette(2)  # NPG 前 2 色 (蓝-红)

tryCatch({
  km_p <- ggsurvplot(fit, data = surv_data,
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE, risk.table.height = 0.25,
                  palette = km_pal,
                  legend.labs = c(paste0(gene, " Low (n=", n_low, ")"),
                                  paste0(gene, " High (n=", n_high, ")")),
                  xlab = "Time", ylab = "Survival Probability",
                  title = paste0("Kaplan-Meier Curve: ", gene),
                  ggtheme = theme_biof3())
  save_grid_biof3(file.path(output_dir, "km_curve"),
                   width = 8, height = 7, expr = print(km_p))
}, error = function(e) cat(sprintf("  KM plot failed: %s\n", e$message)))

# ============================================================
# Cox Regression + Forest plot (ggplot 系)
# ============================================================
cat("[km] Cox regression...\n")

cox_uni <- coxph(Surv(time, event) ~ gene_expr, data = surv_data)
cox_uni_summary <- summary(cox_uni)
hr <- cox_uni_summary$conf.int[1, 1]
hr_lower <- cox_uni_summary$conf.int[1, 3]
hr_upper <- cox_uni_summary$conf.int[1, 4]
cox_p <- cox_uni_summary$coefficients[1, 5]

cat(sprintf("  HR=%.3f (%.3f-%.3f), p=%.2e\n", hr, hr_lower, hr_upper, cox_p))

cox_table <- data.frame(
  Variable = gene,
  HR = round(hr, 3),
  HR_lower = round(hr_lower, 3),
  HR_upper = round(hr_upper, 3),
  P_value = signif(cox_p, 3),
  stringsAsFactors = FALSE
)
write.csv(cox_table, file.path(output_dir, "cox_results.csv"), row.names = FALSE)

tryCatch({
  forest_df <- data.frame(
    variable = gene, hr = hr, lower = hr_lower, upper = hr_upper, p = cox_p
  )
  forest_color <- if (hr > 1) biof3_palette_div(3)[3] else biof3_palette_div(3)[1]  # 红=风险, 蓝=保护
  p <- ggplot(forest_df, aes(x = hr, y = variable)) +
    geom_point(size = 4, color = forest_color) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = forest_color, linewidth = 0.6) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    scale_x_log10() +
    labs(x = "Hazard Ratio (95% CI)", y = NULL,
         title = paste0("Cox Regression: ", gene),
         subtitle = sprintf("HR = %.3f (%.3f-%.3f), p = %.2e", hr, hr_lower, hr_upper, cox_p)) +
    theme_biof3() +
    theme(plot.subtitle = element_text(size = 10, color = "gray40"))
  ggsave_biof3(p, file.path(output_dir, "forest_plot"), width = 7, height = 3.5)
}, error = function(e) cat(sprintf("  Forest plot failed: %s\n", e$message)))

# ============================================================
# Time-dependent ROC (base R, 用 save_grid_biof3)
# ============================================================
cat("[km] Time-dependent ROC...\n")

valid_times <- NULL
valid_labels <- NULL
aucs <- NULL

tryCatch({
  suppressMessages(library(timeROC))

  time_vals <- surv_data$time
  if (max(time_vals, na.rm = TRUE) > 365) {
    roc_times <- c(365, 365*3, 365*5)
    roc_labels <- c("1-year", "3-year", "5-year")
  } else if (max(time_vals, na.rm = TRUE) > 12) {
    roc_times <- c(12, 36, 60)
    roc_labels <- c("1-year", "3-year", "5-year")
  } else {
    roc_times <- c(1, 3, 5)
    roc_labels <- c("1-year", "3-year", "5-year")
  }

  max_time <- max(time_vals) * 0.9
  valid_times <<- roc_times[roc_times <= max_time]
  valid_labels <<- roc_labels[roc_times <= max_time]

  if (length(valid_times) == 0) {
    valid_times <<- c(quantile(time_vals, 0.25), median(time_vals), quantile(time_vals, 0.75))
    valid_labels <<- c("Q1", "Median", "Q3")
  }

  roc_res <- timeROC(T = surv_data$time, delta = surv_data$event,
                     marker = surv_data$gene_expr,
                     cause = 1, times = valid_times, iid = TRUE)
  aucs <<- roc_res$AUC
  cat(sprintf("  AUCs: %s\n", paste(sprintf("%.3f", aucs), collapse = ", ")))

  roc_pal <- c(biof3_palette_div(3)[3], biof3_palette_div(3)[1], biof3_palette()[3])  # 红/蓝/绿
  save_grid_biof3(file.path(output_dir, "roc_curve"),
                   width = 6.5, height = 6, expr = {
    par(mar = c(4.5, 4.5, 2, 1), bty = "l")
    plot(roc_res, time = valid_times[1], col = roc_pal[1], lwd = 2,
         title = paste0("Time-dependent ROC: ", gene))
    if (length(valid_times) >= 2) plot(roc_res, time = valid_times[2], col = roc_pal[2], lwd = 2, add = TRUE)
    if (length(valid_times) >= 3) plot(roc_res, time = valid_times[3], col = roc_pal[3], lwd = 2, add = TRUE)
    abline(0, 1, lty = 2, col = "gray50")
    legend("bottomright",
           legend = paste0(valid_labels[1:length(valid_times)], " (AUC=", sprintf("%.3f", aucs), ")"),
           col = roc_pal[1:length(valid_times)], lwd = 2, bty = "n", cex = 0.95)
  })

  roc_data <- data.frame(Time = valid_labels[1:length(valid_times)], AUC = round(aucs, 3))
  write.csv(roc_data, file.path(output_dir, "roc_auc.csv"), row.names = FALSE)

}, error = function(e) cat(sprintf("  ROC failed: %s\n", e$message)))

# ============================================================
# Expression distribution (ggplot 系)
# ============================================================
tryCatch({
  dist_pal <- biof3_palette(2)  # 与 KM 同色
  p <- ggplot(surv_data, aes(x = gene_expr, fill = group)) +
    geom_histogram(bins = 30, alpha = 0.75, position = "identity") +
    geom_vline(xintercept = cutoff_val, linetype = "dashed", color = "gray30", linewidth = 0.7) +
    scale_fill_manual(values = c(Low = dist_pal[1], High = dist_pal[2])) +
    labs(x = paste0(gene, " Expression"), y = "Count",
         title = paste0("Expression Distribution: ", gene),
         subtitle = paste0("Cutoff = ", round(cutoff_val, 3), " (", cutoff_method, ")")) +
    theme_biof3() +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 10, color = "gray40"))
  ggsave_biof3(p, file.path(output_dir, "expression_dist"), width = 7, height = 5)
}, error = function(e) cat(sprintf("  Distribution plot failed: %s\n", e$message)))

# ============================================================
# 中间数据
# ============================================================
write.csv(surv_data, file.path(output_dir, "survival_data.csv"), row.names = FALSE)

# ============================================================
# Summary
# ============================================================
summary_text <- sprintf(
  "KM Survival Analysis Summary\n\nGene: %s\nCutoff method: %s (value = %.4f)\nTotal samples: %d (High: %d, Low: %d)\n\nLog-rank test p-value: %.2e\n\nCox regression:\n  HR = %.3f (95%% CI: %.3f - %.3f)\n  p-value = %.2e\n\nInterpretation:\n  %s\n",
  gene, cutoff_method, cutoff_val, nrow(surv_data), n_high, n_low,
  logrank_p, hr, hr_lower, hr_upper, cox_p,
  ifelse(hr > 1 & cox_p < 0.05,
         paste0(gene, " high expression is significantly associated with worse prognosis (risk factor)"),
         ifelse(hr < 1 & cox_p < 0.05,
                paste0(gene, " high expression is significantly associated with better prognosis (protective factor)"),
                paste0(gene, " expression is not significantly associated with prognosis")))
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest (双格式)
# ============================================================
files_list <- list(
  list(name = "summary.txt",       type = "text",  label = "分析摘要"),
  list(name = "cox_results.csv",   type = "table", label = "Cox 回归结果"),
  list(name = "survival_data.csv", type = "table", label = "生存数据表")
)

# Plots: PNG + PDF
plot_pairs <- list(
  c("km_curve",        "KM 生存曲线"),
  c("forest_plot",     "Cox 森林图"),
  c("roc_curve",       "时间依赖 ROC"),
  c("expression_dist", "表达分布图")
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
  files_list <- c(files_list, list(list(name = "roc_auc.csv", type = "table", label = "ROC AUC 值")))

files_list <- c(files_list, list(list(name = "report.html", type = "file", label = "解读报告(HTML)")))

manifest <- list(
  files = files_list,
  summary = list(
    gene = gene,
    cutoff_method = cutoff_method,
    cutoff_value = round(cutoff_val, 4),
    n_samples = nrow(surv_data),
    n_high = n_high,
    n_low = n_low,
    logrank_p = signif(logrank_p, 3),
    hr = round(hr, 3),
    hr_ci = paste0(round(hr_lower, 3), "-", round(hr_upper, 3)),
    cox_p = signif(cox_p, 3),
    sci_style = TRUE, dual_format = TRUE
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ============================================================
# HTML 报告
# ============================================================
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  cat("[km] Generating report...\n")
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
  for (img_path in png_files) {
    img_name <- basename(img_path)
    b64 <- base64enc::base64encode(img_path)
    data_uri <- paste0("data:image/png;base64,", b64)
    report_html <- gsub(paste0("{{", img_name, "}}"), data_uri, report_html, fixed = TRUE)
  }

  report_html <- gsub("{{gene}}", gene, report_html, fixed = TRUE)
  report_html <- gsub("{{cutoff_method}}", cutoff_method, report_html, fixed = TRUE)
  report_html <- gsub("{{cutoff_value}}", round(cutoff_val, 4), report_html, fixed = TRUE)
  report_html <- gsub("{{n_samples}}", nrow(surv_data), report_html, fixed = TRUE)
  report_html <- gsub("{{n_high}}", n_high, report_html, fixed = TRUE)
  report_html <- gsub("{{n_low}}", n_low, report_html, fixed = TRUE)
  report_html <- gsub("{{logrank_p}}", sprintf("%.2e", logrank_p), report_html, fixed = TRUE)
  report_html <- gsub("{{hr}}", round(hr, 3), report_html, fixed = TRUE)
  report_html <- gsub("{{hr_lower}}", round(hr_lower, 3), report_html, fixed = TRUE)
  report_html <- gsub("{{hr_upper}}", round(hr_upper, 3), report_html, fixed = TRUE)
  report_html <- gsub("{{cox_p}}", sprintf("%.2e", cox_p), report_html, fixed = TRUE)

  interpretation <- ifelse(hr > 1 & cox_p < 0.05,
    paste0(gene, " 高表达与较差预后显著相关，是潜在的风险因子（HR > 1, p < 0.05）。"),
    ifelse(hr < 1 & cox_p < 0.05,
      paste0(gene, " 高表达与较好预后显著相关，是潜在的保护因子（HR < 1, p < 0.05）。"),
      paste0(gene, " 表达水平与患者预后无显著统计学关联（p = ", sprintf("%.2e", cox_p), "）。")))
  report_html <- gsub("{{interpretation}}", interpretation, report_html, fixed = TRUE)

  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[km] Report generated.\n")
} else {
  cat("[km] Report template not found, skipping.\n")
}

cat(sprintf("[km] Pipeline complete. gene=%s HR=%.3f p=%.2e\n", gene, hr, logrank_p))
