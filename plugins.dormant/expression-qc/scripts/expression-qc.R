#!/usr/bin/env Rscript
# ============================================================
# expression-qc
# 科学子问题: 表达/counts 矩阵是否具有足够的数据质量，可以进入差异表达建模？
#
# 只做 QC 和质量门，不执行差异表达。
# 不根据预期结论删除样本。
# 发现疑似离群样本时输出 warning，并要求确认。
# 不能自动删除用户样本。
# 过滤规则必须写入参数和报告。
# ============================================================

suppressMessages({
  library(jsonlite)
  library(digest)
  library(tools)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a

# ---- 工具函数 ----
sha256_file <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sha256_string <- function(str) {
  digest::digest(str, algo = "sha256", serialize = FALSE)
}

read_matrix_auto <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f, ignore.case = TRUE)) {
    read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# ---- 参数解析 ----
args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 获取脚本自身路径
all_args <- commandArgs(FALSE)
script_path <- NA_character_
file_arg <- grep("^--file=", all_args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
}

cat("[expression-qc] Starting QC\n")

# ---- 运行时检查 ----
runtime_ok <- requireNamespace("jsonlite", quietly = TRUE) &&
              requireNamespace("digest", quietly = TRUE)

# 检查平台
platform <- Sys.info()[["sysname"]]
arch <- Sys.info()[["machine"]]
# R 返回 "Darwin"，统一映射为 "macos"
platform_norm <- tolower(platform)
if (platform_norm == "darwin") platform_norm <- "macos"
platform_id <- paste0(platform_norm, "-", tolower(arch))
supported_platforms <- c("macos-arm64")

if (!(platform_id %in% supported_platforms)) {
  diag <- list(
    plugin = "expression-qc",
    version = "1.0.0",
    status = "blocked",
    failureMode = "platform_not_supported",
    message = paste0("Platform '", platform_id, "' not supported. Supported: ",
                     paste(supported_platforms, collapse = ", ")),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  writeLines(toJSON(diag, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "qc_summary.json"))
  manifest <- list(
    pluginId = "expression-qc",
    pluginVersion = "1.0.0",
    status = "blocked",
    failureMode = "platform_not_supported",
    outputs = list(),
    lineage = list(parents = list(), dagStep = "qc"),
    reproducibility = list(
      rVersion = R.version.string,
      scriptSha256 = if (!is.na(script_path) && file.exists(script_path)) sha256_file(script_path) else NA_character_,
      replayCommand = paste("Rscript", script_path %||% "expression-qc.R", job_dir)
    )
  )
  writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "artifact-manifest.json"))
  cat("[expression-qc] BLOCKED: platform_not_supported\n")
  quit(status = 1)
}

if (!runtime_ok) {
  diag <- list(
    plugin = "expression-qc",
    version = "1.0.0",
    status = "blocked",
    failureMode = "runtime_not_ready",
    message = "Required R packages (jsonlite, digest) not available",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  writeLines(toJSON(diag, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "qc_summary.json"))
  manifest <- list(
    pluginId = "expression-qc",
    pluginVersion = "1.0.0",
    status = "blocked",
    failureMode = "runtime_not_ready",
    outputs = list(),
    lineage = list(parents = list(), dagStep = "qc"),
    reproducibility = list(
      rVersion = R.version.string,
      scriptSha256 = if (!is.na(script_path) && file.exists(script_path)) sha256_file(script_path) else NA_character_,
      replayCommand = paste("Rscript", script_path %||% "expression-qc.R", job_dir)
    )
  )
  writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "artifact-manifest.json"))
  cat("[expression-qc] BLOCKED: runtime_not_ready\n")
  quit(status = 1)
}

# ---- 读取参数 ----
params <- fromJSON(file.path(job_dir, "params.json"))

min_count <- as.numeric(params$min_count %||% 10)
min_samples <- as.integer(params$min_samples %||% 2)
transform_mode <- params$transformation %||% params$transform_mode %||% "none"
# 上游 validated_metadata (可选)
validated_metadata_ref <- params$validated_metadata %||% params$metadata_artifact %||% NULL

# ---- 诊断结构 ----
diagnostics <- list(
  plugin = "expression-qc",
  version = "1.0.0",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  checks = list(),
  status = "valid",
  failureMode = NA_character_,
  warnings = list(),
  summary = list(),
  filterRules = list(
    minCount = min_count,
    minSamples = min_samples,
    transformation = transform_mode
  )
)

add_check <- function(name, passed, severity = "info", detail = "") {
  diagnostics$checks[[length(diagnostics$checks) + 1]] <<- list(
    name = name, passed = passed, severity = severity, detail = detail
  )
}

add_warning <- function(msg) {
  diagnostics$warnings[[length(diagnostics$warnings) + 1]] <<- msg
}

set_blocked <- function(fm, msg) {
  diagnostics$status <<- "blocked"
  diagnostics$failureMode <<- fm
  diagnostics$blockedMessage <<- msg
}

set_warning_state <- function(msg) {
  if (diagnostics$status != "blocked") {
    diagnostics$status <<- "warning"
  }
  add_warning(msg)
}

# ---- 查找输入文件 ----
counts_file <- list.files(job_dir, pattern = "^counts", full.names = TRUE)[1]
metadata_file <- list.files(job_dir, pattern = "^metadata|^coldata|^validated_metadata", full.names = TRUE)[1]

# ---- 检查 1: 输入文件存在 ----
if (is.na(counts_file)) {
  add_check("counts_file_present", FALSE, "error", "counts file not found in job dir")
  set_blocked("missing_required_file_inputs", "counts file is required")
} else {
  add_check("counts_file_present", TRUE, "info", paste("Found:", basename(counts_file)))
}

# ---- 如果有 validated_metadata，记录 lineage ----
parent_artifacts <- list()
if (!is.null(validated_metadata_ref)) {
  parent_artifacts <- list(list(
    artifactId = validated_metadata_ref,
    plugin = "sample-metadata-validator",
    outputId = "validated_metadata"
  ))
  add_check("upstream_artifact_linked", TRUE, "info",
            paste("Linked to upstream:", validated_metadata_ref))
} else {
  add_check("upstream_artifact_linked", FALSE, "warning",
            "No upstream validated_metadata artifact reference provided")
  set_warning_state("No upstream validated_metadata artifact — running standalone QC")
}

# ---- 读取并检查 counts ----
if (diagnostics$status != "blocked") {

  counts_raw <- read_matrix_auto(counts_file)

  # ---- 检查 2: gene × sample 方向 ----
  # 第一列应为 gene ID
  gene_col <- colnames(counts_raw)[1]
  sample_cols <- colnames(counts_raw)[-1]
  n_genes <- nrow(counts_raw)
  n_samples <- length(sample_cols)

  if (n_samples < 2) {
    add_check("gene_sample_orientation", FALSE, "error",
              "Fewer than 2 sample columns — cannot perform QC")
    set_blocked("invalid_input", "counts must have at least 2 sample columns")
  } else {
    add_check("gene_sample_orientation", TRUE, "info",
              paste("gene × sample:", n_genes, "genes ×", n_samples, "samples"))
    diagnostics$summary$nGenes <- n_genes
    diagnostics$summary$nSamples <- n_samples
  }
}

# ---- 数值检查 ----
if (diagnostics$status != "blocked") {

  counts_data <- counts_raw[, -1, drop = FALSE]
  # 尝试转换为数值
  counts_numeric <- suppressWarnings(as.data.frame(lapply(counts_data, function(x) as.numeric(x))))

  # ---- 检查 3: 非数值列 ----
  non_numeric_cols <- c()
  for (col in sample_cols) {
    original <- counts_data[[col]]
    converted <- counts_numeric[[col]]
    na_ratio <- sum(is.na(converted)) / length(converted)
    # 如果原始值不是 NA 但转换后是 NA，说明是非数值
    original_not_na <- !is.na(original) & original != ""
    converted_na <- is.na(converted)
    if (any(original_not_na & converted_na)) {
      non_numeric_cols <- c(non_numeric_cols, col)
    }
  }
  if (length(non_numeric_cols) > 0) {
    add_check("no_non_numeric_columns", FALSE, "error",
              paste("Non-numeric columns:", paste(non_numeric_cols, collapse = ", ")))
    set_blocked("invalid_input",
                paste("counts contains non-numeric data in columns:",
                      paste(non_numeric_cols, collapse = ", ")))
  } else {
    add_check("no_non_numeric_columns", TRUE, "info", "All sample columns are numeric")
  }
}

# ---- NA / Inf 检查 ----
if (diagnostics$status != "blocked") {

  # ---- 检查 4: NA/Inf ----
  na_count <- sum(is.na(counts_numeric))
  inf_count <- sum(is.infinite(as.matrix(counts_numeric)))
  if (na_count > 0 || inf_count > 0) {
    add_check("no_na_inf", FALSE, "error",
              paste("NA:", na_count, "Inf:", inf_count))
    set_blocked("invalid_input",
                paste("counts contains", na_count, "NA and", inf_count, "Inf values"))
  } else {
    add_check("no_na_inf", TRUE, "info", "No NA or Inf values")
  }

  # ---- 检查 5: 负 counts ----
  neg_count <- sum(counts_numeric < 0, na.rm = TRUE)
  if (neg_count > 0) {
    add_check("no_negative_counts", FALSE, "error",
              paste(neg_count, "negative values found"))
    set_blocked("invalid_input",
                paste("counts contains", neg_count, "negative values — counts must be non-negative"))
  } else {
    add_check("no_negative_counts", TRUE, "info", "No negative values")
  }

  # ---- 检查 6: 小数 counts ----
  decimal_count <- sum(counts_numeric != floor(counts_numeric) & counts_numeric > 0, na.rm = TRUE)
  if (decimal_count > 0) {
    add_check("no_decimal_counts", FALSE, "warning",
              paste(decimal_count, "non-integer values found — counts are expected to be integers"))
    set_warning_state(paste("counts contains", decimal_count,
                            "non-integer values — may be TPM/FPKM or pre-normalized data"))
  } else {
    add_check("no_decimal_counts", TRUE, "info", "All values are integers")
  }

  # ---- 检查 7: 重复 gene ID ----
  gene_ids <- as.character(counts_raw[[gene_col]])
  dup_genes <- gene_ids[duplicated(gene_ids)]
  if (length(dup_genes) > 0) {
    add_check("no_duplicate_gene_ids", FALSE, "warning",
              paste(length(unique(dup_genes)), "duplicate gene IDs found"))
    set_warning_state(paste("Duplicate gene IDs:", paste(unique(dup_genes)[1:min(5, length(unique(dup_genes)))], collapse = ", ")))
  } else {
    add_check("no_duplicate_gene_ids", TRUE, "info", "No duplicate gene IDs")
  }
}

# ---- 基因过滤统计 ----
if (diagnostics$status != "blocked") {

  rownames(counts_numeric) <- gene_ids

  # ---- 检查 8: 全零基因 ----
  all_zero <- rowSums(counts_numeric) == 0
  n_all_zero <- sum(all_zero)
  add_check("all_zero_genes", n_all_zero == 0, "warning",
            paste(n_all_zero, "all-zero genes"))
  if (n_all_zero > 0) add_warning(paste(n_all_zero, "all-zero genes detected"))

  # ---- 检查 9: 低表达基因过滤 ----
  # 基因在至少 min_samples 个样本中 count >= min_count 才保留
  pass_filter <- rowSums(counts_numeric >= min_count) >= min_samples
  n_filtered_out <- sum(!pass_filter)
  n_retained <- sum(pass_filter)
  pct_filtered <- round(100 * n_filtered_out / n_genes, 2)

  add_check("low_expression_filter", TRUE, "info",
            paste("Retained:", n_retained, "| Filtered:", n_filtered_out,
                  "(", pct_filtered, "% ) | Rule: count >=", min_count,
                  "in >=", min_samples, "samples"))
  diagnostics$summary$nGenesRetained <- n_retained
  diagnostics$summary$nGenesFiltered <- n_filtered_out
  diagnostics$summary$pctFiltered <- pct_filtered

  # 如果过滤后基因太少
  if (n_retained < 10) {
    add_check("sufficient_genes_after_filter", FALSE, "error",
              paste("Only", n_retained, "genes remain after filtering"))
    set_blocked("invalid_input",
                paste("Insufficient genes after filtering:", n_retained,
                      "— need at least 10"))
  } else {
    add_check("sufficient_genes_after_filter", TRUE, "info",
              paste(n_retained, "genes retained"))
  }
}

# ---- Library size ----
if (diagnostics$status != "blocked") {

  lib_sizes <- colSums(counts_numeric, na.rm = TRUE)
  lib_size_df <- data.frame(
    sample = names(lib_sizes),
    librarySize = as.numeric(lib_sizes),
    stringsAsFactors = FALSE
  )
  write.csv(lib_size_df, file.path(output_dir, "sample_library_size.csv"), row.names = FALSE)

  # ---- 检查 10: Library size 差异 ----
  lib_min <- min(lib_sizes)
  lib_max <- max(lib_sizes)
  lib_ratio <- lib_max / lib_min
  if (lib_ratio > 5) {
    add_check("library_size_uniformity", FALSE, "warning",
              paste("Max/min ratio:", round(lib_ratio, 2),
                    "— large library size variation"))
    set_warning_state(paste("Library size variation (max/min ratio =", round(lib_ratio, 2),
                            ") — consider normalization"))
  } else {
    add_check("library_size_uniformity", TRUE, "info",
              paste("Max/min ratio:", round(lib_ratio, 2)))
  }
  diagnostics$summary$librarySizeMin <- lib_min
  diagnostics$summary$librarySizeMax <- lib_max
  diagnostics$summary$librarySizeRatio <- round(lib_ratio, 2)
}

# ---- 样本相关性 ----
if (diagnostics$status != "blocked" && n_samples >= 2) {

  # 使用 log2(counts + 1) 计算相关性
  log_counts <- log2(counts_numeric + 1)
  cor_mat <- cor(log_counts, method = "pearson")

  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("sample_x", "sample_y", "correlation")
  write.csv(cor_df, file.path(output_dir, "sample_correlation.csv"), row.names = FALSE)

  # ---- 检查 11: 样本相关性 / 离群 ----
  # 计算每个样本与其他样本的平均相关性
  diag(cor_mat) <- NA
  mean_cor <- colMeans(cor_mat, na.rm = TRUE)
  low_cor_samples <- names(mean_cor)[mean_cor < 0.8]

  if (length(low_cor_samples) > 0) {
    add_check("sample_correlation", FALSE, "warning",
              paste("Low correlation samples (<0.8):", paste(low_cor_samples, collapse = ", ")))
    set_warning_state(paste("Potential outlier samples (mean correlation < 0.8):",
                            paste(low_cor_samples, collapse = ", "),
                            "— manual confirmation required"))
  } else {
    add_check("sample_correlation", TRUE, "info",
              paste("All samples have mean correlation > 0.8 (min:", round(min(mean_cor), 3), ")"))
  }

  # 生成相关性热图 PNG (使用基础 R，不依赖 ggplot2)
  png_path <- file.path(output_dir, "sample_correlation.png")
  grDevices::png(png_path, width = 800, height = 800, res = 120)
  par(mar = c(8, 8, 4, 2))
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(cor_mat,
      main = "Sample Correlation (Pearson, log2(counts+1))",
      color = grDevices::colorRampPalette(c("#2563eb", "white", "#dc2626"))(50),
      cluster_rows = TRUE, cluster_cols = TRUE,
      display_numbers = TRUE, number_format = "%.2f",
      fontsize_number = 8)
  } else {
    # fallback: base R heatmap
    heatmap(cor_mat, main = "Sample Correlation (Pearson, log2(counts+1))",
            col = grDevices::colorRampPalette(c("#2563eb", "white", "#dc2626"))(50),
            margins = c(10, 10), scale = "none")
  }
  grDevices::dev.off()

  # ---- PCA ----
  # 使用过滤后的基因做 PCA
  filtered_counts <- counts_numeric[pass_filter, , drop = FALSE]
  log_filtered <- t(log2(filtered_counts + 1))

  # 中心化
  log_filtered_centered <- scale(log_filtered, center = TRUE, scale = FALSE)

  # PCA via svd
  pca_result <- svd(log_filtered_centered)
  pca_coords <- pca_result$u[, 1:min(2, ncol(pca_result$u)), drop = FALSE]
  pca_var <- pca_result$d^2 / sum(pca_result$d^2)

  pca_df <- data.frame(
    sample = rownames(log_filtered),
    PC1 = pca_coords[, 1],
    PC2 = if (ncol(pca_coords) >= 2) pca_coords[, 2] else NA_real_,
    stringsAsFactors = FALSE
  )
  write.csv(pca_df, file.path(output_dir, "pca_coordinates.csv"), row.names = FALSE)

  # PCA 图
  pca_png_path <- file.path(output_dir, "pca.png")
  grDevices::png(pca_png_path, width = 800, height = 600, res = 120)
  par(mar = c(5, 5, 4, 2))
  plot(pca_df$PC1, pca_df$PC2,
       xlab = paste0("PC1 (", round(100 * pca_var[1], 1), "%)"),
       ylab = paste0("PC2 (", round(100 * pca_var[2], 1), "%)"),
       main = "PCA (log2(counts+1), filtered genes)",
       pch = 19, col = "#2563eb", cex = 1.5)
  text(pca_df$PC1, pca_df$PC2, pca_df$sample, pos = 4, cex = 0.7, col = "#4b5563")
  grDevices::dev.off()

  # ---- 检查 12: PCA 离群 ----
  # 计算 PC1+PC2 空间中每个样本到中心的距离
  pc_center <- c(mean(pca_df$PC1), mean(pca_df$PC2))
  pc_dist <- sqrt((pca_df$PC1 - pc_center[1])^2 + (pca_df$PC2 - pc_center[2])^2)
  pc_dist_z <- scale(pc_dist)[, 1]
  outlier_samples <- pca_df$sample[abs(pc_dist_z) > 2]

  if (length(outlier_samples) > 0) {
    add_check("pca_outlier_detection", FALSE, "warning",
              paste("Potential outliers (z > 2):", paste(outlier_samples, collapse = ", ")))
    set_warning_state(paste("PCA potential outliers:", paste(outlier_samples, collapse = ", "),
                            "— manual confirmation required, samples NOT removed"))
  } else {
    add_check("pca_outlier_detection", TRUE, "info", "No PCA outliers detected (z < 2)")
  }
}

# ---- metadata 对齐状态 ----
if (diagnostics$status != "blocked" && !is.na(metadata_file)) {
  metadata_df <- read_matrix_auto(metadata_file)
  meta_samples <- if ("sample" %in% colnames(metadata_df)) as.character(metadata_df$sample) else rownames(metadata_df)
  counts_samples <- sample_cols

  common <- intersect(counts_samples, meta_samples)
  if (length(common) < length(counts_samples)) {
    add_check("metadata_alignment", FALSE, "warning",
              paste("Only", length(common), "of", length(counts_samples),
                    "samples found in metadata"))
    set_warning_state(paste("Not all counts samples found in metadata:",
                            length(counts_samples) - length(common), "missing"))
  } else {
    add_check("metadata_alignment", TRUE, "info",
              "All counts samples present in metadata")
  }
}

# ---- 生成 validated_counts.csv 或引用 ----
if (diagnostics$status != "blocked") {
  # 输出过滤后的 counts (保留所有样本，只过滤基因)
  filtered_counts_out <- counts_raw[pass_filter, , drop = FALSE]
  write.csv(filtered_counts_out, file.path(output_dir, "validated_counts.csv"), row.names = FALSE)
  diagnostics$summary$validatedCountsFile <- "validated_counts.csv"
}

# ---- 基因过滤统计表 ----
if (diagnostics$status != "blocked") {
  gene_filter_stats <- data.frame(
    metric = c("total_genes", "all_zero_genes", "low_expression_filtered",
               "retained_genes", "pct_filtered"),
    value = c(n_genes, n_all_zero, n_filtered_out, n_retained, pct_filtered),
    stringsAsFactors = FALSE
  )
  write.csv(gene_filter_stats, file.path(output_dir, "gene_filter_statistics.csv"), row.names = FALSE)
}

# ============================================================
# 生成输出文件
# ============================================================

# qc_summary.json (always)
writeLines(toJSON(diagnostics, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "qc_summary.json"))

# ============================================================
# HTML 报告
# ============================================================
generate_report <- function() {
  status_color <- switch(diagnostics$status,
    valid = "#16a34a", warning = "#f59e0b", blocked = "#dc2626", "#6b7280")
  status_icon <- switch(diagnostics$status,
    valid = "✓", warning = "⚠", blocked = "✗", "?")

  checks_html <- ""
  for (chk in diagnostics$checks) {
    chk_color <- if (chk$passed) "#16a34a" else switch(chk$severity,
      error = "#dc2626", warning = "#f59e0b", "#6b7280")
    chk_icon <- if (chk$passed) "✓" else switch(chk$severity,
      error = "✗", warning = "⚠", "•")
    checks_html <- paste0(checks_html,
      '<div style="padding:8px 12px;border-left:3px solid ', chk_color,
      ';margin-bottom:4px;background:#f9fafb;">',
      '<span style="color:', chk_color, ';font-weight:bold;">', chk_icon, '</span> ',
      '<strong>', chk$name, '</strong>',
      if (nzchar(chk$detail)) paste0(' — <span style="color:#4b5563;">', chk$detail, '</span>'),
      '</div>')
  }

  warnings_html <- ""
  if (length(diagnostics$warnings) > 0) {
    warnings_html <- '<div style="margin:12px 0;padding:10px;background:#fef3c7;border-radius:6px;"><strong>⚠ Warnings:</strong><ul>'
    for (w in diagnostics$warnings) {
      warnings_html <- paste0(warnings_html, '<li>', w, '</li>')
    }
    warnings_html <- paste0(warnings_html, '</ul></div>')
  }

  filter_rules_html <- paste0(
    '<div style="margin:12px 0;padding:12px;background:#eff6ff;border-radius:6px;">',
    '<h3>Filter Rules</h3>',
    '<table style="border-collapse:collapse;">',
    '<tr><td style="padding:4px;font-weight:bold;">min_count</td><td>', min_count, '</td></tr>',
    '<tr><td style="padding:4px;font-weight:bold;">min_samples</td><td>', min_samples, '</td></tr>',
    '<tr><td style="padding:4px;font-weight:bold;">transformation</td><td>', transform_mode, '</td></tr>',
    '</table></div>'
  )

  summary_html <- ""
  if (diagnostics$status != "blocked") {
    summary_html <- paste0(
      '<div style="margin:12px 0;padding:12px;background:#f0fdf4;border-radius:6px;">',
      '<h3>QC Summary</h3>',
      '<table style="border-collapse:collapse;width:100%;">',
      '<tr><td style="padding:4px;font-weight:bold;">Total Genes</td><td>', diagnostics$summary$nGenes %||% 0, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Samples</td><td>', diagnostics$summary$nSamples %||% 0, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Genes Retained</td><td>', diagnostics$summary$nGenesRetained %||% 0, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Genes Filtered</td><td>', diagnostics$summary$nGenesFiltered %||% 0,
      ' (', diagnostics$summary$pctFiltered %||% 0, '%)</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Library Size Range</td><td>',
      diagnostics$summary$librarySizeMin %||% 0, ' – ', diagnostics$summary$librarySizeMax %||% 0,
      ' (ratio: ', diagnostics$summary$librarySizeRatio %||% 0, ')</td></tr>',
      '</table></div>'
    )
  }

  if (diagnostics$status == "blocked") {
    summary_html <- paste0(
      '<div style="margin:12px 0;padding:12px;background:#fef2f2;border-radius:6px;">',
      '<h3>Blocked</h3><p><strong>Failure Mode:</strong> ', diagnostics$failureMode, '</p>',
      '<p>', diagnostics$blockedMessage %||% "", '</p></div>'
    )
  }

  # 图片引用
  images_html <- ""
  if (file.exists(file.path(output_dir, "sample_correlation.png"))) {
    images_html <- paste0(images_html,
      '<div style="margin:12px 0;"><h3>Sample Correlation</h3>',
      '<img src="sample_correlation.png" style="max-width:100%;border:1px solid #e5e7eb;border-radius:4px;"/></div>')
  }
  if (file.exists(file.path(output_dir, "pca.png"))) {
    images_html <- paste0(images_html,
      '<div style="margin:12px 0;"><h3>PCA</h3>',
      '<img src="pca.png" style="max-width:100%;border:1px solid #e5e7eb;border-radius:4px;"/></div>')
  }

  html <- paste0(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Expression QC Report</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         margin: 0; padding: 20px; background: #ffffff; color: #1f2937; }
  .header { border-bottom: 2px solid #e5e7eb; padding-bottom: 12px; margin-bottom: 16px; }
  .header h1 { font-size: 20px; margin: 0; }
  .header .meta { font-size: 12px; color: #6b7280; margin-top: 4px; }
  .status-badge { display: inline-block; padding: 4px 12px; border-radius: 4px;
                  font-weight: bold; color: white; background: ', status_color, '; }
</style>
</head>
<body>
<div class="header">
  <h1>Expression QC Report</h1>
  <div class="meta">
    Plugin: expression-qc v1.0.0 |
    Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ' |
    <span class="status-badge">', status_icon, ' ', toupper(diagnostics$status), '</span>
  </div>
</div>

<h2>QC Checks</h2>
', checks_html, '
', warnings_html, '
', filter_rules_html, '
', summary_html, '
', images_html, '

<div style="margin-top:20px;padding-top:12px;border-top:1px solid #e5e7eb;font-size:11px;color:#6b7280;">
  Generated by BioF3 expression-qc |
  R: ', R.version.string, ' |
  Status: ', diagnostics$status, '
</div>
</body>
</html>')

  writeLines(html, file.path(output_dir, "report.html"))
}

generate_report()

# ============================================================
# Artifact Manifest
# ============================================================
output_files <- list.files(output_dir, full.names = TRUE)
output_hashes <- list()
for (f in output_files) {
  output_hashes[basename(f)] <- sha256_file(f)
}

# 脚本 SHA-256
if (is.na(script_path) || !file.exists(script_path)) {
  script_sha <- sha256_string(paste(readLines(file.path(job_dir, "run-script.R")), collapse = "\n"))
} else {
  script_sha <- sha256_file(script_path)
}

# 参数 hash
params_json <- toJSON(params, auto_unbox = TRUE)
params_hash <- sha256_string(params_json)

# 输入 SHA-256
input_hashes <- list()
if (!is.na(counts_file) && file.exists(counts_file)) {
  input_hashes[basename(counts_file)] <- sha256_file(counts_file)
}
if (!is.na(metadata_file) && file.exists(metadata_file)) {
  input_hashes[basename(metadata_file)] <- sha256_file(metadata_file)
}

# R 包版本
pkg_versions <- list()
for (p in c("jsonlite", "digest", "tools")) {
  if (requireNamespace(p, quietly = TRUE)) {
    pkg_versions[[p]] <- as.character(packageVersion(p))
  }
}
if (requireNamespace("pheatmap", quietly = TRUE)) {
  pkg_versions[["pheatmap"]] <- as.character(packageVersion("pheatmap"))
}

manifest <- list(
  pluginId = "expression-qc",
  pluginVersion = "1.0.0",
  status = diagnostics$status,
  failureMode = diagnostics$failureMode,
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  inputs = list(
    files = input_hashes,
    params = params,
    paramsHash = params_hash
  ),
  outputs = list(
    files = output_hashes
  ),
  lineage = list(
    parents = parent_artifacts,
    dagStep = "qc",
    downstreamConsumers = c("deseq2")
  ),
  reproducibility = list(
    rVersion = R.version.string,
    rLibPaths = .libPaths(),
    packageVersions = pkg_versions,
    scriptSha256 = script_sha,
    parentScriptSha256 = NA_character_,
    paramsHash = params_hash,
    inputSha256 = input_hashes,
    outputSha256 = output_hashes,
    replayCommand = paste("Rscript", script_path %||% "expression-qc.R", job_dir)
  )
)

writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "artifact-manifest.json"))

# 保存 run-script.R
if (!is.na(script_path) && file.exists(script_path)) {
  snapshot_path <- file.path(job_dir, "run-script.R")
  if (!identical(normalizePath(script_path, mustWork = FALSE), normalizePath(snapshot_path, mustWork = FALSE))) {
    file.copy(script_path, snapshot_path, overwrite = TRUE)
  }
}

cat(sprintf("[expression-qc] Done. Status: %s\n", diagnostics$status))
if (diagnostics$status == "blocked") {
  quit(status = 1)
}
