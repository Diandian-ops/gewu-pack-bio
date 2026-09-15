#!/usr/bin/env Rscript
# ============================================================
# deg-standardizer
# 科学子问题: 差异分析结果是否能够转换成稳定、可供下游插件消费的标准 DEG Artifact？
#
# 标准输出列:
#   gene, gene_id_type, log2FoldChange, pvalue, padj,
#   direction, significant, source_method
#
# 原则:
#   - 不覆盖上游原始结果
#   - 不伪造缺失 p 值或 padj
#   - 缺少 padj 时必须明确阻断或标记 degraded
#   - ranked_genes 必须保留足够基因，不能只保留显著基因
#   - 阈值属于参数和科学决策，必须记录
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

read_table_auto <- function(f) {
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

cat("[deg-standardizer] Starting standardization\n")

# ---- 平台检查 ----
platform <- Sys.info()[["sysname"]]
arch <- Sys.info()[["machine"]]
platform_norm <- tolower(platform)
if (platform_norm == "darwin") platform_norm <- "macos"
platform_id <- paste0(platform_norm, "-", tolower(arch))
supported_platforms <- c("macos-arm64")

# ---- 运行时检查 ----
runtime_ok <- requireNamespace("jsonlite", quietly = TRUE) &&
              requireNamespace("digest", quietly = TRUE)

# ---- 读取参数 ----
params <- fromJSON(file.path(job_dir, "params.json"))

gene_col <- params$gene_col %||% params$gene_id_column %||% "gene"
gene_id_type <- params$gene_id_type %||% "symbol"
effect_col <- params$effect_col %||% params$effect_size_column %||% "log2FoldChange"
pval_col <- params$pval_col %||% params$p_value_column %||% "pvalue"
padj_col <- params$padj_col %||% params$adjusted_p_value_column %||% "padj"
padj_threshold <- as.numeric(params$padj_threshold %||% 0.05)
effect_threshold <- as.numeric(params$effect_threshold %||% params$effect_size_threshold %||% 1)
source_method <- params$source_method %||% "DESeq2"
# 上游 artifact 引用
deg_artifact_ref <- params$deg_artifact %||% params$deg_result %||% NULL

# ---- 诊断结构 ----
diagnostics <- list(
  plugin = "deg-standardizer",
  version = "1.0.0",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  checks = list(),
  status = "valid",
  failureMode = NA_character_,
  warnings = list(),
  summary = list(),
  filterRules = list(
    padjThreshold = padj_threshold,
    effectThreshold = effect_threshold,
    sourceMethod = source_method
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

# ---- 平台检查 ----
if (!(platform_id %in% supported_platforms)) {
  set_blocked("platform_not_supported",
              paste0("Platform '", platform_id, "' not supported"))
}

# ---- 运行时检查 ----
if (diagnostics$status != "blocked" && !runtime_ok) {
  set_blocked("runtime_not_ready", "Required R packages (jsonlite, digest) not available")
}

# ---- 上游 lineage ----
parent_artifacts <- list()
if (!is.null(deg_artifact_ref)) {
  parent_artifacts <- list(list(
    artifactId = deg_artifact_ref,
    plugin = source_method,
    outputId = "deg_results"
  ))
  add_check("upstream_artifact_linked", TRUE, "info",
            paste("Linked to upstream:", deg_artifact_ref))
} else {
  add_check("upstream_artifact_linked", FALSE, "warning",
            "No upstream DEG artifact reference — running standalone")
  set_warning_state("No upstream DEG artifact reference provided")
}

# ---- 查找输入文件 ----
deg_file <- list.files(job_dir, pattern = "^deg|^deseq|^diff|^results", full.names = TRUE)[1]

# ---- 检查 1: 输入文件存在 ----
if (is.na(deg_file)) {
  add_check("deg_file_present", FALSE, "error", "DEG result file not found in job dir")
  set_blocked("missing_required_file_inputs", "DEG result file is required")
} else {
  add_check("deg_file_present", TRUE, "info", paste("Found:", basename(deg_file)))
}

# ---- 读取并检查 DEG 表 ----
if (diagnostics$status != "blocked") {

  deg_raw <- read_table_auto(deg_file)
  n_rows_raw <- nrow(deg_raw)

  # ---- 检查 2: gene 列存在 ----
  if (!(gene_col %in% colnames(deg_raw))) {
    add_check("gene_column_present", FALSE, "error",
              paste("Column '", gene_col, "' not found. Available: ",
                    paste(colnames(deg_raw), collapse = ", ")))
    set_blocked("invalid_input",
                paste("Gene column '", gene_col, "' not found in DEG results"))
  } else {
    add_check("gene_column_present", TRUE, "info", paste("Gene column:", gene_col))
  }
}

# ---- 检查 effect / pvalue / padj 列 ----
if (diagnostics$status != "blocked") {

  # ---- 检查 3: effect size 列 ----
  if (!(effect_col %in% colnames(deg_raw))) {
    add_check("effect_column_present", FALSE, "error",
              paste("Effect size column '", effect_col, "' not found"))
    set_blocked("invalid_input",
                paste("Effect size column '", effect_col, "' not found"))
  } else {
    add_check("effect_column_present", TRUE, "info", paste("Effect column:", effect_col))
  }

  # ---- 检查 4: pvalue 列 ----
  if (!(pval_col %in% colnames(deg_raw))) {
    add_check("pvalue_column_present", FALSE, "error",
              paste("P-value column '", pval_col, "' not found"))
    set_blocked("invalid_input",
                paste("P-value column '", pval_col, "' not found"))
  } else {
    add_check("pvalue_column_present", TRUE, "info", paste("P-value column:", pval_col))
  }

  # ---- 检查 5: padj 列 ----
  has_padj <- padj_col %in% colnames(deg_raw)
  if (!has_padj) {
    add_check("padj_column_present", FALSE, "error",
              paste("Adjusted p-value column '", padj_col, "' not found — cannot determine significance without FDR control"))
    # 缺少 padj 时阻断或标记 degraded
    set_blocked("invalid_input",
                paste("Adjusted p-value column '", padj_col, "' not found. ",
                      "Cannot standardize without FDR-controlled padj. ",
                      "If upstream method does not produce padj, ",
                      "consider using a method that does."))
  } else {
    add_check("padj_column_present", TRUE, "info", paste("Adjusted p-value column:", padj_col))
  }
}

# ---- 标准化 ----
if (diagnostics$status != "blocked") {

  # 提取标准列
  deg_std <- data.frame(
    gene = as.character(deg_raw[[gene_col]]),
    gene_id_type = gene_id_type,
    log2FoldChange = as.numeric(deg_raw[[effect_col]]),
    pvalue = as.numeric(deg_raw[[pval_col]]),
    padj = as.numeric(deg_raw[[padj_col]]),
    stringsAsFactors = FALSE
  )

  # ---- 检查 6: NA 在关键列 ----
  na_gene <- sum(is.na(deg_std$gene) | deg_std$gene == "")
  na_lfc <- sum(is.na(deg_std$log2FoldChange))
  na_pval <- sum(is.na(deg_std$pvalue))
  na_padj <- sum(is.na(deg_std$padj))

  if (na_gene > 0) {
    add_check("no_na_genes", FALSE, "warning",
              paste(na_gene, "rows with missing gene ID"))
    set_warning_state(paste(na_gene, "rows with missing gene ID — will be excluded"))
    deg_std <- deg_std[!is.na(deg_std$gene) & deg_std$gene != "", ]
  } else {
    add_check("no_na_genes", TRUE, "info", "No missing gene IDs")
  }

  if (na_lfc > 0) {
    add_check("no_na_effect", FALSE, "warning",
              paste(na_lfc, "rows with missing effect size"))
    set_warning_state(paste(na_lfc, "rows with missing log2FoldChange"))
  } else {
    add_check("no_na_effect", TRUE, "info", "No missing effect sizes")
  }

  if (na_pval > 0) {
    add_check("no_na_pvalue", FALSE, "warning",
              paste(na_pval, "rows with missing p-value"))
    set_warning_state(paste(na_pval, "rows with missing pvalue"))
  } else {
    add_check("no_na_pvalue", TRUE, "info", "No missing p-values")
  }

  if (na_padj > 0) {
    add_check("no_na_padj", FALSE, "warning",
              paste(na_padj, "rows with missing padj"))
    set_warning_state(paste(na_padj, "rows with missing padj — marked as non-significant"))
    # 不伪造 padj，将 NA padj 标记为非显著
    deg_std$padj[is.na(deg_std$padj)] <- 1
  } else {
    add_check("no_na_padj", TRUE, "info", "No missing padj values")
  }

  # ---- 检查 7: 重复 gene ID ----
  dup_genes <- deg_std$gene[duplicated(deg_std$gene)]
  if (length(dup_genes) > 0) {
    add_check("no_duplicate_genes", FALSE, "warning",
              paste(length(unique(dup_genes)), "duplicate gene IDs"))
    set_warning_state(paste("Duplicate gene IDs in DEG results:", length(unique(dup_genes)), "genes"))
  } else {
    add_check("no_duplicate_genes", TRUE, "info", "No duplicate gene IDs")
  }

  # ---- 计算 direction 和 significant ----
  deg_std$direction <- ifelse(deg_std$log2FoldChange > 0, "up",
                       ifelse(deg_std$log2FoldChange < 0, "down", "neutral"))
  deg_std$significant <- (deg_std$padj < padj_threshold &
                          abs(deg_std$log2FoldChange) >= effect_threshold)
  deg_std$source_method <- source_method

  # ---- 检查 8: 空结果 ----
  n_sig <- sum(deg_std$significant, na.rm = TRUE)
  n_up <- sum(deg_std$significant & deg_std$direction == "up", na.rm = TRUE)
  n_down <- sum(deg_std$significant & deg_std$direction == "down", na.rm = TRUE)
  n_total <- nrow(deg_std)

  diagnostics$summary$totalGenes <- n_total
  diagnostics$summary$significantGenes <- n_sig
  diagnostics$summary$upGenes <- n_up
  diagnostics$summary$downGenes <- n_down

  if (n_total == 0) {
    add_check("non_empty_results", FALSE, "error", "DEG table has 0 rows after cleaning")
    set_blocked("invalid_input", "DEG table is empty after removing missing gene IDs")
  } else {
    add_check("non_empty_results", TRUE, "info",
              paste(n_total, "genes total,", n_sig, "significant"))
  }

  if (n_sig == 0 && n_total > 0) {
    add_check("has_significant_genes", FALSE, "warning",
              "No significant genes found — this is a valid scientific result")
    set_warning_state("No significant genes found — valid empty result, not a failure")
  } else {
    add_check("has_significant_genes", TRUE, "info",
              paste(n_sig, "significant genes (", n_up, "up,", n_down, "down)"))
  }
}

# ============================================================
# 生成输出文件
# ============================================================

if (diagnostics$status != "blocked") {

  # deg_standardized.csv — 完整标准化表
  write.csv(deg_std, file.path(output_dir, "deg_standardized.csv"), row.names = FALSE)

  # deg_up.csv — 显著上调
  deg_up <- deg_std[deg_std$significant & deg_std$direction == "up", , drop = FALSE]
  write.csv(deg_up, file.path(output_dir, "deg_up.csv"), row.names = FALSE)

  # deg_down.csv — 显著下调
  deg_down <- deg_std[deg_std$significant & deg_std$direction == "down", , drop = FALSE]
  write.csv(deg_down, file.path(output_dir, "deg_down.csv"), row.names = FALSE)

  # ranked_genes.csv — 按统计显著性排序的完整基因列表
  # 必须保留足够基因，不能只保留显著基因
  ranked <- deg_std[order(deg_std$padj, -abs(deg_std$log2FoldChange)), , drop = FALSE]
  ranked$rank <- seq_len(nrow(ranked))
  write.csv(ranked, file.path(output_dir, "ranked_genes.csv"), row.names = FALSE)

  # filter_summary.json
  filter_summary <- list(
    padjThreshold = padj_threshold,
    effectThreshold = effect_threshold,
    sourceMethod = source_method,
    totalGenes = n_total,
    significantGenes = n_sig,
    upGenes = n_up,
    downGenes = n_down,
    nonSignificantGenes = n_total - n_sig,
    note = if (n_sig == 0) "No significant genes — valid empty result" else NULL
  )
  writeLines(toJSON(filter_summary, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "filter_summary.json"))

  diagnostics$summary$rankedGenesFile <- "ranked_genes.csv"
  diagnostics$summary$rankedGenesCount <- nrow(ranked)
}

# filter_summary.json (always, even if blocked)
if (diagnostics$status == "blocked") {
  writeLines(toJSON(list(
    status = "blocked",
    failureMode = diagnostics$failureMode,
    message = diagnostics$blockedMessage %||% "",
    padjThreshold = padj_threshold,
    effectThreshold = effect_threshold,
    sourceMethod = source_method
  ), auto_unbox = TRUE, pretty = TRUE),
  file.path(output_dir, "filter_summary.json"))
}

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
    '<h3>Filter Rules (Scientific Decisions)</h3>',
    '<table style="border-collapse:collapse;">',
    '<tr><td style="padding:4px;font-weight:bold;">padj threshold</td><td>', padj_threshold, '</td></tr>',
    '<tr><td style="padding:4px;font-weight:bold;">|log2FC| threshold</td><td>', effect_threshold, '</td></tr>',
    '<tr><td style="padding:4px;font-weight:bold;">source method</td><td>', source_method, '</td></tr>',
    '</table></div>'
  )

  summary_html <- ""
  if (diagnostics$status != "blocked") {
    summary_html <- paste0(
      '<div style="margin:12px 0;padding:12px;background:#f0fdf4;border-radius:6px;">',
      '<h3>Standardization Summary</h3>',
      '<table style="border-collapse:collapse;width:100%;">',
      '<tr><td style="padding:4px;font-weight:bold;">Total Genes</td><td>', diagnostics$summary$totalGenes %||% 0, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Significant</td><td>', diagnostics$summary$significantGenes %||% 0, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Up</td><td>', diagnostics$summary$upGenes %||% 0, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Down</td><td>', diagnostics$summary$downGenes %||% 0, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Ranked Genes (all)</td><td>', diagnostics$summary$rankedGenesCount %||% 0, '</td></tr>',
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

  html <- paste0(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>DEG Standardization Report</title>
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
  <h1>DEG Standardization Report</h1>
  <div class="meta">
    Plugin: deg-standardizer v1.0.0 |
    Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ' |
    <span class="status-badge">', status_icon, ' ', toupper(diagnostics$status), '</span>
  </div>
</div>

<h2>Standardization Checks</h2>
', checks_html, '
', warnings_html, '
', filter_rules_html, '
', summary_html, '

<div style="margin-top:20px;padding-top:12px;border-top:1px solid #e5e7eb;font-size:11px;color:#6b7280;">
  Generated by BioF3 deg-standardizer |
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
if (!is.na(deg_file) && file.exists(deg_file)) {
  input_hashes[basename(deg_file)] <- sha256_file(deg_file)
}

# R 包版本
pkg_versions <- list()
for (p in c("jsonlite", "digest", "tools")) {
  if (requireNamespace(p, quietly = TRUE)) {
    pkg_versions[[p]] <- as.character(packageVersion(p))
  }
}

manifest <- list(
  pluginId = "deg-standardizer",
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
    dagStep = "deg-standardization",
    downstreamConsumers = c("volcano-plot", "pca-explorer", "go-kegg", "gsea")
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
    replayCommand = paste("Rscript", script_path %||% "deg-standardizer.R", job_dir)
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

cat(sprintf("[deg-standardizer] Done. Status: %s\n", diagnostics$status))
if (diagnostics$status == "blocked") {
  quit(status = 1)
}
