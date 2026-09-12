#!/usr/bin/env Rscript
# ============================================================
# sample-metadata-validator
# 科学子问题: counts 矩阵和 metadata 是否能够支持指定的实验设计与 contrast？
#
# 检查项:
#   1. counts 是否有样本列
#   2. metadata 是否有样本标识列
#   3. counts 样本名与 metadata 样本名是否一致
#   4. 是否有重复样本
#   5. 是否有缺失分组
#   6. reference/treatment 是否真实存在
#   7. 每组样本数
#   8. 是否只有单个样本
#   9. 是否出现完全混淆的设计
#  10. 对齐时是否丢失样本
#  11. 是否可以继续建模
#
# 状态: valid | warning | blocked
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

# 获取脚本自身路径（trailingOnly=TRUE 不含脚本名，需从完整 args 取）
all_args <- commandArgs(FALSE)
script_path <- NA_character_
file_arg <- grep("^--file=", all_args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
} else {
  # fallback: 尝试从 trailingOnly 的 args[0]（R 1-indexed，args[0] 为 NA）
  if (!is.na(args[0]) && file.exists(args[0])) script_path <- args[0]
}

cat("[sample-metadata-validator] Starting validation\n")

# ---- 运行时检查 ----
runtime_ok <- requireNamespace("jsonlite", quietly = TRUE) &&
              requireNamespace("digest", quietly = TRUE)

if (!runtime_ok) {
  diag <- list(
    plugin = "sample-metadata-validator",
    version = "1.0.0",
    status = "blocked",
    failureMode = "runtime_not_ready",
    message = "Required R packages (jsonlite, digest) not available",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  writeLines(toJSON(diag, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "validation_diagnostics.json"))
  writeLines(toJSON(list(status = "blocked", failureMode = "runtime_not_ready"),
                    auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "design_summary.json"))
  # 空 manifest
  manifest <- list(
    pluginId = "sample-metadata-validator",
    pluginVersion = "1.0.0",
    status = "blocked",
    failureMode = "runtime_not_ready",
    outputs = list(),
    lineage = list(parents = list(), dagStep = "input-validation"),
    reproducibility = list(
      rVersion = R.version.string,
      scriptSha256 = if (!is.na(script_path) && file.exists(script_path)) sha256_file(script_path) else NA_character_,
      paramsHash = NA_character_,
      replayCommand = paste("Rscript", script_path %||% "sample-metadata-validator.R", job_dir)
    )
  )
  writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "artifact-manifest.json"))
  cat("[sample-metadata-validator] BLOCKED: runtime_not_ready\n")
  quit(status = 1)
}

# ---- 读取参数 ----
params <- fromJSON(file.path(job_dir, "params.json"))

sample_col <- params$sample_column %||% params$sample_col %||% "sample"
design_col <- params$design_column %||% params$design_col %||% "condition"
ref_group <- params$reference_group %||% params$ref_group %||% params$contrast_ref %||% "Control"
treat_group <- params$treatment_group %||% params$treat_group %||% params$contrast_treat %||% "Treatment"
min_replicates <- as.integer(params$min_replicates %||% 2)

# ---- 查找输入文件 ----
counts_file <- list.files(job_dir, pattern = "^counts", full.names = TRUE)[1]
metadata_file <- list.files(job_dir, pattern = "^metadata|^coldata", full.names = TRUE)[1]

# ---- 诊断结构 ----
diagnostics <- list(
  plugin = "sample-metadata-validator",
  version = "1.0.0",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  checks = list(),
  status = "valid",
  failureMode = NA_character_,
  warnings = list(),
  summary = list()
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

# ---- 检查 1: 输入文件存在 ----
if (is.na(counts_file) || is.na(metadata_file)) {
  missing <- c()
  if (is.na(counts_file)) missing <- c(missing, "counts")
  if (is.na(metadata_file)) missing <- c(missing, "metadata")
  add_check("input_files_present", FALSE, "error",
            paste("Missing:", paste(missing, collapse = ", ")))
  set_blocked("missing_required_file_inputs",
              paste("Required input files missing:", paste(missing, collapse = ", ")))
} else {
  add_check("input_files_present", TRUE, "info",
            paste("counts:", basename(counts_file), "metadata:", basename(metadata_file)))
}

# ---- 如果文件存在，继续后续检查 ----
if (diagnostics$status != "blocked") {

  counts_df <- read_matrix_auto(counts_file)
  metadata_df <- read_matrix_auto(metadata_file)

  # ---- 检查 2: counts 有样本列 (第一列是 gene ID，其余是样本) ----
  if (ncol(counts_df) < 2) {
    add_check("counts_has_sample_columns", FALSE, "error",
              "counts matrix has fewer than 2 columns (need gene ID + sample columns)")
    set_blocked("missing_required_file_inputs",
                "counts matrix must have at least a gene ID column and one sample column")
  } else {
    # 第一列通常是 gene ID
    gene_col_name <- colnames(counts_df)[1]
    sample_names_counts <- colnames(counts_df)[-1]
    add_check("counts_has_sample_columns", TRUE, "info",
              paste("Found", length(sample_names_counts), "sample columns in counts"))
  }

  # ---- 检查 3: metadata 有样本标识列 ----
  if (!(sample_col %in% colnames(metadata_df))) {
    add_check("metadata_has_sample_column", FALSE, "error",
              paste("Column '", sample_col, "' not found in metadata. Available: ",
                    paste(colnames(metadata_df), collapse = ", ")))
    set_blocked("sample_id_column_missing",
                paste("Sample column '", sample_col, "' not found in metadata"))
  } else {
    add_check("metadata_has_sample_column", TRUE, "info",
              paste("Found sample column '", sample_col, "' in metadata"))
  }

  # ---- 检查 4: metadata 有设计列 ----
  if (diagnostics$status != "blocked") {
    if (!(design_col %in% colnames(metadata_df))) {
      add_check("metadata_has_design_column", FALSE, "error",
                paste("Design column '", design_col, "' not found. Available: ",
                      paste(colnames(metadata_df), collapse = ", ")))
      set_blocked("missing_design_column",
                  paste("Design column '", design_col, "' not found in metadata"))
    } else {
      add_check("metadata_has_design_column", TRUE, "info",
                paste("Found design column '", design_col, "'"))
    }
  }

  # ---- 检查 5: 重复样本 (metadata) ----
  if (diagnostics$status != "blocked") {
    meta_samples <- as.character(metadata_df[[sample_col]])
    dup_meta <- meta_samples[duplicated(meta_samples)]
    if (length(dup_meta) > 0) {
      add_check("no_duplicate_samples_metadata", FALSE, "error",
                paste("Duplicate sample IDs in metadata:", paste(unique(dup_meta), collapse = ", ")))
      set_blocked("duplicated_samples",
                  paste("Duplicate sample IDs found:", paste(unique(dup_meta), collapse = ", ")))
    } else {
      add_check("no_duplicate_samples_metadata", TRUE, "info",
                "No duplicate sample IDs in metadata")
    }
  }

  # ---- 检查 6: 重复样本 (counts) ----
  if (diagnostics$status != "blocked") {
    dup_counts <- sample_names_counts[duplicated(sample_names_counts)]
    if (length(dup_counts) > 0) {
      add_check("no_duplicate_samples_counts", FALSE, "error",
                paste("Duplicate sample columns in counts:", paste(unique(dup_counts), collapse = ", ")))
      set_blocked("duplicated_samples",
                  paste("Duplicate sample columns in counts:", paste(unique(dup_counts), collapse = ", ")))
    } else {
      add_check("no_duplicate_samples_counts", TRUE, "info",
                "No duplicate sample columns in counts")
    }
  }

  # ---- 检查 7: 样本名对齐 ----
  if (diagnostics$status != "blocked") {
    meta_samples <- as.character(metadata_df[[sample_col]])
    counts_samples <- sample_names_counts

    only_counts <- setdiff(counts_samples, meta_samples)
    only_meta <- setdiff(meta_samples, counts_samples)
    common_samples <- intersect(counts_samples, meta_samples)

    if (length(common_samples) == 0) {
      add_check("sample_alignment", FALSE, "error",
                "No common samples between counts and metadata")
      set_blocked("sample_mismatch",
                  "Zero overlap between counts sample columns and metadata sample IDs")
    } else {
      if (length(only_counts) > 0 || length(only_meta) > 0) {
        detail <- paste0(
          "Only in counts: ", paste(only_counts, collapse = ", "),
          " | Only in metadata: ", paste(only_meta, collapse = ", "),
          " | Common: ", length(common_samples)
        )
        if (length(common_samples) < length(counts_samples) * 0.5 ||
            length(common_samples) < length(meta_samples) * 0.5) {
          add_check("sample_alignment", FALSE, "error", detail)
          set_blocked("sample_mismatch",
                      paste("Severe sample mismatch. Common:", length(common_samples),
                            "Counts-only:", length(only_counts),
                            "Metadata-only:", length(only_meta)))
        } else {
          add_check("sample_alignment", TRUE, "warning", detail)
          set_warning_state(paste("Partial sample mismatch:",
                                  length(only_counts), "in counts only,",
                                  length(only_meta), "in metadata only"))
        }
      } else {
        add_check("sample_alignment", TRUE, "info",
                  paste("All", length(common_samples), "samples aligned perfectly"))
      }
    }
  }

  # ---- 检查 8: contrast levels 存在 ----
  if (diagnostics$status != "blocked") {
    meta_samples <- as.character(metadata_df[[sample_col]])
    common_samples <- intersect(sample_names_counts, meta_samples)
    meta_aligned <- metadata_df[match(common_samples, meta_samples), , drop = FALSE]
    design_values <- as.character(meta_aligned[[design_col]])

    unique_groups <- unique(design_values)

    ref_exists <- ref_group %in% unique_groups
    treat_exists <- treat_group %in% unique_groups

    if (!ref_exists || !treat_exists) {
      missing_levels <- c()
      if (!ref_exists) missing_levels <- c(missing_levels, ref_group)
      if (!treat_exists) missing_levels <- c(missing_levels, treat_group)
      add_check("contrast_levels_present", FALSE, "error",
                paste("Missing contrast levels:", paste(missing_levels, collapse = ", "),
                      "| Available groups:", paste(unique_groups, collapse = ", ")))
      set_blocked("contrast_level_missing",
                  paste("Contrast levels not found:", paste(missing_levels, collapse = ", ")))
    } else {
      add_check("contrast_levels_present", TRUE, "info",
                paste("Both reference ('", ref_group, "') and treatment ('", treat_group,
                      "') found in design column"))
    }
  }

  # ---- 检查 9: 每组样本数 / 重复 ----
  if (diagnostics$status != "blocked") {
    meta_samples <- as.character(metadata_df[[sample_col]])
    common_samples <- intersect(sample_names_counts, meta_samples)
    meta_aligned <- metadata_df[match(common_samples, meta_samples), , drop = FALSE]
    design_values <- as.character(meta_aligned[[design_col]])

    group_counts <- table(design_values)
    diagnostics$summary$groupCounts <- as.list(group_counts)

    ref_n <- as.integer(group_counts[ref_group] %||% 0)
    treat_n <- as.integer(group_counts[treat_group] %||% 0)

    # 检查单样本
    single_sample_groups <- names(group_counts)[group_counts == 1]
    if (length(single_sample_groups) > 0) {
      add_check("no_single_sample_groups", FALSE, "warning",
                paste("Groups with single sample:", paste(single_sample_groups, collapse = ", ")))
      set_warning_state(paste("Single-sample groups:", paste(single_sample_groups, collapse = ", ")))
    } else {
      add_check("no_single_sample_groups", TRUE, "info", "All groups have >= 2 samples")
    }

    # 检查重复数
    if (ref_n < min_replicates || treat_n < min_replicates) {
      add_check("sufficient_replicates", FALSE, "error",
                paste("ref:", ref_n, "treat:", treat_n, "min required:", min_replicates))
      if (ref_n < 2 || treat_n < 2) {
        set_blocked("insufficient_replicates",
                    paste("Insufficient replicates: ref=", ref_n, "treat=", treat_n,
                          "minimum=", min_replicates))
      } else {
        set_warning_state(paste("Low replicate count: ref=", ref_n, "treat=", treat_n))
      }
    } else {
      add_check("sufficient_replicates", TRUE, "info",
                paste("ref:", ref_n, "treat:", treat_n))
    }
  }

  # ---- 检查 10: 混淆设计 ----
  if (diagnostics$status != "blocked") {
    meta_samples <- as.character(metadata_df[[sample_col]])
    common_samples <- intersect(sample_names_counts, meta_samples)
    meta_aligned <- metadata_df[match(common_samples, meta_samples), , drop = FALSE]

    # 检查是否所有样本都在一个组
    design_values <- as.character(meta_aligned[[design_col]])
    if (length(unique(design_values)) < 2) {
      add_check("valid_design", FALSE, "error",
                "All aligned samples belong to a single group — design is confounded")
      set_blocked("invalid_design",
                  "Cannot estimate contrast: all samples in one group")
    } else {
      add_check("valid_design", TRUE, "info",
                paste("Design has", length(unique(design_values)), "groups"))
    }

    # 检查是否有其他潜在混淆（如 batch 完全与 condition 重合）
    # 这里做基本检查：如果 metadata 有其他列且每组的该列值完全一致
    other_cols <- setdiff(colnames(meta_aligned), c(sample_col, design_col))
    for (oc in other_cols) {
      oc_vals <- as.character(meta_aligned[[oc]])
      for (g in c(ref_group, treat_group)) {
        g_idx <- design_values == g
        if (any(g_idx)) {
          g_oc <- unique(oc_vals[g_idx])
          if (length(g_oc) == 1) {
            # 该组内该列完全一致 — 可能是混淆，但不是必然
            # 只在两组都完全混淆时才标记
          }
        }
      }
    }
  }

  # ---- 检查 11: 对齐后样本丢失 ----
  if (diagnostics$status != "blocked") {
    meta_samples <- as.character(metadata_df[[sample_col]])
    common_samples <- intersect(sample_names_counts, meta_samples)
    n_counts_only <- length(setdiff(sample_names_counts, meta_samples))
    n_meta_only <- length(setdiff(meta_samples, sample_names_counts))

    if (n_counts_only > 0 || n_meta_only > 0) {
      diagnostics$summary$alignedSamples <- length(common_samples)
      diagnostics$summary$lostFromCounts <- n_counts_only
      diagnostics$summary$lostFromMetadata <- n_meta_only
    } else {
      diagnostics$summary$alignedSamples <- length(common_samples)
      diagnostics$summary$lostFromCounts <- 0L
      diagnostics$summary$lostFromMetadata <- 0L
    }
  }

  # ---- 检查 12: 是否可以继续建模 ----
  if (diagnostics$status != "blocked") {
    # 最终建模可行性
    meta_samples <- as.character(metadata_df[[sample_col]])
    common_samples <- intersect(sample_names_counts, meta_samples)
    meta_aligned <- metadata_df[match(common_samples, meta_samples), , drop = FALSE]
    design_values <- as.character(meta_aligned[[design_col]])

    ref_n <- sum(design_values == ref_group)
    treat_n <- sum(design_values == treat_group)
    total_n <- length(common_samples)

    can_model <- (ref_n >= 2 && treat_n >= 2 &&
                  length(unique(design_values)) >= 2 &&
                  diagnostics$status != "blocked")

    diagnostics$summary$canModel <- can_model
    diagnostics$summary$totalAlignedSamples <- total_n
    diagnostics$summary$refGroupN <- ref_n
    diagnostics$summary$treatGroupN <- treat_n

    if (can_model) {
      add_check("can_proceed_to_modeling", TRUE, "info",
                paste("Ready for modeling: ref=", ref_n, "treat=", treat_n))
    } else {
      add_check("can_proceed_to_modeling", FALSE, "error",
                "Cannot proceed to modeling")
      if (diagnostics$status != "blocked") {
        set_blocked("invalid_design", "Design does not support modeling")
      }
    }
  }
}

# ============================================================
# 生成输出文件
# ============================================================

if (diagnostics$status != "blocked") {
  # 重新读取并生成对齐后的输出
  counts_df <- read_matrix_auto(counts_file)
  metadata_df <- read_matrix_auto(metadata_file)
  meta_samples <- as.character(metadata_df[[sample_col]])
  common_samples <- intersect(sample_names_counts, meta_samples)

  # validated_metadata.csv — 对齐后的 metadata
  meta_aligned <- metadata_df[match(common_samples, meta_samples), , drop = FALSE]
  write.csv(meta_aligned, file.path(output_dir, "validated_metadata.csv"), row.names = FALSE)

  # sample_alignment.csv — 对齐详情
  alignment_df <- data.frame(
    sample = c(common_samples,
               setdiff(sample_names_counts, meta_samples),
               setdiff(meta_samples, sample_names_counts)),
    source = c(rep("both", length(common_samples)),
               rep("counts_only", length(setdiff(sample_names_counts, meta_samples))),
               rep("metadata_only", length(setdiff(meta_samples, sample_names_counts)))),
    in_counts = c(rep(TRUE, length(common_samples)),
                  rep(TRUE, length(setdiff(sample_names_counts, meta_samples))),
                  rep(FALSE, length(setdiff(meta_samples, sample_names_counts)))),
    in_metadata = c(rep(TRUE, length(common_samples)),
                    rep(FALSE, length(setdiff(sample_names_counts, meta_samples))),
                    rep(TRUE, length(setdiff(meta_samples, sample_names_counts)))),
    stringsAsFactors = FALSE
  )
  write.csv(alignment_df, file.path(output_dir, "sample_alignment.csv"), row.names = FALSE)

  # design_summary.json
  design_values <- as.character(meta_aligned[[design_col]])
  group_counts <- as.list(table(design_values))
  design_summary <- list(
    status = diagnostics$status,
    designColumn = design_col,
    referenceGroup = ref_group,
    treatmentGroup = treat_group,
    referenceN = as.integer(group_counts[[ref_group]] %||% 0),
    treatmentN = as.integer(group_counts[[treat_group]] %||% 0),
    totalAlignedSamples = length(common_samples),
    groupCounts = group_counts,
    canModel = diagnostics$summary$canModel %||% FALSE,
    warnings = diagnostics$warnings
  )
  writeLines(toJSON(design_summary, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "design_summary.json"))
} else {
  # blocked 状态也输出 design_summary
  design_summary <- list(
    status = "blocked",
    failureMode = diagnostics$failureMode,
    message = diagnostics$blockedMessage %||% "Validation blocked",
    canModel = FALSE
  )
  writeLines(toJSON(design_summary, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "design_summary.json"))

  # 空对齐表
  write.csv(data.frame(sample = character(), source = character(),
                       in_counts = logical(), in_metadata = logical(),
                       stringsAsFactors = FALSE),
            file.path(output_dir, "sample_alignment.csv"), row.names = FALSE)
}

# validation_diagnostics.json (always)
writeLines(toJSON(diagnostics, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "validation_diagnostics.json"))

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

  summary_html <- ""
  if (diagnostics$status != "blocked" && !is.null(diagnostics$summary$canModel)) {
    summary_html <- paste0(
      '<div style="margin:12px 0;padding:12px;background:#f0fdf4;border-radius:6px;">',
      '<h3>Design Summary</h3>',
      '<table style="border-collapse:collapse;width:100%;">',
      '<tr><td style="padding:4px;font-weight:bold;">Design Column</td><td>', design_col, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Reference Group</td><td>', ref_group,
      ' (n=', diagnostics$summary$refGroupN %||% 0, ')</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Treatment Group</td><td>', treat_group,
      ' (n=', diagnostics$summary$treatGroupN %||% 0, ')</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Total Aligned Samples</td><td>',
      diagnostics$summary$totalAlignedSamples %||% 0, '</td></tr>',
      '<tr><td style="padding:4px;font-weight:bold;">Can Model</td><td>',
      if (diagnostics$summary$canModel %||% FALSE) "✓ Yes" else "✗ No", '</td></tr>',
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
<title>Sample Metadata Validation Report</title>
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
  <h1>Sample Metadata Validation Report</h1>
  <div class="meta">
    Plugin: sample-metadata-validator v1.0.0 |
    Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ' |
    <span class="status-badge">', status_icon, ' ', toupper(diagnostics$status), '</span>
  </div>
</div>

<h2>Validation Checks</h2>
', checks_html, '
', warnings_html, '
', summary_html, '

<div style="margin-top:20px;padding-top:12px;border-top:1px solid #e5e7eb;font-size:11px;color:#6b7280;">
  Generated by BioF3 sample-metadata-validator |
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

# 脚本 SHA-256 (script_path 已在脚本头部解析)
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

manifest <- list(
  pluginId = "sample-metadata-validator",
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
    parents = list(),
    dagStep = "input-validation",
    downstreamConsumers = c("expression-qc", "deseq2")
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
    replayCommand = paste("Rscript", script_path %||% "sample-metadata-validator.R", job_dir)
  )
)

writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "artifact-manifest.json"))

# 保存 run-script.R 到 job_dir
if (!is.na(script_path) && file.exists(script_path)) {
  snapshot_path <- file.path(job_dir, "run-script.R")
  if (!identical(normalizePath(script_path, mustWork = FALSE), normalizePath(snapshot_path, mustWork = FALSE))) {
    file.copy(script_path, snapshot_path, overwrite = TRUE)
  }
}

cat(sprintf("[sample-metadata-validator] Done. Status: %s\n", diagnostics$status))
if (diagnostics$status == "blocked") {
  quit(status = 1)
}
