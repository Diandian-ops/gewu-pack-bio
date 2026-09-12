#!/usr/bin/env Rscript
# ============================================================
# __TOOL_TITLE__ 工具脚本
# Version targeted: <pkg> <ver>+ (R 4.5.x)
#
# Steps (in-scope):
#   1. <step 1>
#   2. <step 2>
#   3. <step 3>
#   4. 可视化 + 中间数据保存 + HTML 报告
#
# 范围卡: 见 _<TOOL_ID>-scope.md
# 4 条死规则: 禁 theme_classic / 禁 dpi=150 / 禁硬编码颜色 / 禁字号过小
#
# 修订记录:
#   YYYY-MM-DD: 初始版本
# ============================================================

suppressMessages({
  library(ggplot2)
  library(jsonlite)
  # library(<your-pkg>)   # ← 按需加载
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

# ============================================================
# 0. 环境: job_dir / output_dir / params
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
job_dir   <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题 (提供 theme_biof3 / ggsave_biof3 / save_grid_biof3)
source(file.path(job_dir, "_biof3-theme.R"))

cat("[__TOOL_ID__] Starting pipeline\n")

params <- fromJSON(file.path(job_dir, "params.json"))
# param_x <- as.numeric(params$param_x %||% 0.05)
# param_y <- as.character(params$param_y %||% "default")

# ============================================================
# 1. 读取输入
# ============================================================
# input_file <- list.files(job_dir, pattern = "^input|^data", full.names = TRUE)[1]
# if (is.na(input_file)) stop("No input file found")
# dat <- read.csv(input_file, row.names = 1, check.names = FALSE)

# ============================================================
# 2. 核心分析逻辑
# ============================================================
# ... 在此实现 ...

# ============================================================
# 3. 可视化 (双格式 PNG + PDF, 遵守 4 条死规则)
# ============================================================
# --- 3a. ggplot 系图: ggsave_biof3() ---
# p <- ggplot(df, aes(x, y)) +
#   geom_point() +
#   theme_biof3() +
#   scale_color_manual(values = biof3_palette(n))
# ggsave_biof3(p, "result_scatter", width = 6, height = 4.5)
#
# --- 3b. grid 系图 (ComplexHeatmap / pheatmap / base R): save_grid_biof3() ---
# save_grid_biof3("result_heatmap", width = 6, height = 4.5, expr = {
#   pheatmap(mat, color = biof3_palette_div(100))
# })

# ============================================================
# 4. 中间数据保存 (CSV 供下载/展示 + RDS 供后续复用)
# ============================================================
# write.csv(result_df, file.path(output_dir, "result_table.csv"), row.names = FALSE)
# saveRDS(model_obj, file.path(output_dir, "model.rds"))

# ============================================================
# 5. Summary (summary.txt)
# ============================================================
summary_text <- sprintf(
  "__TOOL_TITLE__ Analysis Summary\n\nSamples: %d\nFeatures: %d\n\nResults:\n  <metric 1>: %s\n  <metric 2>: %s\n",
  NA, NA, "TODO", "TODO"
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# 6. Manifest (manifest.json) - 核心清单
# ============================================================
cat("[__TOOL_ID__] Generating manifest...\n")

files_list <- list()

# --- 6a. 图形文件 (双格式, plot=PNG 预览, file=PDF 投稿) ---
plot_pairs <- list(
  # c("filename_no_ext", "中文标签")
  # c("result_scatter", "散点图"),
  # c("result_heatmap", "热图")
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

# --- 6b. 报告 ---
files_list <- c(files_list, list(list(name = "report.html", type = "file", label = "解读报告(HTML)")))

# --- 6c. 数据文件 (table=可预览表格, file=下载, text=纯文本) ---
data_files <- list(
  # c("result_table.csv", "结果表", "table"),
  # c("model.rds",        "模型",    "file"),
  # c("summary.txt",      "分析摘要", "text")
)
for (df in data_files) {
  if (file.exists(file.path(output_dir, df[1]))) {
    files_list <- c(files_list, list(list(name = df[1], type = df[3], label = df[2])))
  }
}

manifest <- list(
  files = files_list,
  summary = list(
    samples   = NA,            # ← 替换为实际值
    features  = NA,
    sci_style = TRUE,
    dual_format = TRUE
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ============================================================
# 7. HTML 报告 (从 template 替换 {{占位符}}, 图片 base64 内嵌)
# ============================================================
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  cat("[__TOOL_ID__] Generating report...\n")
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  # 图片 base64 内嵌
  png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
  for (img_path in png_files) {
    img_name <- basename(img_path)
    b64 <- base64enc::base64encode(img_path)
    data_uri <- paste0("data:image/png;base64,", b64)
    report_html <- gsub(paste0("{{", img_name, "}}"), data_uri, report_html, fixed = TRUE)
  }

  # 文本占位符替换
  # report_html <- gsub("{{samples}}", as.character(n_samples), report_html, fixed = TRUE)

  # 清理未匹配的图片占位符块
  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[__TOOL_ID__] Report generated.\n")
} else {
  cat("[__TOOL_ID__] Report template not found, skipping.\n")
}

cat("[__TOOL_ID__] Pipeline complete.\n")
