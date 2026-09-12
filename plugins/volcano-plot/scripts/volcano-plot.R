#!/usr/bin/env Rscript
# ============================================================
# volcano-plot.R — 火山图绘制插件
#
# 输入: CSV (gene, log2FoldChange|log2FC, pvalue[, padj])
#      lfc_col 参数指定 LFC 列名，默认自动探测：
#      - log2FoldChange（deg-standardizer 输出，canonical）
#      - log2FC（DESeq2 直接输出）
# 输出: 火山图 PNG/PDF + 分类结果表 CSV
# ============================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("用法: Rscript volcano-plot.R <work_dir>")
work_dir <- args[1]
setwd(work_dir)

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

canonicalize_png <- function(file_path, text = NULL) {
  if (!file.exists(file_path)) return(invisible(FALSE))
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("[volcano-plot] R package 'png' is required for deterministic Result Studio artifacts")
  }
  png::writePNG(png::readPNG(file_path), file_path, text = text)
  invisible(TRUE)
}

# 读参数
params <- fromJSON(file.path(work_dir, "params.json"))
deg_path <- params$deg_table
gene_col  <- params$gene_col  %||% params$gene_id_column  %||% "gene"
lfc_col   <- params$lfc_col   %||% params$lfc_column     %||% NULL
pval_col  <- params$pval_col  %||% params$pval_column    %||% "padj"
fc_threshold <- as.numeric(params$fc_threshold %||% 1)
pval_threshold <- as.numeric(params$pval_threshold %||% 1.3)
top_n <- as.integer(params$top_n %||% 10)
# P5 selection-aware presentation input.  This is only a bounded list of
# stable gene IDs to display; it never changes the test statistic, threshold
# or Up/Down/NS classification below.
label_genes_raw <- as.character(params$label_genes %||% "")

# G0-RS-P4 presentation 参数。
# 这些只影响呈现，不参与显著性判定：阈值与分类逻辑仍只由 fc_threshold /
# pval_threshold 决定，改配色或字号不会改变任何科学结论。
plot_title <- params$plot_title %||% "Volcano Plot"
color_up <- params$color_up %||% "#dc2626"
color_down <- params$color_down %||% "#2563eb"
color_ns <- params$color_ns %||% "#94a3b8"
point_size <- as.numeric(params$point_size %||% 1.8)
point_alpha <- as.numeric(params$point_alpha %||% 0.6)
label_size <- as.numeric(params$label_size %||% 3)
base_font_size <- as.numeric(params$base_font_size %||% 12)
legend_position <- params$legend_position %||% "right"
plot_width <- as.numeric(params$plot_width %||% 7)
plot_height <- as.numeric(params$plot_height %||% 5)
export_dpi <- as.numeric(params$export_dpi %||% 300)

# 呈现参数越界时收敛到契约声明的范围，不静默接受非法值也不中断科学计算。
clamp <- function(value, low, high, fallback) {
  if (!is.finite(value)) return(fallback)
  min(max(value, low), high)
}
is_hex_color <- function(value) {
  is.character(value) && length(value) == 1 && grepl("^#[0-9A-Fa-f]{6}$", value)
}
if (!is_hex_color(color_up)) color_up <- "#dc2626"
if (!is_hex_color(color_down)) color_down <- "#2563eb"
if (!is_hex_color(color_ns)) color_ns <- "#94a3b8"
point_size <- clamp(point_size, 0.2, 6, 1.8)
point_alpha <- clamp(point_alpha, 0.05, 1, 0.6)
label_size <- clamp(label_size, 1, 8, 3)
base_font_size <- clamp(base_font_size, 6, 24, 12)
plot_width <- clamp(plot_width, 3, 20, 7)
plot_height <- clamp(plot_height, 3, 20, 5)
export_dpi <- clamp(export_dpi, 72, 1200, 300)
if (!(legend_position %in% c("right", "left", "top", "bottom", "none"))) legend_position <- "right"
plot_title <- if (is.character(plot_title) && length(plot_title) == 1) substr(plot_title, 1, 120) else "Volcano Plot"

# 读数据
deg <- read.csv(deg_path, stringsAsFactors = FALSE)

# 列名探测：canonical = log2FoldChange (deg-standardizer)，fallback = log2FC (DESeq2 direct)
known_lfc_names <- c("log2FoldChange", "log2FC", "logFC", "avg_log2FC", "log2.Ratio")
if (is.null(lfc_col)) {
  found_lfc <- intersect(known_lfc_names, names(deg))
  if (length(found_lfc) == 0) {
    stop(paste0(
      "[volcano-plot] 无法找到 LFC 列。期望列名: ",
      paste(known_lfc_names, collapse = ", "),
      "。实际列名: ",
      paste(names(deg), collapse = ", ")
    ))
  }
  lfc_col <- found_lfc[1]
  cat(sprintf("[volcano-plot] 探测到 LFC 列: '%s'\n", lfc_col))
} else {
  if (!(lfc_col %in% names(deg))) {
    stop(paste0(
      "[volcano-plot] 指定的 LFC 列 '", lfc_col, "' 不存在。实际列名: ",
      paste(names(deg), collapse = ", ")
    ))
  }
}

# gene 列名参数化（默认为 gene）
if (!(gene_col %in% names(deg))) {
  stop(paste0("[volcano-plot] 指定的 gene 列 '", gene_col, "' 不存在。实际列名: ", paste(names(deg), collapse = ", ")))
}
# 重命名基因列
if (gene_col != "gene") {
  names(deg)[names(deg) == gene_col] <- "gene"
}

# 标准化命名：内部统一用 log2FoldChange，输出用 canonical 名称
deg$log2FoldChange <- as.numeric(deg[[lfc_col]])
if (any(is.na(deg$log2FoldChange)) || any(!is.finite(deg$log2FoldChange))) {
  stop(sprintf("LFC 列 '%s' 含缺失、非数值或无限值，拒绝生成火山图", lfc_col))
}

# 显著性列必须由调用方明确选择；默认是多重校正后的 padj。
# 不能在缺列时静默回退到裸 p 值，否则会把假阳性画成论文结论。
if (!(pval_col %in% names(deg))) {
  if (pval_col == "padj" && "adj.P.Val" %in% names(deg)) {
    pval_col <- "adj.P.Val"
  } else {
    stop(sprintf("显著性列 '%s' 不存在；火山图默认要求校正后的 padj", pval_col))
  }
}
parse_probability_column <- function(column_name) {
  values <- suppressWarnings(as.numeric(deg[[column_name]]))
  if (any(is.na(values)) || any(!is.finite(values))) {
    stop(sprintf("显著性列 '%s' 含缺失、非数值或无限值，拒绝生成火山图", column_name))
  }
  if (any(values < 0 | values > 1)) {
    stop(sprintf("显著性列 '%s' 必须位于 0 到 1 之间", column_name))
  }
  values
}
selected_significance <- parse_probability_column(pval_col)
raw_pvalue <- if ("pvalue" %in% names(deg)) parse_probability_column("pvalue") else NULL
adjusted_pvalue <- if ("padj" %in% names(deg)) parse_probability_column("padj") else NULL
deg$neg_log10_significance <- -log10(selected_significance + 1e-300)
deg$sig <- "NS"
deg$sig[deg$log2FoldChange > fc_threshold & deg$neg_log10_significance > pval_threshold] <- "Up"
deg$sig[deg$log2FoldChange < -fc_threshold & deg$neg_log10_significance > pval_threshold] <- "Down"
deg$sig <- factor(deg$sig, levels = c("Up", "NS", "Down"))

# top N 标注
deg$rank_score <- abs(deg$log2FoldChange) * deg$neg_log10_significance
deg$label <- ""
top_idx <- order(deg$rank_score, decreasing = TRUE)[1:min(top_n, nrow(deg))]
deg$label[top_idx] <- deg$gene[top_idx]
label_genes <- unique(trimws(unlist(strsplit(label_genes_raw, ",", fixed = TRUE))))
label_genes <- label_genes[nzchar(label_genes)]
if (length(label_genes) > 0) {
  selected_idx <- which(deg$gene %in% label_genes)
  deg$label[selected_idx] <- deg$gene[selected_idx]
}
rendered_label_genes <- sort(unique(as.character(deg$gene[deg$label != ""])))

# 配色
colors <- c(Up = color_up, NS = color_ns, Down = color_down)

# 绘图
p <- ggplot(deg, aes(x = log2FoldChange, y = neg_log10_significance, color = sig)) +
  geom_point(alpha = point_alpha, size = point_size) +
  geom_text_repel(
    aes(label = label),
    data = subset(deg, label != ""),
    size = label_size, max.overlaps = 20,
    box.padding = 0.4, segment.color = "grey50", seed = 20260729
  ) +
  scale_color_manual(values = colors, name = "") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  labs(
    x = expression(log[2](Fold~Change)),
    y = bquote(-log[10](.(pval_col))),
    title = plot_title
  ) +
  theme_minimal(base_size = base_font_size) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    legend.position = legend_position
  )

# 保存
render_contract <- list(
  renderContractVersion = 1,
  artifact = "volcano",
  renderedLabelGenes = rendered_label_genes,
  selectedLabelGenes = sort(label_genes),
  palette = list(up = color_up, down = color_down, nonSignificant = color_ns),
  pointSize = point_size,
  pointAlpha = point_alpha,
  labelSize = label_size,
  baseFontSize = base_font_size,
  width = plot_width,
  height = plot_height,
  dpi = export_dpi,
  legendPosition = legend_position
)
render_contract_json <- toJSON(render_contract, auto_unbox = TRUE, null = "null")
ggsave("volcano.png", p, width = plot_width, height = plot_height, dpi = export_dpi, bg = "white")
canonicalize_png("volcano.png", text = c(BioF3RenderContract = render_contract_json))
svg("volcano.svg", width = plot_width, height = plot_height, bg = "white")
print(p)
dev.off()
svg_lines <- readLines("volcano.svg", warn = FALSE)
svg_start <- grep("<svg\\b", svg_lines, perl = TRUE)[1]
if (is.na(svg_start)) stop("[volcano-plot] generated SVG lacks an svg root")
render_contract_hex <- paste(sprintf("%02x", as.integer(charToRaw(render_contract_json))), collapse = "")
svg_lines <- append(
  svg_lines,
  sprintf('<metadata id="biof3-render-contract" data-json-hex="%s"/>', render_contract_hex),
  after = svg_start
)
writeLines(svg_lines, "volcano.svg", useBytes = TRUE)
ggsave("volcano.pdf", p, width = plot_width, height = plot_height)

# 分类表（输出 canonical 列名）
classified <- deg[, c("gene", "log2FoldChange", "sig")]
classified$significance <- selected_significance
classified$significance_column <- pval_col
if (!is.null(raw_pvalue)) classified$pvalue <- raw_pvalue
if (!is.null(adjusted_pvalue)) classified$padj <- adjusted_pvalue
write.csv(classified, "classified_table.csv", row.names = FALSE)

summary_counts <- list(
  up = sum(deg$sig == "Up"),
  down = sum(deg$sig == "Down"),
  ns = sum(deg$sig == "NS")
)

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

top_rows <- head(classified[order(-abs(classified$log2FoldChange), classified$significance), ], 12)
top_gene_rows_html <- paste(apply(top_rows, 1, function(row) {
  paste0(
    "<tr>",
    "<td>", html_escape(row[["gene"]]), "</td>",
    "<td>", sprintf("%.3f", as.numeric(row[["log2FoldChange"]])), "</td>",
    "<td>", format(as.numeric(row[["significance"]]), scientific = TRUE, digits = 3), "</td>",
    "<td><span class=\"badge badge-", tolower(html_escape(row[["sig"]])), "\">", html_escape(row[["sig"]]), "</span></td>",
    "</tr>"
  )
}), collapse = "\n")

replace_token <- function(text, token, value) {
  gsub(paste0("{{", token, "}}"), as.character(value), text, fixed = TRUE)
}

report_generated <- FALSE
template_path <- file.path(work_dir, "report-template.html")
if (file.exists(template_path)) {
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
  replacements <- list(
    up_count = summary_counts$up,
    down_count = summary_counts$down,
    ns_count = summary_counts$ns,
    total_count = nrow(deg),
    fc_threshold = fc_threshold,
    pval_threshold = pval_threshold,
    significance_column = html_escape(pval_col),
    top_n = top_n,
    top_gene_rows = top_gene_rows_html
  )
  for (name in names(replacements)) {
    report_html <- replace_token(report_html, name, replacements[[name]])
  }
  writeLines(report_html, "report.html", useBytes = TRUE)
  report_generated <- TRUE
}

# manifest
manifest <- list(
  status = "done",
  outputs = list(
    list(id = "volcano", filename = "volcano.png", type = "plot", label = "火山图"),
    list(id = "volcano_svg", filename = "volcano.svg", type = "plot", label = "火山图 SVG"),
    list(id = "classified_table", filename = "classified_table.csv", type = "table", label = "分类结果表")
  ),
  summary = c(summary_counts, list(
    significance_column = pval_col,
    presentation = list(
      label_genes = label_genes,
      color_up = color_up,
      color_down = color_down,
      point_size = point_size,
      label_size = label_size,
      plot_width = plot_width,
      plot_height = plot_height,
      export_dpi = export_dpi
    )
  )),
  title = "Honeycomb 火山图 Demo 结果",
  files = list(
    list(name = "volcano.png", type = "plot", label = "火山图", title = "火山图", preview = TRUE, download = TRUE),
    list(name = "volcano.svg", type = "plot", label = "火山图 SVG", title = "火山图 SVG", preview = TRUE, download = TRUE),
    list(name = "classified_table.csv", type = "table", label = "分类结果表", title = "分类结果表", preview = TRUE, download = TRUE),
    list(name = "manifest.json", type = "text", label = "结果 Manifest", title = "结果 Manifest", preview = TRUE, download = TRUE),
    list(name = "volcano.pdf", type = "file", label = "火山图 PDF", title = "火山图 PDF", preview = FALSE, download = TRUE)
  ),
  sections = list(
    list(
      id = "report",
      label = "预览报告",
      artifactNames = list("report.html"),
      description = "汇总图表、统计和关键输出的主报告。"
    ),
    list(
      id = "figures",
      label = "图表结果",
      artifactNames = list("volcano.png", "volcano.svg", "volcano.pdf"),
      description = "火山图 PNG、SVG 与可下载 PDF。"
    ),
    list(
      id = "tables",
      label = "结果表格",
      artifactNames = list("classified_table.csv"),
      description = "差异分类结果。"
    ),
    list(
      id = "records",
      label = "运行记录",
      artifactNames = list("manifest.json"),
      description = "结果 manifest 与输出结构记录。"
    )
  )
)
if (report_generated) {
  manifest$outputs <- c(
    list(list(id = "report", filename = "report.html", type = "report", label = "预览报告")),
    manifest$outputs
  )
  manifest$files <- c(
    list(list(name = "report.html", type = "report", label = "预览报告", title = "预览报告", preview = TRUE, download = TRUE)),
    manifest$files
  )
} else {
  manifest$sections <- manifest$sections[2:length(manifest$sections)]
}
write(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), "manifest.json")

cat("DONE\n")
