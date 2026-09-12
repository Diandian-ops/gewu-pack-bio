#!/usr/bin/env Rscript
# ============================================================
# gene-correlation.R — 基因相关性热图插件
#
# 输入: 表达矩阵 CSV (行=基因, 列=样本)
# 输出: 相关性热图 PNG/PDF + 相关系数矩阵 CSV
# ============================================================

suppressMessages({
  library(pheatmap)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

canonicalize_png <- function(file_path, text = NULL) {
  if (!file.exists(file_path)) return(invisible(FALSE))
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("[gene-correlation] R package 'png' is required for deterministic Result Studio artifacts")
  }
  png::writePNG(png::readPNG(file_path), file_path, text = text)
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("用法: Rscript gene-correlation.R <work_dir>")
work_dir <- args[1]
setwd(work_dir)

params <- fromJSON(file.path(work_dir, "params.json"))
expr_path <- params$expr_matrix
top_n <- as.integer(params$top_n %||% 30)
method <- params$method %||% "pearson"
color_low <- params$color_low %||% "#2563eb"
color_mid <- params$color_mid %||% "#ffffff"
color_high <- params$color_high %||% "#dc2626"
cluster_rows <- isTRUE(params$cluster_rows %||% TRUE)
cluster_cols <- isTRUE(params$cluster_cols %||% TRUE)
show_correlation_values <- isTRUE(params$show_correlation_values %||% FALSE)
base_font_size <- as.numeric(params$base_font_size %||% 8)

# 读表达矩阵（行=基因，列=样本）
expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)
if (ncol(expr) < 2) stop("表达矩阵至少需要 2 列（样本）")
if (nrow(expr) < 2) stop("表达矩阵至少需要 2 行（基因）")

# 按方差取 top N
gene_var <- apply(expr, 1, var, na.rm = TRUE)
top_idx <- order(gene_var, decreasing = TRUE)[1:min(top_n, nrow(expr))]
expr_top <- expr[top_idx, , drop = FALSE]

# 相关系数矩阵
corr_mat <- cor(t(expr_top), method = method, use = "pairwise.complete.obs")

# 保存矩阵
write.csv(as.data.frame(corr_mat), "corr_matrix.csv")

# 热图
png("heatmap.png", width = 7, height = 6, units = "in", res = 300, bg = "white")
rendered_heatmap <- pheatmap(
  corr_mat,
  color = colorRampPalette(c(color_low, color_mid, color_high))(100),
  border_color = NA,
  fontsize_row = base_font_size,
  fontsize_col = base_font_size,
  main = paste0("Gene Correlation (", method, ")"),
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  display_numbers = show_correlation_values,
  number_format = "%.2f"
)
dev.off()
rendered_row_order <- if (inherits(rendered_heatmap$tree_row, "hclust")) {
  rownames(corr_mat)[rendered_heatmap$tree_row$order]
} else {
  rownames(corr_mat)
}
rendered_column_order <- if (inherits(rendered_heatmap$tree_col, "hclust")) {
  colnames(corr_mat)[rendered_heatmap$tree_col$order]
} else {
  colnames(corr_mat)
}
render_contract <- list(
  renderContractVersion = 1,
  artifact = "heatmap",
  palette = list(low = color_low, mid = color_mid, high = color_high),
  clusterRows = cluster_rows,
  clusterCols = cluster_cols,
  renderedRowOrder = rendered_row_order,
  renderedColumnOrder = rendered_column_order,
  showCorrelationValues = show_correlation_values,
  displayedCorrelationValueCount = if (show_correlation_values) length(corr_mat) else 0,
  baseFontSize = base_font_size
)
render_contract_json <- toJSON(render_contract, auto_unbox = TRUE, null = "null")
canonicalize_png("heatmap.png", text = c(BioF3RenderContract = render_contract_json))

pdf("heatmap.pdf", width = 7, height = 6)
pheatmap(
  corr_mat,
  color = colorRampPalette(c(color_low, color_mid, color_high))(100),
  border_color = NA,
  fontsize_row = base_font_size,
  fontsize_col = base_font_size,
  main = paste0("Gene Correlation (", method, ")"),
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  display_numbers = show_correlation_values,
  number_format = "%.2f"
)
dev.off()

# manifest
manifest <- list(
  status = "done",
  outputs = list(
    list(id = "heatmap", filename = "heatmap.png", type = "plot", label = "相关性热图"),
    list(id = "corr_matrix", filename = "corr_matrix.csv", type = "table", label = "相关系数矩阵")
  ),
  summary = list(
    genes = nrow(expr_top),
    method = method,
    display = list(
      color_low = color_low,
      color_mid = color_mid,
      color_high = color_high,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      show_correlation_values = show_correlation_values,
      base_font_size = base_font_size
    )
  )
)
write(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), "manifest.json")

cat("DONE\n")
