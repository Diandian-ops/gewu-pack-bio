#!/usr/bin/env Rscript
# ============================================================
# Tool: complexheatmap
# Version targeted: ComplexHeatmap 2.24.1 (R 4.5.1)
# Canonical workflow source:
#   https://jokergoo.github.io/ComplexHeatmap-reference/book/
# Reference:
#   Gu Z. (2022). Complex heatmap visualization. iMeta 1(3):e43.
#
# 这是 BioF3 SCI 视觉风格规范 (theme_biof3 + ggsave_biof3 + save_grid_biof3)
# 的首个落地实例 — 全图双格式输出 (PNG + PDF), 直接发表水平.
#
# Steps (in-scope) — 标准热图工作流 5 步:
#   1. Load + 检查输入
#   2. Z-score 标准化 (可选)
#   3. 行/列聚类
#   4. Annotation (categorical → 离散色 / numeric → 渐变)
#   5. 输出: 主热图 + 未标准化版 + 行/列树 + 相关性矩阵
#
# Out-of-scope:
#   - OncoPrint (留给 maftools)
#   - UpSet plot (独立工具)
#   - density / interactive heatmap
#
# Default params (按官方 + BioF3 SCI 默认):
#   z_score=TRUE
#   cluster_rows=TRUE, cluster_columns=TRUE
#   clustering_distance="pearson", clustering_method="ward.D2"
#   color_scheme="div" (蓝-白-红)
#   top_n_label=20 (行 > 50 时仅标记 top variable genes)
#
# Input (job_dir/):
#   - matrix      CSV/TSV: 行=基因 / 特征, 列=样本
#   - annotation  (可选) CSV/TSV: 行=样本, 列=分组属性
#   - params.json
#
# Output (job_dir/output/):
#   5 plots × 双格式 (PNG + PDF) = 10 files
#   + 4 中间数据 + summary.txt + manifest.json + report.html
# ============================================================

suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(jsonlite)
  library(matrixStats)
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题 (tools.js 自动复制到 jobDir)
source(file.path(job_dir, "_biof3-theme.R"))

cat("[complexheatmap] Starting pipeline\n")
cat("[complexheatmap] R:", R.version.string, "\n")
cat("[complexheatmap] ComplexHeatmap:", as.character(packageVersion("ComplexHeatmap")), "\n")

# ---- Read params ----
params <- fromJSON(file.path(job_dir, "params.json"))
z_score             <- as.logical(params$z_score %||% TRUE)
cluster_rows        <- as.logical(params$cluster_rows %||% TRUE)
cluster_columns     <- as.logical(params$cluster_columns %||% TRUE)
clustering_distance <- params$clustering_distance %||% "pearson"
clustering_method   <- params$clustering_method %||% "ward.D2"
color_scheme        <- params$color_scheme %||% "div"
top_n_label         <- as.integer(params$top_n_label %||% 20)

cat(sprintf("[complexheatmap] params: z_score=%s, cluster_rows=%s, cluster_columns=%s, dist=%s, method=%s, color=%s\n",
            z_score, cluster_rows, cluster_columns, clustering_distance, clustering_method, color_scheme))

# ---- Read matrix ----
read_table_auto <- function(file, row_names_col = 1) {
  first_line <- readLines(file, n = 1, warn = FALSE)
  sep <- if (grepl("\t", first_line) && !grepl(",", first_line)) "\t" else ","
  read.table(file, header = TRUE, sep = sep, row.names = row_names_col,
             check.names = FALSE, stringsAsFactors = FALSE)
}

matrix_files <- list.files(job_dir, pattern = "^matrix$", full.names = TRUE)
if (length(matrix_files) == 0) stop("No matrix file found in job_dir")

cat("[complexheatmap] Loading matrix...\n")
mat_raw <- as.matrix(read_table_auto(matrix_files[1]))
storage.mode(mat_raw) <- "numeric"

cat(sprintf("  matrix: %d rows × %d columns\n", nrow(mat_raw), ncol(mat_raw)))

# ---- Validate dimensions ----
if (nrow(mat_raw) > 5000 || ncol(mat_raw) > 500) {
  stop(sprintf("Matrix too large: %d rows × %d cols. Max: 5000 × 500 (pre-filter to high-variance rows / subset samples).",
               nrow(mat_raw), ncol(mat_raw)))
}

# ---- Filter zero-variance rows (causes NaN in z-score / pearson dist) ----
row_vars <- rowVars(mat_raw, na.rm = TRUE)
keep_rows <- !is.na(row_vars) & row_vars > 0
n_removed_zero_var <- sum(!keep_rows)
if (n_removed_zero_var > 0) {
  cat(sprintf("  Removed %d zero-variance rows\n", n_removed_zero_var))
  mat_raw <- mat_raw[keep_rows, ]
}

# ---- Read annotation (optional) ----
anno_files <- list.files(job_dir, pattern = "^annotation$", full.names = TRUE)
top_annotation <- NULL
anno_df_clean <- NULL
if (length(anno_files) > 0) {
  cat("[complexheatmap] Loading annotation...\n")
  anno_raw <- read_table_auto(anno_files[1])

  # Match samples
  common_samples <- intersect(colnames(mat_raw), rownames(anno_raw))
  cat(sprintf("  annotation: %d cols, matched %d / %d samples\n",
              ncol(anno_raw), length(common_samples), ncol(mat_raw)))
  if (length(common_samples) < ncol(mat_raw) * 0.5) {
    cat("  WARN: < 50% samples have annotation, dropping annotation\n")
  } else if (length(common_samples) >= 1) {
    # Keep only matched samples in BOTH mat and anno
    mat_raw <- mat_raw[, common_samples, drop = FALSE]
    anno_raw <- anno_raw[common_samples, , drop = FALSE]

    # Hard-cap to 6 annotation columns (avoid visual explosion)
    if (ncol(anno_raw) > 6) {
      cat(sprintf("  WARN: %d annotation columns > 6, keeping first 6\n", ncol(anno_raw)))
      anno_raw <- anno_raw[, 1:6, drop = FALSE]
    }
    anno_df_clean <- anno_raw

    # Build annotation: numeric → continuous; character/factor → discrete
    anno_args <- list()
    anno_colors <- list()
    pal_cycle <- biof3_palette(20)
    color_idx <- 1
    for (col in colnames(anno_raw)) {
      val <- anno_raw[[col]]
      if (is.numeric(val)) {
        anno_args[[col]] <- val
        # Continuous color ramp via colorRamp2 (function form, ComplexHeatmap-compatible)
        rng <- range(val, na.rm = TRUE)
        mid <- mean(rng)
        anno_colors[[col]] <- circlize::colorRamp2(c(rng[1], mid, rng[2]),
                                                     biof3_palette_div(3))
      } else {
        val <- as.character(val)
        anno_args[[col]] <- val
        # Discrete colors (cycle palette) — must be NAMED vector
        levels_unique <- unique(val[!is.na(val)])
        named_colors <- pal_cycle[((color_idx - 1):(color_idx - 2 + length(levels_unique))) %% 20 + 1]
        anno_colors[[col]] <- setNames(named_colors[seq_along(levels_unique)], levels_unique)
        color_idx <- color_idx + length(levels_unique)
      }
    }
    top_annotation <- do.call(HeatmapAnnotation,
                               c(anno_args, list(col = anno_colors,
                                                  annotation_name_side = "left",
                                                  annotation_name_gp = gpar(fontsize = 9))))
  }
}

# ---- Z-score normalization ----
if (z_score) {
  cat("[complexheatmap] Z-score normalizing rows...\n")
  mat <- t(scale(t(mat_raw)))
  # Z-score may produce NaN for rows that had near-zero variance after rounding
  mat[is.nan(mat)] <- 0
} else {
  mat <- mat_raw
}

# Save z-score matrix
mat_out <- data.frame(feature = rownames(mat), mat, check.names = FALSE)
write.csv(mat_out, file.path(output_dir, "matrix_zscore.csv"), row.names = FALSE)
write.csv(data.frame(feature = rownames(mat_raw), mat_raw, check.names = FALSE),
          file.path(output_dir, "matrix_unscaled.csv"), row.names = FALSE)

# ---- Decide whether to show row names ----
n_rows <- nrow(mat)
show_row_names_flag <- n_rows <= 50

# Mark top N variable genes (for big matrix when row names are hidden)
top_var_genes <- character(0)
if (!show_row_names_flag && n_rows > top_n_label) {
  row_vars_z <- rowVars(mat, na.rm = TRUE)
  top_idx <- order(row_vars_z, decreasing = TRUE)[1:top_n_label]
  top_var_genes <- rownames(mat)[top_idx]
  cat(sprintf("  Marking top %d variable genes (matrix has %d rows)\n", top_n_label, n_rows))
}

# ---- Color function ----
col_fun <- biof3_heat_colors(mat, kind = if (color_scheme == "div") "div" else "seq")

# ---- Build main heatmap ----
ht_main_args <- list(
  matrix = mat,
  col = col_fun,
  cluster_rows = cluster_rows,
  cluster_columns = cluster_columns,
  clustering_distance_rows = clustering_distance,
  clustering_distance_columns = clustering_distance,
  clustering_method_rows = clustering_method,
  clustering_method_columns = clustering_method,
  show_row_names = show_row_names_flag,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = if (z_score) "Z-score" else "Value",
    legend_height = unit(3, "cm"),
    title_position = "leftcenter-rot"
  ),
  top_annotation = top_annotation
)

# Add right annotation marking top variable genes
if (length(top_var_genes) > 0) {
  ht_main_args$right_annotation <- rowAnnotation(
    foo = anno_mark(at = match(top_var_genes, rownames(mat)),
                     labels = top_var_genes,
                     labels_gp = gpar(fontsize = 8))
  )
}

ht_main <- do.call(Heatmap, ht_main_args)

# ---- Save main heatmap (PNG + PDF) ----
heat_w <- max(8, ncol(mat) * 0.08 + 4)
heat_h <- max(6, min(20, nrow(mat) * 0.05 + 4))

cat(sprintf("[complexheatmap] Saving main heatmap %g × %g in (PNG + PDF)...\n", heat_w, heat_h))

# Save main heatmap (this draw call also fixes the row/column order)
save_grid_biof3(file.path(output_dir, "heatmap_main"),
                 width = heat_w, height = heat_h, expr = {
  ht_drawn <- draw(ht_main, merge_legend = TRUE,
                    heatmap_legend_side = "right",
                    annotation_legend_side = "right")
  # Capture order on first draw (within the open device); store in parent for later use
  if (cluster_rows) {
    ro <- row_order(ht_drawn); if (is.list(ro)) ro <- unlist(ro)
    .GlobalEnv$.row_order <- ro
  }
  if (cluster_columns) {
    co <- column_order(ht_drawn); if (is.list(co)) co <- unlist(co)
    .GlobalEnv$.col_order <- co
  }
})

# ---- Save reordered matrix (for users to inspect / re-plot) ----
final_row_order <- if (exists(".row_order", .GlobalEnv) && length(.GlobalEnv$.row_order) == nrow(mat)) .GlobalEnv$.row_order else seq_len(nrow(mat))
final_col_order <- if (exists(".col_order", .GlobalEnv) && length(.GlobalEnv$.col_order) == ncol(mat)) .GlobalEnv$.col_order else seq_len(ncol(mat))
mat_clustered <- mat[final_row_order, final_col_order]
mat_clustered_out <- data.frame(feature = rownames(mat_clustered), mat_clustered, check.names = FALSE)
write.csv(mat_clustered_out, file.path(output_dir, "matrix_clustered.csv"), row.names = FALSE)

# ---- Save unscaled heatmap (no z-score, raw values) ----
if (z_score) {
  # Make a separate unscaled version
  col_fun_raw <- biof3_heat_colors(mat_raw, kind = "seq")
  ht_unscaled <- Heatmap(mat_raw,
                          col = col_fun_raw,
                          cluster_rows = cluster_rows,
                          cluster_columns = cluster_columns,
                          clustering_distance_rows = clustering_distance,
                          clustering_distance_columns = clustering_distance,
                          clustering_method_rows = clustering_method,
                          clustering_method_columns = clustering_method,
                          show_row_names = show_row_names_flag,
                          show_column_names = TRUE,
                          row_names_gp = gpar(fontsize = 8),
                          column_names_gp = gpar(fontsize = 8),
                          column_names_rot = 45,
                          heatmap_legend_param = list(title = "Value", legend_height = unit(3, "cm")),
                          top_annotation = top_annotation)
  save_grid_biof3(file.path(output_dir, "heatmap_unscaled"),
                   width = heat_w, height = heat_h, expr = {
    draw(ht_unscaled, merge_legend = TRUE)
  })
} else {
  # If z_score=FALSE, the main IS unscaled — produce a copy for consistency
  save_grid_biof3(file.path(output_dir, "heatmap_unscaled"),
                   width = heat_w, height = heat_h, expr = {
    draw(ht_main, merge_legend = TRUE)
  })
}

# ---- Sample correlation heatmap (QC) ----
cat("[complexheatmap] Sample correlation matrix...\n")
sample_cor <- cor(mat_raw, method = "pearson", use = "p")
ht_cor <- Heatmap(sample_cor,
                   col = colorRamp2(c(min(sample_cor, 0.5), 0.75, 1), biof3_palette_div(3)),
                   cluster_rows = TRUE, cluster_columns = TRUE,
                   clustering_method_rows = "ward.D2",
                   clustering_method_columns = "ward.D2",
                   show_row_names = ncol(mat_raw) <= 60,
                   show_column_names = ncol(mat_raw) <= 60,
                   row_names_gp = gpar(fontsize = 7),
                   column_names_gp = gpar(fontsize = 7),
                   column_names_rot = 45,
                   heatmap_legend_param = list(title = "Pearson r"),
                   top_annotation = top_annotation)
save_grid_biof3(file.path(output_dir, "correlation_heatmap"),
                 width = max(7, ncol(mat_raw) * 0.08 + 3),
                 height = max(6, ncol(mat_raw) * 0.08 + 3), expr = {
  draw(ht_cor, merge_legend = TRUE)
})

# ---- Standalone dendrograms (debug / supplementary) ----
if (cluster_rows && nrow(mat) <= 200) {  # only draw if rows ≤ 200 (else illegible)
  row_dist <- as.dist(1 - cor(t(mat), method = clustering_distance))
  if (clustering_distance == "euclidean") row_dist <- dist(mat)
  row_hc <- hclust(row_dist, method = clustering_method)
  save_grid_biof3(file.path(output_dir, "dendrogram_rows"),
                   width = 12, height = max(4, nrow(mat) * 0.04 + 2), expr = {
    par(mar = c(0, 4, 2, 1), cex = 0.5)
    plot(as.dendrogram(row_hc), horiz = FALSE,
         main = "Row (feature) clustering")
  })
}
if (cluster_columns) {
  col_dist <- as.dist(1 - cor(mat, method = clustering_distance))
  if (clustering_distance == "euclidean") col_dist <- dist(t(mat))
  col_hc <- hclust(col_dist, method = clustering_method)
  save_grid_biof3(file.path(output_dir, "dendrogram_columns"),
                   width = max(8, ncol(mat) * 0.1 + 2), height = 5, expr = {
    par(mar = c(8, 4, 2, 1), cex = 0.6)
    plot(as.dendrogram(col_hc), horiz = FALSE,
         main = "Column (sample) clustering")
  })
}

# ---- Summary ----
clustering_summary <- sprintf(
  "ComplexHeatmap Pipeline Summary\n\nInput:\n  Matrix: %d rows × %d cols (after dropping %d zero-variance rows)\n  Annotation: %s\n\nProcessing:\n  Z-score normalized: %s\n  Row names shown: %s (n_rows=%d, threshold=50)\n  Top variable genes marked: %d (threshold=%d when row names hidden)\n\nClustering:\n  Cluster rows: %s\n  Cluster columns: %s\n  Distance: %s\n  Linkage method: %s\n\nColor scheme: %s\n\nOutputs:\n  heatmap_main.{png,pdf}    Main heatmap (z-score + annotation + top markers)\n  heatmap_unscaled.{png,pdf} Unscaled heatmap\n  correlation_heatmap.{png,pdf} Sample correlation matrix\n  dendrogram_rows.{png,pdf}  Row dendrogram (if rows ≤ 200)\n  dendrogram_columns.{png,pdf} Column dendrogram\n  matrix_zscore.csv         Z-score matrix\n  matrix_clustered.csv      Reordered after clustering\n  matrix_unscaled.csv       Original matrix (post-filter)\n",
  nrow(mat_raw), ncol(mat_raw), n_removed_zero_var,
  if (is.null(anno_df_clean)) "none" else sprintf("%d cols", ncol(anno_df_clean)),
  z_score, show_row_names_flag, nrow(mat),
  length(top_var_genes), top_n_label,
  cluster_rows, cluster_columns, clustering_distance, clustering_method,
  color_scheme
)
writeLines(clustering_summary, file.path(output_dir, "summary.txt"))
writeLines(clustering_summary, file.path(output_dir, "clustering_summary.txt"))
cat(clustering_summary)

# ---- Manifest ----
manifest <- list(
  files = list(
    list(name = "heatmap_main.png",         type = "plot",  label = "主热图 (PNG)"),
    list(name = "heatmap_main.pdf",         type = "file",  label = "主热图 (PDF, 矢量)"),
    list(name = "heatmap_unscaled.png",     type = "plot",  label = "未标准化热图 (PNG)"),
    list(name = "heatmap_unscaled.pdf",     type = "file",  label = "未标准化热图 (PDF)"),
    list(name = "correlation_heatmap.png",  type = "plot",  label = "样本相关性 (PNG)"),
    list(name = "correlation_heatmap.pdf",  type = "file",  label = "样本相关性 (PDF)"),
    list(name = "dendrogram_columns.png",   type = "plot",  label = "列聚类树 (PNG)"),
    list(name = "dendrogram_columns.pdf",   type = "file",  label = "列聚类树 (PDF)"),
    list(name = "dendrogram_rows.png",      type = "plot",  label = "行聚类树 (PNG)"),
    list(name = "dendrogram_rows.pdf",      type = "file",  label = "行聚类树 (PDF)"),
    list(name = "matrix_zscore.csv",        type = "table", label = "Z-score 矩阵"),
    list(name = "matrix_unscaled.csv",      type = "table", label = "原始矩阵"),
    list(name = "matrix_clustered.csv",     type = "table", label = "聚类排序后矩阵"),
    list(name = "summary.txt",              type = "text",  label = "聚类参数摘要"),
    list(name = "report.html",              type = "file",  label = "解读报告 HTML")
  ),
  summary = list(
    n_rows_input = nrow(mat_raw) + n_removed_zero_var,
    n_rows_final = nrow(mat),
    n_cols = ncol(mat),
    z_score = z_score,
    n_clusters_modes = ifelse(cluster_rows && cluster_columns, "both",
                       ifelse(cluster_rows, "rows", ifelse(cluster_columns, "columns", "none"))),
    has_annotation = !is.null(anno_df_clean)
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ---- HTML report ----
cat("[complexheatmap] Generating HTML report...\n")
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  png_files <- c("heatmap_main.png", "heatmap_unscaled.png", "correlation_heatmap.png",
                 "dendrogram_columns.png", "dendrogram_rows.png")
  for (png in png_files) {
    img_path <- file.path(output_dir, png)
    placeholder <- paste0("{{", png, "}}")
    if (file.exists(img_path)) {
      b64 <- base64enc::base64encode(img_path)
      data_uri <- paste0("data:image/png;base64,", b64)
      report_html <- gsub(placeholder, data_uri, report_html, fixed = TRUE)
    }
  }

  report_html <- gsub("{{n_rows_input}}",  as.character(nrow(mat_raw) + n_removed_zero_var), report_html, fixed = TRUE)
  report_html <- gsub("{{n_rows_final}}",  as.character(nrow(mat)),                          report_html, fixed = TRUE)
  report_html <- gsub("{{n_cols}}",        as.character(ncol(mat)),                          report_html, fixed = TRUE)
  report_html <- gsub("{{z_score}}",       if (z_score) "TRUE" else "FALSE",                 report_html, fixed = TRUE)
  report_html <- gsub("{{cluster_rows}}",  if (cluster_rows) "TRUE" else "FALSE",            report_html, fixed = TRUE)
  report_html <- gsub("{{cluster_columns}}", if (cluster_columns) "TRUE" else "FALSE",       report_html, fixed = TRUE)
  report_html <- gsub("{{clustering_distance}}", clustering_distance,                         report_html, fixed = TRUE)
  report_html <- gsub("{{clustering_method}}",   clustering_method,                           report_html, fixed = TRUE)
  report_html <- gsub("{{color_scheme}}", color_scheme,                                       report_html, fixed = TRUE)

  # Remove fig blocks for any unreplaced {{*.png}}
  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)

  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[complexheatmap] Report generated.\n")
} else {
  cat("[complexheatmap] WARN: report template not found.\n")
}

cat("[complexheatmap] Pipeline complete.\n")
