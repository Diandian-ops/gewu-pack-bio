#!/usr/bin/env Rscript
# ============================================================
# DESeq2 差异分析工具脚本
# Version targeted: DESeq2 1.46+ (R 4.5.x)
# Canonical workflow source:
#   https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#
# Steps (in-scope) — Bioconductor 主 vignette 主流程:
#   1. read counts + coldata
#   2. low-count filter
#   3. DESeq() 拟合
#   4. results() + lfcShrink(apeglm)
#   5. plotPCA / plotDispEsts / plotMA / 火山图 / 热图 / sample distance / library size
#   6. 中间数据保存 (vst_matrix / coldata_info / dispersion_data)
#   7. HTML 报告 + manifest
#
# 2026-05-25: 升级到 BioF3 SCI 视觉风格规范
#   - 全图 PNG + PDF 双格式 (300dpi PNG, Cairo PDF 矢量)
#   - ggplot2 系图: theme_biof3() + biof3_palette() + ggsave_biof3()
#   - grid 系图 (pheatmap, plotMA, plotDispEsts): save_grid_biof3()
#   - 配色: biof3_palette_div() (蓝-白-红, 火山/热图二极用)
# ============================================================

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(jsonlite)
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

canonicalize_png <- function(file_path) {
  if (!file.exists(file_path)) return(invisible(FALSE))
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("[deseq2] R package 'png' is required for deterministic Result Studio artifacts")
  }
  png::writePNG(png::readPNG(file_path), file_path)
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题 (tools.js 自动复制到 jobDir)
source(file.path(job_dir, "_biof3-theme.R"))

cat("[deseq2] Starting pipeline (SCI style v2026-05-25)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
normalize_design_col <- function(x) {
  x <- trimws(as.character(x %||% "condition"))
  x <- sub("^~\\s*", "", x)
  if (!nzchar(x)) "condition" else x
}

design_col <- normalize_design_col(params$design %||% params$design_col %||% "condition")
ref <- params$contrast_ref %||% params$ref %||% "Control"
treat <- params$contrast_treat %||% params$treat %||% "Treatment"
padj_cut <- as.numeric(params$padj_cutoff %||% params$padj_threshold %||% 0.05)
lfc_cut <- as.numeric(params$lfc_cutoff %||% params$lfc_threshold %||% 1)
min_count <- as.numeric(params$min_count %||% 10)

# Read data
read_matrix <- function(f) {
  if (grepl("\\.tsv$|\\.txt$", f)) read.delim(f, row.names = 1, check.names = FALSE)
  else read.csv(f, row.names = 1, check.names = FALSE)
}

# Preflight checks — fail fast before expensive computation
counts_file <- list.files(job_dir, pattern = "^counts", full.names = TRUE)[1]
coldata_file <- list.files(job_dir, pattern = "^coldata", full.names = TRUE)[1]

if (is.na(counts_file) || !file.exists(counts_file)) {
  stop("[deseq2] counts 文件未找到。期望 job 目录包含以 'counts' 开头的 CSV 文件。")
}
if (is.na(coldata_file) || !file.exists(coldata_file)) {
  stop("[deseq2] coldata 文件未找到。期望 job 目录包含以 'coldata' 开头的 CSV 文件。")
}

counts_raw <- as.matrix(read_matrix(counts_file))
coldata_raw <- read_matrix(coldata_file)

# Preflight 1: sample name alignment
coldata_samples <- rownames(coldata_raw)
counts_samples <- colnames(counts_raw)
aligned <- intersect(counts_samples, coldata_samples)
if (length(aligned) == 0) {
  stop(paste0(
    "[deseq2] 样本名对齐失败。counts 列名 (", length(counts_samples), "): ",
    paste(head(counts_samples, 5), collapse=", "),
    if(length(counts_samples) > 5) "..." else "",
    "。coldata 行名 (", length(coldata_samples), "): ",
    paste(head(coldata_samples, 5), collapse=", "),
    if(length(coldata_samples) > 5) "..." else "",
    ". 无交集，无法对齐。"
  ))
}
cat(sprintf("[deseq2] Preflight: %d/%d counts 列名与 coldata 行名对齐\n", length(aligned), ncol(counts_raw)))

# Preflight 2: design column must exist
if (!design_col %in% colnames(coldata_raw)) {
  available <- paste(colnames(coldata_raw), collapse=", ")
  stop(paste0(
    "[deseq2] design 列 '", design_col, "' 不存在于 coldata。",
    "可用的列名: ", available,
    "。请在 Job 参数中设置正确的 design 列名，或确保 coldata 包含 '", design_col, "' 列。"
  ))
}

# Preflight 3: contrast levels must exist in design column
design_levels <- unique(as.character(coldata_raw[[design_col]]))
if (!(ref %in% design_levels)) {
  stop(paste0(
    "[deseq2] 对照组 '", ref, "' 不在 design 列 '", design_col, "' 中。",
    "实际 levels: ", paste(design_levels, collapse=", ")
  ))
}
if (!(treat %in% design_levels)) {
  stop(paste0(
    "[deseq2] 实验组 '", treat, "' 不在 design 列 '", design_col, "' 中。",
    "实际 levels: ", paste(design_levels, collapse=", ")
  ))
}

# Preflight 4: minimum replicates
n_ref <- sum(coldata_raw[[design_col]] == ref)
n_treat <- sum(coldata_raw[[design_col]] == treat)
if (n_ref < 2 || n_treat < 2) {
  stop(paste0(
    "[deseq2] 每个条件至少需要 2 个生物学重复。当前: 对照组(", ref, ")=", n_ref,
    ", 实验组(", treat, ")=", n_treat, "。"
  ))
}

# Preflight 5: integer counts
if (any(round(counts_raw) < 0)) {
  stop("[deseq2] counts 矩阵包含负值，不是合法的 counts 数据。")
}

cat(sprintf("[deseq2] Preflight passed: %d samples (%s=%d, %s=%d), design='%s'\n",
  ncol(counts_raw), ref, n_ref, treat, n_treat, design_col))

# Apply sample alignment
counts <- counts_raw[, aligned, drop = FALSE]
coldata <- coldata_raw[aligned, , drop = FALSE]

# Ensure counts are integers
counts <- round(counts)
storage.mode(counts) <- "integer"

coldata[[design_col]] <- factor(coldata[[design_col]], levels = c(ref, treat))

# Filter low-expression genes
keep <- rowSums(counts >= min_count) >= ncol(counts) * 0.3
counts <- counts[keep, ]

cat(sprintf("[deseq2] Samples: %d, Genes after filter: %d\n", ncol(counts), nrow(counts)))

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                               design = as.formula(paste0("~ ", design_col)))
dds <- DESeq(dds)
res <- results(dds, contrast = c(design_col, treat, ref))
res_shrink <- lfcShrink(dds, coef = paste0(design_col, "_", treat, "_vs_", ref), type = "apeglm", quiet = TRUE)

# Results table
res_df <- as.data.frame(res_shrink)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(output_dir, "deg_table.csv"), row.names = FALSE)

# Significant genes
sig <- res_df[!is.na(res_df$padj) & res_df$padj < padj_cut & abs(res_df$log2FoldChange) > lfc_cut, ]
n_up <- sum(sig$log2FoldChange > 0)
n_down <- sum(sig$log2FoldChange < 0)

cat(sprintf("[deseq2] Significant: %d up, %d down\n", n_up, n_down))

# ============================================================
# Step 5a: Volcano plot (ggplot2 系)
# ============================================================
cat("[deseq2] Step 5a: volcano...\n")
res_df$sig <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cut & res_df$log2FoldChange > lfc_cut, "Up",
              ifelse(!is.na(res_df$padj) & res_df$padj < padj_cut & res_df$log2FoldChange < -lfc_cut, "Down", "NS"))
top <- head(res_df[res_df$sig != "NS", ], 15)
volcano_pal <- biof3_palette_div(3)  # 蓝-白-红 → 用 Down/NS/Up 三档

p_vol <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_text_repel(data = top, aes(label = gene), size = 3, max.overlaps = 15, color = "black", seed = 20260729) +
  scale_color_manual(values = c(Up = volcano_pal[3], Down = volcano_pal[1], NS = "gray70")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "gray40") +
  labs(title = sprintf("Volcano (%d up, %d down)", n_up, n_down),
       x = expression(log[2]~"Fold Change"), y = expression(-log[10]~italic(P)[adj])) +
  theme_biof3()
ggsave_biof3(p_vol, file.path(output_dir, "volcano"), width = 6.5, height = 5)

# ============================================================
# Step 5b: PCA (ggplot2 系)
# ============================================================
cat("[deseq2] Step 5b: PCA...\n")
# vst() 默认从最多 1000 个基因拟合趋势；小型验证集或靶向 panel 少于该数量时，
# 使用同一 DESeq2 官方变换的直接实现，避免把合法的小矩阵误判为执行异常。
vsd <- if (nrow(dds) < 1000) {
  varianceStabilizingTransformation(dds, blind = TRUE)
} else {
  vst(dds, blind = TRUE)
}
pca_data <- plotPCA(vsd, intgroup = design_col, returnData = TRUE)
pca_var <- round(100 * attr(pca_data, "percentVar"))
pca_pal <- biof3_palette(2)  # 2 个组用 NPG 前 2 色

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = pca_pal) +
  labs(title = "Principal Component Analysis",
       x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"),
       color = design_col) +
  theme_biof3()
ggsave_biof3(p_pca, file.path(output_dir, "pca"), width = 6, height = 4.5)

# ============================================================
# Step 5c: Library Size barplot (ggplot2 系)
# ============================================================
cat("[deseq2] Step 5c: library_size...\n")
lib_df <- data.frame(sample = colnames(counts), millions = colSums(counts) / 1e6, group = coldata[[design_col]])
lib_pal <- biof3_palette(2)

p_lib <- ggplot(lib_df, aes(sample, millions, fill = group)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = lib_pal) +
  labs(title = "Library Size", x = NULL, y = "Reads (millions)", fill = design_col) +
  theme_biof3() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
ggsave_biof3(p_lib, file.path(output_dir, "library_size"), width = 6, height = 4.5)

# ============================================================
# Step 5d: Dispersion plot (base R / grid 系)
# ============================================================
cat("[deseq2] Step 5d: dispersion...\n")
save_grid_biof3(file.path(output_dir, "dispersion"), width = 6.5, height = 5, expr = {
  par(mar = c(4.5, 4.5, 2, 1), bty = "l")
  plotDispEsts(dds, main = "Dispersion Estimates",
               cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)
})

# ============================================================
# Step 5e: MA plot (base R / grid 系)
# ============================================================
cat("[deseq2] Step 5e: MA plot...\n")
save_grid_biof3(file.path(output_dir, "ma_plot"), width = 6.5, height = 5, expr = {
  par(mar = c(4.5, 4.5, 2, 1), bty = "l")
  plotMA(res_shrink, main = "MA Plot (apeglm-shrunken)",
         ylim = c(-5, 5),
         cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)
})

# ============================================================
# Step 5f: Sample Distance heatmap (pheatmap / grid 系)
# ============================================================
cat("[deseq2] Step 5f: sample_distance...\n")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste0(colnames(counts), " (", coldata[[design_col]], ")")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)

dist_colors <- biof3_palette_seq(50, option = "ylorbr")  # 浅黄→深棕，单极
save_grid_biof3(file.path(output_dir, "sample_distance"), width = 7, height = 6, expr = {
  pheatmap(sampleDistMatrix,
           main = "Sample Distance (VST Euclidean)",
           color = dist_colors,
           border_color = NA,
           fontsize = 9,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists)
})

# ============================================================
# Step 5g: Heatmap (top 50 DEGs, pheatmap / grid 系)
# ============================================================
cat("[deseq2] Step 5g: heatmap top 50...\n")
top50 <- head(rownames(res_df[res_df$sig != "NS", ]), 50)
if (length(top50) > 2) {
  mat <- assay(vsd)[top50, ]
  mat_z <- t(scale(t(mat)))
  anno_col <- data.frame(Group = coldata[[design_col]], row.names = colnames(mat))
  group_pal <- biof3_palette(2)
  names(group_pal) <- levels(coldata[[design_col]])
  anno_colors <- list(Group = group_pal)

  heat_pal <- biof3_palette_div(100)  # 蓝-白-红，DEG z-score 标准

  save_grid_biof3(file.path(output_dir, "heatmap"), width = 7, height = 8.5, expr = {
    pheatmap(mat_z,
             annotation_col = anno_col,
             annotation_colors = anno_colors,
             show_colnames = FALSE,
             clustering_method = "ward.D2",
             fontsize = 9, fontsize_row = 7,
             color = heat_pal,
             border_color = NA,
             main = "Top 50 DEGs (VST z-score)")
  })
}

# ============================================================
# Step 6: Summary + 上下调列表
# ============================================================
summary_text <- sprintf(
  "DESeq2 Differential Expression Analysis\n\nSamples: %d (%s: %d, %s: %d)\nGenes tested: %d (after filtering)\nSignificant (padj < %g, |log2FC| > %g): %d\n  Up-regulated: %d\n  Down-regulated: %d\n",
  ncol(counts), treat, sum(coldata[[design_col]] == treat),
  ref, sum(coldata[[design_col]] == ref),
  nrow(counts), padj_cut, lfc_cut, nrow(sig), n_up, n_down
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

up_genes <- sig[sig$log2FoldChange > 0, ]
down_genes <- sig[sig$log2FoldChange < 0, ]
write.csv(up_genes, file.path(output_dir, "up_genes.csv"), row.names = FALSE)
write.csv(down_genes, file.path(output_dir, "down_genes.csv"), row.names = FALSE)

# ============================================================
# Step 7: 中间数据 (确保在线绘图可复现)
# ============================================================
write.csv(as.data.frame(assay(vsd)), file.path(output_dir, "vst_matrix.csv"))
write.csv(as.data.frame(coldata), file.path(output_dir, "coldata_info.csv"))
disp_df <- data.frame(gene = rownames(counts),
                      baseMean = mcols(dds)$baseMean,
                      dispGeneEst = mcols(dds)$dispGeneEst,
                      dispFit = mcols(dds)$dispFit,
                      dispersion = mcols(dds)$dispersion)
write.csv(disp_df, file.path(output_dir, "dispersion_data.csv"), row.names = FALSE)
cat("[deseq2] Intermediate data saved.\n")

# The formal HTML report embeds every generated PNG. R's PNG devices may emit
# different compressed byte streams for identical pixels, so normalize the
# complete PNG set before either plots or the report enter Artifact Truth.
for (png_path in list.files(output_dir, pattern = "\\.png$", full.names = TRUE)) {
  canonicalize_png(png_path)
}

# ============================================================
# Manifest (含 PNG + PDF 双格式)
# ============================================================
manifest <- list(
  files = list(
    list(name = "deg_table.csv",        type = "table", label = "差异基因表"),
    list(name = "up_genes.csv",         type = "table", label = "上调基因"),
    list(name = "down_genes.csv",       type = "table", label = "下调基因"),
    list(name = "volcano.png",          type = "plot",  label = "火山图 (PNG)"),
    list(name = "volcano.pdf",          type = "file",  label = "火山图 (PDF)"),
    list(name = "pca.png",              type = "plot",  label = "PCA 图 (PNG)"),
    list(name = "pca.pdf",              type = "file",  label = "PCA 图 (PDF)"),
    list(name = "heatmap.png",          type = "plot",  label = "Top50 热图 (PNG)"),
    list(name = "heatmap.pdf",          type = "file",  label = "Top50 热图 (PDF)"),
    list(name = "library_size.png",     type = "plot",  label = "文库大小 (PNG)"),
    list(name = "library_size.pdf",     type = "file",  label = "文库大小 (PDF)"),
    list(name = "dispersion.png",       type = "plot",  label = "离散度估计 (PNG)"),
    list(name = "dispersion.pdf",       type = "file",  label = "离散度估计 (PDF)"),
    list(name = "ma_plot.png",          type = "plot",  label = "MA 图 (PNG)"),
    list(name = "ma_plot.pdf",          type = "file",  label = "MA 图 (PDF)"),
    list(name = "sample_distance.png",  type = "plot",  label = "样本距离 (PNG)"),
    list(name = "sample_distance.pdf",  type = "file",  label = "样本距离 (PDF)"),
    list(name = "vst_matrix.csv",       type = "table", label = "VST 矩阵"),
    list(name = "coldata_info.csv",     type = "table", label = "样本信息"),
    list(name = "dispersion_data.csv",  type = "table", label = "离散度数据"),
    list(name = "summary.txt",          type = "text",  label = "分析摘要"),
    list(name = "report.html",          type = "file",  label = "解读报告 HTML")
  ),
  summary = list(
    total_genes = nrow(counts),
    sig_genes = nrow(sig),
    up = n_up,
    down = n_down,
    sci_style = TRUE,
    dual_format = TRUE
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ============================================================
# HTML 报告 (从模板 base64 内嵌)
# ============================================================
template_path <- file.path(job_dir, "report-template.html")
if (!file.exists(template_path)) {
  # 平台安全兜底：从宿主注入的内置插件资源目录解析（main.js 启动 r-server 时设置
  # BIOF3_BUILT_IN_PLUGINS_DIR=resources/built-in-plugins），不再硬编码 POSIX 路径
  # （/opt/... 在 Windows 与 macOS app bundle 上均不存在）。
  built_in_dir <- Sys.getenv("BIOF3_BUILT_IN_PLUGINS_DIR", "")
  if (nzchar(built_in_dir)) {
    candidate <- file.path(built_in_dir, "deseq2", "deseq2-report-template.html")
    if (file.exists(candidate)) template_path <- candidate
  }
}
if (file.exists(template_path)) {
  report_html <- readLines(template_path, warn = FALSE)
  report_html <- paste(report_html, collapse = "\n")
  # 填充摘要占位符（模板中的 {{...}} token）
  report_html <- gsub("{{TOTAL_GENES}}", as.character(nrow(counts)), report_html, fixed = TRUE)
  report_html <- gsub("{{SIG_GENES}}", as.character(nrow(sig)), report_html, fixed = TRUE)
  report_html <- gsub("{{UP}}", as.character(n_up), report_html, fixed = TRUE)
  report_html <- gsub("{{DOWN}}", as.character(n_down), report_html, fixed = TRUE)
  report_html <- gsub("{{DESIGN}}", design_col, report_html, fixed = TRUE)
  report_html <- gsub("{{REF}}", ref, report_html, fixed = TRUE)
  report_html <- gsub("{{TREAT}}", treat, report_html, fixed = TRUE)
  report_html <- gsub("{{PADJ}}", formatC(padj_cut, format = "g"), report_html, fixed = TRUE)
  report_html <- gsub("{{LFC}}", formatC(lfc_cut, format = "g"), report_html, fixed = TRUE)
  # 填充差异基因结果表（{{DEG_TABLE}} token）——sig 为已算好的显著差异基因框
  deg_rows <- ""
  if (nrow(sig) > 0) {
    sig_sorted <- sig[order(sig$padj, decreasing = FALSE), ]
    max_rows <- min(nrow(sig_sorted), 2000)
    for (i in seq_len(max_rows)) {
      g <- rownames(sig_sorted)[i]
      bm <- formatC(sig_sorted$baseMean[i], format = "g", digits = 4)
      lfc <- formatC(sig_sorted$log2FoldChange[i], format = "f", digits = 3)
      pj <- formatC(sig_sorted$padj[i], format = "g", digits = 3)
      if (sig_sorted$log2FoldChange[i] > 0) {
        reg <- "上调"; reg_cls <- "up"
      } else {
        reg <- "下调"; reg_cls <- "down"
      }
      deg_rows <- paste0(deg_rows,
        sprintf('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class="%s">%s</td></tr>\n',
                g, bm, lfc, pj, reg_cls, reg))
    }
    if (nrow(sig_sorted) > max_rows) {
      deg_rows <- paste0(deg_rows,
        sprintf('<tr><td colspan="5" class="more">… 仅显示前 %d 个（共 %d 个），完整列表见 deg_table.csv</td></tr>\n',
                max_rows, nrow(sig_sorted)))
    }
  } else {
    deg_rows <- '<tr><td colspan="5" class="more">未检测到显著差异基因</td></tr>\n'
  }
  report_html <- gsub("{{DEG_TABLE}}", deg_rows, report_html, fixed = TRUE)
  img_files <- c("library_size.png", "dispersion.png", "pca.png",
                 "sample_distance.png", "volcano.png", "ma_plot.png", "heatmap_top50.png")
  for (img in img_files) {
    img_path <- file.path(output_dir, sub("heatmap_top50.png", "heatmap.png", img))
    if (file.exists(img_path)) {
      b64 <- base64enc::base64encode(img_path)
      data_uri <- paste0("data:image/png;base64,", b64)
      report_html <- gsub(paste0("https://biof3.com/api/r/tools/demo/deseq2/", img), data_uri, report_html, fixed = TRUE)
      report_html <- gsub(paste0("https://biof3.com/api/platform/tools/demo/deseq2/", img), data_uri, report_html, fixed = TRUE)
    }
  }
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[deseq2] Report generated.\n")
} else {
  cat("[deseq2] Report template not found, skipping.\n")
}

cat("[deseq2] Pipeline complete.\n")
