#!/usr/bin/env Rscript
# ============================================================
# 免疫浸润评估工具（ssGSEA / Bindea 28 immune cells）
#
# 输入（job_dir/）：
#   - tpm_matrix 文件 (CSV/TSV: 行=基因Symbol, 列=样本, 值=TPM)
#   - group_file（可选，两列：sample, group）
#   - params.json: {method}
#
# 输出（job_dir/output/）：
#   - 免疫细胞比例矩阵 + 堆叠条形图 + 相关性热图 + 组间箱线图 + report.html
# ============================================================

suppressMessages({
  library(GSVA)
  library(ggplot2)
  library(jsonlite)
  library(reshape2)
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题 (W3-2 升级 2026-05-25)
source(file.path(job_dir, "_biof3-theme.R"))

params <- fromJSON(file.path(job_dir, "params.json"))
method <- params$method %||% "ssgsea"

cat("Parameters loaded.\n")
cat("  Method:", method, "\n")

# ============================================================
# Bindea immune cell gene sets (curated)
# ============================================================
immune_gene_sets <- list(
  "aDCs" = c("CD1A","CD1B","CD1E","CCL13","CCL17","CCL22","CXCR4","NRP1","FLT3","ITGAX"),
  "B cells" = c("CD19","MS4A1","CD79A","CD79B","BLK","FCRL2","TNFRSF17","TCL1A","SPIB","PNOC"),
  "CD8 T cells" = c("CD8A","CD8B","GZMK","GZMA","GZMB","PRF1","CXCR3","IFNG","TBX21","EOMES"),
  "DCs" = c("CCL13","CCL17","CCL22","CD209","HSD11B1","CCR7","LAMP3","FCER1A","CD1C","CLEC10A"),
  "iDCs" = c("CCR6","CD1A","CD1B","CD1E","CXCL13","FLT3","NRP1","CXCR4","MRC1","CD209"),
  "Macrophages" = c("CD68","CD163","CSF1R","MSR1","MRC1","MARCO","VSIG4","SIGLEC1","CD14","FCGR1A"),
  "Mast cells" = c("CPA3","TPSAB1","TPSB2","MS4A2","HDC","CTSG","KIT","CMA1","GATA2","IL1RL1"),
  "Neutrophils" = c("CEACAM8","FCGR3B","CSF3R","CXCR1","CXCR2","FPR1","FPR2","SIGLEC5","S100A12","MNDA"),
  "NK cells" = c("NCR1","NCR3","KLRD1","KLRF1","KLRC2","KIR2DL4","CD160","NKG7","GNLY","GZMB"),
  "pDCs" = c("CLEC4C","IL3RA","NRP1","JCHAIN","LILRA4","GZMB","SERPINF1","ITM2C","IRF7","TCF4"),
  "T cells" = c("CD3D","CD3E","CD3G","CD2","CD28","LCK","TRAT1","CD5","ITK","TRBC1"),
  "T helper cells" = c("CD4","IL2RA","FOXP3","CTLA4","ICOS","CD40LG","IL7R","CCR4","GATA3","RORC"),
  "Tcm" = c("CCR7","SELL","IL7R","CD27","CD28","LEF1","TCF7","CD44","LMNA","STAT5A"),
  "Tem" = c("GZMK","GZMA","CCL5","NKG7","CST7","KLRG1","CX3CR1","FGFBP2","FCGR3A","GNLY"),
  "Tfh" = c("CXCR5","ICOS","BCL6","SH2D1A","CD200","PDCD1","IL21","MAF","ASCL2","BTLA"),
  "Th1 cells" = c("TBX21","IFNG","IL12RB2","STAT4","CXCR3","CCR5","IL2","TNF","LTA","FASLG"),
  "Th2 cells" = c("GATA3","IL4","IL5","IL13","IL10","CCR3","CCR4","PTGDR2","IL1RL1","STAT6"),
  "Th17 cells" = c("RORC","IL17A","IL17F","IL22","IL23R","CCR6","IL26","AHR","BATF","STAT3"),
  "Treg" = c("FOXP3","IL2RA","CTLA4","TNFRSF18","IKZF2","LRRC32","TIGIT","TNFRSF4","ENTPD1","IL10"),
  "Cytotoxic cells" = c("GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","GNLY","NKG7","KLRK1","KLRD1"),
  "Eosinophils" = c("CCR3","SIGLEC8","IL5RA","ALOX15","EPX","PRG2","CLC","RNASE2","RNASE3","IL4"),
  "Monocytes" = c("CD14","FCGR1A","CSF1R","VCAN","S100A8","S100A9","FCN1","LYZ","MNDA","CD68")
)

# ============================================================
# Read input files
# ============================================================
tpm_file <- list.files(job_dir, pattern = "^tpm|^expression|^expr|^fpkm",
                       full.names = TRUE, ignore.case = TRUE)[1]
if (is.na(tpm_file)) {
  all_files <- list.files(job_dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
  all_files <- all_files[!grepl("params|group|clinical|coldata", basename(all_files), ignore.case = TRUE)]
  tpm_file <- all_files[1]
}

group_file <- list.files(job_dir, pattern = "group|clinical|coldata|pheno",
                         full.names = TRUE, ignore.case = TRUE)[1]

if (is.na(tpm_file)) stop("No TPM expression file found")

cat("Reading TPM matrix:", basename(tpm_file), "\n")
tpm <- tryCatch(
  read.csv(tpm_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) read.delim(tpm_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
)
tpm <- as.matrix(tpm)
cat(sprintf("  TPM matrix: %d genes x %d samples\n", nrow(tpm), ncol(tpm)))

# Read group info if available
group_info <- NULL
if (!is.na(group_file)) {
  cat("Reading group file:", basename(group_file), "\n")
  group_info <- tryCatch(
    read.csv(group_file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) read.delim(group_file, stringsAsFactors = FALSE, check.names = FALSE)
  )
  if (ncol(group_info) >= 2) {
    colnames(group_info)[1:2] <- c("sample", "group")
    cat(sprintf("  Groups: %s\n", paste(unique(group_info$group), collapse = ", ")))
  } else {
    group_info <- NULL
  }
}

# ============================================================
# Run ssGSEA
# ============================================================
cat("Running ssGSEA (GSVA)...\n")

filtered_sets <- lapply(immune_gene_sets, function(genes) intersect(genes, rownames(tpm)))
set_sizes <- sapply(filtered_sets, length)
cat(sprintf("  Gene set coverage: %d/%d sets have >= 3 genes\n", sum(set_sizes >= 3), length(filtered_sets)))
filtered_sets <- filtered_sets[set_sizes >= 3]

if (length(filtered_sets) < 5) {
  stop("Too few immune gene sets have matching genes. Please ensure input uses gene symbols (e.g., CD8A, FOXP3).")
}

gsva_res <- tryCatch({
  # GSVA >= 2.0 new API
  param <- ssgseaParam(tpm, filtered_sets)
  gsva(param, verbose = FALSE)
}, error = function(e) {
  # Fallback to old API (GSVA < 2.0)
  gsva(tpm, filtered_sets, method = "ssgsea", verbose = FALSE)
})
cat(sprintf("  ssGSEA complete: %d cell types x %d samples\n", nrow(gsva_res), ncol(gsva_res)))

# Normalize to 0-1
gsva_norm <- t(apply(gsva_res, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-10)))

write.csv(as.data.frame(gsva_res), file.path(output_dir, "immune_scores_raw.csv"))
write.csv(as.data.frame(gsva_norm), file.path(output_dir, "immune_scores_normalized.csv"))
write.csv(as.data.frame(gsva_res), file.path(output_dir, "immune_matrix.csv"))

# ============================================================
# Plot 1: Stacked barplot
# ============================================================
cat("Generating plots...\n")

tryCatch({
  bar_df <- melt(gsva_norm)
  colnames(bar_df) <- c("CellType", "Sample", "Score")
  bar_df$Sample <- factor(bar_df$Sample, levels = colnames(gsva_norm))

  n_cells_bar <- length(unique(bar_df$CellType))
  bar_pal <- biof3_palette(n_cells_bar)  # 自动 ramp 到 N 个免疫细胞类型

  p <- ggplot(bar_df, aes(x = Sample, y = Score, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill", width = 0.9) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = bar_pal) +
    labs(x = NULL, y = "Relative Proportion",
         title = "Immune Cell Composition (ssGSEA)", fill = "Cell Type") +
    theme_biof3() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
          legend.text = element_text(size = 9),
          legend.key.size = unit(0.4, "cm"))

  ggsave_biof3(p, file.path(output_dir, "barplot_composition"),
                width = max(8, ncol(gsva_norm) * 0.15 + 4), height = 7)
}, error = function(e) cat("  Barplot failed:", e$message, "\n"))

# ============================================================
# Plot 2: Correlation heatmap
# ============================================================
tryCatch({
  cor_mat <- cor(t(gsva_res), method = "spearman")
  cor_df <- melt(cor_mat)
  colnames(cor_df) <- c("CellType1", "CellType2", "Correlation")

  div_pal <- biof3_palette_div(3)  # 蓝-白-红
  p <- ggplot(cor_df, aes(x = CellType1, y = CellType2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = div_pal[1], mid = "white", high = div_pal[3],
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Immune Cell Type Correlation (Spearman)", x = NULL, y = NULL) +
    theme_biof3() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9))

  ggsave_biof3(p, file.path(output_dir, "correlation_heatmap"), width = 10, height = 9)
}, error = function(e) cat("  Correlation heatmap failed:", e$message, "\n"))

# ============================================================
# Plot 3: Group comparison boxplot
# ============================================================
if (!is.null(group_info)) {
  tryCatch({
    common <- intersect(colnames(gsva_res), group_info$sample)
    if (length(common) >= 10) {
      box_df <- melt(gsva_res[, common])
      colnames(box_df) <- c("CellType", "Sample", "Score")
      box_df$Group <- group_info$group[match(box_df$Sample, group_info$sample)]
      box_df <- box_df[!is.na(box_df$Group), ]

      cell_var <- apply(gsva_res[, common], 1, var)
      top_cells <- names(sort(cell_var, decreasing = TRUE))[1:min(12, nrow(gsva_res))]
      box_df <- box_df[box_df$CellType %in% top_cells, ]

      n_groups <- length(unique(box_df$Group))
      group_pal <- biof3_palette(max(2, n_groups))

      p <- ggplot(box_df, aes(x = CellType, y = Score, fill = Group)) +
        geom_boxplot(outlier.size = 0.5, width = 0.7) +
        scale_fill_manual(values = group_pal) +
        labs(x = NULL, y = "ssGSEA Score",
             title = "Immune Infiltration: Group Comparison") +
        theme_biof3() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
              legend.position = "top")

      ggsave_biof3(p, file.path(output_dir, "boxplot_group"), width = 10, height = 6)
    }
  }, error = function(e) cat("  Group boxplot failed:", e$message, "\n"))
}

# ============================================================
# Plot 4: Lollipop (mean infiltration level)
# ============================================================
tryCatch({
  mean_scores <- sort(rowMeans(gsva_res), decreasing = TRUE)
  lollipop_df <- data.frame(
    CellType = factor(names(mean_scores), levels = rev(names(mean_scores))),
    MeanScore = mean_scores
  )

  lolli_color <- biof3_palette()[3]  # NPG 第 3 色 (NPG 绿)
  p <- ggplot(lollipop_df, aes(x = MeanScore, y = CellType)) +
    geom_segment(aes(x = 0, xend = MeanScore, y = CellType, yend = CellType),
                 color = "gray90", linewidth = 0.8) +
    geom_point(size = 3, color = lolli_color) +
    labs(x = "Mean ssGSEA Score", y = NULL, title = "Immune Cell Infiltration Level") +
    theme_biof3()

  ggsave_biof3(p, file.path(output_dir, "lollipop_scores"),
                width = 8, height = max(5, nrow(gsva_res) * 0.3))
}, error = function(e) cat("  Lollipop failed:", e$message, "\n"))

# ============================================================
# Summary + Manifest + Report
# ============================================================
n_cells <- nrow(gsva_res)
n_samples <- ncol(gsva_res)
top3 <- names(sort(rowMeans(gsva_res), decreasing = TRUE))[1:min(3, n_cells)]

summary_text <- sprintf(
  "Immune Infiltration Analysis (ssGSEA)\n\nMethod: %s\nSamples: %d\nImmune cell types: %d\n\nTop 3 infiltrating:\n  1. %s\n  2. %s\n  3. %s\n\nGene sets: Bindea et al. (2013)\n",
  method, n_samples, n_cells, top3[1], top3[2], top3[3])
writeLines(summary_text, file.path(output_dir, "summary.txt"))

files_list <- list(
  list(name = "summary.txt", type = "text", label = "\u5206\u6790\u6458\u8981"),
  list(name = "immune_scores_normalized.csv", type = "table", label = "\u514d\u75ab\u7ec6\u80de\u6bd4\u4f8b"),
  list(name = "immune_scores_raw.csv", type = "table", label = "ssGSEA \u539f\u59cb\u5206\u6570")
)
plot_pairs <- list(
  c("barplot_composition",  "堆叠条形图"),
  c("correlation_heatmap",  "相关性热图"),
  c("boxplot_group",        "组间箱线图"),
  c("lollipop_scores",      "浸润水平图")
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
files_list <- c(files_list, list(list(name = "report.html", type = "file", label = "\u89e3\u8bfb\u62a5\u544a")))

manifest <- list(files = files_list, summary = list(
  method = method, n_samples = n_samples, n_cell_types = n_cells,
  top_cells = paste(top3, collapse = ", "), has_group = !is.null(group_info)))
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), file.path(output_dir, "manifest.json"))

# Report
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
  report_html <- gsub("{{method}}", method, report_html, fixed = TRUE)
  report_html <- gsub("{{n_samples}}", n_samples, report_html, fixed = TRUE)
  report_html <- gsub("{{n_cell_types}}", n_cells, report_html, fixed = TRUE)
  report_html <- gsub("{{top_cells}}", paste(top3, collapse = ", "), report_html, fixed = TRUE)
  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("Report generated.\n")
}

cat("\n========================================\n")
cat("Immune infiltration complete!\n")
cat(sprintf("  %d cell types x %d samples | Top: %s\n", n_cells, n_samples, paste(top3, collapse = ", ")))
cat("========================================\n")
