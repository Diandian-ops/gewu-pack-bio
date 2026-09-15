#!/usr/bin/env Rscript
# ============================================================
# GSEA 基因集富集分析工具脚本
# Validated: clusterProfiler 4.18.4 / enrichplot 1.30.5 /
# ggplot2 4.0.3 (R 4.5.2, Bioconductor 3.22)
#
# Steps (in-scope):
#   1. read ranked gene list (gene + log2FC, 或 DESeq2 result)
#   2. SYMBOL → ENTREZID
#   3. gseGO BP / gseKEGG / Hallmark (msigdbr)
#   4. dotplot / waterfall / running score / ridge / emapplot / cnetplot / treeplot
#   5. 中间数据 (gsea_res.rds) + HTML 报告
#
# 2026-05-25: 升级到 BioF3 SCI 视觉风格规范
#   - dotplot / waterfall (ggplot 系): + theme_biof3() + ggsave_biof3()
#   - running_score / ridgeplot / cnetplot / emapplot / treeplot (enrichplot 复合): save_grid_biof3()
#   - waterfall 配色: biof3_palette_div(2 极) 红激活 / 蓝抑制
# ============================================================

suppressMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(jsonlite)
  library(DOSE)
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

find_biof3_kegg_snapshot <- function(organism) {
  candidates <- file.path(.libPaths(), ".biof3-data", paste0("kegg-", organism, ".rds"))
  matches <- candidates[file.exists(candidates)]
  if (length(matches) == 0) return(NA_character_)
  matches[1]
}

run_kegg_gsea <- function(gene_list, organism, min_size, max_size) {
  snapshot_path <- find_biof3_kegg_snapshot(organism)
  if (!is.na(snapshot_path)) {
    offline <- tryCatch({
      snapshot <- readRDS(snapshot_path)
      GSEA(
        geneList = gene_list,
        TERM2GENE = snapshot$term2gene,
        TERM2NAME = snapshot$term2name,
        minGSSize = min_size,
        maxGSSize = max_size,
        pvalueCutoff = 0.25,
        verbose = FALSE
      )
    }, error = function(e) {
      cat(sprintf("  Offline KEGG snapshot failed: %s\n", e$message))
      NULL
    })
    if (!is.null(offline)) {
      cat(sprintf("  KEGG source: offline snapshot %s\n", basename(snapshot_path)))
      return(offline)
    }
  }

  cat("  KEGG offline snapshot unavailable; falling back to KEGG REST\n")
  gseKEGG(
    geneList = gene_list,
    organism = organism,
    minGSSize = min_size,
    maxGSSize = max_size,
    pvalueCutoff = 0.25,
    verbose = FALSE,
    use_internal_data = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题
source(file.path(job_dir, "_biof3-theme.R"))

cat("[gsea] Starting pipeline (SCI style v2026-05-25)\n")

params <- fromJSON(file.path(job_dir, "params.json"))
organism <- params$organism %||% "hsa"
gene_set <- params$gene_set %||% "GO_BP"
min_size <- as.integer(params$min_size %||% 15)
max_size <- as.integer(params$max_size %||% 500)

cat(sprintf("  organism=%s, gene_set=%s, size %d-%d\n",
            organism, gene_set, min_size, max_size))

if (organism == "mmu") {
  suppressMessages(library(org.Mm.eg.db))
  orgdb <- org.Mm.eg.db
  kegg_org <- "mmu"
} else {
  suppressMessages(library(org.Hs.eg.db))
  orgdb <- org.Hs.eg.db
  kegg_org <- "hsa"
}

# ============================================================
# Read input
# ============================================================
input_file <- list.files(job_dir, pattern = "^ranked|^genes|^deg_table|^input",
                         full.names = TRUE)[1]
if (is.na(input_file)) {
  input_file <- list.files(job_dir, pattern = "\\.(csv|tsv|txt)$",
                           full.names = TRUE)[1]
}
if (is.na(input_file)) stop("No input gene file found")

cat(sprintf("[gsea] Reading: %s\n", basename(input_file)))

raw <- tryCatch(
  read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)
)

if (all(c("log2FoldChange", "padj") %in% colnames(raw)) ||
    all(c("logFC", "adj.P.Val") %in% colnames(raw))) {
  cat("[gsea] Detected DESeq2/limma table\n")
  if ("logFC" %in% colnames(raw)) {
    colnames(raw)[colnames(raw) == "logFC"] <- "log2FoldChange"
  }
  gene_col <- NULL
  for (cn in c("gene", "Gene", "SYMBOL", "gene_name", "gene_id")) {
    if (cn %in% colnames(raw)) { gene_col <- cn; break }
  }
  if (is.null(gene_col)) {
    if (is.character(raw[[1]])) gene_col <- colnames(raw)[1]
    else { raw$gene <- rownames(raw); gene_col <- "gene" }
  }
  gene_list <- raw$log2FoldChange
  names(gene_list) <- as.character(raw[[gene_col]])
} else {
  cat("[gsea] Detected ranked list\n")
  if (ncol(raw) >= 2) {
    gene_list <- as.numeric(raw[[2]])
    names(gene_list) <- as.character(raw[[1]])
  } else {
    stop("Input must have ≥2 columns: gene_symbol + log2FC/score")
  }
}

gene_list <- gene_list[!is.na(gene_list) & !is.na(names(gene_list)) & nchar(names(gene_list)) > 0]
gene_list <- sort(gene_list, decreasing = TRUE)
cat(sprintf("  Ranked: %d genes\n", length(gene_list)))

if (length(gene_list) < 100) stop("Too few genes (need >= 100)")

# ============================================================
# Gene ID conversion
# ============================================================
if (all(grepl("^[0-9]+$", head(names(gene_list), 20)))) {
  cat("  IDs are ENTREZID\n")
  entrez_list <- gene_list
} else {
  id_map <- tryCatch(
    bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb),
    error = function(e) {
      stop(sprintf(
        "Gene ID mapping failed: none of the input values are valid %s SYMBOLs",
        if (organism == "mmu") "mouse" else "human"
      ), call. = FALSE)
    }
  )
  mapped <- gene_list[names(gene_list) %in% id_map$SYMBOL]
  names(mapped) <- id_map$ENTREZID[match(names(mapped), id_map$SYMBOL)]
  mapped <- mapped[!duplicated(names(mapped))]
  entrez_list <- sort(mapped, decreasing = TRUE)
  cat(sprintf("  Converted: %d → %d ENTREZID\n", length(gene_list), length(entrez_list)))
  if (length(entrez_list) < 100) {
    stop(sprintf(
      "Gene ID mapping failed: only %d of %d input SYMBOLs mapped to ENTREZID (need >= 100)",
      length(entrez_list), length(gene_list)
    ), call. = FALSE)
  }
  write.csv(id_map, file.path(output_dir, "gene_id_mapping.csv"), row.names = FALSE)
}

# ============================================================
# Run GSEA
# ============================================================
cat(sprintf("[gsea] Running GSEA (%s)...\n", gene_set))
gsea_res <- NULL

if (gene_set == "GO_BP") {
  gsea_res <- tryCatch(
    gseGO(geneList = entrez_list, OrgDb = orgdb, ont = "BP",
           minGSSize = min_size, maxGSSize = max_size,
           pvalueCutoff = 0.25, verbose = FALSE),
    error = function(e) { cat(sprintf("  gseGO failed: %s\n", e$message)); NULL })
} else if (gene_set == "KEGG") {
  gsea_res <- tryCatch(
    run_kegg_gsea(entrez_list, kegg_org, min_size, max_size),
    error = function(e) { cat(sprintf("  gseKEGG failed: %s\n", e$message)); NULL })
} else if (gene_set == "Hallmark") {
  tryCatch({
    suppressMessages(library(msigdbr))
    species_name <- if (organism == "mmu") "Mus musculus" else "Homo sapiens"
    hallmark <- msigdbr(species = species_name, category = "H")
    t2g <- data.frame(term = hallmark$gs_name, gene = as.character(hallmark$entrez_gene))
    gsea_res <- GSEA(geneList = entrez_list, TERM2GENE = t2g,
                     minGSSize = min_size, maxGSSize = max_size,
                     pvalueCutoff = 0.25, verbose = FALSE)
  }, error = function(e) cat(sprintf("  Hallmark GSEA failed: %s\n", e$message)))
}

if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
  writeLines("No significantly enriched gene sets found.", file.path(output_dir, "summary.txt"))
  manifest <- list(files = list(list(name = "summary.txt", type = "text", label = "分析摘要")),
                   summary = list(n_enriched = 0, gene_set = gene_set))
  writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), file.path(output_dir, "manifest.json"))
  cat("[gsea] No enriched gene sets. Exit OK.\n")
  quit(save = "no", status = 0)
}

n_enriched <- nrow(as.data.frame(gsea_res))
cat(sprintf("  Enriched: %d gene sets\n", n_enriched))

if (gene_set == "GO_BP") {
  gsea_res <- setReadable(gsea_res, OrgDb = orgdb, keyType = "ENTREZID")
}

# ============================================================
# Plots — SCI 风格双格式
# ============================================================
cat("[gsea] Generating plots...\n")

# ---- 1. Dotplot (ggplot 系) ----
tryCatch({
  p <- dotplot(gsea_res, showCategory = min(20, n_enriched)) +
    labs(title = paste("GSEA", gene_set, "Dotplot")) +
    theme_biof3(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
  ggsave_biof3(p, file.path(output_dir, "gsea_dotplot"), width = 8, height = 7)
}, error = function(e) cat(sprintf("  dotplot skipped: %s\n", e$message)))

# ---- 2. NES Waterfall (ggplot 系，自建) ----
tryCatch({
  gsea_df <- as.data.frame(gsea_res)
  top_up <- head(gsea_df[gsea_df$NES > 0, ], 10)
  top_down <- tail(gsea_df[gsea_df$NES < 0, ], 10)
  plot_df <- rbind(top_up, top_down)
  plot_df$Direction <- ifelse(plot_df$NES > 0, "Activated", "Suppressed")
  plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description[order(plot_df$NES)])

  pal <- biof3_palette_div(3)  # 蓝-白-红 → 取 1（蓝）和 3（红）
  p <- ggplot(plot_df, aes(x = NES, y = Description, fill = Direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(Activated = pal[3], Suppressed = pal[1])) +
    geom_vline(xintercept = 0, color = "gray30", linewidth = 0.4) +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL, title = "GSEA NES Waterfall") +
    theme_biof3() +
    theme(legend.position = "top",
          axis.text.y = element_text(size = 9))
  ggsave_biof3(p, file.path(output_dir, "gsea_waterfall"), width = 9, height = 7)
}, error = function(e) cat(sprintf("  waterfall skipped: %s\n", e$message)))

# ---- 3. Running Score (top 3，gseaplot2 复合对象，用 save_grid_biof3) ----
tryCatch({
  top_ids <- head(gsea_res@result$ID, 3)
  for (i in seq_along(top_ids)) {
    pp <- gseaplot2(gsea_res, geneSetID = top_ids[i],
                     title = gsea_res@result$Description[i])
    save_grid_biof3(file.path(output_dir, paste0("running_score_", i)),
                     width = 8, height = 6, expr = print(pp))
  }
}, error = function(e) cat(sprintf("  running score skipped: %s\n", e$message)))

# ---- 4. Ridge Plot (ggplot 系) ----
tryCatch({
  p <- ridgeplot(gsea_res, showCategory = min(15, n_enriched)) +
    labs(title = paste("GSEA", gene_set, "Ridge Plot")) +
    theme_biof3(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
  ggsave_biof3(p, file.path(output_dir, "gsea_ridge"), width = 9, height = 8)
}, error = function(e) cat(sprintf("  ridge skipped: %s\n", e$message)))

# ---- 5. Emapplot (enrichplot ggraph 复合, 用 save_grid_biof3) ----
gsea_sim <- NULL
if (n_enriched >= 5) {
  tryCatch({
    gsea_sim <- pairwise_termsim(gsea_res)
    p <- emapplot(gsea_sim, showCategory = min(30, n_enriched)) +
      ggtitle(paste("GSEA", gene_set, "Enrichment Map"))
    save_grid_biof3(file.path(output_dir, "gsea_emapplot"), width = 10, height = 9, expr = print(p))
  }, error = function(e) cat(sprintf("  emapplot skipped: %s\n", e$message)))
}

# ---- 6. Cnetplot ----
tryCatch({
  p <- cnetplot(gsea_res, showCategory = min(5, n_enriched), foldChange = entrez_list) +
    ggtitle(paste("GSEA", gene_set, "Gene-Concept Network"))
  save_grid_biof3(file.path(output_dir, "gsea_cnetplot"), width = 10, height = 8, expr = print(p))
}, error = function(e) cat(sprintf("  cnetplot skipped: %s\n", e$message)))

# ---- 7. Treeplot ----
if (n_enriched >= 10) {
  tryCatch({
    if (is.null(gsea_sim)) gsea_sim <- pairwise_termsim(gsea_res)
    p <- treeplot(gsea_sim, showCategory = min(30, n_enriched)) +
      ggtitle(paste("GSEA", gene_set, "Tree"))
    save_grid_biof3(file.path(output_dir, "gsea_treeplot"), width = 12, height = 8, expr = print(p))
  }, error = function(e) cat(sprintf("  treeplot skipped: %s\n", e$message)))
}

# ============================================================
# 中间数据
# ============================================================
cat("[gsea] Saving intermediate data...\n")
write.csv(as.data.frame(gsea_res), file.path(output_dir, "gsea_results.csv"), row.names = FALSE)
saveRDS(gsea_res, file.path(output_dir, "gsea_res.rds"))
write.csv(data.frame(gene = names(gene_list), log2FC = gene_list),
          file.path(output_dir, "ranked_gene_list.csv"), row.names = FALSE)

# ============================================================
# Summary
# ============================================================
n_activated <- sum(as.data.frame(gsea_res)$NES > 0)
n_suppressed <- sum(as.data.frame(gsea_res)$NES < 0)

summary_text <- sprintf(
  "GSEA Analysis Summary\n\nInput genes: %d\nOrganism: %s\nGene set: %s\nSize range: %d - %d\n\nResults:\n  Enriched gene sets: %d\n  Activated (NES > 0): %d\n  Suppressed (NES < 0): %d\n",
  length(gene_list), organism, gene_set, min_size, max_size,
  n_enriched, n_activated, n_suppressed
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# Manifest (双格式 + RDS)
# ============================================================
files_list <- list(
  list(name = "summary.txt",       type = "text",  label = "分析摘要"),
  list(name = "gsea_results.csv",  type = "table", label = "GSEA 结果表")
)
if (file.exists(file.path(output_dir, "gene_id_mapping.csv")))
  files_list <- c(files_list, list(list(name = "gene_id_mapping.csv", type = "table", label = "基因 ID 转换")))

# Plots: PNG + PDF 双格式
plot_pairs <- list(
  c("gsea_dotplot",     "GSEA 气泡图"),
  c("gsea_waterfall",   "NES Waterfall"),
  c("running_score_1",  "Running Score #1"),
  c("running_score_2",  "Running Score #2"),
  c("running_score_3",  "Running Score #3"),
  c("gsea_ridge",       "Ridge Plot"),
  c("gsea_emapplot",    "富集网络图"),
  c("gsea_cnetplot",    "基因网络图"),
  c("gsea_treeplot",    "富集树图")
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

# Intermediate data
for (extra in list(
  c("gsea_res.rds",        "file",  "GSEA RDS"),
  c("ranked_gene_list.csv", "table", "排序基因列表")
)) {
  if (file.exists(file.path(output_dir, extra[1]))) {
    files_list <- c(files_list, list(list(name = extra[1], type = extra[2], label = extra[3])))
  }
}

files_list <- c(files_list, list(list(name = "report.html", type = "file", label = "解读报告(HTML)")))

manifest <- list(
  files = files_list,
  summary = list(
    input_genes = length(gene_list),
    n_enriched = n_enriched,
    n_activated = n_activated,
    n_suppressed = n_suppressed,
    gene_set = gene_set,
    organism = organism,
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
  cat("[gsea] Generating report...\n")
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
  for (img_path in png_files) {
    img_name <- basename(img_path)
    b64 <- base64enc::base64encode(img_path)
    data_uri <- paste0("data:image/png;base64,", b64)
    report_html <- gsub(paste0("{{", img_name, "}}"), data_uri, report_html, fixed = TRUE)
  }

  report_html <- gsub("{{input_genes}}", as.character(length(gene_list)), report_html, fixed = TRUE)
  report_html <- gsub("{{n_enriched}}", as.character(n_enriched), report_html, fixed = TRUE)
  report_html <- gsub("{{n_activated}}", as.character(n_activated), report_html, fixed = TRUE)
  report_html <- gsub("{{n_suppressed}}", as.character(n_suppressed), report_html, fixed = TRUE)
  report_html <- gsub("{{gene_set}}", gene_set, report_html, fixed = TRUE)
  report_html <- gsub("{{organism}}", organism, report_html, fixed = TRUE)

  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[gsea] Report generated.\n")
} else {
  cat("[gsea] Report template not found, skipping.\n")
}

cat(sprintf("[gsea] Pipeline complete. Enriched=%d (up=%d, down=%d)\n",
            n_enriched, n_activated, n_suppressed))
