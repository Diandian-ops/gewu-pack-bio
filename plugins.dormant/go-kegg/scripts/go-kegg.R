#!/usr/bin/env Rscript
# ============================================================
# GO/KEGG 富集分析工具脚本
# Validated: clusterProfiler 4.18.4 / enrichplot 1.30.5 /
# ggplot2 4.0.3 (R 4.5.2, Bioconductor 3.22)
# enrichplot 1.28.4 曾与 ggplot2 4.0+ 不兼容；当前组合已完成 dotplot + ggsave 实测。
#
# Steps (in-scope):
#   1. read input (gene list 或 DESeq2 result)
#   2. gene ID → ENTREZID
#   3. enrichGO (BP / MF / CC) + enrichKEGG
#   4. compareCluster (Up vs Down，如有 DESeq2 输入)
#   5. dotplot / barplot / cnetplot / treeplot / heatplot / emapplot
#   6. 中间数据保存 (ego_*.rds + *_data.csv) + HTML 报告
#
# 2026-05-25: 升级到 BioF3 SCI 视觉风格规范
#   - dotplot / barplot / compareCluster: + theme_biof3() + ggsave_biof3() 双格式
#   - cnetplot / treeplot / heatplot / emapplot: save_grid_biof3() 包装 (这些是 ggraph
#     复合对象，主题套上易破布局，保留 enrichplot 默认 + 双格式输出)
#   - 配色: 保留 enrichplot 默认 viridis/red gradient (已合规且与 p 值/q 值映射对齐)
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

canonicalize_png <- function(file_path) {
  if (!file.exists(file_path)) return(invisible(FALSE))
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("[go-kegg] R package 'png' is required for deterministic Result Studio artifacts")
  }
  png::writePNG(png::readPNG(file_path), file_path)
  invisible(TRUE)
}

find_biof3_kegg_snapshot <- function(organism) {
  candidates <- file.path(.libPaths(), ".biof3-data", paste0("kegg-", organism, ".rds"))
  matches <- candidates[file.exists(candidates)]
  if (length(matches) == 0) return(NA_character_)
  matches[1]
}

run_kegg_ora <- function(gene_ids, organism, pvalue_cutoff, qvalue_cutoff) {
  snapshot_path <- find_biof3_kegg_snapshot(organism)
  if (!is.na(snapshot_path)) {
    offline <- tryCatch({
      snapshot <- readRDS(snapshot_path)
      enricher(
        gene = gene_ids,
        TERM2GENE = snapshot$term2gene,
        TERM2NAME = snapshot$term2name,
        pAdjustMethod = "BH",
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff
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
  enrichKEGG(
    gene = gene_ids,
    organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    use_internal_data = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题
source(file.path(job_dir, "_biof3-theme.R"))

cat("[go-kegg] Starting pipeline (SCI style v2026-05-25)\n")

# Read params
params <- fromJSON(file.path(job_dir, "params.json"))
species <- params$species %||% "human"
analysis_type <- params$analysis_type %||% "all"
gene_id_type <- params$gene_id_type %||% "SYMBOL"
padj_cut <- as.numeric(params$padj_cutoff %||% 0.05)
lfc_cut <- as.numeric(params$lfc_cutoff %||% 1)
enrich_pval <- as.numeric(params$enrich_pvalue %||% 0.05)
enrich_qval <- as.numeric(params$enrich_qvalue %||% 0.2)
top_n <- as.integer(params$top_n %||% 20)

cat(sprintf("  species=%s, type=%s, id=%s, top_n=%d\n",
            species, analysis_type, gene_id_type, top_n))

# Load species-specific OrgDb
if (species == "mouse") {
  suppressMessages(library(org.Mm.eg.db))
  orgdb <- org.Mm.eg.db
  kegg_org <- "mmu"
} else {
  suppressMessages(library(org.Hs.eg.db))
  orgdb <- org.Hs.eg.db
  kegg_org <- "hsa"
}

# ============================================================
# Read input: support gene list OR DESeq2 result table
# ============================================================
input_file <- list.files(job_dir, pattern = "^genes|^deg_table|^input",
                         full.names = TRUE)[1]
if (is.na(input_file)) {
  input_file <- list.files(job_dir, pattern = "\\.(csv|tsv|txt)$",
                           full.names = TRUE)[1]
}
if (is.na(input_file)) stop("No input gene file found")

cat(sprintf("[go-kegg] Reading input: %s\n", basename(input_file)))

raw <- tryCatch(
  read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)
)

is_deseq_table <- all(c("log2FoldChange", "padj") %in% colnames(raw)) ||
                  all(c("logFC", "adj.P.Val") %in% colnames(raw))

up_genes <- NULL
down_genes <- NULL

if (is_deseq_table) {
  cat("[go-kegg] Detected DESeq2/limma result table\n")
  if ("logFC" %in% colnames(raw)) {
    colnames(raw)[colnames(raw) == "logFC"] <- "log2FoldChange"
    colnames(raw)[colnames(raw) == "adj.P.Val"] <- "padj"
  }
  gene_col <- NULL
  for (cn in c("gene", "Gene", "SYMBOL", "gene_name", "gene_id")) {
    if (cn %in% colnames(raw)) { gene_col <- cn; break }
  }
  if (is.null(gene_col)) {
    if (is.character(raw[[1]]) || is.factor(raw[[1]])) gene_col <- colnames(raw)[1]
    else { raw$gene <- rownames(raw); gene_col <- "gene" }
  }
  raw$gene <- as.character(raw[[gene_col]])
  raw <- raw[!is.na(raw$padj), ]

  sig <- raw[raw$padj < padj_cut & abs(raw$log2FoldChange) > lfc_cut, ]
  up_genes <- sig$gene[sig$log2FoldChange > 0]
  down_genes <- sig$gene[sig$log2FoldChange < 0]
  gene_list <- sig$gene
  cat(sprintf("  Significant: %d (up=%d, down=%d)\n",
              length(gene_list), length(up_genes), length(down_genes)))
} else {
  cat("[go-kegg] Detected simple gene list\n")
  gene_list <- if (ncol(raw) == 1) as.character(raw[[1]]) else as.character(raw[[1]])
  gene_list <- unique(gene_list[nchar(gene_list) > 0 & !is.na(gene_list)])
  cat(sprintf("  Input genes: %d\n", length(gene_list)))
}

if (length(gene_list) < 3) stop("Too few genes (need >= 3)")

# ============================================================
# Gene ID conversion
# ============================================================
detect_id_type <- function(genes, orgdb) {
  sample_genes <- head(genes, 20)
  if (any(grepl("^ENS[A-Z]*G", sample_genes))) return("ENSEMBL")
  if (all(grepl("^[0-9]+$", sample_genes))) return("ENTREZID")
  return("SYMBOL")
}

if (gene_id_type == "auto") {
  gene_id_type <- detect_id_type(gene_list, orgdb)
  cat(sprintf("  Auto-detected ID: %s\n", gene_id_type))
}

if (gene_id_type != "ENTREZID") {
  id_map <- bitr(gene_list, fromType = gene_id_type, toType = "ENTREZID", OrgDb = orgdb)
  entrez_ids <- unique(id_map$ENTREZID)
  cat(sprintf("  Converted: %d → %d ENTREZID\n", length(gene_list), length(entrez_ids)))
  write.csv(id_map, file.path(output_dir, "gene_annotation.csv"), row.names = FALSE)
} else {
  entrez_ids <- unique(gene_list)
  id_map <- data.frame(ENTREZID = entrez_ids, stringsAsFactors = FALSE)
}

up_entrez <- NULL
down_entrez <- NULL
up_readable <- NULL
down_readable <- NULL
if (!is.null(up_genes)) {
  up_entrez <- character(0)
  up_readable <- character(0)
}
if (!is.null(down_genes)) {
  down_entrez <- character(0)
  down_readable <- character(0)
}
if (!is.null(up_genes) && length(up_genes) > 0) {
  up_map <- bitr(up_genes, fromType = gene_id_type, toType = "ENTREZID", OrgDb = orgdb)
  up_entrez <- unique(up_map$ENTREZID)
  up_readable <- if (gene_id_type == "SYMBOL") unique(up_genes) else {
    unique(bitr(up_genes, fromType = gene_id_type, toType = "SYMBOL", OrgDb = orgdb)$SYMBOL)
  }
}
if (!is.null(down_genes) && length(down_genes) > 0) {
  down_map <- bitr(down_genes, fromType = gene_id_type, toType = "ENTREZID", OrgDb = orgdb)
  down_entrez <- unique(down_map$ENTREZID)
  down_readable <- if (gene_id_type == "SYMBOL") unique(down_genes) else {
    unique(bitr(down_genes, fromType = gene_id_type, toType = "SYMBOL", OrgDb = orgdb)$SYMBOL)
  }
}

write.csv(data.frame(gene = gene_list), file.path(output_dir, "gene_list.csv"), row.names = FALSE)
if (!is.null(up_genes)) write.csv(data.frame(gene = up_genes), file.path(output_dir, "up_genes.csv"), row.names = FALSE)
if (!is.null(down_genes)) write.csv(data.frame(gene = down_genes), file.path(output_dir, "down_genes.csv"), row.names = FALSE)

# ============================================================
# GO & KEGG Enrichment
# ============================================================
run_go <- function(ont) {
  tryCatch({
    enrichGO(gene = entrez_ids, OrgDb = orgdb, keyType = "ENTREZID",
             ont = ont, pAdjustMethod = "BH",
             pvalueCutoff = enrich_pval, qvalueCutoff = enrich_qval, readable = TRUE)
  }, error = function(e) {
    cat(sprintf("  GO-%s failed: %s\n", ont, e$message)); NULL
  })
}

# Direction belongs to the signed DEG genes overlapping each term, not to the
# enrichment p-value itself. Preserve the rank statistics and add an auditable
# overlap summary for Result Studio Compare.
annotate_enrichment_direction <- function(table, up_ids, down_ids) {
  result <- as.data.frame(table)
  if (!("geneID" %in% colnames(result))) {
    stop("Enrichment result is missing the declared geneID overlap column")
  }
  if (is.null(up_ids) || is.null(down_ids)) {
    result$up_gene_count <- NA_integer_
    result$down_gene_count <- NA_integer_
    result$direction <- "unavailable_without_signed_deg"
    return(result)
  }
  overlaps <- strsplit(as.character(result$geneID), "/", fixed = TRUE)
  result$up_gene_count <- vapply(overlaps, function(ids) sum(ids %in% up_ids), integer(1))
  result$down_gene_count <- vapply(overlaps, function(ids) sum(ids %in% down_ids), integer(1))
  result$direction <- ifelse(
    result$up_gene_count > result$down_gene_count,
    "up",
    ifelse(
      result$down_gene_count > result$up_gene_count,
      "down",
      ifelse(result$up_gene_count + result$down_gene_count > 0, "mixed", "unresolved")
    )
  )
  result
}

ego_bp <- NULL; ego_mf <- NULL; ego_cc <- NULL
if (analysis_type %in% c("all", "GO-BP", "GO")) {
  cat("[go-kegg] enrichGO BP...\n")
  ego_bp <- run_go("BP")
  if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    write.csv(annotate_enrichment_direction(ego_bp, up_readable, down_readable),
              file.path(output_dir, "go_bp_results.csv"), row.names = FALSE)
    cat(sprintf("  GO-BP: %d terms\n", nrow(as.data.frame(ego_bp))))
  }
}
if (analysis_type %in% c("all", "GO-MF", "GO")) {
  cat("[go-kegg] enrichGO MF...\n")
  ego_mf <- run_go("MF")
  if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
    write.csv(annotate_enrichment_direction(ego_mf, up_readable, down_readable),
              file.path(output_dir, "go_mf_results.csv"), row.names = FALSE)
    cat(sprintf("  GO-MF: %d terms\n", nrow(as.data.frame(ego_mf))))
  }
}
if (analysis_type %in% c("all", "GO-CC", "GO")) {
  cat("[go-kegg] enrichGO CC...\n")
  ego_cc <- run_go("CC")
  if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0) {
    write.csv(annotate_enrichment_direction(ego_cc, up_readable, down_readable),
              file.path(output_dir, "go_cc_results.csv"), row.names = FALSE)
    cat(sprintf("  GO-CC: %d terms\n", nrow(as.data.frame(ego_cc))))
  }
}

kegg_res <- NULL
if (analysis_type %in% c("all", "KEGG")) {
  cat("[go-kegg] enrichKEGG...\n")
  tryCatch({
    kegg_res <- run_kegg_ora(entrez_ids, kegg_org, enrich_pval, enrich_qval)
    if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
      write.csv(annotate_enrichment_direction(kegg_res, up_entrez, down_entrez),
                file.path(output_dir, "kegg_results.csv"), row.names = FALSE)
      cat(sprintf("  KEGG: %d pathways\n", nrow(as.data.frame(kegg_res))))
    }
  }, error = function(e) cat(sprintf("  KEGG failed: %s\n", e$message)))
}

# ============================================================
# Plots — SCI 风格双格式
# ============================================================
cat("[go-kegg] Generating plots...\n")
# enrichplot network layouts consume R's RNG. Bind the seed before any plot so
# formal reports that embed cnetplot/emapplot remain byte-identical on replay.
set.seed(20260729)

# ---- 1-2: GO-BP dotplot + barplot (ggplot 系) ----
if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  n_show <- min(top_n, nrow(as.data.frame(ego_bp)))
  p1 <- dotplot(ego_bp, showCategory = n_show) +
    labs(title = "GO Biological Process") +
    theme_biof3(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
  ggsave_biof3(p1, file.path(output_dir, "go_bp_dotplot"), width = 8, height = 7)

  p2 <- barplot(ego_bp, showCategory = n_show) +
    labs(title = "GO-BP Top Enriched Terms") +
    theme_biof3(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
  ggsave_biof3(p2, file.path(output_dir, "go_bp_barplot"), width = 8, height = 6)

  # ---- 3-6: cnetplot / treeplot / heatplot / emapplot (ggraph 复合对象, 用 save_grid_biof3) ----
  tryCatch({
    p3 <- cnetplot(ego_bp, showCategory = min(5, nrow(as.data.frame(ego_bp))),
                   foldChange = NULL) +
      ggtitle("GO-BP Gene-Concept Network")
    save_grid_biof3(file.path(output_dir, "go_bp_cnetplot"), width = 10, height = 8, expr = print(p3))
  }, error = function(e) cat(sprintf("  cnetplot skipped: %s\n", e$message)))

  ego_bp_sim <- NULL
  tryCatch({
    ego_bp_sim <- pairwise_termsim(ego_bp)
    p4 <- treeplot(ego_bp_sim, showCategory = min(30, nrow(as.data.frame(ego_bp)))) +
      ggtitle("GO-BP Enrichment Tree")
    save_grid_biof3(file.path(output_dir, "go_bp_treeplot"), width = 12, height = 8, expr = print(p4))
  }, error = function(e) cat(sprintf("  treeplot skipped: %s\n", e$message)))

  tryCatch({
    p5 <- heatplot(ego_bp, showCategory = min(15, nrow(as.data.frame(ego_bp)))) +
      ggtitle("GO-BP Enrichment Heatmap")
    save_grid_biof3(file.path(output_dir, "go_bp_heatplot"), width = 12, height = 6, expr = print(p5))
  }, error = function(e) cat(sprintf("  heatplot skipped: %s\n", e$message)))

  if (!is.null(ego_bp_sim) && nrow(as.data.frame(ego_bp)) >= 5) {
    tryCatch({
      p6 <- emapplot(ego_bp_sim, showCategory = min(30, nrow(as.data.frame(ego_bp)))) +
        ggtitle("GO-BP Enrichment Map")
      save_grid_biof3(file.path(output_dir, "go_bp_emapplot"), width = 10, height = 9, expr = print(p6))
    }, error = function(e) cat(sprintf("  emapplot skipped: %s\n", e$message)))
  }
}

# ---- 7: GO-MF dotplot ----
if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
  p <- dotplot(ego_mf, showCategory = min(top_n, nrow(as.data.frame(ego_mf)))) +
    labs(title = "GO Molecular Function") +
    theme_biof3(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
  ggsave_biof3(p, file.path(output_dir, "go_mf_dotplot"), width = 8, height = 7)
}

# ---- 8: GO-CC dotplot ----
if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0) {
  p <- dotplot(ego_cc, showCategory = min(top_n, nrow(as.data.frame(ego_cc)))) +
    labs(title = "GO Cellular Component") +
    theme_biof3(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
  ggsave_biof3(p, file.path(output_dir, "go_cc_dotplot"), width = 8, height = 7)
}

# ---- 9-10: KEGG dotplot + barplot ----
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  n_show <- min(top_n, nrow(as.data.frame(kegg_res)))
  p1 <- dotplot(kegg_res, showCategory = n_show) +
    labs(title = "KEGG Pathway Enrichment") +
    theme_biof3(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
  ggsave_biof3(p1, file.path(output_dir, "kegg_dotplot"), width = 8, height = 7)

  p2 <- barplot(kegg_res, showCategory = n_show) +
    labs(title = "KEGG Top Pathways") +
    theme_biof3(base_size = 10) +
    theme(axis.text.y = element_text(size = 9))
  ggsave_biof3(p2, file.path(output_dir, "kegg_barplot"), width = 8, height = 6)
}

# ---- 11: compareCluster (Up vs Down) ----
if (!is.null(up_entrez) && !is.null(down_entrez) &&
    length(up_entrez) >= 3 && length(down_entrez) >= 3) {
  cat("[go-kegg] compareCluster Up vs Down...\n")
  tryCatch({
    gene_clusters <- list(Up = up_entrez, Down = down_entrez)
    cc_go <- compareCluster(gene_clusters, fun = "enrichGO",
                            OrgDb = orgdb, ont = "BP",
                            pvalueCutoff = enrich_pval, readable = TRUE)
    if (!is.null(cc_go) && nrow(as.data.frame(cc_go)) > 0) {
      p <- dotplot(cc_go, showCategory = min(10, nrow(as.data.frame(cc_go)))) +
        labs(title = "GO-BP: Up vs Down Regulated") +
        theme_biof3(base_size = 10) +
        theme(axis.text.y = element_text(size = 9))
      ggsave_biof3(p, file.path(output_dir, "compare_up_down"), width = 9, height = 8)
      write.csv(as.data.frame(cc_go), file.path(output_dir, "compare_cluster_results.csv"), row.names = FALSE)
    }
  }, error = function(e) cat(sprintf("  compareCluster skipped: %s\n", e$message)))
}

# The formal HTML report embeds every generated PNG. Canonical encoding removes
# compression metadata/stream drift from the complete report dependency set.
for (png_path in list.files(output_dir, pattern = "\\.png$", full.names = TRUE)) {
  canonicalize_png(png_path)
}

# ============================================================
# Summary
# ============================================================
n_bp <- if (!is.null(ego_bp)) nrow(as.data.frame(ego_bp)) else 0
n_mf <- if (!is.null(ego_mf)) nrow(as.data.frame(ego_mf)) else 0
n_cc <- if (!is.null(ego_cc)) nrow(as.data.frame(ego_cc)) else 0
n_kegg <- if (!is.null(kegg_res)) nrow(as.data.frame(kegg_res)) else 0

summary_text <- sprintf(
  "GO/KEGG Enrichment Analysis Summary\n\nInput genes: %d\nSpecies: %s\nGene ID type: %s\n\nResults:\n  GO-BP: %d enriched terms\n  GO-MF: %d enriched terms\n  GO-CC: %d enriched terms\n  KEGG: %d enriched pathways\n\nCutoffs:\n  Enrichment p-value: %g\n  Enrichment q-value: %g\n",
  length(gene_list), species, gene_id_type, n_bp, n_mf, n_cc, n_kegg,
  enrich_pval, enrich_qval
)
if (!is.null(up_genes)) {
  summary_text <- paste0(summary_text,
    sprintf("\nDESeq2 input:\n  Up-regulated: %d\n  Down-regulated: %d\n  padj cutoff: %g\n  |log2FC| cutoff: %g\n",
            length(up_genes), length(down_genes), padj_cut, lfc_cut))
}
writeLines(summary_text, file.path(output_dir, "summary.txt"))

# ============================================================
# 中间数据 (RDS for online replotting)
# ============================================================
cat("[go-kegg] Saving intermediate data...\n")
if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  write.csv(annotate_enrichment_direction(ego_bp, up_readable, down_readable), file.path(output_dir, "go_bp_data.csv"), row.names = FALSE)
  saveRDS(ego_bp, file.path(output_dir, "ego_bp.rds"))
}
if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
  write.csv(annotate_enrichment_direction(ego_mf, up_readable, down_readable), file.path(output_dir, "go_mf_data.csv"), row.names = FALSE)
  saveRDS(ego_mf, file.path(output_dir, "ego_mf.rds"))
}
if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0) {
  write.csv(annotate_enrichment_direction(ego_cc, up_readable, down_readable), file.path(output_dir, "go_cc_data.csv"), row.names = FALSE)
  saveRDS(ego_cc, file.path(output_dir, "ego_cc.rds"))
}
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  write.csv(annotate_enrichment_direction(kegg_res, up_entrez, down_entrez), file.path(output_dir, "kegg_data.csv"), row.names = FALSE)
  saveRDS(kegg_res, file.path(output_dir, "kegg_res.rds"))
}

# ============================================================
# Manifest (含 PNG + PDF 双格式)
# ============================================================
files_list <- list(
  list(name = "summary.txt", type = "text", label = "分析摘要")
)

# tables
for (csv_label in list(
  c("go_bp_results.csv", "GO-BP 结果表"),
  c("go_mf_results.csv", "GO-MF 结果表"),
  c("go_cc_results.csv", "GO-CC 结果表"),
  c("kegg_results.csv",  "KEGG 结果表"),
  c("gene_annotation.csv", "基因 ID 转换"),
  c("compare_cluster_results.csv", "上下调对比")
)) {
  if (file.exists(file.path(output_dir, csv_label[1]))) {
    files_list <- c(files_list, list(list(name = csv_label[1], type = "table", label = csv_label[2])))
  }
}

# plots — PNG + PDF 双格式
plot_pairs <- list(
  c("go_bp_dotplot",   "GO-BP 气泡图"),
  c("go_bp_barplot",   "GO-BP 条形图"),
  c("go_bp_cnetplot",  "GO-BP 基因网络"),
  c("go_bp_treeplot",  "GO-BP 富集树图"),
  c("go_bp_heatplot",  "GO-BP 热图"),
  c("go_bp_emapplot",  "GO-BP 富集网络"),
  c("go_mf_dotplot",   "GO-MF 气泡图"),
  c("go_cc_dotplot",   "GO-CC 气泡图"),
  c("kegg_dotplot",    "KEGG 气泡图"),
  c("kegg_barplot",    "KEGG 条形图"),
  c("compare_up_down", "上下调对比")
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

files_list <- c(files_list, list(list(name = "report.html", type = "file", label = "解读报告(HTML)")))

# Add RDS intermediate data files (for online replotting)
for (rds_label in list(
  c("ego_bp.rds", "GO-BP RDS"),
  c("ego_mf.rds", "GO-MF RDS"),
  c("ego_cc.rds", "GO-CC RDS"),
  c("kegg_res.rds", "KEGG RDS"),
  c("go_bp_data.csv", "GO-BP 数据"),
  c("go_mf_data.csv", "GO-MF 数据"),
  c("go_cc_data.csv", "GO-CC 数据"),
  c("kegg_data.csv",  "KEGG 数据"),
  c("gene_list.csv",  "基因列表"),
  c("up_genes.csv",   "上调基因"),
  c("down_genes.csv", "下调基因")
)) {
  if (file.exists(file.path(output_dir, rds_label[1]))) {
    typ <- if (endsWith(rds_label[1], ".rds")) "file" else "table"
    files_list <- c(files_list, list(list(name = rds_label[1], type = typ, label = rds_label[2])))
  }
}

manifest <- list(
  files = files_list,
  summary = list(
    input_genes = length(gene_list),
    go_bp = n_bp, go_mf = n_mf, go_cc = n_cc, kegg = n_kegg,
    species = species,
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
  cat("[go-kegg] Generating report...\n")
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
  for (img_path in png_files) {
    img_name <- basename(img_path)
    b64 <- base64enc::base64encode(img_path)
    data_uri <- paste0("data:image/png;base64,", b64)
    report_html <- gsub(paste0("{{", img_name, "}}"), data_uri, report_html, fixed = TRUE)
  }

  report_html <- gsub("{{input_genes}}", as.character(length(gene_list)), report_html, fixed = TRUE)
  report_html <- gsub("{{species}}", species, report_html, fixed = TRUE)
  report_html <- gsub("{{n_bp}}", as.character(n_bp), report_html, fixed = TRUE)
  report_html <- gsub("{{n_mf}}", as.character(n_mf), report_html, fixed = TRUE)
  report_html <- gsub("{{n_cc}}", as.character(n_cc), report_html, fixed = TRUE)
  report_html <- gsub("{{n_kegg}}", as.character(n_kegg), report_html, fixed = TRUE)

  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[go-kegg] Report generated.\n")
} else {
  cat("[go-kegg] Report template not found, skipping.\n")
}

cat(sprintf("[go-kegg] Pipeline complete. GO-BP=%d MF=%d CC=%d KEGG=%d\n",
            n_bp, n_mf, n_cc, n_kegg))
