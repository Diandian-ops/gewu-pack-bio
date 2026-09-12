#!/usr/bin/env Rscript
# ============================================================
# Tool: cellchat
# Version targeted: CellChat 2.2.0 (R 4.5.1)
# Canonical workflow source:
#   https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
# Reference:
#   Jin S et al. (2025). CellChat for systematic analysis of cell-cell communication
#   from single-cell transcriptomics. Nature Protocols.
#
# 第二个使用 BioF3 SCI 视觉风格规范的工具 (theme_biof3 + save_grid_biof3 双格式).
#
# Steps (in-scope) — Single-Dataset Workflow:
#   1. Load expression + meta
#   2. createCellChat → 选 LR database
#   3. subsetData → identifyOverExpressedGenes → identifyOverExpressedInteractions
#   4. computeCommunProb → filterCommunication
#   5. computeCommunProbPathway → aggregateNet
#   6. 可视化: 网络图 + 热图 + 气泡 + 角色 + top pathway chord/violin
#
# Out-of-scope:
#   - 多数据集对比 (mergeCellChat)
#   - Spatial CellChat
#   - NMF pattern analysis
#   - 单 pathway 切片热图 (留给 RunModal)
#
# Default params:
#   species=human, signaling_db="Secreted Signaling"
#   min_cells_per_group=10, top_pathways_to_visualize=5
#
# Input (job_dir/):
#   - expression  CSV/TSV: 行=基因, 列=细胞 (log-normalized, NOT raw)
#   - meta        CSV/TSV: cell_id + cell_type
#   - params.json
#
# Output (job_dir/output/):
#   8 plots × PNG + PDF = 16 files
#   + 4 中间数据 + cellchat.rds + report.html + manifest.json
# ============================================================

suppressMessages({
  library(CellChat)
  library(circlize)
  library(jsonlite)
  library(dplyr)
  library(ggplot2)
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(job_dir, "_biof3-theme.R"))

cat("[cellchat] Starting pipeline\n")
cat("[cellchat] CellChat:", as.character(packageVersion("CellChat")), "\n")

# ---- Read params ----
params <- fromJSON(file.path(job_dir, "params.json"))
species              <- params$species %||% "human"
signaling_db_choice  <- params$signaling_db %||% "Secreted Signaling"
min_cells_per_group  <- as.integer(params$min_cells_per_group %||% 10)
top_n_pathways       <- as.integer(params$top_pathways_to_visualize %||% 5)

cat(sprintf("[cellchat] params: species=%s, db=%s, min_cells=%d, top_n=%d\n",
            species, signaling_db_choice, min_cells_per_group, top_n_pathways))

# ---- Read inputs ----
read_table_auto <- function(file, row_names_col = 1) {
  first_line <- readLines(file, n = 1, warn = FALSE)
  sep <- if (grepl("\t", first_line) && !grepl(",", first_line)) "\t" else ","
  read.table(file, header = TRUE, sep = sep, row.names = row_names_col,
             check.names = FALSE, stringsAsFactors = FALSE)
}

expr_files <- list.files(job_dir, pattern = "^expression$", full.names = TRUE)
meta_files <- list.files(job_dir, pattern = "^meta$", full.names = TRUE)
if (length(expr_files) == 0) stop("No expression file found")
if (length(meta_files) == 0) stop("No meta file found")

cat("[cellchat] Loading expression + meta...\n")
expr <- as.matrix(read_table_auto(expr_files[1]))
storage.mode(expr) <- "numeric"

# Meta — could have cell_id as first col or as rownames
meta_raw <- read.table(meta_files[1], header = TRUE,
                        sep = if (grepl("\t", readLines(meta_files[1], 1))) "\t" else ",",
                        check.names = FALSE, stringsAsFactors = FALSE)
if ("cell_id" %in% colnames(meta_raw)) {
  rownames(meta_raw) <- meta_raw$cell_id
  meta_raw$cell_id <- NULL
} else {
  # Assume first column is cell_id
  rownames(meta_raw) <- meta_raw[, 1]
  meta_raw <- meta_raw[, -1, drop = FALSE]
}

if (!"cell_type" %in% colnames(meta_raw)) {
  stop("meta must have a 'cell_type' column")
}

cat(sprintf("  expression: %d genes × %d cells\n", nrow(expr), ncol(expr)))
cat(sprintf("  meta: %d rows, cell_type distribution:\n", nrow(meta_raw)))
print(table(meta_raw$cell_type))

# Match cells
common <- intersect(colnames(expr), rownames(meta_raw))
if (length(common) < ncol(expr) * 0.8) {
  stop(sprintf("Only %d / %d cells matched between expression and meta", length(common), ncol(expr)))
}
expr <- expr[, common, drop = FALSE]
meta <- meta_raw[common, , drop = FALSE]

# Validate
if (ncol(expr) > 5000) {
  stop(sprintf("Too many cells (%d) — desktop limit is 5000. Subsample first.", ncol(expr)))
}
n_types <- length(unique(meta$cell_type))
if (n_types < 3) {
  stop(sprintf("Need at least 3 cell types for CellChat (got %d)", n_types))
}

# ---- Step 1: Create CellChat object ----
cat("[cellchat] Step 1: createCellChat...\n")
cellchat <- createCellChat(object = expr, meta = meta, group.by = "cell_type")

# ---- Step 2: Choose database ----
cat("[cellchat] Step 2: Loading LR database...\n")
db <- if (species == "mouse") CellChatDB.mouse else CellChatDB.human

if (signaling_db_choice != "all") {
  db_use <- subsetDB(db, search = signaling_db_choice)
  cat(sprintf("  Subset to '%s': %d interactions\n", signaling_db_choice, nrow(db_use$interaction)))
} else {
  db_use <- db
  cat(sprintf("  Using all subsets: %d interactions\n", nrow(db_use$interaction)))
}
cellchat@DB <- db_use

# ---- Step 3: Preprocessing ----
cat("[cellchat] Step 3: Preprocessing...\n")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cat(sprintf("  Highly variable LR pairs: %d\n", nrow(cellchat@LR$LRsig)))

if (nrow(cellchat@LR$LRsig) < 5) {
  warning(sprintf("Only %d highly variable LR pairs — results may be sparse. Try expanding the signaling_db or relaxing thresholds.",
                  nrow(cellchat@LR$LRsig)))
}

# ---- Step 4: Inference ----
cat("[cellchat] Step 4: computeCommunProb (this can take 1-3 min)...\n")
t_start <- Sys.time()
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = min_cells_per_group)
elapsed <- as.numeric(Sys.time() - t_start)
cat(sprintf("  Done in %.1fs\n", elapsed))

# ---- Step 5: Pathway-level aggregation ----
cat("[cellchat] Step 5: Pathway aggregation...\n")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

n_pathways <- length(cellchat@netP$pathways)
cat(sprintf("  Found %d significant pathways\n", n_pathways))
if (n_pathways == 0) {
  stop("No significant pathways detected — check input data quality (cell types coverage / gene depth)")
}
cat(sprintf("  Top pathways: %s\n", paste(head(cellchat@netP$pathways, 5), collapse = ", ")))

# ============================================================
# Visualizations
# ============================================================
cat("[cellchat] Generating plots...\n")
cell_types <- levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

# ---- Plot 1+2: Circle network (count + weight) ----
save_grid_biof3(file.path(output_dir, "network_circle_count"),
                 width = 7, height = 7, expr = {
  par(xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Number of interactions")
})

save_grid_biof3(file.path(output_dir, "network_circle_weight"),
                 width = 7, height = 7, expr = {
  par(xpd = TRUE)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Interaction strength")
})

# ---- Plot 3: Communication heatmap ----
save_grid_biof3(file.path(output_dir, "network_heatmap"),
                 width = 8, height = 6, expr = {
  print(netVisual_heatmap(cellchat, measure = "weight",
                           color.heatmap = "Reds",
                           title.name = "Interaction strength heatmap"))
})

# ---- Plot 4: Top LR bubble ----
# Use all sources × all targets (CellChat's bubble plot)
ggsave_biof3(
  netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL,
                    remove.isolate = FALSE, return.data = FALSE) +
    theme_biof3(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9)),
  file.path(output_dir, "bubble_top_LR"),
  width = max(8, length(cell_types) * 1.2),
  height = max(8, 12)
)

# ---- Plot 5: Signaling role scatter ----
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

save_grid_biof3(file.path(output_dir, "signaling_role_scatter"),
                 width = 7, height = 6, expr = {
  print(netAnalysis_signalingRole_scatter(cellchat) +
          theme_biof3(base_size = 11))
})

# ---- Plot 6: Top pathway aggregate (chord) ----
top_pathway <- cellchat@netP$pathways[1]
save_grid_biof3(file.path(output_dir, "pathway_aggregate_top1"),
                 width = 7, height = 7, expr = {
  par(xpd = TRUE)
  netVisual_aggregate(cellchat, signaling = top_pathway, layout = "chord")
})

# ---- Plot 7: Top pathway violin ----
top_pathways_violin <- head(cellchat@netP$pathways, top_n_pathways)
save_grid_biof3(file.path(output_dir, "pathway_violin_top"),
                 width = max(8, length(top_pathways_violin) * 2),
                 height = 6, expr = {
  print(plotGeneExpression(cellchat, signaling = top_pathway, enriched.only = FALSE) +
          theme_biof3(base_size = 10))
})

# ---- Plot 8: Pathway count per cell type ----
# Compute outgoing / incoming pathway counts per cell type
n_pw <- length(cellchat@netP$pathways)
out_count <- in_count <- numeric(length(cell_types))
names(out_count) <- names(in_count) <- cell_types
for (pw in cellchat@netP$pathways) {
  prob <- cellchat@netP$prob[, , pw]
  for (ct in cell_types) {
    if (any(prob[ct, ] > 0)) out_count[ct] <- out_count[ct] + 1
    if (any(prob[, ct] > 0)) in_count[ct]  <- in_count[ct]  + 1
  }
}
pw_df <- data.frame(
  cell_type = rep(cell_types, 2),
  count     = c(out_count, in_count),
  direction = factor(rep(c("Outgoing", "Incoming"), each = length(cell_types)),
                     levels = c("Outgoing", "Incoming"))
)
p_pw <- ggplot(pw_df, aes(x = reorder(cell_type, -count), y = count, fill = direction)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  scale_fill_manual(values = c(Outgoing = biof3_palette_div(3)[3], Incoming = biof3_palette_div(3)[1])) +
  labs(x = NULL, y = "Number of pathways", fill = NULL,
       title = sprintf("Pathway count per cell type (%d total pathways)", n_pw)) +
  theme_biof3() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave_biof3(p_pw, file.path(output_dir, "pathway_count_barplot"),
              width = max(8, length(cell_types) * 0.8 + 2), height = 5)

# ============================================================
# Intermediate data
# ============================================================
cat("[cellchat] Saving intermediate data...\n")

# Significant LR pairs
df_net <- subsetCommunication(cellchat)
write.csv(df_net, file.path(output_dir, "lr_significant.csv"), row.names = FALSE)

# Pathway significance ranking
pathway_strength <- sapply(cellchat@netP$pathways, function(pw) {
  sum(cellchat@netP$prob[, , pw])
})
pw_rank_df <- data.frame(
  pathway = names(pathway_strength),
  total_strength = pathway_strength
)
pw_rank_df <- pw_rank_df[order(-pw_rank_df$total_strength), ]
write.csv(pw_rank_df, file.path(output_dir, "pathway_significance.csv"), row.names = FALSE)

# Signaling role data
role_in  <- rowSums(cellchat@net$weight)
role_out <- colSums(cellchat@net$weight)
role_df <- data.frame(
  cell_type = names(role_out),
  outgoing_strength = role_out,
  incoming_strength = role_in,
  outgoing_pathway_count = out_count[names(role_out)],
  incoming_pathway_count = in_count[names(role_in)]
)
write.csv(role_df, file.path(output_dir, "signaling_role_data.csv"), row.names = FALSE)

# Save full CellChat object
saveRDS(cellchat, file.path(output_dir, "cellchat.rds"))

# ============================================================
# Summary
# ============================================================
summary_text <- sprintf(
"CellChat Single-Dataset Analysis Summary

Input:
  Expression: %d genes × %d cells
  Cell types: %d (%s)
  Database: %s (%s, %d interactions)

Analysis:
  Highly variable LR pairs: %d
  Significant pathways: %d
  Top 5 pathways: %s

Communication network:
  Total significant LR pairs: %d
  Cell types with most outgoing: %s (%d pathways)
  Cell types with most incoming: %s (%d pathways)

Runtime: %.1f sec\n",
  nrow(expr), ncol(expr),
  n_types, paste(cell_types, collapse = ", "),
  species, signaling_db_choice, nrow(db_use$interaction),
  nrow(cellchat@LR$LRsig),
  n_pathways,
  paste(head(cellchat@netP$pathways, 5), collapse = ", "),
  nrow(df_net),
  names(which.max(out_count)), max(out_count),
  names(which.max(in_count)), max(in_count),
  elapsed
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))
cat(summary_text)

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = list(
    list(name = "network_circle_count.png",   type = "plot",  label = "通讯网络（互作数）PNG"),
    list(name = "network_circle_count.pdf",   type = "file",  label = "通讯网络（互作数）PDF"),
    list(name = "network_circle_weight.png",  type = "plot",  label = "通讯网络（互作强度）PNG"),
    list(name = "network_circle_weight.pdf",  type = "file",  label = "通讯网络（互作强度）PDF"),
    list(name = "network_heatmap.png",        type = "plot",  label = "通讯热图 PNG"),
    list(name = "network_heatmap.pdf",        type = "file",  label = "通讯热图 PDF"),
    list(name = "bubble_top_LR.png",          type = "plot",  label = "Top LR 气泡 PNG"),
    list(name = "bubble_top_LR.pdf",          type = "file",  label = "Top LR 气泡 PDF"),
    list(name = "signaling_role_scatter.png", type = "plot",  label = "信号角色 PNG"),
    list(name = "signaling_role_scatter.pdf", type = "file",  label = "信号角色 PDF"),
    list(name = "pathway_aggregate_top1.png", type = "plot",  label = "Top pathway chord PNG"),
    list(name = "pathway_aggregate_top1.pdf", type = "file",  label = "Top pathway chord PDF"),
    list(name = "pathway_violin_top.png",     type = "plot",  label = "Top pathway violin PNG"),
    list(name = "pathway_violin_top.pdf",     type = "file",  label = "Top pathway violin PDF"),
    list(name = "pathway_count_barplot.png",  type = "plot",  label = "通路数 barplot PNG"),
    list(name = "pathway_count_barplot.pdf",  type = "file",  label = "通路数 barplot PDF"),
    list(name = "lr_significant.csv",         type = "table", label = "显著 LR 对"),
    list(name = "pathway_significance.csv",   type = "table", label = "通路显著性排序"),
    list(name = "signaling_role_data.csv",    type = "table", label = "信号角色数据"),
    list(name = "cellchat.rds",               type = "file",  label = "完整 CellChat 对象"),
    list(name = "summary.txt",                type = "text",  label = "分析摘要"),
    list(name = "report.html",                type = "file",  label = "解读报告 HTML")
  ),
  summary = list(
    n_genes = nrow(expr),
    n_cells = ncol(expr),
    n_cell_types = n_types,
    n_lr_pairs = nrow(cellchat@LR$LRsig),
    n_significant_pathways = n_pathways,
    top_pathways = head(cellchat@netP$pathways, 5),
    top_outgoing = names(which.max(out_count)),
    top_incoming = names(which.max(in_count)),
    runtime_sec = round(elapsed, 1)
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ============================================================
# HTML report
# ============================================================
cat("[cellchat] Generating HTML report...\n")
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  png_files <- c("network_circle_count.png", "network_circle_weight.png",
                 "network_heatmap.png", "bubble_top_LR.png",
                 "signaling_role_scatter.png", "pathway_aggregate_top1.png",
                 "pathway_violin_top.png", "pathway_count_barplot.png")
  for (png in png_files) {
    img_path <- file.path(output_dir, png)
    placeholder <- paste0("{{", png, "}}")
    if (file.exists(img_path)) {
      b64 <- base64enc::base64encode(img_path)
      data_uri <- paste0("data:image/png;base64,", b64)
      report_html <- gsub(placeholder, data_uri, report_html, fixed = TRUE)
    }
  }

  report_html <- gsub("{{n_genes}}",         as.character(nrow(expr)),                         report_html, fixed = TRUE)
  report_html <- gsub("{{n_cells}}",         as.character(ncol(expr)),                         report_html, fixed = TRUE)
  report_html <- gsub("{{n_cell_types}}",    as.character(n_types),                            report_html, fixed = TRUE)
  report_html <- gsub("{{species}}",         species,                                          report_html, fixed = TRUE)
  report_html <- gsub("{{signaling_db}}",    signaling_db_choice,                              report_html, fixed = TRUE)
  report_html <- gsub("{{n_lr_pairs}}",      as.character(nrow(cellchat@LR$LRsig)),            report_html, fixed = TRUE)
  report_html <- gsub("{{n_pathways}}",      as.character(n_pathways),                         report_html, fixed = TRUE)
  report_html <- gsub("{{top_pathways}}",    paste(head(cellchat@netP$pathways, 5), collapse = ", "), report_html, fixed = TRUE)
  report_html <- gsub("{{top_outgoing}}",    names(which.max(out_count)),                      report_html, fixed = TRUE)
  report_html <- gsub("{{top_incoming}}",    names(which.max(in_count)),                       report_html, fixed = TRUE)
  report_html <- gsub("{{n_significant_LR}}", as.character(nrow(df_net)),                      report_html, fixed = TRUE)

  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[cellchat] Report generated.\n")
}

cat("[cellchat] Pipeline complete.\n")
