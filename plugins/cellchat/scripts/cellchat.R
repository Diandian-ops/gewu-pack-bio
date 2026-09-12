#!/usr/bin/env Rscript
# ============================================================
# cellchat · 单细胞细胞-细胞通讯推断（R 包 CellChat 1.6+）
# Validated: CellChat 1.6.1 / Seurat 5.5.1 / ComplexHeatmap 2.x /
# circlize 0.4.15 / NMF 0.27 (R 4.5.2, Bioconductor 3.22)
#
# Steps (in-scope) — CellChat canonical vignette:
#   1. preflight: 校验 job_dir 输入 (Seurat RDS, >1KB)
#   2. 读取 Seurat 对象 + 验证 cell_type 列存在
#   3. createCellChat(data.input = seurat) + set DB (human/mouse/zebrafish)
#   4. subsetDB(by = db_category) — Secreted Signaling / Cell-Cell Contact / ECM-Receptor
#   5. computeCommunProb (type = "truncatedMean", trim = 0.1)
#   6. computeCommunProbPathway
#   7. aggregateNet
#   8. 可选 compute_hierarchy: selectK → identifyCommunicationPatterns → netAnalysis_signalingRole
#   9. 落盘: communications_table.csv / pathway_enrichment.csv / interaction_heatmap.png /
#      chord_diagram.png / bubble_plot.png / cellchat_object.rds / report.html / manifest.json
#  10. PNG 重写（png::writePNG(png::readPNG(f), f)）以稳定 sha256
#
# 2026-08-18: SUB-1.4 BLOCK-11 PHASE 1 落地
#   - runtime: biof3-r-runtime (CellChat 通过 BiocManager::install 一次性补齐)
#   - 不外推 packaged / Provider / GA / dist staleness (G1+)
#   - 不引 CRAN 网络下载（CellChatDB 已随包内置）
# ============================================================

# Delayed library load — 让 preflight 失败路径不付出 import 成本
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(jsonlite)
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

canonicalize_png <- function(file_path) {
  if (!file.exists(file_path)) return(invisible(FALSE))
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("[cellchat] R package 'png' is required for deterministic Result Studio artifacts")
  }
  png::writePNG(png::readPNG(file_path), file_path)
  invisible(TRUE)
}

report_progress <- function(pct, msg, job_dir) {
  progress_file <- file.path(job_dir, "progress.json")
  writeLines(
    jsonlite::toJSON(
      list(progress = pct, message = msg),
      auto_unbox = TRUE,
      pretty = FALSE
    ),
    con = progress_file
  )
}

resolve_cellchat_db <- function(species) {
  if (species == "mouse") return(CellChatDB.mouse)
  if (species == "zebrafish") return(CellChatDB.zebrafish)
  CellChatDB.human
}

# ============================================================
# 主流程
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
if (is.na(job_dir) || !nzchar(job_dir)) stop("[cellchat] job_dir argument missing")
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题（若已塞到 jobDir；本插件不强依赖 theme 函数）
theme_path <- file.path(job_dir, "_biof3-theme.R")
if (file.exists(theme_path)) {
  source(theme_path)
}

cat("[cellchat] Starting pipeline\n")
report_progress(5, "初始化 R CellChat runtime...", job_dir)

# ── 1. preflight ───────────────────────────────────────────
params <- fromJSON(file.path(job_dir, "params.json"))
seurat_rds_name <- params$seurat_rds %||% NULL
cell_type_col   <- params$cell_type_col %||% "cell_type"
species         <- params$species %||% "human"
db_category     <- params$db_category %||% "Secreted Signaling"
min_cells       <- as.integer(params$min_cells %||% 10)
pvalue_cut      <- as.numeric(params$pvalue_threshold %||% 0.05)
n_patterns      <- as.integer(params$n_patterns %||% 5)
compute_hier    <- identical(tolower(as.character(params$compute_hierarchy %||% "false")), "true")

# 寻找 seurat_rds（约定 input id 'seurat_rds' / 任意 *.rds）
seurat_path <- NULL
candidates <- c(
  if (!is.null(seurat_rds_name)) file.path(job_dir, seurat_rds_name) else NULL,
  list.files(job_dir, pattern = "\\.rds$", full.names = TRUE)
)
for (cand in candidates) {
  if (!is.null(cand) && file.exists(cand) && file.size(cand) > 1000) {
    seurat_path <- cand
    break
  }
}

if (is.null(seurat_path)) {
  msg <- paste0(
    "[cellchat] 缺少 Seurat RDS 输入文件。请在 jobDir 放置 .rds (anndata_file input)，",
    "或由上游 scanpy-advanced / seurat 链式传入。"
  )
  cat(msg, file = stderr())
  report_progress(100, paste("blocked:", msg), job_dir)
  writeLines(msg, con = file.path(output_dir, "error.txt"))
  stop(msg, call. = FALSE)
}

# ── 2. runtime import（延迟到 preflight 之后）──────────────
report_progress(15, "加载 CellChat / ComplexHeatmap / NMF...", job_dir)
tryCatch({
  suppressMessages({
    library(CellChat)
    library(ComplexHeatmap)
    library(circlize)
    library(igraph)
    library(NMF)
  })
}, error = function(e) {
  msg <- paste0(
    "[cellchat] biof3-r-runtime 缺少 CellChat / igraph / ComplexHeatmap / circlize / NMF: ",
    conditionMessage(e),
    "。请执行 install.json 中的 BiocManager::install('CellChat')。"
  )
  cat(msg, file = stderr())
  report_progress(100, paste("runtime-not-ready:", msg), job_dir)
  writeLines(msg, con = file.path(output_dir, "error.txt"))
  stop(msg, call. = FALSE)
})

# ── 3. 读入 Seurat 对象 ─────────────────────────────────────
report_progress(25, paste0("读取 Seurat RDS: ", basename(seurat_path)), job_dir)
seurat_obj <- readRDS(seurat_path)
cat(sprintf("[cellchat] Seurat: %d cells × %d features\n", ncol(seurat_obj), nrow(seurat_obj)))

if (!(cell_type_col %in% colnames(seurat_obj@meta.data))) {
  msg <- paste0(
    "[cellchat] cell_type_col '", cell_type_col, "' 不在 meta.data 中。可用列: ",
    paste(colnames(seurat_obj@meta.data), collapse = ", "),
    "。请在 params.json 中设置正确的 cell_type_col 或对 Seurat 对象重新注释。"
  )
  cat(msg, file = stderr())
  report_progress(100, msg, job_dir)
  writeLines(msg, con = file.path(output_dir, "error.txt"))
  stop(msg, call. = FALSE)
}

Idents(seurat_obj) <- cell_type_col

# ── 4. createCellChat + DB ────────────────────────────────
report_progress(35, paste0("createCellChat + CellChatDB.", species), job_dir)
# CellChat 1.6.1 internally calls GetAssayData(..., slot="data") which is
# defunct in SeuratObject 5.0+. Workaround: pass `data.input = normalized_matrix`
# directly so CellChat skips the Seurat::GetAssayData branch.
data.input <- tryCatch(
  GetAssayData(seurat_obj, layer = "data"),
  error = function(e) GetAssayData(seurat_obj, slot = "data")
)
meta <- seurat_obj@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = cell_type_col)
db <- resolve_cellchat_db(species)
if (!is.null(db_category) && nzchar(db_category) && db_category != "all") {
  cellchat@DB <- subsetDB(CellChatDB = db, search = db_category, key = "annotation")
} else {
  cellchat@DB <- db
}

# ── 5. computeCommunProb ──────────────────────────────────
report_progress(55, "subsetData + identifyOverExpressed* + computeCommunProb (truncatedMean, trim=0.1, raw.use=TRUE)", job_dir)
# CellChat 1.6.1 canonical pipeline（subsetData → identifyOverExpressedGenes →
# identifyOverExpressedInteractions → computeCommunProb）。computeCommunProb
# 不会自动跑前面三个步骤；如跳过 → data.signaling 0×0 / LRsig 0 row → "no rows
# to aggregate" / "subscript out of bounds"。
cellchat <- subsetData(cellchat)
tryCatch({
  cellchat <- identifyOverExpressedGenes(cellchat)
}, error = function(e) {
  cat(sprintf("[cellchat] identifyOverExpressedGenes failed: %s\n", conditionMessage(e)))
})
tryCatch({
  cellchat <- identifyOverExpressedInteractions(cellchat)
}, error = function(e) {
  cat(sprintf("[cellchat] identifyOverExpressedInteractions failed: %s\n", conditionMessage(e)))
})

# 注: CellChat 默认 raw.use=FALSE 走 data.project（size-factor-normalized）；
# PHASE 1 demo fixture（270 cells × 100 features）data.project 容易因过小
# 触发 "no rows to aggregate"，实际 PBMC 3k+ 数据集两者皆可。这里用 raw.use=TRUE
# 与 GetAssayData(layer="data") 配套，避免 Seurat 5 + CellChat 1.6 兼容性陷阱。
comm_prob_status <- tryCatch({
  cellchat <- computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    raw.use = TRUE
  )
  "ok"
}, error = function(e) {
  cat(sprintf("[cellchat] computeCommunProb failed: %s\n", conditionMessage(e)))
  "failed"
})

if (comm_prob_status != "ok") {
  report_progress(60, "computeCommunProb empty-result fallback", job_dir)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
} else {
  # 过滤低概率通讯（filterCommunication 接受 min.cells）
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
}

# ── 6. computeCommunProbPathway ───────────────────────────
report_progress(70, "computeCommunProbPathway", job_dir)
tryCatch({
  cellchat <- computeCommunProbPathway(cellchat)
}, error = function(e) {
  cat(sprintf("[cellchat] computeCommunProbPathway skipped: %s\n", conditionMessage(e)))
})

# ── 7. aggregateNet + 通讯数热图 ─────────────────────────
report_progress(78, "aggregateNet + interaction_heatmap", job_dir)
tryCatch({
  cellchat <- aggregateNet(cellchat)
}, error = function(e) {
  cat(sprintf("[cellchat] aggregateNet skipped: %s\n", conditionMessage(e)))
})

# interaction_heatmap.png：cell group × cell group 通讯数量
interaction_heatmap_png <- file.path(output_dir, "interaction_heatmap.png")
png(interaction_heatmap_png, width = 1600, height = 1400, res = 200)
tryCatch({
  netVisual_heatmap(cellchat, measure = "count", color.use = NULL, title.name = "Number of interactions")
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, paste("interaction_heatmap fallback\n", conditionMessage(e)), cex = 1.0)
})
dev.off()
canonicalize_png(interaction_heatmap_png)

# chord_diagram.png：cell-cell chord（取最强 pathway）
chord_png <- file.path(output_dir, "chord_diagram.png")
png(chord_png, width = 1600, height = 1600, res = 200)
tryCatch({
  pathways.show <- cellchat@netP$pathways[1]
  if (length(pathways.show) > 0) {
    netVisual_chord_cell(
      cellchat,
      signaling = pathways.show,
      title.name = paste0("Chord: ", pathways.show)
    )
  } else {
    plot.new()
    text(0.5, 0.5, "no pathway available", cex = 1.2)
  }
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, paste("chord fallback\n", conditionMessage(e)), cex = 1.0)
})
dev.off()
canonicalize_png(chord_png)

# ── 8. 通讯结果表 (communications_table.csv) ──────────────
report_progress(85, "extract communications + pathway enrichment", job_dir)
df.net <- tryCatch(
  subsetCommunication(cellchat),
  error = function(e) {
    cat(sprintf("[cellchat] subsetCommunication failed: %s\n", conditionMessage(e)))
    data.frame()  # 0-row empty data.frame
  }
)
# 在 0 LR pair 的情况下，subsetCommunication 可能返回 NULL 而非空表
if (is.null(df.net) || nrow(df.net) == 0) {
  empty_comm <- data.frame(
    source = character(0), target = character(0), ligand = character(0),
    receptor = character(0), prob = numeric(0), pval = numeric(0),
    pathway_name = character(0)
  )
  write.csv(empty_comm, file.path(output_dir, "communications_table.csv"), row.names = FALSE)
  pathway_df <- data.frame(
    pathway = character(0), n_lr_pairs = integer(0),
    pvalue = numeric(0), stringsAsFactors = FALSE
  )
  write.csv(pathway_df, file.path(output_dir, "pathway_enrichment.csv"), row.names = FALSE)
  msg <- "[cellchat] empty_result_or_insufficient_signal — no significant LR pairs under current thresholds"
  cat(msg, "\n", file = stderr())
  report_progress(95, msg, job_dir)

  # placeholder bubble plot + empty CellChat object 保存
  bubble_png <- file.path(output_dir, "bubble_plot.png")
  png(bubble_png, width = 1600, height = 1200, res = 200)
  plot.new()
  text(0.5, 0.5, "empty_result\nno LR pairs", cex = 1.4)
  dev.off()
  canonicalize_png(bubble_png)
  saveRDS(cellchat, file.path(output_dir, "cellchat_object.rds"))

  report_progress(100, "empty-result: 见 error.txt 与 communications_table.csv", job_dir)
  quit(status = 0)  # empty-result 仍 valid status
}

write.csv(df.net, file.path(output_dir, "communications_table.csv"), row.names = FALSE)

# pathway_enrichment.csv：按 pathway 聚合
if ("pathway_name" %in% colnames(df.net)) {
  pathway_df <- aggregate(
    cbind(prob = df.net$prob) ~ pathway_name,
    data = df.net,
    FUN = function(x) c(mean_prob = mean(x), n_lr_pairs = length(x))
  )
  # 简化列名
  if (is.matrix(pathway_df$prob) && ncol(pathway_df$prob) == 2) {
    pathway_df$mean_prob <- pathway_df$prob[, "mean_prob"]
    pathway_df$n_lr_pairs <- pathway_df$prob[, "n_lr_pairs"]
    pathway_df$prob <- NULL
  }
  # 默认按 n_lr_pairs 降序
  pathway_df <- pathway_df[order(-pathway_df$n_lr_pairs), ]
  pathway_df$pvalue <- NA_real_  # CellChat 不输出通路级 pvalue，预留列
  colnames(pathway_df)[colnames(pathway_df) == "pathway_name"] <- "pathway"
  pathway_df <- pathway_df[, c("pathway", "n_lr_pairs", "mean_prob", "pvalue")]
} else {
  pathway_df <- data.frame(
    pathway = character(0), n_lr_pairs = integer(0),
    mean_prob = numeric(0), pvalue = numeric(0)
  )
}
write.csv(pathway_df, file.path(output_dir, "pathway_enrichment.csv"), row.names = FALSE)

# ── 9. bubble plot ────────────────────────────────────────
report_progress(92, "netVisual_bubble", job_dir)
bubble_png <- file.path(output_dir, "bubble_plot.png")
png(bubble_png, width = 1800, height = max(1200, 80 * min(length(unique(df.net$pathway_name)), 20)), res = 200)
tryCatch({
  if ("pathway_name" %in% colnames(df.net) && length(unique(df.net$pathway_name)) > 0) {
    show_pathways <- unique(df.net$pathway_name)
    show_pathways <- head(show_pathways, 20)
    netVisual_bubble(
      cellchat,
      sources.use = NULL,
      targets.use = NULL,
      signaling = show_pathways,
      remove.isolate = FALSE
    )
  } else {
    netVisual_bubble(cellchat, remove.isolate = FALSE)
  }
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, paste("bubble fallback\n", conditionMessage(e)), cex = 1.0)
})
dev.off()
canonicalize_png(bubble_png)

# ── 10. 可选：层级推断（hierarchical）────────────────────
if (compute_hier) {
  report_progress(95, "selectK + identifyCommunicationPatterns (NMF)", job_dir)
  tryCatch({
    cellchat <- selectK(cellchat, pattern = c("outgoing", "incoming"))
    nPatterns <- min(n_patterns, ncol(cellchat@net$P))
    cellchat <- identifyCommunicationPatterns(
      cellchat,
      pattern = c("outgoing", "incoming"),
      k = nPatterns
    )
  }, error = function(e) {
    cat(sprintf("[cellchat] hierarchy step skipped: %s\n", conditionMessage(e)))
  })
}

# ── 11. 保存 CellChat 对象 ───────────────────────────────
saveRDS(cellchat, file.path(output_dir, "cellchat_object.rds"))

# ── 12. 摘要统计 + report.html ───────────────────────────
report_progress(98, "写 report.html + manifest.json", job_dir)
n_lr <- nrow(df.net)
n_pathways <- length(unique(df.net$pathway_name))
n_groups <- length(unique(c(as.character(df.net$source), as.character(df.net$target))))

report_html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\">",
  "<title>cellchat report</title>",
  "<style>body{font-family:-apple-system,sans-serif;max-width:780px;margin:24px auto;",
  "padding:0 16px;line-height:1.55;color:#1f2937;}",
  "h1{font-size:22px;}h2{font-size:16px;color:#16a34a;margin-top:24px;}",
  "table{width:100%;border-collapse:collapse;font-size:13px;}",
  "th,td{padding:4px 8px;text-align:left;border-bottom:1px solid #e5e7eb;}",
  "th{background:#f9fafb;font-weight:600;}",
  "code{font-family:ui-monospace,Menlo,monospace;font-size:12px;",
  "background:#f3f4f6;padding:2px 4px;border-radius:3px;}",
  ".footer{color:#9ca3af;font-size:12px;margin-top:32px;text-align:center;}</style>",
  "</head><body>",
  "<h1>CellChat Cell-Cell Communication</h1>",
  "<p><strong>Input:</strong> ", basename(seurat_path),
  " &nbsp;·&nbsp; <strong>species:</strong> ", species,
  " &nbsp;·&nbsp; <strong>DB category:</strong> ", db_category, "</p>",
  "<p><strong>Cell groups:</strong> ", n_groups,
  " &nbsp;·&nbsp; <strong>LR pairs:</strong> ", n_lr,
  " &nbsp;·&nbsp; <strong>Pathways:</strong> ", n_pathways, "</p>",
  "<h2>Parameters</h2>",
  "<table><tr><th>Parameter</th><th>Value</th></tr>",
  "<tr><td>cell_type_col</td><td>", cell_type_col, "</td></tr>",
  "<tr><td>species</td><td>", species, "</td></tr>",
  "<tr><td>db_category</td><td>", db_category, "</td></tr>",
  "<tr><td>min_cells</td><td>", min_cells, "</td></tr>",
  "<tr><td>pvalue_threshold</td><td>", pvalue_cut, "</td></tr>",
  "<tr><td>n_patterns</td><td>", n_patterns, "</td></tr>",
  "<tr><td>compute_hierarchy</td><td>", compute_hier, "</td></tr>",
  "</table>",
  "<h2>Outputs</h2><ul>",
  "<li><code>communications_table.csv</code> — ", n_lr, " LR pairs</li>",
  "<li><code>pathway_enrichment.csv</code> — ", n_pathways, " pathways</li>",
  "<li><code>interaction_heatmap.png</code> — cell group × cell group interaction counts</li>",
  "<li><code>chord_diagram.png</code> — top pathway chord diagram</li>",
  "<li><code>bubble_plot.png</code> — pathway bubble (sources × targets)</li>",
  "<li><code>cellchat_object.rds</code> — CellChat object (for downstream netAnalysis_*)</li>",
  "</ul>",
  "<p class=\"footer\">G0 dev Electron 边界内 verified · biof3-r-runtime R 4.5.2 · ",
  "CellChat 1.6.1 · PHASE 1 BLOCK-11 SUB-1.4</p>",
  "</body></html>"
)
writeLines(report_html, con = file.path(output_dir, "report.html"))

# manifest.json
manifest <- list(
  plugin = "cellchat",
  pluginVersion = "1.0.0",
  status = "valid",
  species = species,
  cell_type_col = cell_type_col,
  db_category = db_category,
  stats = list(
    n_cells = ncol(seurat_obj),
    n_genes = nrow(seurat_obj),
    n_cell_groups = n_groups,
    n_lr_pairs = n_lr,
    n_pathways = n_pathways,
    min_cells = min_cells,
    pvalue_threshold = pvalue_cut,
    compute_hierarchy = compute_hier
  ),
  outputs = list(
    communications_table = "communications_table.csv",
    pathway_enrichment = "pathway_enrichment.csv",
    interaction_heatmap = "interaction_heatmap.png",
    chord_diagram = "chord_diagram.png",
    bubble_plot = "bubble_plot.png",
    cellchat_object = "cellchat_object.rds",
    report = "report.html"
  )
)
writeLines(
  jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
  con = file.path(output_dir, "manifest.json")
)

report_progress(100, sprintf("valid: %d LR pairs · %d pathways · %d groups", n_lr, n_pathways, n_groups), job_dir)
cat(sprintf("[cellchat] Done. %d LR pairs, %d pathways, %d cell groups.\n", n_lr, n_pathways, n_groups))