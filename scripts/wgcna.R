#!/usr/bin/env Rscript
# ============================================================
# Tool: wgcna
# Validated: WGCNA 1.74 / flashClust 1.1-4 / impute 1.84.0
# (R 4.5.2, Bioconductor 3.22)
# Canonical workflow source:
#   https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#   Tutorial I — Network construction and module detection
#
# Steps (in-scope) — Tutorial I 7 步:
#   1. Load + QC          → goodSamplesGenes() + sample dendrogram
#   2. Soft threshold     → pickSoftThreshold() 自动选 power
#   3. Network            → blockwiseModules(networkType="signed")
#   4. Module dendrogram  → plotDendroAndColors() 动态合并
#   5. Module-trait corr  → moduleEigengenes() + cor() + corPvalueStudent()
#   6. Hub genes          → intramodularConnectivity() top kIM
#   7. Gene-trait sig     → 最显著模块 GS vs MM 散点
#
# Out-of-scope:
#   - Tutorial III consensus network (多数据集)
#   - module preservation analysis
#   - bicor robust correlation (默认 pearson)
#   - GO 富集 (用 BioF3 GO/KEGG 工具下游)
#   - CytoScape 网络导出
#
# Default params (官方 / Tutorial I 推荐):
#   power=auto (R²≥0.85 最小)
#   networkType="signed", TOMType="signed"
#   minModuleSize=30, mergeCutHeight=0.25
#   maxBlockSize=5000 (BioF3 限定)
#   randomSeed=54321 (官方)
#
# Input (job_dir/):
#   - expression  CSV/TSV: 行=基因, 列=样本
#   - traits      CSV/TSV: 行=样本, 含数值型 traits 列
#   - params.json: {min_module_size, merge_cut_height, network_type, soft_power}
#
# Output (job_dir/output/):
#   8 plots + 6 tables + network.rds + report.html + manifest.json + summary.txt
# ============================================================

suppressMessages({
  library(WGCNA)
  library(ggplot2)
  library(jsonlite)
  library(dplyr)
})

# WGCNA settings (官方推荐)
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 2)  # keep desktop analysis responsive

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题 (W3-4 升级 2026-05-25)
source(file.path(job_dir, "_biof3-theme.R"))

cat("[wgcna] Starting WGCNA pipeline\n")
cat("[wgcna] R:", R.version.string, "\n")
cat("[wgcna] WGCNA:", as.character(packageVersion("WGCNA")), "\n")

# ---- Read params ----
params <- fromJSON(file.path(job_dir, "params.json"))
min_module_size  <- as.integer(params$min_module_size %||% 30)
merge_cut_height <- as.numeric(params$merge_cut_height %||% 0.25)
network_type     <- params$network_type %||% "signed"
user_power       <- if (!is.null(params$soft_power) && params$soft_power != "") as.integer(params$soft_power) else NULL
random_seed      <- 54321
max_block_size   <- 5000

cat(sprintf("[wgcna] params: networkType=%s, minModuleSize=%d, mergeCutHeight=%g, user_power=%s\n",
            network_type, min_module_size, merge_cut_height,
            if (is.null(user_power)) "auto" else as.character(user_power)))

# ---- Read input files (auto-detect format) ----
read_table_auto <- function(file, row_names_col = 1) {
  first_line <- readLines(file, n = 1, warn = FALSE)
  sep <- if (grepl("\t", first_line) && !grepl(",", first_line)) "\t" else ","
  read.table(file, header = TRUE, sep = sep, row.names = row_names_col,
             check.names = FALSE, stringsAsFactors = FALSE)
}

# NOTE: file inputs are copied into job_dir as "<inputId><ext>" (e.g. expression.csv)
# by both the HTTP executor and the MCP registry, so match by unanchored prefix
# (mirrors deseq2 "^counts" / go-kegg "^genes") rather than an exact "^expression$".
expr_files <- list.files(job_dir, pattern = "^expression", full.names = TRUE)
traits_files <- list.files(job_dir, pattern = "^traits", full.names = TRUE)

if (length(expr_files) == 0) stop("No expression file found in job_dir")
if (length(traits_files) == 0) stop("No traits file found in job_dir")

cat("[wgcna] Step 1: Loading data...\n")
expr_raw <- read_table_auto(expr_files[1])
traits_raw <- read_table_auto(traits_files[1])

cat(sprintf("  expression: %d genes × %d samples\n", nrow(expr_raw), ncol(expr_raw)))
cat(sprintf("  traits: %d samples × %d traits\n", nrow(traits_raw), ncol(traits_raw)))

# ---- Validate dimensions ----
if (ncol(expr_raw) > 500 || nrow(expr_raw) > 10000) {
  stop(sprintf("Dataset too large: %d genes × %d samples exceeds the 10000 × 500 limit. Pre-filter to high-variance genes.",
               nrow(expr_raw), ncol(expr_raw)))
}

# Match samples between expr and traits
common_samples <- intersect(colnames(expr_raw), rownames(traits_raw))
if (length(common_samples) < 15) {
  stop(sprintf("Only %d samples shared between expression and traits — need at least 15 for WGCNA.", length(common_samples)))
}
cat(sprintf("  matched samples: %d\n", length(common_samples)))

datExpr <- t(expr_raw[, common_samples, drop = FALSE])  # rows = samples, cols = genes
datTraits <- traits_raw[common_samples, , drop = FALSE]

# Ensure traits are numeric
non_numeric <- sapply(datTraits, function(x) !is.numeric(x))
if (any(non_numeric)) {
  bad <- names(datTraits)[non_numeric]
  stop(sprintf("Traits must be numeric. Non-numeric columns: %s. Encode categorical traits as 0/1 first.",
               paste(bad, collapse = ", ")))
}

# ---- Step 1 QC: goodSamplesGenes ----
gsg <- goodSamplesGenes(datExpr, verbose = 0)
n_bad_genes <- sum(!gsg$goodGenes)
n_bad_samples <- sum(!gsg$goodSamples)
if (n_bad_genes > 0 || n_bad_samples > 0) {
  cat(sprintf("  QC: removed %d bad genes + %d bad samples\n", n_bad_genes, n_bad_samples))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  datTraits <- datTraits[gsg$goodSamples, , drop = FALSE]
}

# Sample dendrogram (visual check for outliers)
sample_tree <- hclust(dist(datExpr), method = "average")
save_grid_biof3(file.path(output_dir, "sample_dendrogram"),
                 width = 9, height = 5, expr = {
  par(cex = 0.6, mar = c(0, 4, 2, 0))
  plot(sample_tree, main = "Sample clustering — check for outliers",
       xlab = "", sub = "", cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.2)
  abline(h = mean(sample_tree$height) + 3 * sd(sample_tree$height),
         col = biof3_palette_div(3)[3], lty = 2)
})

# Save QC summary table
qc_summary <- data.frame(
  metric = c("input_genes", "input_samples", "removed_genes_qc", "removed_samples_qc",
             "final_genes", "final_samples"),
  value = c(nrow(expr_raw), length(common_samples), n_bad_genes, n_bad_samples,
            ncol(datExpr), nrow(datExpr))
)
write.csv(qc_summary, file.path(output_dir, "qc_summary.csv"), row.names = FALSE)

# ============================================================
# Step 2: pickSoftThreshold — find the power that achieves
# scale-free topology (R² ≥ 0.85)
# ============================================================
cat("[wgcna] Step 2: pickSoftThreshold...\n")
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = network_type,
                         verbose = 0, blockSize = max_block_size)

# Save sft results
sft_results <- sft$fitIndices
write.csv(sft_results, file.path(output_dir, "sft_results.csv"), row.names = FALSE)

# Decide power
auto_power <- sft$powerEstimate
if (!is.null(user_power)) {
  chosen_power <- user_power
  power_source <- sprintf("user-specified (%d)", user_power)
} else if (!is.na(auto_power)) {
  chosen_power <- auto_power
  power_source <- sprintf("auto (R² ≥ 0.85, power=%d)", auto_power)
} else {
  # fallback: take the power with maximum R² fit
  best <- which.max(-sign(sft_results$slope) * sft_results$SFT.R.sq)
  chosen_power <- sft_results$Power[best]
  power_source <- sprintf("fallback max-R² (R²=%.2f, power=%d)",
                          sft_results$SFT.R.sq[best], chosen_power)
  cat(sprintf("  WARN: no power achieves R² ≥ 0.85, using %s\n", power_source))
}
cat(sprintf("  Chosen power: %s\n", power_source))

# Plot soft threshold (2-panel: SFT R² vs power, mean connectivity vs power)
save_grid_biof3(file.path(output_dir, "soft_threshold"),
                 width = 11, height = 4.5, expr = {
  par(mfrow = c(1, 2), bty = "l")
  cex1 <- 0.9
  red_c  <- biof3_palette_div(3)[3]
  blue_c <- biof3_palette_div(3)[1]
  plot(sft_results$Power, -sign(sft_results$slope) * sft_results$SFT.R.sq,
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R²",
       type = "n", main = "Scale independence")
  text(sft_results$Power, -sign(sft_results$slope) * sft_results$SFT.R.sq,
       labels = powers, cex = cex1, col = red_c)
  abline(h = 0.85, col = red_c, lty = 2)
  abline(v = chosen_power, col = blue_c, lty = 2)

  plot(sft_results$Power, sft_results$mean.k.,
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = "Mean connectivity")
  text(sft_results$Power, sft_results$mean.k., labels = powers, cex = cex1, col = red_c)
  abline(v = chosen_power, col = blue_c, lty = 2)
})

# ============================================================
# Step 3-4: Network construction + module detection
# ============================================================
cat("[wgcna] Step 3-4: blockwiseModules...\n")
TOM_type <- if (network_type == "signed") "signed" else "unsigned"
net <- blockwiseModules(
  datExpr,
  power            = chosen_power,
  networkType      = network_type,
  TOMType          = TOM_type,
  minModuleSize    = min_module_size,
  reassignThreshold = 0,
  mergeCutHeight   = merge_cut_height,
  numericLabels    = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs         = FALSE,
  randomSeed       = random_seed,
  maxBlockSize     = max_block_size,
  verbose          = 0
)

# Convert numeric labels to colors
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
n_modules <- length(unique(moduleColors)) - as.integer("grey" %in% moduleColors)
cat(sprintf("  Found %d modules (excluding grey)\n", n_modules))

if (n_modules < 2) {
  stop(sprintf("Only %d non-grey module detected. Try lowering min_module_size or merge_cut_height.", n_modules))
}

# Save module assignment
module_assignment <- data.frame(
  gene = colnames(datExpr),
  module_label = moduleLabels,
  module_color = moduleColors
)
write.csv(module_assignment, file.path(output_dir, "module_assignment.csv"), row.names = FALSE)

# Plot gene dendrogram + module colors
save_grid_biof3(file.path(output_dir, "gene_dendrogram"),
                 width = 9, height = 5, expr = {
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
})

# Module size barplot
mod_sizes <- as.data.frame(table(moduleColors))
mod_sizes <- mod_sizes[order(-mod_sizes$Freq), ]
mod_sizes$moduleColors <- factor(mod_sizes$moduleColors, levels = mod_sizes$moduleColors)

p_mod_size <- ggplot(mod_sizes, aes(x = moduleColors, y = Freq, fill = moduleColors)) +
  geom_col() +
  scale_fill_identity() +
  labs(x = "Module", y = "Number of genes",
       title = sprintf("Module sizes (%d modules + grey)", n_modules)) +
  theme_biof3() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave_biof3(p_mod_size, file.path(output_dir, "module_size_barplot"),
              width = max(8, nrow(mod_sizes) * 0.4), height = 5)

# ============================================================
# Step 5: Module-trait correlation
# ============================================================
cat("[wgcna] Step 5: module-trait correlation...\n")
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Save MEs
MEs_out <- data.frame(sample = rownames(MEs), MEs, check.names = FALSE)
write.csv(MEs_out, file.path(output_dir, "module_eigengenes.csv"), row.names = FALSE)

# Correlation
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Save correlation table
mtc_df <- data.frame(module = rownames(moduleTraitCor), moduleTraitCor, check.names = FALSE)
mtp_df <- data.frame(module = rownames(moduleTraitPvalue), moduleTraitPvalue, check.names = FALSE)
write.csv(mtc_df, file.path(output_dir, "module_trait_corr.csv"), row.names = FALSE)
write.csv(mtp_df, file.path(output_dir, "module_trait_pvalue.csv"), row.names = FALSE)

# Heatmap with cell labels (r value above, p value below in parens)
text_matrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
dim(text_matrix) <- dim(moduleTraitCor)

save_grid_biof3(file.path(output_dir, "module_trait_heatmap"),
                 width = max(7, ncol(datTraits) * 0.7 + 2),
                 height = max(5, nrow(MEs) * 0.25 + 2), expr = {
  par(mar = c(6, 9, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(datTraits),
                 yLabels = colnames(MEs),
                 ySymbols = colnames(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = text_matrix,
                 setStdMargins = FALSE,
                 cex.text = 0.55,
                 zlim = c(-1, 1),
                 main = "Module-Trait correlation (r above, p below)")
})

# Module eigengene barplot (one panel per top-correlated module + most-related trait)
# Use the top 4 modules by max |r| across traits
mod_max_r <- apply(abs(moduleTraitCor), 1, max)
top_mods <- names(sort(mod_max_r, decreasing = TRUE))[1:min(4, length(mod_max_r))]
top_mods <- top_mods[top_mods != "MEgrey"]

if (length(top_mods) > 0) {
  me_long <- do.call(rbind, lapply(top_mods, function(m) {
    data.frame(sample = rownames(MEs), ME = MEs[[m]], module = m)
  }))
  me_long$module <- factor(me_long$module, levels = top_mods)
  p_me <- ggplot(me_long, aes(x = sample, y = ME, fill = module)) +
    geom_col() +
    facet_wrap(~ module, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = setNames(gsub("ME", "", top_mods), top_mods)) +
    labs(x = "Sample", y = "Module Eigengene", title = "Top 4 Module Eigengenes per sample") +
    theme_biof3() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none")
  ggsave_biof3(p_me, file.path(output_dir, "module_eigengenes"),
                width = 12, height = 2 * length(top_mods) + 1)
}

# ============================================================
# Step 6: Hub genes (top intramodular connectivity per module)
# ============================================================
cat("[wgcna] Step 6: hub genes...\n")
adj <- adjacency(datExpr, power = chosen_power, type = network_type)
kIM <- intramodularConnectivity(adj, moduleColors)

hub_df <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors,
  kTotal = kIM$kTotal,
  kWithin = kIM$kWithin,
  kOut = kIM$kOut
)

# Top 10 per module by kWithin (exclude grey)
hub_top <- hub_df %>%
  filter(module != "grey") %>%
  group_by(module) %>%
  slice_max(order_by = kWithin, n = 10) %>%
  ungroup()
write.csv(hub_top, file.path(output_dir, "hub_genes.csv"), row.names = FALSE)

# Plot top hub genes for the top modules (2x2 grid)
hub_plot_modules <- gsub("ME", "", top_mods)
hub_subset <- hub_top %>% filter(module %in% hub_plot_modules)

if (nrow(hub_subset) > 0) {
  p_hub <- ggplot(hub_subset, aes(x = reorder(gene, kWithin), y = kWithin, fill = module)) +
    geom_col() +
    facet_wrap(~ module, scales = "free", ncol = 2) +
    coord_flip() +
    scale_fill_manual(values = setNames(hub_plot_modules, hub_plot_modules)) +
    labs(x = NULL, y = "Intramodular connectivity (kWithin)",
         title = "Top 10 hub genes per top-correlated module") +
    theme_biof3() +
    theme(legend.position = "none")
  ggsave_biof3(p_hub, file.path(output_dir, "hub_gene_table"),
                width = 11, height = 7)
}

# ============================================================
# Step 7: Gene-trait significance (GS) vs Module Membership (MM)
# Use the trait with maximum |r| against the most-correlated module
# ============================================================
cat("[wgcna] Step 7: GS vs MM scatter...\n")

# Find the highest |r| cell in moduleTraitCor (module × trait)
mtc_abs <- abs(moduleTraitCor)
max_idx <- which(mtc_abs == max(mtc_abs[!is.na(mtc_abs)]), arr.ind = TRUE)[1, ]
top_module <- rownames(moduleTraitCor)[max_idx[1]]
top_trait <- colnames(moduleTraitCor)[max_idx[2]]
cat(sprintf("  Top association: %s ~ %s (r=%.3f, p=%.2g)\n",
            top_module, top_trait,
            moduleTraitCor[max_idx[1], max_idx[2]],
            moduleTraitPvalue[max_idx[1], max_idx[2]]))

# MM = correlation of each gene with its module's ME
# GS = correlation of each gene with the trait
trait_vec <- datTraits[[top_trait]]
GS <- as.numeric(cor(datExpr, trait_vec, use = "p"))
GS_p <- corPvalueStudent(GS, nrow(datExpr))

top_module_color <- gsub("ME", "", top_module)
genes_in_top <- moduleColors == top_module_color
MM_top <- as.numeric(cor(datExpr[, genes_in_top, drop = FALSE], MEs[[top_module]], use = "p"))

gs_mm_data <- data.frame(
  gene = colnames(datExpr)[genes_in_top],
  MM = MM_top,
  GS = GS[genes_in_top],
  GS_pvalue = GS_p[genes_in_top]
)
write.csv(gs_mm_data, file.path(output_dir, "gs_mm_data.csv"), row.names = FALSE)

p_scatter <- ggplot(gs_mm_data, aes(x = abs(MM), y = abs(GS))) +
  geom_point(color = top_module_color, alpha = 0.65, size = 1.4) +
  geom_smooth(method = "lm", color = "gray40", se = FALSE, linewidth = 0.6) +
  labs(x = sprintf("Module Membership in '%s' module", top_module_color),
       y = sprintf("Gene Significance for %s", top_trait),
       title = sprintf("GS vs MM: %s module ~ %s\n(cor=%.2f, n=%d)",
                       top_module_color, top_trait,
                       cor(abs(MM_top), abs(GS[genes_in_top]), use = "p"),
                       sum(genes_in_top))) +
  theme_biof3()
ggsave_biof3(p_scatter, file.path(output_dir, "gs_mm_scatter"),
              width = 7, height = 6)

# ============================================================
# Save full network as RDS
# ============================================================
saveRDS(list(net = net, moduleColors = moduleColors, MEs = MEs,
             datTraits = datTraits, chosen_power = chosen_power,
             top_module = top_module, top_trait = top_trait),
        file.path(output_dir, "network.rds"))

# ============================================================
# Summary
# ============================================================
sft_r2 <- sft_results$SFT.R.sq[sft_results$Power == chosen_power]
top_r <- moduleTraitCor[max_idx[1], max_idx[2]]
top_p <- moduleTraitPvalue[max_idx[1], max_idx[2]]

summary_text <- sprintf(
  "WGCNA Pipeline Summary\n\nDataset:\n  Input: %d genes × %d samples\n  After QC: %d genes × %d samples\n  Traits: %d numeric\n\nNetwork:\n  Type: %s\n  Soft power: %d (%s, R²=%.3f)\n  TOM type: %s\n\nModules:\n  Total: %d (+ grey)\n  Min size: %d, merge height: %.2f\n\nTop Module-Trait association:\n  %s ~ %s\n  r = %.3f, p = %.2g\n\nHub genes:\n  Top 10 per module saved to hub_genes.csv\n",
  nrow(expr_raw), length(common_samples),
  ncol(datExpr), nrow(datExpr),
  ncol(datTraits),
  network_type, chosen_power, power_source, sft_r2,
  TOM_type,
  n_modules, min_module_size, merge_cut_height,
  top_module, top_trait, top_r, top_p
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))
cat(summary_text)

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = c(list(
    list(name = "module_assignment.csv",   type = "table", label = "基因 → 模块"),
    list(name = "module_trait_corr.csv",   type = "table", label = "模块-性状 r"),
    list(name = "module_trait_pvalue.csv", type = "table", label = "模块-性状 p"),
    list(name = "module_eigengenes.csv",   type = "table", label = "模块 ME 矩阵"),
    list(name = "hub_genes.csv",           type = "table", label = "Hub 基因"),
    list(name = "sft_results.csv",         type = "table", label = "pickSoftThreshold"),
    list(name = "gs_mm_data.csv",          type = "table", label = "GS / MM 数据"),
    list(name = "qc_summary.csv",          type = "table", label = "QC 概况"),
    list(name = "network.rds",             type = "file",  label = "完整网络 RDS"),
    list(name = "summary.txt",             type = "text",  label = "分析摘要"),
    list(name = "report.html",             type = "file",  label = "解读报告 HTML")
  ), unlist(lapply(list(
    c("sample_dendrogram",   "样本聚类树"),
    c("soft_threshold",      "软阈值选择"),
    c("gene_dendrogram",     "基因聚类 + 模块"),
    c("module_size_barplot", "模块大小"),
    c("module_trait_heatmap","模块-性状热图"),
    c("module_eigengenes",   "Top 模块 Eigengene"),
    c("hub_gene_table",      "Hub 基因图"),
    c("gs_mm_scatter",       "GS vs MM 散点")
  ), function(pp) {
    out <- list()
    png_path <- file.path(output_dir, paste0(pp[1], ".png"))
    pdf_path <- file.path(output_dir, paste0(pp[1], ".pdf"))
    if (file.exists(png_path)) out <- c(out, list(list(name = paste0(pp[1], ".png"), type = "plot", label = paste0(pp[2], " (PNG)"))))
    if (file.exists(pdf_path)) out <- c(out, list(list(name = paste0(pp[1], ".pdf"), type = "file", label = paste0(pp[2], " (PDF)"))))
    out
  }), recursive = FALSE)),
  summary = list(
    n_genes_input = nrow(expr_raw),
    n_genes_final = ncol(datExpr),
    n_samples_final = nrow(datExpr),
    n_traits = ncol(datTraits),
    chosen_power = chosen_power,
    sft_r_squared = sft_r2,
    n_modules = n_modules,
    top_module = top_module,
    top_trait = top_trait,
    top_correlation = signif(top_r, 3),
    top_pvalue = signif(top_p, 3),
    sci_style = TRUE,
    dual_format = TRUE
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ============================================================
# Generate HTML report
# ============================================================
cat("[wgcna] Generating HTML report...\n")
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  png_files <- c("sample_dendrogram.png", "soft_threshold.png", "gene_dendrogram.png",
                 "module_size_barplot.png", "module_trait_heatmap.png",
                 "module_eigengenes.png", "hub_gene_table.png", "gs_mm_scatter.png")
  for (png in png_files) {
    img_path <- file.path(output_dir, png)
    placeholder <- paste0("{{", png, "}}")
    if (file.exists(img_path)) {
      b64 <- base64enc::base64encode(img_path)
      data_uri <- paste0("data:image/png;base64,", b64)
      report_html <- gsub(placeholder, data_uri, report_html, fixed = TRUE)
    }
  }

  # Replace {{variable}} stat placeholders
  report_html <- gsub("{{n_genes_input}}",   as.character(nrow(expr_raw)),         report_html, fixed = TRUE)
  report_html <- gsub("{{n_genes_final}}",   as.character(ncol(datExpr)),          report_html, fixed = TRUE)
  report_html <- gsub("{{n_samples}}",       as.character(nrow(datExpr)),          report_html, fixed = TRUE)
  report_html <- gsub("{{n_traits}}",        as.character(ncol(datTraits)),        report_html, fixed = TRUE)
  report_html <- gsub("{{chosen_power}}",    as.character(chosen_power),           report_html, fixed = TRUE)
  report_html <- gsub("{{sft_r_squared}}",   sprintf("%.3f", sft_r2),              report_html, fixed = TRUE)
  report_html <- gsub("{{network_type}}",    network_type,                          report_html, fixed = TRUE)
  report_html <- gsub("{{n_modules}}",       as.character(n_modules),              report_html, fixed = TRUE)
  report_html <- gsub("{{top_module}}",      gsub("ME", "", top_module),           report_html, fixed = TRUE)
  report_html <- gsub("{{top_trait}}",       top_trait,                             report_html, fixed = TRUE)
  report_html <- gsub("{{top_correlation}}", sprintf("%.3f", top_r),               report_html, fixed = TRUE)
  report_html <- gsub("{{top_pvalue}}",      sprintf("%.2g", top_p),               report_html, fixed = TRUE)

  # Remove fig blocks for any unreplaced {{*.png}}
  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)

  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[wgcna] Report generated.\n")
} else {
  cat("[wgcna] WARN: report template not found.\n")
}

cat("[wgcna] Pipeline complete.\n")
