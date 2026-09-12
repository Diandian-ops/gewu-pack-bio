#!/usr/bin/env Rscript
# ============================================================
# Tool: seurat-standard
# Version targeted: Seurat 5.4.0 (R 4.5.1)
# Canonical workflow source:
#   https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#
# Steps (in-scope):
#   1. Load → 2. QC + filter → 3. Normalize → 4. HVGs → 5. Scale →
#   6. PCA + ElbowPlot → 7. Cluster (FindNeighbors+FindClusters) →
#   8. UMAP → 9. FindAllMarkers → DoHeatmap + VlnPlot
#
# Out-of-scope (留给后续工具):
#   - Integration / Harmony / SCTransform
#   - WNN / spatial / multi-modal
#   - Pseudotime / RNA velocity
#   - Auto cell-type annotation (SingleR etc.)
#   - Doublet detection
#
# Default params (from official PBMC 3k tutorial, no deviation):
#   min.cells=3 / min.features=200 / nFeature_max=2500 / percent_mt_max=5
#   normalization=LogNormalize, scale.factor=10000
#   selection.method=vst, nfeatures=2000
#   dims=1:10, resolution=0.5, only.pos=TRUE
#
# Inputs (job_dir/):
#   - matrix (.csv/.tsv/.h5/.mtx) or 10X三件套 (matrix.mtx[.gz] + barcodes.tsv[.gz] + features.tsv[.gz])
#   - params.json: {species, min_features, max_features, max_mt, nfeatures, dims, resolution, top_n_markers}
#
# Outputs (job_dir/output/):
#   Plots (10): qc_violin, qc_scatter, hvg_plot, pca_dimplot, pca_dim_loadings,
#               elbow_plot, umap_clusters, markers_heatmap, markers_violin, feature_plot
#   Tables (4): qc_metrics, hvg_table, umap_coords, markers_table
#   Intermediate: pca_stdev.csv, pca_loadings.csv, seurat_obj.rds
#   Report: report.html, manifest.json, summary.txt
# ============================================================

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(jsonlite)
  library(Matrix)
})

# %||% helper (some R versions have it, some don't)
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ★ Load BioF3 SCI 主题 (W3-3 升级 2026-05-25)
source(file.path(job_dir, "_biof3-theme.R"))

cat("[seurat] Starting seurat-standard pipeline\n")
cat("[seurat] R:", R.version.string, "\n")
cat("[seurat] Seurat:", as.character(packageVersion("Seurat")), "\n")

# ---- Read params ----
params <- fromJSON(file.path(job_dir, "params.json"))
species        <- params$species %||% "human"
min_features   <- as.integer(params$min_features %||% 200)
max_features   <- as.integer(params$max_features %||% 2500)
max_mt         <- as.numeric(params$max_mt %||% 5)
nfeatures      <- as.integer(params$nfeatures %||% 2000)
dims_val       <- as.integer(params$dims %||% 10)
resolution     <- as.numeric(params$resolution %||% 0.5)
top_n_markers  <- as.integer(params$top_n_markers %||% 10)

mt_pattern <- if (species == "mouse") "^mt-" else "^MT-"
dims_use   <- 1:dims_val

cat(sprintf("[seurat] params: species=%s, min/max_features=%d/%d, max_mt=%g, nfeatures=%d, dims=1:%d, resolution=%g, top_n=%d\n",
            species, min_features, max_features, max_mt, nfeatures, dims_val, resolution, top_n_markers))

# ============================================================
# Step 1. Load data — auto-detect format
# ============================================================
cat("[seurat] Step 1: Loading data...\n")

# Look for input file(s)
input_files <- list.files(job_dir, pattern = "^matrix", full.names = TRUE)
input_files <- input_files[!grepl("\\.json$", input_files)]

# Sanity: detect 10X-style triple-file (matrix.mtx + barcodes.tsv + features.tsv)
all_files <- list.files(job_dir, full.names = TRUE)
has_mtx       <- any(grepl("matrix\\.mtx(\\.gz)?$", all_files))
has_barcodes  <- any(grepl("barcodes\\.tsv(\\.gz)?$", all_files))
has_features  <- any(grepl("(features|genes)\\.tsv(\\.gz)?$", all_files))

counts <- NULL

if (has_mtx && has_barcodes && has_features) {
  # 10X三件套 — Read10X expects a directory containing the three files
  cat("[seurat]   Detected 10X三件套 (.mtx + barcodes + features)\n")
  counts <- Read10X(data.dir = job_dir)
} else {
  input_file <- input_files[1]
  if (is.null(input_file) || is.na(input_file)) {
    stop("No input matrix file found in job_dir")
  }
  ext <- tolower(tools::file_ext(input_file))
  cat(sprintf("[seurat]   Reading %s (ext=%s)\n", basename(input_file), ext))

  # If no extension (file was uploaded via API as just "matrix"), peek at content
  if (ext == "") {
    # Read first line to detect format
    first_line <- readLines(input_file, n = 1, warn = FALSE)
    if (grepl("^%%MatrixMarket", first_line)) {
      ext <- "mtx"
    } else if (grepl(",", first_line) || grepl("\t", first_line)) {
      # Looks like a delimited table — use tab if tab-separated, else comma
      ext <- if (grepl("\t", first_line) && !grepl(",", first_line)) "tsv" else "csv"
    } else {
      stop(sprintf("Cannot auto-detect file format. First line: %.80s", first_line))
    }
    cat(sprintf("[seurat]   Auto-detected format: %s\n", ext))
  }

  if (ext == "h5") {
    counts <- Read10X_h5(input_file)
    if (is.list(counts)) {
      # CITE-seq style — take Gene Expression assay only (per out-of-scope)
      cat("[seurat]   Multi-modal H5 detected, using only Gene Expression assay\n")
      counts <- counts[["Gene Expression"]]
    }
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- if (ext == "csv") "," else "\t"
    counts <- read.table(input_file, header = TRUE, sep = sep, row.names = 1, check.names = FALSE)
    counts <- as(as.matrix(counts), "CsparseMatrix")
  } else if (ext == "mtx") {
    counts <- Matrix::readMM(input_file)
  } else {
    stop(sprintf("Unsupported file extension: %s", ext))
  }
}

cat(sprintf("[seurat]   Loaded matrix: %d genes × %d cells\n", nrow(counts), ncol(counts)))

# Sanity check on dataset size to avoid exhausting desktop memory.
if (ncol(counts) > 50000) {
  stop(sprintf("Dataset too large: %d cells exceeds the 50,000 cell limit for this online tool. Please filter or downsample first.", ncol(counts)))
}

# ============================================================
# Step 2. CreateSeuratObject + QC metrics
# ============================================================
cat("[seurat] Step 2: Creating Seurat object + computing QC metrics...\n")

obj <- CreateSeuratObject(counts = counts, project = "biof3",
                           min.cells = 3, min.features = min_features)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

# QC violin (3-panel) — patchwork composition, use & for theme broadcast
p_qc_violin <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       ncol = 3, pt.size = 0.1) &
  theme_biof3(base_size = 10) &
  theme(plot.title = element_text(size = 11), legend.position = "none")
ggsave_biof3(p_qc_violin, file.path(output_dir, "qc_violin"),
              width = 11, height = 4)

# QC scatter (patchwork)
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  NoLegend() + ggtitle("Count vs % MT") + theme_biof3(base_size = 10)
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  NoLegend() + ggtitle("Count vs nFeature") + theme_biof3(base_size = 10)
p_qc_scatter <- plot1 + plot2
ggsave_biof3(p_qc_scatter, file.path(output_dir, "qc_scatter"),
              width = 10, height = 4)

# Save QC metrics table (for reproducible re-plotting)
qc_metrics <- data.frame(
  cell        = rownames(obj@meta.data),
  nFeature_RNA = obj$nFeature_RNA,
  nCount_RNA   = obj$nCount_RNA,
  percent.mt   = obj$percent.mt
)
write.csv(qc_metrics, file.path(output_dir, "qc_metrics.csv"), row.names = FALSE)

n_before_qc <- ncol(obj)

# ============================================================
# Step 3. Filter cells (apply QC thresholds)
# ============================================================
cat(sprintf("[seurat] Step 3: Filtering cells (nFeature: %d-%d, percent.mt < %g%%)...\n",
            min_features, max_features, max_mt))

obj <- subset(obj, subset = nFeature_RNA > min_features &
                            nFeature_RNA < max_features &
                            percent.mt < max_mt)
n_after_qc <- ncol(obj)
cat(sprintf("[seurat]   Cells: %d → %d after QC\n", n_before_qc, n_after_qc))

if (n_after_qc < 50) {
  stop(sprintf("Too few cells (%d) remaining after QC. Try relaxing thresholds.", n_after_qc))
}

# ============================================================
# Step 4. Normalize + Step 5. Find variable features
# ============================================================
cat("[seurat] Step 4-5: Normalize + FindVariableFeatures...\n")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)

# HVG plot — labeled top 10
top10 <- head(VariableFeatures(obj), 10)
p_hvg_base <- VariableFeaturePlot(obj)
p_hvg <- LabelPoints(plot = p_hvg_base, points = top10, repel = TRUE) +
  theme_biof3(base_size = 10) +
  theme(legend.position = "bottom")
ggsave_biof3(p_hvg, file.path(output_dir, "hvg_plot"), width = 8, height = 6)

# Save HVG table
hvg_df <- HVFInfo(obj)
hvg_df$gene <- rownames(hvg_df)
hvg_df$is_top_variable <- hvg_df$gene %in% VariableFeatures(obj)
hvg_df <- hvg_df[order(-hvg_df$variance.standardized), ]
write.csv(hvg_df, file.path(output_dir, "hvg_table.csv"), row.names = FALSE)

# ============================================================
# Step 6. Scale data + Step 7. RunPCA
# ============================================================
cat("[seurat] Step 6: ScaleData (all genes)...\n")
all_genes <- rownames(obj)
obj <- ScaleData(obj, features = all_genes, verbose = FALSE)

cat("[seurat] Step 7: RunPCA + ElbowPlot...\n")
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)

# PCA DimPlot
p_pca <- DimPlot(obj, reduction = "pca") + NoLegend() +
  ggtitle("PCA") + theme_biof3(base_size = 10)
ggsave_biof3(p_pca, file.path(output_dir, "pca_dimplot"), width = 6, height = 5)

# PCA top loadings
p_loadings <- VizDimLoadings(obj, dims = 1:2, reduction = "pca", balanced = TRUE) &
  theme_biof3(base_size = 9)
ggsave_biof3(p_loadings, file.path(output_dir, "pca_dim_loadings"),
              width = 8, height = 6)

# Elbow plot
p_elbow <- ElbowPlot(obj, ndims = max(20, dims_val + 5)) +
  theme_biof3(base_size = 10)
ggsave_biof3(p_elbow, file.path(output_dir, "elbow_plot"), width = 6, height = 4)

# Save PCA stdev for re-plotting
pca_stdev <- data.frame(
  PC = paste0("PC_", 1:length(Stdev(obj, reduction = "pca"))),
  stdev = Stdev(obj, reduction = "pca")
)
write.csv(pca_stdev, file.path(output_dir, "pca_stdev.csv"), row.names = FALSE)

# Save PC1-2 loadings
loadings_mat <- Loadings(obj[["pca"]])[, 1:2, drop = FALSE]
loadings_df <- data.frame(
  gene  = rownames(loadings_mat),
  PC_1  = loadings_mat[, 1],
  PC_2  = loadings_mat[, 2]
)
loadings_df <- loadings_df[order(-abs(loadings_df$PC_1)), ]
write.csv(head(loadings_df, 50), file.path(output_dir, "pca_loadings.csv"), row.names = FALSE)

# ============================================================
# Step 8. FindNeighbors + FindClusters + Step 9. RunUMAP
# ============================================================
cat(sprintf("[seurat] Step 8-9: Cluster (dims=1:%d, res=%g) + UMAP...\n", dims_val, resolution))
obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)

n_clusters <- length(levels(Idents(obj)))
cat(sprintf("[seurat]   Found %d clusters\n", n_clusters))

# UMAP colored by cluster
cluster_pal <- biof3_palette(n_clusters)
p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.5,
                   cols = cluster_pal) +
  NoLegend() + ggtitle(sprintf("UMAP — %d clusters", n_clusters)) +
  theme_biof3(base_size = 10)
ggsave_biof3(p_umap, file.path(output_dir, "umap_clusters"), width = 7, height = 6)

# Save UMAP coords + cluster
umap_coords <- data.frame(
  cell    = colnames(obj),
  UMAP_1  = Embeddings(obj, "umap")[, 1],
  UMAP_2  = Embeddings(obj, "umap")[, 2],
  cluster = as.character(Idents(obj))
)
write.csv(umap_coords, file.path(output_dir, "umap_coords.csv"), row.names = FALSE)

# ============================================================
# Step 10. FindAllMarkers + visualize top markers
# ============================================================
cat("[seurat] Step 10: FindAllMarkers...\n")
markers <- FindAllMarkers(obj, only.pos = TRUE, verbose = FALSE)

# Save full markers table
write.csv(markers, file.path(output_dir, "markers_table.csv"), row.names = FALSE)

# Top N markers per cluster (by avg_log2FC)
top_markers <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = top_n_markers) %>%
  dplyr::ungroup()

# Markers heatmap (top N per cluster, on scaled data)
# DoHeatmap 内部布局复杂, 不套 theme_biof3 (会破坏 cluster bar 排列), 但走双格式输出
top_marker_genes <- unique(top_markers$gene)
if (length(top_marker_genes) > 1) {
  p_heat <- DoHeatmap(obj, features = top_marker_genes,
                       size = 3, raster = TRUE) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 9))
  # Heatmap height grows with number of marker genes
  heat_h <- max(8, min(20, length(top_marker_genes) * 0.12))
  ggsave_biof3(p_heat, file.path(output_dir, "markers_heatmap"),
                width = 11, height = heat_h, limitsize = FALSE)
}

# Top 6 markers (one per top cluster) for VlnPlot + FeaturePlot
top6 <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 1) %>%
  dplyr::ungroup() %>%
  head(6) %>%
  dplyr::pull(gene)
top6 <- unique(top6)

if (length(top6) >= 1) {
  # VlnPlot 是 patchwork, 用 & 广播 theme
  p_vln <- VlnPlot(obj, features = top6, ncol = 3, pt.size = 0) &
    theme_biof3(base_size = 9) &
    theme(legend.position = "none")
  ggsave_biof3(p_vln, file.path(output_dir, "markers_violin"),
                width = 11, height = 7)

  # FeaturePlot 是 patchwork, 同样用 &
  p_feat <- FeaturePlot(obj, features = top6, ncol = 3,
                         reduction = "umap", order = TRUE) &
    theme_biof3(base_size = 9)
  ggsave_biof3(p_feat, file.path(output_dir, "feature_plot"),
                width = 11, height = 7)
}

# ============================================================
# Save the full Seurat object (RDS) for reproducible online plotting
# ============================================================
cat("[seurat] Saving Seurat object (.rds)...\n")
saveRDS(obj, file.path(output_dir, "seurat_obj.rds"))

# ============================================================
# Summary
# ============================================================
summary_text <- sprintf(
  "Seurat Standard Pipeline Summary\n\nDataset:\n  Input cells: %d (after gene/feature pre-filter)\n  After QC filter: %d cells\n    nFeature_RNA: %d - %d, percent.mt < %g%%\n  Genes: %d\n\nWorkflow:\n  Normalization: LogNormalize (scale.factor=10000)\n  Variable features: %d (vst)\n  PCA dims used: 1-%d\n  Cluster resolution: %g\n  Clusters found: %d\n\nMarkers:\n  Total markers (only.pos=TRUE): %d\n  Top %d per cluster shown in heatmap\n",
  n_before_qc, n_after_qc,
  min_features, max_features, max_mt,
  nrow(obj),
  nfeatures, dims_val, resolution, n_clusters,
  nrow(markers), top_n_markers
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))
cat(summary_text)

# ============================================================
# Manifest
# ============================================================
manifest <- list(
  files = c(list(
    list(name = "qc_metrics.csv",      type = "table", label = "QC 指标表"),
    list(name = "hvg_table.csv",       type = "table", label = "高变基因表"),
    list(name = "umap_coords.csv",     type = "table", label = "UMAP 坐标 + cluster"),
    list(name = "markers_table.csv",   type = "table", label = "Marker 基因表"),
    list(name = "pca_stdev.csv",       type = "table", label = "PCA stdev"),
    list(name = "pca_loadings.csv",    type = "table", label = "PCA top loadings"),
    list(name = "seurat_obj.rds",      type = "file",  label = "Seurat 对象 (.rds)"),
    list(name = "summary.txt",         type = "text",  label = "分析摘要"),
    list(name = "report.html",         type = "file",  label = "解读报告 (HTML)")
  ), unlist(lapply(list(
    c("qc_violin",        "QC 三联小提琴"),
    c("qc_scatter",       "QC scatter"),
    c("hvg_plot",         "高变基因"),
    c("pca_dimplot",      "PCA DimPlot"),
    c("pca_dim_loadings", "PCA Top loadings"),
    c("elbow_plot",       "Elbow Plot"),
    c("umap_clusters",    "UMAP 聚类"),
    c("markers_heatmap",  "Marker 热图"),
    c("markers_violin",   "Top markers 小提琴"),
    c("feature_plot",     "Top markers FeaturePlot")
  ), function(pp) {
    out <- list()
    png_path <- file.path(output_dir, paste0(pp[1], ".png"))
    pdf_path <- file.path(output_dir, paste0(pp[1], ".pdf"))
    if (file.exists(png_path)) out <- c(out, list(list(name = paste0(pp[1], ".png"), type = "plot", label = paste0(pp[2], " (PNG)"))))
    if (file.exists(pdf_path)) out <- c(out, list(list(name = paste0(pp[1], ".pdf"), type = "file", label = paste0(pp[2], " (PDF)"))))
    out
  }), recursive = FALSE)),
  summary = list(
    cells_input    = n_before_qc,
    cells_after_qc = n_after_qc,
    genes          = nrow(obj),
    n_variable     = nfeatures,
    dims_used      = dims_val,
    resolution     = resolution,
    n_clusters     = n_clusters,
    n_markers      = nrow(markers),
    sci_style      = TRUE,
    dual_format    = TRUE
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ============================================================
# Generate HTML interpretation report
# ============================================================
cat("[seurat] Generating HTML report...\n")
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  # Replace {{filename.png}} placeholders with base64 data URIs
  png_files <- c("qc_violin.png", "qc_scatter.png", "hvg_plot.png",
                 "pca_dimplot.png", "pca_dim_loadings.png", "elbow_plot.png",
                 "umap_clusters.png", "markers_heatmap.png",
                 "markers_violin.png", "feature_plot.png")
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
  report_html <- gsub("{{cells_input}}",    as.character(n_before_qc),  report_html, fixed = TRUE)
  report_html <- gsub("{{cells_after_qc}}", as.character(n_after_qc),   report_html, fixed = TRUE)
  report_html <- gsub("{{n_genes}}",        as.character(nrow(obj)),    report_html, fixed = TRUE)
  report_html <- gsub("{{n_variable}}",     as.character(nfeatures),    report_html, fixed = TRUE)
  report_html <- gsub("{{dims_used}}",      as.character(dims_val),     report_html, fixed = TRUE)
  report_html <- gsub("{{resolution}}",     as.character(resolution),   report_html, fixed = TRUE)
  report_html <- gsub("{{n_clusters}}",     as.character(n_clusters),   report_html, fixed = TRUE)
  report_html <- gsub("{{n_markers}}",      as.character(nrow(markers)),report_html, fixed = TRUE)
  report_html <- gsub("{{species}}",        species,                    report_html, fixed = TRUE)

  # Remove fig blocks for any unreplaced {{*.png}} (shouldn't happen but safety net)
  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)

  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[seurat] Report generated.\n")
} else {
  cat("[seurat] Warning: report template not found, skipping report.html.\n")
}

cat("[seurat] Pipeline complete.\n")
