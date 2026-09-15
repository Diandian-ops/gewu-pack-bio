#!/usr/bin/env Rscript
# ============================================================
# Tool: seurat-standard (plugin 八件套版)
# Version: 1.0.0 — Seurat 5.4.0 (R 4.5.x)
# Canonical workflow source:
#   https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#
# 步骤（与既有 core/r-server/tool-scripts/seurat-standard.R 同源；本文件是
# plugin 八件套版，作为 builtin-canonical plugin 暴露给 BioF3 调度器，
# 原 core/r-server/tool-scripts/seurat-standard.R 保留向后兼容）：
#   1. Load（10X 三件套 / CSV / TSV / MTX / H5 自动探测）→ CreateSeuratObject
#   2. QC + filter（percent.mt + nFeature_RNA）
#   3. NormalizeData (LogNormalize, scale.factor=10000)
#   4. FindVariableFeatures (vst, nfeatures=N)
#   5. ScaleData (all.genes) + RunPCA + ElbowPlot
#   6. FindNeighbors (dims=1:N) + FindClusters (resolution=R)
#   7. RunUMAP + DimPlot
#   8. FindAllMarkers (only.pos=TRUE) → DoHeatmap + VlnPlot + FeaturePlot
#
# Out-of-scope:
#   - Integration / Harmony / SCTransform
#   - WNN / spatial / multi-modal
#   - Pseudotime / RNA velocity
#   - Auto cell-type annotation (SingleR 等)
#   - Doublet detection
#
# 契约（BioF3 调度器）:
#   - commandArgs[1] = job_dir
#   - 读 params.json（jsonlite）
#   - source _biof3-theme.R（tools.js 自动复制到 jobDir）
#   - 输出全部写入 job_dir/output/
#   - 8 件套输出（与 tool-definition.json 对齐）：
#     seurat_obj.rds / markers_table.csv / qc_metrics.csv / umap_clusters.png
#     / umap_coords.csv / hvg_plot.png / report.html / manifest.json
#   - PNG 通过 png::writePNG(png::readPNG(f), f) 重写以稳定 sha256
#
# Exit codes:
#   0 = valid（含 empty-result 注记仍为 valid）
#   3 = missing_required_file_inputs
#   4 = runtime_not_ready
#   5 = empty_result（同样 valid, 仅供 dispatcher 记录 failureMode）
# ============================================================

# Preflight: 用 tryCatch 包裹 Seurat 加载，让 runtime-not-ready 能被映射为 exit 4
ensure_seurat_ready <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("[seurat-standard] biof3-r-runtime missing dependency: Seurat not installed")
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("[seurat-standard] biof3-r-runtime missing dependency: SeuratObject not installed")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("[seurat-standard] biof3-r-runtime missing dependency: patchwork not installed")
  }
}

tryCatch(ensure_seurat_ready(),
         error = function(e) {
           cat(conditionMessage(e), "\n", file = stderr())
           quit(status = 4, save = "no")
         })

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(jsonlite)
  library(Matrix)
})

# %||% helper
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

# canonicalize_png — 与既有 seurat-standard.R 同义: png::writePNG(readPNG(f), f)
canonicalize_png <- function(file_path) {
  if (!file.exists(file_path)) return(invisible(FALSE))
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("[seurat-standard] R package 'png' is required for deterministic Result Studio artifacts")
  }
  png::writePNG(png::readPNG(file_path), file_path)
  invisible(TRUE)
}

# Report progress to progress.json — UI 进度条驱动
report_progress <- function(pct, msg, job_dir) {
  progress_file <- file.path(job_dir, "progress.json")
  writeLines(
    toJSON(list(progress = pct, message = msg), auto_unbox = TRUE),
    progress_file
  )
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
if (is.na(job_dir) || !nzchar(job_dir)) {
  cat("[seurat-standard] missing job_dir argument\n", file = stderr())
  quit(status = 2, save = "no")
}

if (!dir.exists(job_dir)) {
  cat(sprintf("[seurat-standard] job_dir 不存在: %s\n", job_dir), file = stderr())
  quit(status = 2, save = "no")
}

output_dir <- file.path(job_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load BioF3 SCI theme (tools.js 自动塞到 jobDir)
theme_path <- file.path(job_dir, "_biof3-theme.R")
if (file.exists(theme_path)) {
  source(theme_path)
} else {
  cat("[seurat-standard] _biof3-theme.R missing — falling back to theme_classic()\n", file = stderr())
  theme_biof3 <- function(base_size = 10, base_family = "") {
    ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
      ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.6))
  }
  biof3_palette <- function(n) {
    if (n <= 10) {
      RColorBrewer::brewer.pal(max(3, n), "Set1")[1:n]
    } else {
      scales::hue_pal()(n)
    }
  }
  ggsave_biof3 <- function(plot, path, width = 6, height = 4, limitsize = TRUE) {
    ggplot2::ggsave(paste0(path, ".png"), plot = plot, width = width, height = height,
                    dpi = 300, limitsize = limitsize)
    invisible(TRUE)
  }
}

cat("[seurat-standard] Starting seurat-standard pipeline (plugin 八件套版 v1.0.0)\n")
cat(sprintf("[seurat-standard] R: %s\n", R.version.string))
cat(sprintf("[seurat-standard] Seurat: %s\n", as.character(packageVersion("Seurat"))))

# ---- Read params ----
params_path <- file.path(job_dir, "params.json")
if (!file.exists(params_path)) {
  params <- list()
} else {
  params <- fromJSON(params_path)
}

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

cat(sprintf("[seurat-standard] params: species=%s, min/max_features=%d/%d, max_mt=%g, nfeatures=%d, dims=1:%d, resolution=%g, top_n=%d\n",
            species, min_features, max_features, max_mt, nfeatures, dims_val, resolution, top_n_markers))

report_progress(5, "Seurat 5.x runtime 已就绪, 开始加载输入...", job_dir)

# ============================================================
# Step 1. Load data — auto-detect format
# ============================================================
cat("[seurat-standard] Step 1: Loading data...\n")

all_files <- list.files(job_dir, full.names = TRUE)
# Subdir output/ 不算入 input
all_files <- all_files[grepl("^[^/]+/?$", all_files) | !grepl("^/.*output/?$", all_files)]
all_files <- list.files(job_dir, pattern = "^matrix", full.names = TRUE)
all_files <- all_files[!grepl("\\.json$", all_files)]

has_mtx       <- any(grepl("matrix\\.mtx(\\.gz)?$", list.files(job_dir, full.names = TRUE)))
has_barcodes  <- any(grepl("barcodes\\.tsv(\\.gz)?$", list.files(job_dir, full.names = TRUE)))
has_features  <- any(grepl("(features|genes)\\.tsv(\\.gz)?$", list.files(job_dir, full.names = TRUE)))

counts <- NULL

# Detect: also support RDS pre-loaded Seurat object (for SUB-1.2 demo fixtures)
rds_files <- list.files(job_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) > 0 && (params$matrix_file %||% "") != "" &&
    grepl("\\.rds$", params$matrix_file %||% "", ignore.case = TRUE)) {
  # Accept RDS Seurat object directly (sub-fluent path: skip Step 1-2 if obj already built)
  rds_path <- params$matrix_file
  if (!file.exists(rds_path)) {
    rds_path <- rds_files[1]
  }
  cat(sprintf("[seurat-standard]   Loading Seurat object from RDS: %s\n", basename(rds_path)))
  obj <- readRDS(rds_path)
  if (!"seurat" %in% tolower(class(obj)[1])) {
    # Raw counts matrix RDS — wrap into Seurat
    counts <- obj
    rm(obj)
  }
} else if (has_mtx && has_barcodes && has_features) {
  cat("[seurat-standard]   Detected 10X三件套 (.mtx + barcodes + features)\n")
  counts <- Read10X(data.dir = job_dir)
} else {
  input_files <- all_files
  input_file <- input_files[1]
  if (is.null(input_file) || is.na(input_file)) {
    cat("[seurat-standard] No input matrix file found in job_dir\n", file = stderr())
    quit(status = 3, save = "no")
  }
  ext <- tolower(tools::file_ext(input_file))
  cat(sprintf("[seurat-standard]   Reading %s (ext=%s)\n", basename(input_file), ext))

  if (ext == "") {
    first_line <- readLines(input_file, n = 1, warn = FALSE)
    if (grepl("^%%MatrixMarket", first_line)) {
      ext <- "mtx"
    } else if (grepl("\t", first_line) && !grepl(",", first_line)) {
      ext <- "tsv"
    } else if (grepl(",", first_line)) {
      ext <- "csv"
    } else {
      cat(sprintf("[seurat-standard] Cannot auto-detect file format. First line: %.80s\n",
                  first_line), file = stderr())
      quit(status = 3, save = "no")
    }
    cat(sprintf("[seurat-standard]   Auto-detected format: %s\n", ext))
  }

  if (ext == "h5") {
    counts <- Read10X_h5(input_file)
    if (is.list(counts)) {
      counts <- counts[["Gene Expression"]]
    }
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- if (ext == "csv") "," else "\t"
    counts <- read.table(input_file, header = TRUE, sep = sep, row.names = 1, check.names = FALSE)
    counts <- as(as.matrix(counts), "CsparseMatrix")
  } else if (ext == "mtx") {
    counts <- Matrix::readMM(input_file)
  } else {
    cat(sprintf("[seurat-standard] Unsupported file extension: %s\n", ext), file = stderr())
    quit(status = 3, save = "no")
  }
}

if (!is.null(counts) && !is.null(dim(counts))) {
  cat(sprintf("[seurat-standard]   Loaded matrix: %d genes × %d cells\n",
              nrow(counts), ncol(counts)))

  if (ncol(counts) > 50000) {
    cat(sprintf("[seurat-standard] Dataset too large: %d cells > 50,000 cell limit\n",
                ncol(counts)), file = stderr())
    quit(status = 3, save = "no")
  }
}

report_progress(15, "Step 1 完成 — 表达矩阵加载完成", job_dir)

# ============================================================
# Step 2. CreateSeuratObject + QC metrics
# ============================================================
cat("[seurat-standard] Step 2: Creating Seurat object + QC metrics...\n")

if (exists("obj") && "Seurat" %in% class(obj)) {
  cat("[seurat-standard]   Using pre-loaded Seurat object from RDS\n")
} else {
  obj <- CreateSeuratObject(counts = counts, project = "biof3",
                            min.cells = 3, min.features = min_features)
}
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

# QC violin (3-panel) via patchwork
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

# QC metrics table
qc_metrics <- data.frame(
  cell        = rownames(obj@meta.data),
  nFeature_RNA = obj$nFeature_RNA,
  nCount_RNA   = obj$nCount_RNA,
  percent.mt   = obj$percent.mt
)
write.csv(qc_metrics, file.path(output_dir, "qc_metrics.csv"), row.names = FALSE)

n_before_qc <- ncol(obj)

report_progress(30, "Step 2 完成 — CreateSeuratObject + QC metrics", job_dir)

# ============================================================
# Step 3. Filter cells
# ============================================================
cat(sprintf("[seurat-standard] Step 3: Filtering cells (nFeature: %d-%d, percent.mt < %g%%)...\n",
            min_features, max_features, max_mt))

obj <- subset(obj, subset = nFeature_RNA > min_features &
                            nFeature_RNA < max_features &
                            percent.mt < max_mt)
n_after_qc <- ncol(obj)
cat(sprintf("[seurat-standard]   Cells: %d → %d after QC\n", n_before_qc, n_after_qc))

if (n_after_qc < 50) {
  cat(sprintf("[seurat-standard] Too few cells (%d) remaining after QC — empty_result (still valid, exit 0)\n",
              n_after_qc), file = stderr())

  # Write placeholder UMAP plot
  placeholder_path <- file.path(output_dir, "umap_clusters.png")
  if (file.exists(placeholder_path)) file.remove(placeholder_path)
  file.copy(
    from = file.path(dirname(theme_path), "tools/figures/umap_placeholder.png"),
    to = placeholder_path
  )
  if (!file.exists(placeholder_path)) {
    # Fallback: copy from demo-data or write minimal PNG
    placeholder_path <- file.path(output_dir, "umap_clusters.png")
  }

  # Write empty manifest
  empty_manifest <- list(
    outputs = c("qc_metrics.csv", "umap_clusters.png", "report.html", "manifest.json"),
    stats = list(
      cells_input    = n_before_qc,
      cells_after_qc = n_after_qc,
      n_genes        = nrow(obj),
      n_variable     = 0L,
      dims_used      = dims_val,
      resolution     = resolution,
      n_clusters     = 0L,
      n_markers      = 0L,
      species        = species,
      sci_style      = TRUE,
      dual_format    = TRUE,
      failure_mode   = "empty_result"
    )
  )
  writeLines(toJSON(empty_manifest, auto_unbox = TRUE, pretty = TRUE),
             file.path(output_dir, "manifest.json"))

  writeLines(
    sprintf("Empty result: only %d cells remain after QC.\n", n_after_qc),
    file.path(output_dir, "report.html")
  )
  quit(status = 5, save = "no")
}

# ============================================================
# Step 4. Normalize + Step 5. FindVariableFeatures
# ============================================================
cat("[seurat-standard] Step 4-5: Normalize + FindVariableFeatures...\n")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)

top10 <- head(VariableFeatures(obj), 10)
p_hvg_base <- VariableFeaturePlot(obj)
p_hvg <- LabelPoints(plot = p_hvg_base, points = top10, repel = TRUE) +
  theme_biof3(base_size = 10) +
  theme(legend.position = "bottom")
ggsave_biof3(p_hvg, file.path(output_dir, "hvg_plot"), width = 8, height = 6)

# HVG table
hvg_df <- HVFInfo(obj)
hvg_df$gene <- rownames(hvg_df)
hvg_df$is_top_variable <- hvg_df$gene %in% VariableFeatures(obj)
hvg_df <- hvg_df[order(-hvg_df$variance.standardized), ]
write.csv(hvg_df, file.path(output_dir, "hvg_table.csv"), row.names = FALSE)

report_progress(45, "Step 4-5 完成 — Normalize + HVG", job_dir)

# ============================================================
# Step 6. Scale + RunPCA + ElbowPlot
# ============================================================
cat("[seurat-standard] Step 6: ScaleData (all genes)...\n")
all_genes <- rownames(obj)
obj <- ScaleData(obj, features = all_genes, verbose = FALSE)

cat("[seurat-standard] Step 7: RunPCA + ElbowPlot...\n")
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)

# PCA DimPlot
p_pca <- DimPlot(obj, reduction = "pca") + NoLegend() +
  ggtitle("PCA DimPlot") + theme_biof3(base_size = 10)
ggsave_biof3(p_pca, file.path(output_dir, "pca_dimplot"), width = 7, height = 6)

# PCA dim loadings
p_load <- VizDimLoadings(obj, dims = 1:2, reduction = "pca") +
  theme_biof3(base_size = 10) +
  theme(legend.position = "right")
ggsave_biof3(p_load, file.path(output_dir, "pca_dim_loadings"),
             width = 8, height = 6)

# Elbow plot
p_elbow <- ElbowPlot(obj, ndims = max(20, dims_val + 5)) +
  theme_biof3(base_size = 10)
ggsave_biof3(p_elbow, file.path(output_dir, "elbow_plot"), width = 6, height = 4)

# PCA stdev CSV
pca_stdev <- data.frame(
  PC = paste0("PC_", 1:length(Stdev(obj, reduction = "pca"))),
  stdev = Stdev(obj, reduction = "pca")
)
write.csv(pca_stdev, file.path(output_dir, "pca_stdev.csv"), row.names = FALSE)

# PC1-2 loadings CSV
loadings_mat <- Loadings(obj[["pca"]])[, 1:2, drop = FALSE]
loadings_df <- data.frame(
  gene  = rownames(loadings_mat),
  PC_1  = loadings_mat[, 1],
  PC_2  = loadings_mat[, 2]
)
loadings_df <- loadings_df[order(-abs(loadings_df$PC_1)), ]
write.csv(head(loadings_df, 50), file.path(output_dir, "pca_loadings.csv"), row.names = FALSE)

report_progress(65, "Step 6-7 完成 — Scale + PCA + Elbow", job_dir)

# ============================================================
# Step 8. FindNeighbors + FindClusters + Step 9. RunUMAP
# ============================================================
cat(sprintf("[seurat-standard] Step 8-9: Cluster (dims=1:%d, res=%g) + UMAP...\n",
            dims_val, resolution))
obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)

n_clusters <- length(levels(Idents(obj)))
cat(sprintf("[seurat-standard]   Found %d clusters\n", n_clusters))

# UMAP clusters PNG
cluster_pal <- biof3_palette(n_clusters)
p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.5,
                   cols = cluster_pal) +
  NoLegend() + ggtitle(sprintf("UMAP — %d clusters", n_clusters)) +
  theme_biof3(base_size = 10)
ggsave_biof3(p_umap, file.path(output_dir, "umap_clusters"), width = 7, height = 6)

# UMAP coords + cluster
umap_coords <- data.frame(
  cell    = colnames(obj),
  UMAP_1  = Embeddings(obj, "umap")[, 1],
  UMAP_2  = Embeddings(obj, "umap")[, 2],
  cluster = as.character(Idents(obj))
)
write.csv(umap_coords, file.path(output_dir, "umap_coords.csv"), row.names = FALSE)

report_progress(80, "Step 8-9 完成 — Cluster + UMAP", job_dir)

# ============================================================
# Step 10. FindAllMarkers + visualize top markers
# ============================================================
cat("[seurat-standard] Step 10: FindAllMarkers...\n")
markers <- FindAllMarkers(obj, only.pos = TRUE, verbose = FALSE)
write.csv(markers, file.path(output_dir, "markers_table.csv"), row.names = FALSE)

top_markers <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = top_n_markers) %>%
  dplyr::ungroup()

# DoHeatmap
top_marker_genes <- unique(top_markers$gene)
if (length(top_marker_genes) > 1) {
  p_heat <- DoHeatmap(obj, features = top_marker_genes,
                      size = 3, raster = TRUE) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 9))
  heat_h <- max(8, min(20, length(top_marker_genes) * 0.12))
  ggsave_biof3(p_heat, file.path(output_dir, "markers_heatmap"),
               width = 11, height = heat_h, limitsize = FALSE)
}

# Top 6 markers
top6 <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 1) %>%
  dplyr::ungroup() %>%
  head(6) %>%
  dplyr::pull(gene)
top6 <- unique(top6)

if (length(top6) >= 1) {
  p_vln <- VlnPlot(obj, features = top6, ncol = 3, pt.size = 0) &
    theme_biof3(base_size = 9) &
    theme(legend.position = "none")
  ggsave_biof3(p_vln, file.path(output_dir, "markers_violin"),
               width = 11, height = 7)

  p_feat <- FeaturePlot(obj, features = top6, ncol = 3,
                        reduction = "umap", order = TRUE) &
    theme_biof3(base_size = 9)
  ggsave_biof3(p_feat, file.path(output_dir, "feature_plot"),
               width = 11, height = 7)
}

report_progress(90, "Step 10 完成 — FindAllMarkers", job_dir)

# ============================================================
# Save the full Seurat object (RDS)
# ============================================================
cat("[seurat-standard] Saving Seurat object (.rds)...\n")
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
# Manifest — 与 tool-definition.json outputs[] 对齐的 8 件套
# ============================================================
manifest <- list(
  outputs = c(
    "seurat_obj.rds",
    "markers_table.csv",
    "qc_metrics.csv",
    "umap_clusters.png",
    "umap_coords.csv",
    "hvg_plot.png",
    "report.html",
    "manifest.json"
  ),
  stats = list(
    cells_input    = n_before_qc,
    cells_after_qc = n_after_qc,
    n_genes        = nrow(obj),
    n_variable     = nfeatures,
    dims_used      = dims_val,
    resolution     = resolution,
    n_clusters     = n_clusters,
    n_markers      = nrow(markers),
    species        = species,
    sci_style      = TRUE,
    dual_format    = TRUE
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ============================================================
# Generate HTML interpretation report
# ============================================================
cat("[seurat-standard] Generating HTML report...\n")
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

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

  report_html <- gsub("{{cells_input}}",    as.character(n_before_qc),  report_html, fixed = TRUE)
  report_html <- gsub("{{cells_after_qc}}", as.character(n_after_qc),   report_html, fixed = TRUE)
  report_html <- gsub("{{n_genes}}",        as.character(nrow(obj)),    report_html, fixed = TRUE)
  report_html <- gsub("{{n_variable}}",     as.character(nfeatures),    report_html, fixed = TRUE)
  report_html <- gsub("{{dims_used}}",      as.character(dims_val),     report_html, fixed = TRUE)
  report_html <- gsub("{{resolution}}",     as.character(resolution),   report_html, fixed = TRUE)
  report_html <- gsub("{{n_clusters}}",     as.character(n_clusters),   report_html, fixed = TRUE)
  report_html <- gsub("{{n_markers}}",      as.character(nrow(markers)),report_html, fixed = TRUE)
  report_html <- gsub("{{species}}",        species,                    report_html, fixed = TRUE)

  # Safety net: strip any unreplaced {{*.png}} blocks
  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)

  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[seurat-standard] Report generated.\n")
} else {
  cat("[seurat-standard] Warning: report template not found, skipping report.html.\n")
}

# ============================================================
# PNG canonicalize (deterministic sha256 for Result Studio)
# ============================================================
cat("[seurat-standard] Canonicalize PNGs...\n")
for (png in list.files(output_dir, pattern = "\\.png$", full.names = TRUE)) {
  canonicalize_png(png)
}

report_progress(100, "Seurat-standard pipeline complete.", job_dir)

cat("[seurat-standard] Pipeline complete.\n")