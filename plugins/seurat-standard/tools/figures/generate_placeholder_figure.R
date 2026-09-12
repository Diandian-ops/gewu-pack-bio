#!/usr/bin/env Rscript
# ============================================================
# generate_placeholder_figure.R — 生成 seurat-standard placeholder UMAP PNG
#
# 输出: resources/built-in-plugins/seurat-standard/tools/figures/umap_placeholder.png
#
# 设计：
#   - 256×64 灰色条带（与 scanpy-advanced 占位图风格一致）
#   - png::writePNG(png::readPNG(f), f) 重写以稳定 sha256（与 scanpy Pillow 等价）
#   - 确定性: 不依赖 random / 时间戳
#
# 运行:
#   /path/to/biof3-r-runtime/bin/Rscript generate_placeholder_figure.R
# ============================================================

if (!requireNamespace("png", quietly = TRUE)) {
  stop("[placeholder figure] R package 'png' is required")
}

# Build a minimal valid PNG by hand via png::writePNG
# Output dimensions
w <- 256
h <- 64

# Create a 4-channel RGBA matrix: h × w × 4
# Use a grey gradient with subtle banding (deterministic)
mat <- array(0L, dim = c(h, w, 4L))
for (y in seq_len(h)) {
  for (x in seq_len(w)) {
    intensity <- (x %% 64) * 3 + 128
    r <- (intensity + ((x %/% 8) %% 16)) %% 256
    g <- r
    b <- r
    a <- 255L
    mat[y, x, 1] <- r
    mat[y, x, 2] <- g
    mat[y, x, 3] <- b
    mat[y, x, 4] <- a
  }
}

script_dir <- tryCatch({
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) ".")
out_path <- file.path(script_dir, "umap_placeholder.png")

# Write first time
png::writePNG(mat / 255, target = out_path)

# Re-read + re-write to canonicalize (R png::writePNG(png::readPNG(f), f) equivalent)
buf <- png::readPNG(out_path)
png::writePNG(buf, target = out_path)

sha <- digest::digest(file = out_path, algo = "sha256")
cat(sprintf("[placeholder figure] wrote %s (%d bytes, sha256-prefix=%s)\n",
            out_path,
            file.info(out_path)$size,
            substr(sha, 1, 12)))