# ============================================================
# BioF3 SCI-style 共享绘图主题
# Version: 1.0 (2026-05-24)
#
# 用法（脚本已复制到本地 job 目录）:
#   source(file.path(job_dir, "_biof3-theme.R"))
#
# 提供的工具:
#   theme_biof3()        - ggplot2 标准主题, 替换 theme_classic() / theme_minimal()
#   biof3_palette()      - 分类离散调色板 (NPG-based, 10 色)
#   biof3_palette_div()  - 二极配色 (蓝-白-红, 用于热图/火山)
#   biof3_palette_seq()  - 渐变配色 (Viridis, 用于连续值)
#   ggsave_biof3()       - 标准化的图片保存 (PNG + PDF 双格式, 300dpi)
#   save_grid_biof3()    - 给 grid 系图 (ComplexHeatmap, pheatmap, base R) 用的双格式保存
#
# 设计原则:
#   - 默认 6×4.5 inch / 300dpi → 单栏论文图 (Nature column 88mm = 3.5in, double = 7.2in)
#   - 字体最小 9pt (axis), 标题 13pt face=bold
#   - panel.border 全边框 (sci 标准, 不是只有 axis line)
#   - 移除背景网格 (major=轻灰, minor=无), 不要花哨
#   - 同时出 PDF (Cairo + 矢量) 给用户直接发表用
# ============================================================

# Avoid double-source side effects
if (exists(".biof3_theme_loaded")) {
  invisible(NULL)
} else {
  .biof3_theme_loaded <- TRUE

  # Load only once
  suppressMessages({
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for theme_biof3")
    library(ggplot2)
  })

  # ---- 主题：theme_biof3() ----
  # 替换 theme_classic / theme_minimal, 全工具默认.
  # 单参数 base_size: 全局字体基准, 默认 11 (一个紧凑但能看清的尺寸).
  theme_biof3 <- function(base_size = 11, base_family = "") {
    theme_classic(base_size = base_size, base_family = base_family) +
      theme(
        # 边框 (SCI 风格: 全边框, 不是 L 型)
        panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
        axis.line          = element_blank(),  # by panel.border instead

        # 文本
        axis.text          = element_text(color = "black", size = base_size - 1),
        axis.title         = element_text(color = "black", size = base_size + 1, face = "bold"),
        plot.title         = element_text(color = "black", size = base_size + 2, face = "bold",
                                           hjust = 0, margin = margin(b = 6)),
        plot.subtitle      = element_text(color = "gray30", size = base_size - 1,
                                           margin = margin(b = 8)),
        plot.caption       = element_text(color = "gray40", size = base_size - 2,
                                           hjust = 1, margin = margin(t = 8)),

        # 网格 (移除花哨, 保留 minor 主网格的浅淡)
        panel.grid.major   = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.background   = element_rect(fill = "white", color = NA),

        # 图例
        legend.background  = element_rect(fill = "white", color = NA),
        legend.key         = element_rect(fill = "white", color = NA),
        legend.text        = element_text(size = base_size - 1, color = "black"),
        legend.title       = element_text(size = base_size, face = "bold", color = "black"),
        legend.key.size    = unit(0.4, "cm"),

        # 边距 (整体紧凑)
        plot.margin        = margin(t = 10, r = 12, b = 8, l = 8, unit = "pt"),

        # 分面
        strip.background   = element_rect(fill = "gray95", color = "black", linewidth = 0.6),
        strip.text         = element_text(color = "black", size = base_size, face = "bold")
      )
  }

  # ---- 调色板 ----
  # 离散分类: 优先 ggsci NPG (Nature Publishing Group), 退化到 hand-picked 色
  biof3_palette <- function(n = 10) {
    # NPG 10 色板 (从 ggsci::pal_npg("nrc")(10) 抠出, 防止依赖 ggsci 必装)
    npg <- c(
      "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
      "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
      "#7E6148", "#B09C85"
    )
    if (n <= 10) {
      npg[seq_len(n)]
    } else {
      colorRampPalette(npg)(n)
    }
  }

  # 二极: 蓝-白-红, 适合热图 / 火山图 sig 标记
  biof3_palette_div <- function(n = 100) {
    colorRampPalette(c("#2563eb", "white", "#dc2626"))(n)
  }

  # 连续渐变: viridis-magma 风格 (色盲友好)
  biof3_palette_seq <- function(n = 100, option = "viridis") {
    # 不直接依赖 viridisLite, 用手抠的 viridis 锚点
    if (option == "viridis") {
      colorRampPalette(c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"))(n)
    } else {
      # ylorbr  
      colorRampPalette(c("#FFFFE5", "#FED98E", "#FE9929", "#D95F0E", "#993404"))(n)
    }
  }

  # ---- ggsave_biof3: 标准化图片保存 (PNG + PDF) ----
  # 用法: ggsave_biof3(plot, "filename_no_ext", width=6, height=4.5)
  # 不需要写 .png 后缀, 会自动出双格式
  ggsave_biof3 <- function(plot, filename_base, width = 6, height = 4.5,
                            dpi = 300, units = "in", limitsize = TRUE,
                            png_only = FALSE, pdf_only = FALSE) {
    # filename_base 可能已带 .png / .pdf 扩展名 (兼容旧代码), 自动剥
    base_no_ext <- sub("\\.(png|pdf|svg|jpg|jpeg)$", "", filename_base, ignore.case = TRUE)

    if (!pdf_only) {
      ggsave(paste0(base_no_ext, ".png"), plot,
             width = width, height = height, dpi = dpi, units = units,
             limitsize = limitsize, bg = "white")
    }

    if (!png_only) {
      # PDF: 用 Cairo (向量化, 字体内嵌, 无栅格化伪影). 容器有 Cairo 1.7+
      tryCatch({
        ggsave(paste0(base_no_ext, ".pdf"), plot,
               width = width, height = height, units = units,
               device = grDevices::cairo_pdf,
               limitsize = limitsize, bg = "white")
      }, error = function(e) {
        # Cairo 失败的 fallback (主要环境差异时):
        ggsave(paste0(base_no_ext, ".pdf"), plot,
               width = width, height = height, units = units,
               limitsize = limitsize, bg = "white")
      })
    }
    invisible(NULL)
  }

  # ---- save_grid_biof3: 给 grid 系 (ComplexHeatmap / pheatmap / base R) 双格式保存 ----
  # 用法:
  #   save_grid_biof3("filename_no_ext", width=6, height=4.5, expr={
  #     pheatmap(mat, ...)  或  draw(Heatmap(mat, ...))  或  plotMA(...)
  #   })
  # 第一个参数是文件名 (不带扩展名), expr 是绘图代码 block.
  save_grid_biof3 <- function(filename_base, width = 6, height = 4.5,
                               dpi = 300, units = "in", expr,
                               png_only = FALSE, pdf_only = FALSE) {
    base_no_ext <- sub("\\.(png|pdf|svg|jpg|jpeg)$", "", filename_base, ignore.case = TRUE)

    # Capture expr UNEVALUATED so we can run it twice (once per device).
    # Without substitute(), R's promise mechanism evaluates expr once and caches
    # the result (typically NULL since plotting is for side effects), making the
    # second device write a blank page.
    expr_q <- substitute(expr)
    parent_env <- parent.frame()

    # Convert in-to-px for png() (which uses px not in)
    px_per_in <- dpi
    width_px  <- if (units == "in") width * px_per_in else width
    height_px <- if (units == "in") height * px_per_in else height

    if (!pdf_only) {
      png(paste0(base_no_ext, ".png"),
          width = width_px, height = height_px, res = dpi, bg = "white")
      try(eval(expr_q, envir = parent_env), silent = FALSE)
      dev.off()
    }

    if (!png_only) {
      tryCatch({
        grDevices::cairo_pdf(paste0(base_no_ext, ".pdf"),
                              width = width, height = height, bg = "white")
        try(eval(expr_q, envir = parent_env), silent = FALSE)
        dev.off()
      }, error = function(e) {
        pdf(paste0(base_no_ext, ".pdf"),
            width = width, height = height, bg = "white")
        try(eval(expr_q, envir = parent_env), silent = FALSE)
        dev.off()
      })
    }
    invisible(NULL)
  }

  # ---- 数据可视化辅助 ----
  # 给 ComplexHeatmap 的统一颜色函数: 自动从矩阵值算 [-max(|x|), max(|x|)] 范围
  biof3_heat_colors <- function(mat, kind = c("div", "seq")) {
    kind <- match.arg(kind)
    if (kind == "div") {
      m <- max(abs(mat), na.rm = TRUE)
      circlize::colorRamp2(c(-m, 0, m), c("#2563eb", "white", "#dc2626"))
    } else {
      r <- range(mat, na.rm = TRUE)
      mid <- mean(r)
      circlize::colorRamp2(c(r[1], mid, r[2]),
                            c("#FFFFE5", "#FE9929", "#993404"))
    }
  }

  cat("[biof3-theme] v1.0 loaded (theme_biof3, biof3_palette, ggsave_biof3, save_grid_biof3)\n")
}
