#!/usr/bin/env Rscript
# ============================================================
# Tool: maftools-mutation
# Version targeted: maftools 2.24.0 (R 4.5.1)
# Canonical workflow source:
#   https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html
# Reference:
#   Mayakonda A et al. (2018). Maftools: efficient and comprehensive analysis of
#   somatic variants in cancer. Genome Research 28:1747-1756.
#
# 第三个使用 BioF3 SCI 视觉风格规范的工具.
# 注意: maftools 全是 base R 绘图, 不能套 theme_biof3(); 全部用 save_grid_biof3()
# 把 base R plot 块包装成 PNG + PDF 双格式输出.
#
# Steps (in-scope) — Tutorial 主流程 7 步:
#   1. Load MAF
#   2. plotmafSummary (6 panel)
#   3. oncoplot (top N)
#   4. plotTiTv
#   5. lollipopPlot (top N genes)
#   6. somaticInteractions
#   7. extractSignatures + compareSignatures (COSMIC v3)
#
# Out-of-scope:
#   - tcgaCompare (联网下载, 不稳)
#   - drugInteractions (DGIdb API)
#   - OncogenicPathways
#   - mafSurvival (用 KM 工具)
#
# Default params:
#   genome=hg19, top_n_oncoplot=20, lollipop_top_n=2, n_signatures=3
#
# Input (job_dir/):
#   - maf_file    .maf / .maf.gz / .txt
#   - clinical    (可选) CSV/TSV
#   - params.json
#
# Output (job_dir/output/):
#   8 plots × PNG + PDF = 16 files + 7 中间数据 + maf.rds + report.html + manifest.json
# ============================================================

suppressMessages({
  library(maftools)
  library(jsonlite)
  library(data.table)
  # NOTE: must explicitly attach NMF here. extractSignatures() internally calls
  # NMF::nmf() which dispatches via match.fun("brunet"). Without library(NMF),
  # the algorithm registry is not initialized -> error "'what' must be a function
  # or character string". Confirmed via _test_mafsig3.R against the bundled runtime.
  library(NMF)
})

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 0 && !is.list(a))) b else a
}

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
output_dir <- file.path(job_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(job_dir, "_biof3-theme.R"))

cat("[maftools] Starting pipeline\n")
cat("[maftools] maftools:", as.character(packageVersion("maftools")), "\n")

# ---- Read params ----
params <- fromJSON(file.path(job_dir, "params.json"))
genome           <- params$genome %||% "hg19"
top_n_oncoplot   <- as.integer(params$top_n_oncoplot %||% 20)
lollipop_top_n   <- as.integer(params$lollipop_top_n %||% 2)
n_signatures     <- as.integer(params$n_signatures %||% 3)

cat(sprintf("[maftools] params: genome=%s, top_n_oncoplot=%d, lollipop_top_n=%d, n_sig=%d\n",
            genome, top_n_oncoplot, lollipop_top_n, n_signatures))

# ---- Read inputs ----
maf_files <- list.files(job_dir, pattern = "^maf_file", full.names = TRUE)
if (length(maf_files) == 0) stop("No MAF file found in job_dir (expected name 'maf_file')")
maf_path <- maf_files[1]

# Detect if file is gzipped (auto by content if no extension)
first_bytes <- readBin(maf_path, "raw", n = 2)
is_gz <- length(first_bytes) >= 2 && first_bytes[1] == as.raw(0x1f) && first_bytes[2] == as.raw(0x8b)
cat(sprintf("[maftools] MAF file: %s (gzipped: %s)\n", basename(maf_path), is_gz))

# If gzipped, decompress to a sibling *.maf so data.table::fread (used by read.maf)
# can read it. read.maf relies on extension to choose reader; without .gz / .maf
# it tries plain text and fails on the binary header.
if (is_gz) {
  decompressed_path <- file.path(job_dir, "maf_file.decompressed.maf")
  con_in <- gzfile(maf_path, "rb")
  con_out <- file(decompressed_path, "wb")
  while (length(chunk <- readBin(con_in, "raw", n = 1048576)) > 0) {
    writeBin(chunk, con_out)
  }
  close(con_in); close(con_out)
  maf_path <- decompressed_path
  cat("[maftools] decompressed →", basename(maf_path), "\n")
}

# Clinical (optional)
clin_files <- list.files(job_dir, pattern = "^clinical$", full.names = TRUE)
clin_path <- if (length(clin_files) > 0) clin_files[1] else NULL

# ---- Step 1: read.maf ----
cat("[maftools] Step 1: read.maf...\n")
laml <- if (!is.null(clin_path)) {
  read.maf(maf = maf_path, clinicalData = clin_path, verbose = FALSE)
} else {
  read.maf(maf = maf_path, verbose = FALSE)
}

n_samples <- as.integer(getSampleSummary(laml)[, .N])
n_genes   <- length(unique(getGeneSummary(laml)[[1]]))
cat(sprintf("  Samples: %d, mutated genes: %d\n", n_samples, n_genes))

if (n_samples < 5) {
  stop(sprintf("Too few samples (%d) — need at least 5 for meaningful summary", n_samples))
}

# ---- Step 2: plotmafSummary ----
cat("[maftools] Step 2: plotmafSummary...\n")
save_grid_biof3(file.path(output_dir, "maf_summary"),
                 width = 12, height = 7, expr = {
  plotmafSummary(maf = laml, rmOutlier = TRUE,
                 addStat = "median", dashboard = TRUE, titvRaw = FALSE)
})

# ---- Step 3: oncoplot ----
cat("[maftools] Step 3: oncoplot...\n")
clin_features <- NULL
if (!is.null(clin_path)) {
  clin_data <- getClinicalData(laml)
  # Pick non-id non-numeric columns as annotation features
  candidate <- setdiff(colnames(clin_data), c("Tumor_Sample_Barcode"))
  for (col in candidate) {
    val <- clin_data[[col]]
    if (is.character(val) || is.factor(val)) {
      n_unique <- length(unique(val[!is.na(val)]))
      if (n_unique >= 2 && n_unique <= 10) {
        clin_features <- c(clin_features, col)
        if (length(clin_features) >= 3) break  # cap at 3 annotation tracks
      }
    }
  }
}

if (!is.null(clin_features)) cat(sprintf("  Using clinical features: %s\n", paste(clin_features, collapse = ", ")))

save_grid_biof3(file.path(output_dir, "oncoplot"),
                 width = max(10, top_n_oncoplot * 0.5 + 6),
                 height = max(8, top_n_oncoplot * 0.3 + 4), expr = {
  oncoplot(maf = laml, top = top_n_oncoplot,
           clinicalFeatures = clin_features,
           sortByAnnotation = !is.null(clin_features),
           draw_titv = TRUE)
})

# ---- Step 4: plotTiTv ----
cat("[maftools] Step 4: plotTiTv...\n")
laml.titv <- titv(maf = laml, plot = FALSE, useSyn = TRUE)
save_grid_biof3(file.path(output_dir, "titv"),
                 width = 9, height = 6, expr = {
  plotTiTv(res = laml.titv)
})

# Save TiTv data
fwrite(as.data.table(laml.titv$fraction.contribution),
       file.path(output_dir, "titv_data.csv"))

# ---- Step 5: Lollipop for top N genes ----
cat("[maftools] Step 5: lollipopPlot for top genes...\n")
gene_summary <- getGeneSummary(laml)
top_genes <- head(gene_summary$Hugo_Symbol, lollipop_top_n)
cat(sprintf("  Top %d genes: %s\n", lollipop_top_n, paste(top_genes, collapse = ", ")))

lollipop_files <- character(0)
for (i in seq_along(top_genes)) {
  g <- top_genes[i]
  cat(sprintf("    [%d/%d] %s ... ", i, length(top_genes), g))
  out_file <- file.path(output_dir, sprintf("lollipop_top%d", i))
  ok <- tryCatch({
    save_grid_biof3(out_file, width = 9, height = 5, expr = {
      lollipopPlot(maf = laml, gene = g,
                    AACol = "Protein_Change",
                    showMutationRate = TRUE,
                    labelPos = "all")
    })
    cat("OK\n")
    lollipop_files <- c(lollipop_files, paste0(basename(out_file), c(".png", ".pdf")))
    TRUE
  }, error = function(e) {
    cat(sprintf("FAILED (%s)\n", conditionMessage(e)))
    FALSE
  })
}

# ---- Step 6: somaticInteractions ----
cat("[maftools] Step 6: somaticInteractions...\n")
inter_result <- NULL
save_grid_biof3(file.path(output_dir, "somatic_interactions"),
                 width = 9, height = 7, expr = {
  inter <- somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
  assign("inter_result", inter, envir = .GlobalEnv)
})

if (!is.null(inter_result)) {
  fwrite(as.data.table(inter_result), file.path(output_dir, "interaction_matrix.csv"))
}

# ---- Step 7: Mutational Signatures ----
cat("[maftools] Step 7: Mutational Signatures...\n")
sig_ok <- FALSE
sig_decompose <- NULL

bsg_pkg <- if (genome == "hg38") "BSgenome.Hsapiens.UCSC.hg38" else "BSgenome.Hsapiens.UCSC.hg19"

if (requireNamespace(bsg_pkg, quietly = TRUE)) {
  tnm <- tryCatch({
    suppressMessages(suppressWarnings(
      trinucleotideMatrix(maf = laml, ref_genome = bsg_pkg,
                          prefix = "chr", add = TRUE, useSyn = FALSE)
    ))
  }, error = function(e) {
    cat(sprintf("  trinucleotideMatrix FAILED: %s\n", conditionMessage(e)))
    NULL
  })

  if (!is.null(tnm)) {
    cat(sprintf("  Trinucleotide matrix: %d samples × 96 contexts\n", nrow(tnm$nmf_matrix)))
    sig_decompose <- tryCatch({
      suppressMessages(suppressWarnings(
        # NOTE: maftools 2.24 default is `parallel=4`, which sends `paste0("P", 4) = "P4"`
        # to NMF::nmf's `.opt` arg. NMF 0.28 rejects this with `'what' must be a function`.
        # Pass parallel=NULL explicitly to disable the .opt branch.
        extractSignatures(mat = tnm, n = n_signatures,
                          plotBestFitRes = FALSE,
                          pConstant = 0.1,
                          parallel = NULL)
      ))
    }, error = function(e) {
      cat(sprintf("  extractSignatures (n=%d) FAILED: %s\n", n_signatures, conditionMessage(e)))
      cat("  Retrying with n=2 (smaller signature count)...\n")
      tryCatch({
        suppressMessages(suppressWarnings(
          extractSignatures(mat = tnm, n = 2,
                            plotBestFitRes = FALSE,
                            pConstant = 0.1,
                            parallel = NULL)
        ))
      }, error = function(e2) {
        cat(sprintf("  Retry FAILED: %s\n", conditionMessage(e2)))
        NULL
      })
    })

    if (!is.null(sig_decompose)) {
      sig_ok <- TRUE
      cat(sprintf("  Extracted %d signatures\n", n_signatures))

      # Plot signatures
      save_grid_biof3(file.path(output_dir, "mutation_signature"),
                       width = 11, height = max(4, n_signatures * 2.5), expr = {
        plotSignatures(nmfRes = sig_decompose, contributions = TRUE,
                        title_size = 1.0, sig_db = "SBS")
      })

      # Compare to COSMIC v3
      cosmic_compare <- tryCatch({
        suppressMessages(suppressWarnings(
          compareSignatures(nmfRes = sig_decompose, sig_db = "SBS")
        ))
      }, error = function(e) NULL)

      if (!is.null(cosmic_compare)) {
        save_grid_biof3(file.path(output_dir, "signature_cosmic_compare"),
                         width = max(10, n_signatures * 3 + 4),
                         height = 6, expr = {
          # Heatmap of cosine similarity vs COSMIC
          mat <- cosmic_compare$cosine_similarities
          if (!is.null(mat) && nrow(mat) > 0) {
            par(mar = c(8, 5, 3, 2))
            top_cos <- sort(apply(mat, 2, max), decreasing = TRUE)[1:min(15, ncol(mat))]
            mat_top <- mat[, names(top_cos), drop = FALSE]
            image(t(mat_top), col = biof3_palette_seq(50, "ylorbr"),
                  axes = FALSE, main = "Cosine similarity to COSMIC v3 SBS")
            axis(1, at = seq(0, 1, length.out = ncol(mat_top)),
                 labels = colnames(mat_top), las = 2, cex.axis = 0.8)
            axis(2, at = seq(0, 1, length.out = nrow(mat_top)),
                 labels = rownames(mat_top), las = 1, cex.axis = 0.9)
            box()
          }
        })
      }

      # Save signature contributions
      if (!is.null(sig_decompose$contributions)) {
        contrib <- as.data.frame(sig_decompose$contributions)
        contrib$sample <- rownames(contrib)
        contrib <- contrib[, c("sample", setdiff(colnames(contrib), "sample"))]
        fwrite(contrib, file.path(output_dir, "signature_contributions.csv"))
      }
    }
  }
} else {
  cat(sprintf("  %s not installed, skipping signature analysis\n", bsg_pkg))
}

# ============================================================
# Intermediate data
# ============================================================
fwrite(as.data.table(getGeneSummary(laml)),
       file.path(output_dir, "gene_summary.csv"))
fwrite(as.data.table(getSampleSummary(laml)),
       file.path(output_dir, "sample_summary.csv"))
if (!is.null(clin_path)) {
  fwrite(as.data.table(getClinicalData(laml)),
         file.path(output_dir, "clinical_data.csv"))
}
saveRDS(laml, file.path(output_dir, "maf.rds"))

# ============================================================
# Summary
# ============================================================
top10 <- head(gene_summary$Hugo_Symbol, 10)
top10_freq <- head(gene_summary$AlteredSamples, 10)

summary_text <- sprintf(
"maftools Mutation Analysis Summary

Input:
  Samples: %d
  Mutated genes: %d
  Total mutations: %d
  Genome: %s
  Clinical data: %s

Top 10 most-mutated genes:
%s

Mutational signatures:
  Extracted: %s
  Number of signatures: %d
",
  n_samples,
  n_genes,
  sum(getSampleSummary(laml)$total),
  genome,
  if (!is.null(clin_path)) "yes" else "no",
  paste(sprintf("  %d. %s (%d / %d samples)",
                seq_along(top10), top10, top10_freq, n_samples),
        collapse = "\n"),
  if (sig_ok) "yes" else "no (skipped or failed)",
  if (sig_ok) n_signatures else 0
)
writeLines(summary_text, file.path(output_dir, "summary.txt"))
cat(summary_text)

# ============================================================
# Manifest
# ============================================================
manifest_files <- list(
  list(name = "maf_summary.png",          type = "plot",  label = "MAF 概览 (PNG)"),
  list(name = "maf_summary.pdf",          type = "file",  label = "MAF 概览 (PDF)"),
  list(name = "oncoplot.png",             type = "plot",  label = "Oncoplot (PNG)"),
  list(name = "oncoplot.pdf",             type = "file",  label = "Oncoplot (PDF)"),
  list(name = "titv.png",                 type = "plot",  label = "TiTv (PNG)"),
  list(name = "titv.pdf",                 type = "file",  label = "TiTv (PDF)"),
  list(name = "somatic_interactions.png", type = "plot",  label = "互斥/共现 (PNG)"),
  list(name = "somatic_interactions.pdf", type = "file",  label = "互斥/共现 (PDF)"),
  list(name = "gene_summary.csv",         type = "table", label = "基因汇总"),
  list(name = "sample_summary.csv",       type = "table", label = "样本汇总"),
  list(name = "titv_data.csv",            type = "table", label = "TiTv 矩阵"),
  list(name = "interaction_matrix.csv",   type = "table", label = "互斥共现 p 值"),
  list(name = "maf.rds",                  type = "file",  label = "MAF 对象"),
  list(name = "summary.txt",              type = "text",  label = "分析摘要"),
  list(name = "report.html",              type = "file",  label = "解读报告 HTML")
)
# Add lollipop entries
for (i in seq_along(top_genes)) {
  manifest_files <- c(manifest_files, list(
    list(name = sprintf("lollipop_top%d.png", i), type = "plot", label = sprintf("Lollipop %s (PNG)", top_genes[i])),
    list(name = sprintf("lollipop_top%d.pdf", i), type = "file", label = sprintf("Lollipop %s (PDF)", top_genes[i]))
  ))
}
# Add signature entries (if generated)
if (sig_ok) {
  manifest_files <- c(manifest_files, list(
    list(name = "mutation_signature.png",       type = "plot",  label = "Signature 谱 (PNG)"),
    list(name = "mutation_signature.pdf",       type = "file",  label = "Signature 谱 (PDF)"),
    list(name = "signature_cosmic_compare.png", type = "plot",  label = "COSMIC 比对 (PNG)"),
    list(name = "signature_cosmic_compare.pdf", type = "file",  label = "COSMIC 比对 (PDF)"),
    list(name = "signature_contributions.csv",  type = "table", label = "Signature 贡献度")
  ))
}

manifest <- list(
  files = manifest_files,
  summary = list(
    n_samples       = n_samples,
    n_mutated_genes = n_genes,
    n_total_mutations = sum(getSampleSummary(laml)$total),
    genome          = genome,
    has_clinical    = !is.null(clin_path),
    top_genes       = head(top_genes, 10),
    signatures_extracted = sig_ok,
    n_signatures    = if (sig_ok) n_signatures else 0
  )
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "manifest.json"))

# ============================================================
# HTML report
# ============================================================
cat("[maftools] Generating HTML report...\n")
template_path <- file.path(job_dir, "report-template.html")
if (file.exists(template_path)) {
  report_html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  png_files <- c("maf_summary.png", "oncoplot.png", "titv.png",
                 "somatic_interactions.png",
                 "mutation_signature.png", "signature_cosmic_compare.png")
  for (i in seq_along(top_genes)) {
    png_files <- c(png_files, sprintf("lollipop_top%d.png", i))
  }

  for (png in png_files) {
    img_path <- file.path(output_dir, png)
    placeholder <- paste0("{{", png, "}}")
    if (file.exists(img_path)) {
      b64 <- base64enc::base64encode(img_path)
      data_uri <- paste0("data:image/png;base64,", b64)
      report_html <- gsub(placeholder, data_uri, report_html, fixed = TRUE)
    }
  }

  report_html <- gsub("{{n_samples}}",       as.character(n_samples),     report_html, fixed = TRUE)
  report_html <- gsub("{{n_mutated_genes}}", as.character(n_genes),       report_html, fixed = TRUE)
  report_html <- gsub("{{n_total_mutations}}", as.character(sum(getSampleSummary(laml)$total)), report_html, fixed = TRUE)
  report_html <- gsub("{{genome}}",          genome,                      report_html, fixed = TRUE)
  report_html <- gsub("{{top1_gene}}",       if (length(top_genes) >= 1) top_genes[1] else "N/A", report_html, fixed = TRUE)
  report_html <- gsub("{{top2_gene}}",       if (length(top_genes) >= 2) top_genes[2] else "N/A", report_html, fixed = TRUE)
  report_html <- gsub("{{n_signatures}}",    if (sig_ok) as.character(n_signatures) else "0", report_html, fixed = TRUE)
  report_html <- gsub("{{signatures_status}}", if (sig_ok) "成功提取" else "跳过/失败",   report_html, fixed = TRUE)
  report_html <- gsub("{{has_clinical}}",    if (!is.null(clin_path)) "有" else "无", report_html, fixed = TRUE)

  # Remove fig blocks for any unreplaced {{*.png}}
  report_html <- gsub('<div class="fig">\\s*<img src="\\{\\{[^}]+\\}\\}" [^>]*>\\s*<div class="fig-caption">[^<]*</div>\\s*</div>',
                      '', report_html, perl = TRUE)
  writeLines(report_html, file.path(output_dir, "report.html"))
  cat("[maftools] Report generated.\n")
}

cat("[maftools] Pipeline complete.\n")
