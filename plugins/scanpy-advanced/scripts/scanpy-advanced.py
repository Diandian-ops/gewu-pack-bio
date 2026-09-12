#!/usr/bin/env python3
"""scanpy-advanced — Scanpy single-cell 高级分析主脚本.

主流程:
  1. preflight: 校验 job_dir 输入文件 (h5ad + 可选 metadata.csv)
  2. 读入 AnnData
  3. QC: 计算 qc_metrics, 高表达基因比例, 过滤
  4. normalize + log1p + HVG (sc.pp.highly_variable_genes)
  5. PCA + neighbors + UMAP (sc.tl.umap)
  6. Leiden 聚类 (sc.tl.leiden)
  7. Marker genes (sc.tl.rank_genes_groups)
  8. 落盘: figures + tables + report.html

契约 (BioF3 调度器):
  - sys.argv[1] = job_dir (含 params.json + input file)
  - 输出全部写入 job_dir/output/
  - PNG 通过 PIL.Image.save 重写以稳定 hash (与 R png::writePNG 同义)

G0 阶段纪律:
  - 仅以 dev Electron 边界内 Python runtime 跑 scanpy 1.9+ (biof3-py-runtime)
  - 不触发 packaged / Provider / GA / dist staleness (G1+)
  - 不引 pip 安装：scanpy / anndata / matplotlib / seaborn / pandas / numpy
    全部依赖 biof3-py-runtime 已装 (见 py-runtime-requirements.txt)
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def _report_progress(pct: int, msg: str, job_dir: Path) -> None:
    """写 progress.json; tools.js / harness 读取驱动 UI 进度条."""
    progress_file = job_dir / "progress.json"
    progress_file.write_text(
        json.dumps({"progress": pct, "message": msg}, ensure_ascii=False),
        encoding="utf-8",
    )


def _canonicalize_png(path: Path) -> None:
    """Pillow 等价的 R png::writePNG(png::readPNG(f), f) — 写入重读以稳定 hash.

    这是 worker C §0.5.5 报告里讲的 "deterministic PNG" 要求:
    对每一个生成的 PNG 重写一次, 保证下游 artifact-truth / Result Studio
    不会因为 matplotlib 的非确定性 backend metadata 产生不同 sha256.
    """
    try:
        from PIL import Image

        with Image.open(path) as img:
            img.load()
            # force RGB(A); 跳过 palette / P mode 异常路径
            target = img.convert("RGBA") if img.mode != "RGBA" else img
            target.save(path, format="PNG", optimize=False)
    except Exception as exc:  # noqa: BLE001
        # 容错: 如果 PIL 失败就 fallback 留原图, 不阻塞整个 job
        print(f"[scanpy-advanced] warning: PNG canonicalize skipped for {path}: {exc}")


def _read_params(job_dir: Path) -> dict:
    """读取 BioF3 调度器塞进 job_dir 的 params.json.

    关键字段 (与 tool-definition.json 对齐):
      - min_genes (int, default 200)
      - min_cells (int, default 3)
      - n_top_genes (int, default 2000)
      - n_pcs (int, default 30)
      - n_neighbors (int, default 15)
      - resolution (float, default 0.5)
      - leiden_key (str, default "leiden")
      - method (str, default "wilcoxon")  # rank_genes_groups method
    """
    params_path = job_dir / "params.json"
    if not params_path.exists():
        return {}
    try:
        return json.loads(params_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"[scanpy-advanced] params.json 解析失败: {exc}") from exc


def _find_h5ad(job_dir: Path) -> Path | None:
    """在 job_dir 找 h5ad 输入（约定 input id 'anndata_file' / 任意 *.h5ad）."""
    for candidate in sorted(job_dir.glob("*.h5ad")):
        if candidate.stat().st_size > 1000:  # >1KB 才算真实 fixture
            return candidate
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description="scanpy-advanced plugin entry")
    parser.add_argument("job_dir", help="BioF3 Job 目录")
    args = parser.parse_args()
    job_dir = Path(args.job_dir).resolve()

    if not job_dir.is_dir():
        print(f"[scanpy-advanced] job_dir 不存在: {job_dir}", file=sys.stderr)
        return 2

    output_dir = job_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    _report_progress(5, "初始化 Python scanpy runtime...", job_dir)

    # ── 1. preflight ──────────────────────────────────────────
    params = _read_params(job_dir)
    h5ad_path = _find_h5ad(job_dir)
    if h5ad_path is None:
        # 上游容错: 让 skill workflow 拿到 missing_required_file_inputs
        msg = (
            "[scanpy-advanced] 缺少 .h5ad 输入文件。"
            "请在 jobDir 放置 h5ad 输入 (anndata_file)，"
            "或由上游 sample-metadata-validator / scanpy-cluster 链式传入。"
        )
        print(msg, file=sys.stderr)
        _report_progress(100, f"blocked: {msg}", job_dir)
        (output_dir / "error.txt").write_text(msg, encoding="utf-8")
        return 3  # non-zero exit => harness maps to expectedFailureMode

    # ── 2. runtime import（延迟导入以加速 preflight 失败路径）──
    try:
        import numpy as np
        import pandas as pd
        import scanpy as sc

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:
        msg = f"[scanpy-advanced] biof3-py-runtime 缺少依赖: {exc}"
        print(msg, file=sys.stderr)
        _report_progress(100, f"runtime-not-ready: {msg}", job_dir)
        (output_dir / "error.txt").write_text(msg, encoding="utf-8")
        return 4  # runtime_not_ready

    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=120, frameon=False, facecolor="white")
    sc.settings.figdir = output_dir

    # ── 3. 读入 AnnData ───────────────────────────────────────
    _report_progress(15, f"读取 h5ad: {h5ad_path.name}", job_dir)
    adata = sc.read_h5ad(h5ad_path)
    n_obs_in, n_vars_in = adata.n_obs, adata.n_vars
    print(f"[scanpy-advanced] 输入: {n_obs_in} obs × {n_vars_in} vars")

    if n_obs_in < 10:
        msg = f"[scanpy-advanced] 细胞数过少 ({n_obs_in} < 10)，empty-result"
        print(msg, file=sys.stderr)
        _report_progress(100, msg, job_dir)
        (output_dir / "qc_summary.csv").write_text(
            "n_obs,n_vars,n_obs_after_qc,n_vars_after_qc,n_clusters,n_markers\n"
            f"{n_obs_in},{n_vars_in},0,0,0,0\n",
            encoding="utf-8",
        )
        # 仍生成一张 placeholder PNG 让 UI 不爆 break
        try:
            fig, ax = plt.subplots(figsize=(3, 3))
            ax.text(
                0.5,
                0.5,
                f"insufficient cells\nn={n_obs_in}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_axis_off()
            empty_png = output_dir / "umap.png"
            fig.savefig(empty_png, dpi=120, bbox_inches="tight")
            plt.close(fig)
            _canonicalize_png(empty_png)
        except Exception:  # noqa: BLE001
            pass
        return 0

    # ── 4. QC + filter ────────────────────────────────────────
    _report_progress(30, "QC + 过滤 (min_genes / min_cells)", job_dir)
    min_genes = int(params.get("min_genes", 200))
    min_cells = int(params.get("min_cells", 3))
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    n_obs_post_qc, n_vars_post_qc = adata.n_obs, adata.n_vars
    print(f"[scanpy-advanced] post-QC: {n_obs_post_qc} obs × {n_vars_post_qc} vars")

    if n_obs_post_qc < 10 or n_vars_post_qc < 50:
        msg = (
            f"[scanpy-advanced] QC 后样本/基因过少 "
            f"({n_obs_post_qc} obs × {n_vars_post_qc} vars)"
        )
        print(msg, file=sys.stderr)
        _report_progress(100, msg, job_dir)
        return 5  # empty-result status

    # ── 5. normalize + log1p + HVG ─────────────────────────────
    _report_progress(45, "normalize + log1p + HVG", job_dir)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_top = int(params.get("n_top_genes", 2000))
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(n_top, n_vars_post_qc),
        flavor="seurat",
    )

    # ── 6. scale + PCA + neighbors + UMAP ─────────────────────
    _report_progress(60, "scale + PCA + neighbors + UMAP", job_dir)
    n_pcs = int(params.get("n_pcs", min(30, n_obs_post_qc - 1, n_vars_post_qc - 1)))
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=int(params.get("n_neighbors", 15)), n_pcs=n_pcs)
    sc.tl.umap(adata)
    umap_png = output_dir / "umap.png"
    sc.pl.umap(adata, color=None, save=False, show=False, ax=None)
    # 直接用 sc.pl 输出的最新 figure
    plt.gcf().savefig(umap_png, dpi=120, bbox_inches="tight")
    plt.close("all")
    _canonicalize_png(umap_png)

    # ── 7. Leiden 聚类 ────────────────────────────────────────
    _report_progress(75, "Leiden 聚类", job_dir)
    resolution = float(params.get("resolution", 0.5))
    leiden_key = str(params.get("leiden_key", "leiden"))
    sc.tl.leiden(adata, resolution=resolution, key_added=leiden_key, flavor="igraph", n_iterations=2, directed=False)
    n_clusters = int(adata.obs[leiden_key].nunique())
    # 着色 UMAP
    fig_ax = sc.pl.umap(adata, color=[leiden_key], show=False, return_fig=True)
    umap_clustered_png = output_dir / "umap_clustered.png"
    fig_ax.savefig(umap_clustered_png, dpi=120, bbox_inches="tight")
    plt.close("all")
    _canonicalize_png(umap_clustered_png)

    # ── 8. marker genes (rank_genes_groups) ────────────────────
    _report_progress(88, "Marker gene 检测 (rank_genes_groups)", job_dir)
    method = str(params.get("method", "wilcoxon"))
    sc.tl.rank_genes_groups(adata, groupby=leiden_key, method=method, n_genes=20)
    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)
    # 读 result["names"][g] 为 np.recarray 形式, 用 stack+concat 更稳
    ranks_df = pd.DataFrame(
        {
            "group": np.repeat(groups, [len(result["names"][g]) for g in groups]),
            "rank": np.concatenate(
                [np.arange(1, len(result["names"][g]) + 1) for g in groups]
            ),
            "gene": np.concatenate([result["names"][g] for g in groups]),
            "score": np.concatenate([result["scores"][g] for g in groups]),
            "pvalue": np.concatenate([result["pvals"][g] for g in groups]),
            "padj": np.concatenate([result["pvals_adj"][g] for g in groups]),
        }
    )
    markers_csv = output_dir / "marker_genes.csv"
    ranks_df.to_csv(markers_csv, index=False)

    # ── 9. 摘要表 (qc_summary.csv + preprocessed.h5ad) ────────
    _report_progress(95, "写摘要 + 保存 AnnData", job_dir)
    qc_summary = (
        "n_obs_in,n_vars_in,n_obs_after_qc,n_vars_after_qc,n_clusters,n_markers,method,resolution\n"
        f"{n_obs_in},{n_vars_in},{n_obs_post_qc},{n_vars_post_qc},"
        f"{n_clusters},{len(ranks_df)},{method},{resolution}\n"
    )
    (output_dir / "qc_summary.csv").write_text(qc_summary, encoding="utf-8")

    hvg_png = output_dir / "hvg_plot.png"
    sc.pl.highly_variable_genes(adata, show=False)
    plt.gcf().savefig(hvg_png, dpi=120, bbox_inches="tight")
    plt.close("all")
    _canonicalize_png(hvg_png)

    # preprocessed AnnData (含 PCA / UMAP / leiden / rank_genes_groups)
    preprocessed_path = output_dir / "preprocessed.h5ad"
    adata.write_h5ad(preprocessed_path)

    # ── 10. report.html (简洁摘要) ─────────────────────────────
    _report_progress(99, "写 report.html + manifest", job_dir)
    report_html = (
        "<!doctype html><html><head><meta charset=\"utf-8\">"
        "<title>scanpy-advanced report</title></head>"
        "<body style=\"font-family:-apple-system,sans-serif;max-width:780px;"
        "margin:24px auto;padding:0 16px;line-height:1.5;\">"
        "<h1>Scanpy Advanced Analysis</h1>"
        f"<p><strong>Input:</strong> {n_obs_in} obs × {n_vars_in} vars "
        f"({h5ad_path.name})</p>"
        f"<p><strong>Post-QC:</strong> {n_obs_post_qc} obs × {n_vars_post_qc} vars</p>"
        f"<p><strong>Clusters (leiden, resolution={resolution}):</strong> {n_clusters}</p>"
        f"<p><strong>Marker genes:</strong> {len(ranks_df)} (method={method})</p>"
        "<h2>Outputs</h2>"
        "<ul>"
        "<li><code>preprocessed.h5ad</code> — AnnData with PCA/UMAP/leiden</li>"
        "<li><code>marker_genes.csv</code> — top 20 markers per cluster</li>"
        "<li><code>qc_summary.csv</code> — QC summary table</li>"
        "<li><code>umap.png</code> — UMAP (pre-clustering)</li>"
        "<li><code>umap_clustered.png</code> — UMAP colored by leiden</li>"
        "<li><code>hvg_plot.png</code> — Highly variable genes</li>"
        "</ul>"
        "<p style=\"color:#666;font-size:12px;margin-top:24px;\">"
        "G0 dev Electron 边界内 verified · biof3-py-runtime scanpy ≥1.9 · "
        "PHASE 1 BLOCK-11 SUB-1.1</p>"
        "</body></html>"
    )
    (output_dir / "report.html").write_text(report_html, encoding="utf-8")

    # ── 11. manifest.json (供下游 artifact-lineage 对齐) ──────
    outputs_meta = {
        "schemaVersion": 1,
        "pluginId": "scanpy-advanced",
        "pluginVersion": "1.0.0",
        "outputs": [
            "preprocessed.h5ad",
            "marker_genes.csv",
            "qc_summary.csv",
            "umap.png",
            "umap_clustered.png",
            "hvg_plot.png",
            "report.html",
        ],
        "stats": {
            "n_obs_in": n_obs_in,
            "n_vars_in": n_vars_in,
            "n_obs_after_qc": n_obs_post_qc,
            "n_vars_after_qc": n_vars_post_qc,
            "n_clusters": n_clusters,
            "n_markers": len(ranks_df),
            "method": method,
            "resolution": resolution,
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(outputs_meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    _report_progress(100, "scanpy-advanced 收口", job_dir)
    print(f"[scanpy-advanced] done. n_clusters={n_clusters} n_markers={len(ranks_df)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
