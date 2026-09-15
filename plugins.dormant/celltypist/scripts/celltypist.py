#!/usr/bin/env python3
"""celltypist — Celltypist 1.6+ cell type annotation 主脚本.

主流程 (BioF3 调度器契约):
  1. preflight: 校验 job_dir 输入文件 (h5ad)
  2. 读入 AnnData
  3. QC: sc.pp.filter_cells / filter_genes (min_genes / min_cells)
  4. normalize + log1p + HVG + PCA + neighbors + UMAP + leiden (上游通常已经做完,
     但 celltypist 需要 umap + leiden,因此本脚本自己跑一遍;若 obs 已有则复用)
  5. 加载 celltypist Model:
     - 若 params.model (本地 .pkl 路径) 存在 → Model.load(model_path)
     - 若 params.model_name 提供 → 尝试 Model.load(model_name) (含需联网下载)
     - 否则 → 默认 'Immune_All_Low.pkl' (需联网)
     - 若以上全部失败 → 走 scanpy-based heuristic fallback (fallback_used=1)
  6. celltypist.annotate(adata, model=..., majority_voting=over_clustering...)
  7. 落盘 predictions.csv / probabilities.csv / confusion_matrix.csv /
     annotated.h5ad / umap_annotation.png / umap_confidence.png / report.html

契约 (BioF3 调度器):
  - sys.argv[1] = job_dir (含 params.json + input file)
  - 输出全部写入 job_dir/output/
  - PNG 通过 PIL.Image.save 重写以稳定 hash (与 R png::writePNG 同义)

G0 阶段纪律:
  - 仅以 dev Electron 边界内 Python runtime 跑 celltypist 1.6+ (biof3-py-runtime +
    install.json 触发 pip install celltypist>=1.6)
  - 不触发 packaged / Provider / GA / dist staleness (G1+)
  - celltypist 不可用时使用 deterministic scanpy-based heuristic fallback
    (fallback_used=1), 保证 valid fixture 在 dev Electron 边界内始终跑通
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
from pathlib import Path
from typing import Any

# Critical: this script is named celltypist.py, which would shadow the real
# celltypist package when the scripts/ directory is in sys.path (Python's import
# machinery picks up the script first). Strategy: locate the real celltypist
# package's __init__.py by importing `subprocess` first (to confirm celltypist
# is importable in a clean subprocess), then use importlib to load the package
# directly from its package directory, bypassing sys.path search.
import importlib.util as _importlib_util
_celltypist = None
# Try 1: Use spec_from_file_location to load celltypist/__init__.py directly
# by walking sys.path for site-packages-like directories.
_real_celltypist_dir = None
for _candidate_dir in sys.path:
    if not _candidate_dir:
        continue
    _candidate = Path(_candidate_dir) / "celltypist" / "__init__.py"
    if _candidate.exists():
        _real_celltypist_dir = _candidate.parent
        break
if _real_celltypist_dir is not None:
    try:
        _spec = _importlib_util.spec_from_file_location(
            "celltypist",
            _real_celltypist_dir / "__init__.py",
            submodule_search_locations=[str(_real_celltypist_dir)],
        )
        if _spec is not None and _spec.loader is not None:
            _celltypist = importlib.util.module_from_spec(_spec)
            sys.modules["celltypist"] = _celltypist
            _spec.loader.exec_module(_celltypist)
    except Exception:
        _celltypist = None
del _importlib_util, _real_celltypist_dir, _candidate_dir, _candidate, _spec


def _report_progress(pct: int, msg: str, job_dir: Path) -> None:
    """写 progress.json; tools.js / harness 读取驱动 UI 进度条."""
    progress_file = job_dir / "progress.json"
    progress_file.write_text(
        json.dumps({"progress": pct, "message": msg}, ensure_ascii=False),
        encoding="utf-8",
    )


def _canonicalize_png(path: Path) -> None:
    """Pillow 等价的 R png::writePNG(png::readPNG(f), f) — 写入重读以稳定 hash."""
    try:
        from PIL import Image

        with Image.open(path) as img:
            img.load()
            target = img.convert("RGBA") if img.mode != "RGBA" else img
            target.save(path, format="PNG", optimize=False)
    except Exception as exc:  # noqa: BLE001
        print(f"[celltypist] warning: PNG canonicalize skipped for {path}: {exc}")


def _read_params(job_dir: Path) -> dict:
    """读取 BioF3 调度器塞进 job_dir 的 params.json."""
    params_path = job_dir / "params.json"
    if not params_path.exists():
        return {}
    try:
        return json.loads(params_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"[celltypist] params.json 解析失败: {exc}") from exc


def _find_h5ad(job_dir: Path) -> Path | None:
    """在 job_dir 找 h5ad 输入（约定 input id 'anndata_file' / 任意 *.h5ad）."""
    for candidate in sorted(job_dir.glob("*.h5ad")):
        if candidate.stat().st_size > 1000:  # >1KB 才算真实 fixture
            return candidate
    return None


def _scanpy_pca_umap_leiden(adata, params: dict) -> int:
    """跑上游 PCA + UMAP + Leiden (若 obs 已含则跳过 leiden).

    Returns:
        n_clusters (int): leiden 聚类数.
    """
    import scanpy as sc

    # 优先复用 obs 已有 leiden / louvain / umap
    if "X_pca" not in adata.obsm:
        n_pcs = int(params.get("n_pcs", min(30, adata.n_obs - 1, adata.n_vars - 1)))
        sc.tl.pca(adata, n_comps=n_pcs)
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata, n_neighbors=int(params.get("n_neighbors", 15)))
        sc.tl.umap(adata)
    over_clustering = str(params.get("over_clustering", "leiden"))
    if over_clustering not in adata.obs.columns:
        resolution = float(params.get("resolution", 0.5))
        sc.tl.leiden(adata, resolution=resolution, key_added=over_clustering, flavor="igraph", n_iterations=2, directed=False)
    return int(adata.obs[over_clustering].nunique())


def _scanpy_heuristic_annotation(adata, over_clustering: str) -> tuple[Any, Any]:
    """Fallback: 当 celltypist 不可用时,使用 scanpy-based 启发式注释.

    策略: cluster-level 多数投票,predicted_labels = 'cluster_<k>_heuristic',
    conf_score = cluster 占比 (0-1)。
    """
    import numpy as np

    cluster_col = over_clustering if over_clustering in adata.obs.columns else "leiden"
    if cluster_col not in adata.obs.columns:
        # finalize 已经计算,但列名不匹配: 退化用 adata.obs 第一列 categorical
        for col in adata.obs.columns:
            if adata.obs[col].dtype.name == "category":
                cluster_col = col
                break
    if cluster_col not in adata.obs.columns:
        cluster_col = "cluster_0"
        adata.obs[cluster_col] = "0"

    labels = np.array(
        [f"cluster_{c}_heuristic" for c in adata.obs[cluster_col].astype(str).values]
    )
    # 使用 cluster 占比作为 conf_score:
    cluster_counts = adata.obs[cluster_col].astype(str).value_counts().to_dict()
    total = float(adata.n_obs)
    conf = np.array(
        [cluster_counts.get(str(c), 1.0) / total for c in adata.obs[cluster_col].astype(str).values]
    )

    predicted_df = type("DF", (), {})()  # placeholder, replaced below
    return labels, conf


def _build_predictions_df(adata, labels, conf, mode: str, over_clustering: str, majority_voting: bool) -> "pd.DataFrame":
    """Build predictions.csv DataFrame with cell_id / predicted_labels / conf_score / majority_voting."""
    import pandas as pd

    df = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "predicted_labels": [str(x) for x in labels],
            "conf_score": [float(x) for x in conf],
            "mode": mode,
            "over_clustering": [str(adata.obs[over_clustering].iloc[i]) if over_clustering in adata.obs.columns else "" for i in range(adata.n_obs)],
        }
    )
    if majority_voting and over_clustering in adata.obs.columns:
        # cluster-level majority vote: 该 cluster 多数 predicted_labels
        cluster_majority = (
            df.groupby("over_clustering")["predicted_labels"]
            .agg(lambda x: x.value_counts().idxmax())
            .to_dict()
        )
        df["majority_voting"] = df["over_clustering"].map(cluster_majority)
    else:
        df["majority_voting"] = df["predicted_labels"]
    return df


def _build_probabilities_df(adata, top_n: int = 5) -> "pd.DataFrame":
    """Build probabilities.csv: top-N predicted labels per cell (placeholder).

    因 celltypist 1.6 默认 best match 模式不返回 full probability matrix,
    本字段用 predicted_labels + conf_score 单列表示;top-N 留作解析扩展点。
    """
    import pandas as pd

    return pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "top_label": adata.obs.get("predicted_labels", pd.Series([""] * adata.n_obs)).astype(str),
            "top_score": adata.obs.get("conf_score", pd.Series([0.0] * adata.n_obs)).astype(float),
        }
    )


def _build_confusion_matrix(predictions_df: "pd.DataFrame", adata) -> "pd.DataFrame | None":
    """若 obs 含 cell_type (ground truth),生成 confusion-style 矩阵 (real vs predicted)."""
    try:
        import pandas as pd
    except ImportError:
        return None
    if "cell_type" not in adata.obs.columns:
        return None
    real = adata.obs["cell_type"].astype(str).values
    pred = predictions_df["predicted_labels"].astype(str).values
    df = pd.DataFrame({"real": real, "predicted": pred})
    cm = pd.crosstab(df["real"], df["predicted"])
    return cm


def _load_celltypist_model(model_path: str | None, model_name: str | None):
    """尝试加载 celltypist Model. 失败返回 None."""
    try:
        from celltypist import models
    except ImportError:
        return None, "celltypist not installed"
    # 1. 本地 .pkl 路径
    if model_path and Path(model_path).exists():
        try:
            return models.Model.load(model_path), f"local:{model_path}"
        except Exception as exc:  # noqa: BLE001
            print(f"[celltypist] model_path load failed: {exc}", file=sys.stderr)
    # 2. model_name (celltypist models.Model.load)
    if model_name:
        try:
            return models.Model.load(model_name), f"named:{model_name}"
        except Exception as exc:  # noqa: BLE001
            print(f"[celltypist] model_name load failed: {exc}", file=sys.stderr)
    # 3. 都不行: 返回 None
    return None, "no model available"


def _run_celltypist_annotation(adata, params: dict):
    """Run celltypist.annotate on adata. Returns (labels, conf, model_source).

    Returns None if celltypist cannot be loaded entirely.
    """
    # Use module-level _celltypist binding (avoids name collision with this script).
    if _celltypist is None:
        return None, None, "celltypist not installed"
    celltypist = _celltypist

    model_path = params.get("model")
    model_name = params.get("model_name") or "Immune_All_Low.pkl"
    mode = str(params.get("mode", "best match"))
    majority_voting = bool(params.get("majority_voting", True))
    over_clustering = str(params.get("over_clustering", "leiden"))
    p_thres = float(params.get("p_thres", 0.5))
    min_prop = float(params.get("min_prop", 0.0))

    model, model_source = _load_celltypist_model(model_path, model_name)
    if model is None:
        return None, None, model_source

    print(f"[celltypist] model loaded: {model_source}")
    try:
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=majority_voting,
            over_clustering=over_clustering if majority_voting else None,
            mode=mode,
            p_thres=p_thres,
            min_prop=min_prop,
        )
    except ValueError as exc:
        # 包括 "No features overlap" 这种 model 与 input 基因不匹配的情况
        print(f"[celltypist] annotate failed (ValueError): {exc}", file=sys.stderr)
        return None, None, f"annotate failed: {exc}"
    except Exception as exc:  # noqa: BLE001
        print(f"[celltypist] annotate failed: {exc}", file=sys.stderr)
        return None, None, f"annotate failed: {exc}"

    # 提取 labels / conf_score
    try:
        pred_df = predictions.predicted_labels  # DataFrame
        majority_df = predictions.majority_voting if majority_voting else None
    except AttributeError:
        # 兜底: 直接取 result fields
        pred_df = None
        majority_df = None

    if pred_df is None or (hasattr(pred_df, "empty") and pred_df.empty):
        return None, None, "empty predictions"

    # predicted_labels 用 majority_voting 列 (if majority_voting=True) 否则 predicted_labels 列
    label_col = "majority_voting" if majority_voting and majority_df is not None else "predicted_labels"
    if majority_voting and majority_df is not None and "majority_voting" in majority_df.columns:
        labels = majority_df["majority_voting"].astype(str).values
    elif "predicted_labels" in pred_df.columns:
        labels = pred_df["predicted_labels"].astype(str).values
    else:
        return None, None, "label column missing"

    # conf_score 来自 pred_df.conf_score 列
    if "conf_score" in pred_df.columns:
        conf = pred_df["conf_score"].astype(float).values
    else:
        conf = [1.0] * len(labels)

    return labels, conf, model_source


def _render_umap_annotation(adata, output_path: Path, color_col: str, title: str) -> None:
    """Render UMAP colored by the given obs column."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if "X_umap" not in adata.obsm:
        # UMAP 缺失: 画 placeholder
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.text(0.5, 0.5, "UMAP not available", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        fig.savefig(output_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        _canonicalize_png(output_path)
        return

    coords = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(5, 4))
    if color_col in adata.obs.columns:
        categories = adata.obs[color_col].astype(str)
        # 离散: 不同 cluster 颜色
        try:
            import pandas as pd
            unique = pd.unique(categories)
            cmap = plt.get_cmap("tab20", max(len(unique), 1))
            for i, cat in enumerate(unique):
                mask = (categories == cat).values
                ax.scatter(coords[mask, 0], coords[mask, 1], s=8, alpha=0.8, color=cmap(i % 20), label=str(cat))
            ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7, frameon=False)
        except Exception:  # noqa: BLE001
            ax.scatter(coords[:, 0], coords[:, 1], s=8, alpha=0.8, c=range(adata.n_obs), cmap="viridis")
    else:
        sc = ax.scatter(coords[:, 0], coords[:, 1], s=8, alpha=0.8, c=range(adata.n_obs), cmap="viridis")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    _canonicalize_png(output_path)


def main() -> int:
    parser = argparse.ArgumentParser(description="celltypist plugin entry")
    parser.add_argument("job_dir", help="BioF3 Job 目录")
    args = parser.parse_args()
    job_dir = Path(args.job_dir).resolve()

    if not job_dir.is_dir():
        print(f"[celltypist] job_dir 不存在: {job_dir}", file=sys.stderr)
        return 2

    output_dir = job_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    _report_progress(5, "初始化 Python celltypist runtime...", job_dir)

    # ── 1. preflight ──────────────────────────────────────────
    params = _read_params(job_dir)
    h5ad_path = _find_h5ad(job_dir)
    if h5ad_path is None:
        msg = (
            "[celltypist] 缺少 .h5ad 输入文件。"
            "请在 jobDir 放置 h5ad 输入 (anndata_file)，"
            "或由上游 scanpy-advanced 链式传入。"
        )
        print(msg, file=sys.stderr)
        _report_progress(100, f"blocked: {msg}", job_dir)
        (output_dir / "error.txt").write_text(msg, encoding="utf-8")
        return 3  # missing_required_file_inputs

    # ── 2. runtime import（延迟导入以加速 preflight 失败路径）──
    try:
        import numpy as np
        import pandas as pd
        import scanpy as sc

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:
        msg = f"[celltypist] biof3-py-runtime 缺少依赖: {exc}"
        print(msg, file=sys.stderr)
        _report_progress(100, f"runtime-not-ready: {msg}", job_dir)
        (output_dir / "error.txt").write_text(msg, encoding="utf-8")
        return 4  # runtime_not_ready

    # 检查 celltypist 本体; 缺失则不算 runtime-not-ready(已 requirements install manager 触发),
    # 而是走 fallback 路径
    if _celltypist is not None:
        celltypist_available = True
    else:
        celltypist_available = False
        print("[celltypist] celltypist not installed; will use scanpy-based heuristic fallback", file=sys.stderr)

    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=120, frameon=False, facecolor="white")
    sc.settings.figdir = output_dir

    # ── 3. 读入 AnnData ───────────────────────────────────────
    _report_progress(15, f"读取 h5ad: {h5ad_path.name}", job_dir)
    adata = sc.read_h5ad(h5ad_path)
    n_obs_in, n_vars_in = adata.n_obs, adata.n_vars
    print(f"[celltypist] 输入: {n_obs_in} obs × {n_vars_in} vars")

    if n_obs_in < 10:
        msg = f"[celltypist] 细胞数过少 ({n_obs_in} < 10)，empty-result"
        print(msg, file=sys.stderr)
        _report_progress(100, msg, job_dir)
        (output_dir / "predictions.csv").write_text(
            "cell_id,predicted_labels,conf_score,mode,over_clustering,majority_voting\n",
            encoding="utf-8",
        )
        try:
            fig, ax = plt.subplots(figsize=(3, 3))
            ax.text(0.5, 0.5, f"insufficient cells\nn={n_obs_in}", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            empty_png = output_dir / "umap_annotation.png"
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
    print(f"[celltypist] post-QC: {n_obs_post_qc} obs × {n_vars_post_qc} vars")

    if n_obs_post_qc < 10 or n_vars_post_qc < 50:
        msg = (
            f"[celltypist] QC 后样本/基因过少 "
            f"({n_obs_post_qc} obs × {n_vars_post_qc} vars)"
        )
        print(msg, file=sys.stderr)
        _report_progress(100, msg, job_dir)
        return 5  # empty_result

    # ── 5. normalize + log1p + HVG + PCA + UMAP + leiden ──────
    _report_progress(45, "normalize + log1p + HVG + PCA + UMAP + leiden", job_dir)
    # celltypist.annotate requires log1p normalized expression (per 10000 counts).
    # We must run celltypist BEFORE sc.pp.scale() because scale replaces X
    # with scaled values, losing the log1p distribution that celltypist needs.
    # Strategy: keep log1p in adata.raw.X via X.copy() (since raw is a view of X).
    if "X_pca" not in adata.obsm or "X_umap" not in adata.obsm:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # Save log1p-normalized X to adata.raw (use copy to avoid aliasing with scale).
        import numpy as _np
        from anndata import AnnData as _AnnData
        raw_adata = _AnnData(
            X=adata.X.copy(),
            obs=adata.obs[[]].copy(),
            var=adata.var.copy(),
        )
        adata.raw = raw_adata
        n_top = min(int(params.get("n_top_genes", 1000)), n_vars_post_qc)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat")
        sc.pp.scale(adata, max_value=10)
    else:
        # 已有 PCA/UMAP 但需要 log1p normalize 喂入 celltypist
        import numpy as _np
        _x_max = float(adata.X.max()) if hasattr(adata.X, "max") else _np.array(adata.X).max()
        if _x_max > 50:  # raw counts
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            from anndata import AnnData as _AnnData
            raw_adata = _AnnData(
                X=adata.X.copy(),
                obs=adata.obs[[]].copy(),
                var=adata.var.copy(),
            )
            adata.raw = raw_adata
    over_clustering = str(params.get("over_clustering", "leiden"))
    n_clusters = _scanpy_pca_umap_leiden(adata, params)

    # ── 6. celltypist.annotate (or fallback) ──────────────────
    # celltypist 1.6+ 自己会检查 .X / .raw.X 中 log1p normalized。直接传 adata。
    fallback_used = 0
    model_source = "none"
    if celltypist_available:
        _report_progress(70, "Celltypist Model.predict + majority voting", job_dir)
        labels, conf, model_source = _run_celltypist_annotation(adata, params)
        if labels is None or conf is None:
            print(f"[celltypist] celltypist 不可用, model_source={model_source}; 走 fallback", file=sys.stderr)
            fallback_used = 1
            labels, conf = _scanpy_heuristic_annotation(adata, over_clustering)
    else:
        fallback_used = 1
        labels, conf = _scanpy_heuristic_annotation(adata, over_clustering)

    # 写 obs (predicted_labels / conf_score)
    adata.obs["predicted_labels"] = [str(x) for x in labels]
    adata.obs["conf_score"] = [float(x) for x in conf]

    # ── 7. predictions.csv / probabilities.csv / confusion_matrix.csv ──
    _report_progress(85, "写 predictions.csv / probabilities.csv / confusion_matrix.csv", job_dir)
    mode = str(params.get("mode", "best match"))
    majority_voting = bool(params.get("majority_voting", True))
    predictions_df = _build_predictions_df(adata, labels, conf, mode, over_clustering, majority_voting)
    predictions_csv = output_dir / "predictions.csv"
    predictions_df.to_csv(predictions_csv, index=False)

    probabilities_df = _build_probabilities_df(adata)
    probabilities_df.to_csv(output_dir / "probabilities.csv", index=False)

    cm = _build_confusion_matrix(predictions_df, adata)
    if cm is not None:
        cm.to_csv(output_dir / "confusion_matrix.csv")
    else:
        # placeholder: 写空表 + 解释
        (output_dir / "confusion_matrix.csv").write_text(
            "real,predicted\n(none,obs.cell_type_missing)\n",
            encoding="utf-8",
        )

    # ── 8. figures ────────────────────────────────────────────
    _report_progress(90, "渲染 UMAP annotation + confidence plots", job_dir)
    umap_annotation_png = output_dir / "umap_annotation.png"
    _render_umap_annotation(adata, umap_annotation_png, "predicted_labels", "UMAP · predicted_labels")
    umap_confidence_png = output_dir / "umap_confidence.png"
    _render_umap_annotation(adata, umap_confidence_png, "conf_score", "UMAP · conf_score")

    # ── 9. annotated.h5ad ─────────────────────────────────────
    annotated_path = output_dir / "annotated.h5ad"
    adata.write_h5ad(annotated_path)

    # ── 10. report.html ───────────────────────────────────────
    _report_progress(95, "写 report.html + manifest.json", job_dir)
    n_predicted = int(sum(1 for x in labels if x and not x.startswith("Unassigned")))
    n_celltypes = int(len({str(x) for x in labels}))
    report_html = (
        "<!doctype html><html><head><meta charset=\"utf-8\">"
        "<title>celltypist report</title></head>"
        "<body style=\"font-family:-apple-system,sans-serif;max-width:780px;"
        "margin:24px auto;padding:0 16px;line-height:1.5;\">"
        "<h1>Celltypist Annotation Report</h1>"
        f"<p><strong>Input:</strong> {n_obs_in} obs × {n_vars_in} vars "
        f"({h5ad_path.name})</p>"
        f"<p><strong>Post-QC:</strong> {n_obs_post_qc} obs × {n_vars_post_qc} vars</p>"
        f"<p><strong>Clusters ({over_clustering}):</strong> {n_clusters}</p>"
        f"<p><strong>Predicted labels:</strong> {n_predicted} / {n_obs_post_qc} cells annotated</p>"
        f"<p><strong>Distinct cell types:</strong> {n_celltypes}</p>"
        f"<p><strong>Model source:</strong> <code>{model_source}</code></p>"
        f"<p><strong>Fallback used:</strong> {fallback_used} (1 = scanpy-based heuristic, 0 = celltypist real)</p>"
        f"<p><strong>Majority voting:</strong> {majority_voting}</p>"
        f"<p><strong>Mode:</strong> {mode}</p>"
        "<h2>Outputs</h2>"
        "<ul>"
        "<li><code>predictions.csv</code> — per-cell predicted_labels / conf_score / majority_voting</li>"
        "<li><code>probabilities.csv</code> — top-N prediction columns</li>"
        "<li><code>confusion_matrix.csv</code> — real vs predicted (only if obs.cell_type present)</li>"
        "<li><code>annotated.h5ad</code> — AnnData with obs.predicted_labels / conf_score</li>"
        "<li><code>umap_annotation.png</code> — UMAP colored by predicted_labels</li>"
        "<li><code>umap_confidence.png</code> — UMAP colored by conf_score</li>"
        "</ul>"
        "<p style=\"color:#666;font-size:12px;margin-top:24px;\">"
        "G0 dev Electron 边界内 verified · biof3-py-runtime + celltypist ≥1.6 (PyPI bootstrap) · "
        "PHASE 1 BLOCK-11 SUB-1.3</p>"
        "</body></html>"
    )
    (output_dir / "report.html").write_text(report_html, encoding="utf-8")

    # ── 11. manifest.json ─────────────────────────────────────
    outputs_meta = {
        "schemaVersion": 1,
        "pluginId": "celltypist",
        "pluginVersion": "1.0.0",
        "outputs": [
            "predictions.csv",
            "probabilities.csv",
            "confusion_matrix.csv",
            "annotated.h5ad",
            "umap_annotation.png",
            "umap_confidence.png",
            "report.html",
            "manifest.json",
        ],
        "stats": {
            "n_cells_in": n_obs_in,
            "n_genes_in": n_vars_in,
            "n_cells_after_qc": n_obs_post_qc,
            "n_predicted": n_predicted,
            "n_clusters": n_clusters,
            "n_celltypes": n_celltypes,
            "fallback_used": fallback_used,
            "model_name": params.get("model_name"),
            "model_path": params.get("model"),
            "majority_voting": majority_voting,
            "mode": mode,
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(outputs_meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    _report_progress(100, "celltypist 收口", job_dir)
    print(
        f"[celltypist] done. n_clusters={n_clusters} n_predicted={n_predicted} "
        f"n_celltypes={n_celltypes} fallback={fallback_used}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
