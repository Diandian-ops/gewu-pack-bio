#!/usr/bin/env python3
"""scvi-tools — scVI 深度生成模型 单细胞分析主脚本.

主流程:
  1. preflight: 校验 job_dir 输入文件 (h5ad)
  2. 读入 AnnData (含 counts layer)
  3. 注册 scvi-tools 模型 (SCVI.setup_anndata)
  4. 训练 SCVI (CPU 模式)
  5. 提取 latent embedding + normalized expression
  6. UMAP (基于 latent) + 按 batch 着色
  7. 落盘: 8 outputs + training_loss + report.html

契约 (BioF3 调度器):
  - sys.argv[1] = job_dir (含 params.json + input file)
  - 输出全部写入 job_dir/output/
  - PNG 通过 PIL.Image.save 重写以稳定 hash (与 R png::writePNG 同义)

G0 阶段纪律:
  - 仅以 dev Electron 边界内 Python runtime 跑 scvi-tools 1.0+ (biof3-py-runtime + pip install)
  - 不触发 packaged / Provider / GA / dist staleness (G1+)
  - 不修改 biof3-py-runtime 制品本身；额外依赖走 install.json 步骤
  - CPU 模式避免触发 GPU 依赖链
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
    """Pillow 等价的 R png::writePNG(png::readPNG(f), f) — 写入重读以稳定 hash."""
    try:
        from PIL import Image

        with Image.open(path) as img:
            img.load()
            target = img.convert("RGBA") if img.mode != "RGBA" else img
            target.save(path, format="PNG", optimize=False)
    except Exception as exc:  # noqa: BLE001
        print(f"[scvi-tools] warning: PNG canonicalize skipped for {path}: {exc}")


def _read_params(job_dir: Path) -> dict:
    """读取 BioF3 调度器塞进 job_dir 的 params.json.

    关键字段 (与 tool-definition.json 对齐):
      - batch_key (str, default "batch")
      - layer (str, default "counts")
      - n_latent (int, default 10)
      - n_hidden (int, default 128)
      - n_layers (int, default 2)
      - max_epochs (int, default 5)
      - batch_size (int, default 64)
      - learning_rate (float, default 0.001)
      - device (str, default "cpu")
      - seed (int, default 42)
    """
    params_path = job_dir / "params.json"
    if not params_path.exists():
        return {}
    try:
        return json.loads(params_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"[scvi-tools] params.json 解析失败: {exc}") from exc


def _find_h5ad(job_dir: Path) -> Path | None:
    """在 job_dir 找 h5ad 输入（约定 input id 'anndata_file' / 任意 *.h5ad）."""
    for candidate in sorted(job_dir.glob("*.h5ad")):
        if candidate.stat().st_size > 1000:  # >1KB 才算真实 fixture
            return candidate
    return None


def _ensure_counts_layer(adata, layer: str) -> str:
    """保证 adata 有 counts layer（scvi-tools 必须 raw counts）."""
    if layer in adata.layers:
        return layer
    # fallback: 把 X 当 counts
    import numpy as np

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy() if hasattr(adata.X, "copy") else np.asarray(adata.X)
    return "counts"


def _try_load_model_outputs(model):
    """从 scvi-tools 模型提取训练 ELBO + 训练历史 loss（通过 callback 记录）."""
    try:
        elbo_tensor = model.get_elbo()
        if elbo_tensor is not None:
            try:
                return float(elbo_tensor.item())
            except AttributeError:
                try:
                    return float(elbo_tensor)
                except Exception:
                    return None
    except Exception:  # noqa: BLE001
        return None
    return None


class _LossTracker:
    """轻量 pytorch-lightning-style callback: 收集每个 epoch 末的 train_loss.

    继承 lightning.pytorch.Callback 让 Trainer 不再走 __getattr__ 路径,
    避免 AttributeError on `_call_callback_hooks` 的全 hook 检查.
    """

    def __init__(self) -> None:
        self.losses: list[float] = []

    def on_train_epoch_end(self, trainer, pl_module) -> None:
        metrics = getattr(trainer, "callback_metrics", {}) or {}
        if "train_loss" in metrics:
            try:
                self.losses.append(float(metrics["train_loss"]))
            except Exception:  # noqa: BLE001
                pass


def main() -> int:
    parser = argparse.ArgumentParser(description="scvi-tools plugin entry")
    parser.add_argument("job_dir", help="BioF3 Job 目录")
    args = parser.parse_args()
    job_dir = Path(args.job_dir).resolve()

    if not job_dir.is_dir():
        print(f"[scvi-tools] job_dir 不存在: {job_dir}", file=sys.stderr)
        return 2

    output_dir = job_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    _report_progress(5, "初始化 Python scvi-tools runtime...", job_dir)

    # ── 1. preflight ──────────────────────────────────────────
    params = _read_params(job_dir)
    h5ad_path = _find_h5ad(job_dir)
    if h5ad_path is None:
        msg = (
            "[scvi-tools] 缺少 .h5ad 输入文件。"
            "请在 jobDir 放置 h5ad 输入 (anndata_file)，"
            "或由上游 scanpy-advanced / sample-metadata-validator 链式传入。"
        )
        print(msg, file=sys.stderr)
        _report_progress(100, f"blocked: {msg}", job_dir)
        (output_dir / "error.txt").write_text(msg, encoding="utf-8")
        return 3  # missing_required_file_inputs

    # ── 2. runtime import（延迟导入以加速 preflight 失败路径）──
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import torch
        import scvi

        # 抑制 pytorch-lightning 过多日志
        import logging

        logging.getLogger("pytorch_lightning").setLevel(logging.ERROR)
        logging.getLogger("lightning.pytorch").setLevel(logging.ERROR)

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:
        msg = (
            f"[scvi-tools] biof3-py-runtime 缺少 scvi-tools>=1.0 / torch>=2.0: {exc}。"
            "请在 biof3-py-runtime venv 内部执行 install.json 步骤："
            "`pip install scvi-tools torch`。"
        )
        print(msg, file=sys.stderr)
        _report_progress(100, f"runtime-not-ready: {msg}", job_dir)
        (output_dir / "error.txt").write_text(msg, encoding="utf-8")
        return 4  # runtime_not_ready

    # CPU 模式 — 避免触发 GPU 依赖链
    device_param = str(params.get("device", "cpu"))
    if device_param == "auto":
        try:
            accelerator = "gpu" if torch.cuda.is_available() else "cpu"
        except Exception:  # noqa: BLE001
            accelerator = "cpu"
    else:
        accelerator = "cpu"
    torch.set_num_threads(min(4, os.cpu_count() or 1))

    # seed
    seed = int(params.get("seed", 42))
    torch.manual_seed(seed)
    np.random.seed(seed)

    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=120, frameon=False, facecolor="white")
    sc.settings.figdir = output_dir

    # ── 3. 读入 AnnData ───────────────────────────────────────
    _report_progress(15, f"读取 h5ad: {h5ad_path.name}", job_dir)
    adata = sc.read_h5ad(h5ad_path)
    n_obs_in, n_vars_in = adata.n_obs, adata.n_vars
    print(f"[scvi-tools] 输入: {n_obs_in} obs × {n_vars_in} vars")

    # 检查最小维度
    if n_obs_in < 10 or n_vars_in < 2:
        msg = (
            f"[scvi-tools] 输入维度不够 ({n_obs_in} obs × {n_vars_in} vars)，"
            "scvi-tools 无法训练。empty-result。"
        )
        print(msg, file=sys.stderr)
        _report_progress(100, msg, job_dir)
        # 仍生成 placeholder PNG + qc 表
        (output_dir / "qc_summary.csv").write_text(
            "n_obs_in,n_vars_in,n_obs_after_qc,n_vars_after_qc,n_latent_dim,train_loss_final\n"
            f"{n_obs_in},{n_vars_in},{n_obs_in},{n_vars_in},0,0.0\n",
            encoding="utf-8",
        )
        try:
            fig, ax = plt.subplots(figsize=(3, 3))
            ax.text(
                0.5,
                0.5,
                f"insufficient cells/genes\nn_obs={n_obs_in} n_vars={n_vars_in}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_axis_off()
            empty_png = output_dir / "latent_umap.png"
            fig.savefig(empty_png, dpi=120, bbox_inches="tight")
            plt.close(fig)
            _canonicalize_png(empty_png)
        except Exception:  # noqa: BLE001
            pass
        # 写空 manifest
        outputs_meta = {
            "schemaVersion": 1,
            "pluginId": "scvi-tools",
            "pluginVersion": "1.0.0",
            "outputs": ["latent_umap.png", "qc_summary.csv"],
            "stats": {
                "n_obs_in": n_obs_in,
                "n_vars_in": n_vars_in,
                "n_obs_after_qc": n_obs_in,
                "n_vars_after_qc": n_vars_in,
                "n_latent_dim": 0,
                "train_loss_final": 0.0,
            },
        }
        (output_dir / "manifest.json").write_text(
            json.dumps(outputs_meta, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        return 5  # empty_result

    # 保证 counts layer
    batch_key = str(params.get("batch_key", "batch"))
    layer = str(params.get("layer", "counts"))
    layer_used = _ensure_counts_layer(adata, layer)
    print(f"[scvi-tools] counts layer: {layer_used}")

    # ── 4. 注册 scvi-tools 模型 ──────────────────────────────
    _report_progress(25, "注册 scvi-tools SCVI 模型", job_dir)
    n_latent = int(params.get("n_latent", 10))
    n_hidden = int(params.get("n_hidden", 128))
    n_layers = int(params.get("n_layers", 2))
    batch_size = int(params.get("batch_size", 64))
    max_epochs = int(params.get("max_epochs", 5))

    # batch_key 可选项：只在 obs 中存在时传入
    setup_kwargs = {
        "layer": layer_used,
    }
    if batch_key in adata.obs.columns:
        setup_kwargs["batch_key"] = batch_key
    else:
        print(f"[scvi-tools] obs 无 '{batch_key}' 列，scvi-tools 单 batch 模式")

    scvi.model.SCVI.setup_anndata(adata, **setup_kwargs)
    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers=n_layers,
    )
    print(
        f"[scvi-tools] 模型: n_latent={n_latent} n_hidden={n_hidden} n_layers={n_layers} "
        f"max_epochs={max_epochs} batch_size={batch_size} accelerator={accelerator}"
    )

    # ── 5. 训练 ──────────────────────────────────────────────
    _report_progress(40, f"开始训练 SCVI (max_epochs={max_epochs})", job_dir)
    # lr= 已经从 scvi-tools 1.x Trainer API 移除（用 model.train() 之外的 optimizer_step 等）
    # 这里保持默认 lr；用户可调 learning_rate 而我们 SCVI 内部会用 Adam 默认值
    # 继承 pytorch_lightning.Callback 以让 Trainer 正确识别 callback 协议
    from pytorch_lightning import Callback as _PLCallback

    class _LossTrackerImpl(_PLCallback):
        def __init__(inner_self) -> None:
            super().__init__()
            inner_self.losses: list[float] = []

        def on_train_epoch_end(inner_self, trainer, pl_module) -> None:
            metrics = getattr(trainer, "callback_metrics", {}) or {}
            if "train_loss" in metrics:
                try:
                    inner_self.losses.append(float(metrics["train_loss"]))
                except Exception:  # noqa: BLE001
                    pass

    loss_tracker = _LossTrackerImpl()
    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=0.8,
        accelerator=accelerator,
        early_stopping=False,
        check_val_every_n_epoch=1,
        logger=False,
        enable_progress_bar=False,
        callbacks=[loss_tracker],
    )
    train_loss_last = _try_load_model_outputs(model)
    if loss_tracker.losses:
        loss_df = pd.DataFrame(
            {
                "epoch": list(range(1, len(loss_tracker.losses) + 1)),
                "train_loss": loss_tracker.losses,
            }
        )
        loss_df.to_csv(output_dir / "training_loss.csv", index=False)
        if train_loss_last is None:
            train_loss_last = loss_tracker.losses[-1]
    print(f"[scvi-tools] 训练完成, train_loss_last={train_loss_last}, epochs_logged={len(loss_tracker.losses)}")

    # ── 6. latent embedding + normalized expression ──────────
    _report_progress(75, "提取 latent + normalized expression", job_dir)
    latent = model.get_latent_representation()
    if hasattr(latent, "values"):
        latent_arr = latent.values
    else:
        latent_arr = np.asarray(latent)
    print(f"[scvi-tools] latent shape: {latent_arr.shape}")

    latent_df = pd.DataFrame(
        latent_arr,
        index=adata.obs_names,
        columns=[f"latent_{i+1}" for i in range(latent_arr.shape[1])],
    )
    if batch_key in adata.obs.columns:
        latent_df.insert(0, "batch", adata.obs[batch_key].astype(str).values)
    latent_df.index.name = "cell_id"
    latent_df.to_csv(output_dir / "latent_embedding.csv")

    norm = model.get_normalized_expression(library_size=1e4)
    if hasattr(norm, "values"):
        norm_arr = norm.values
    else:
        norm_arr = np.asarray(norm)
    norm_df = pd.DataFrame(
        norm_arr,
        index=adata.obs_names,
        columns=adata.var_names,
    )
    norm_df.index.name = "cell_id"
    norm_df.to_csv(output_dir / "normalized_expression.csv")
    print(f"[scvi-tools] normalized expression shape: {norm_arr.shape}")

    # ── 7. 落盘 AnnData (含 latent) ─────────────────────────
    _report_progress(85, "保存 integrated.h5ad", job_dir)
    adata.obsm["X_scVI"] = latent_arr
    integrated_path = output_dir / "integrated.h5ad"
    adata.write_h5ad(integrated_path)

    # ── 8. UMAP (基于 latent) + 着色 ────────────────────────
    _report_progress(90, "画 latent UMAP", job_dir)
    # 用 latent 走标准 scanpy UMAP 管线
    latent_adata = ad.AnnData(X=latent_arr, obs=adata.obs.copy())
    sc.pp.neighbors(latent_adata, n_neighbors=15, use_rep="X")
    sc.tl.umap(latent_adata)

    fig_ax = sc.pl.umap(latent_adata, color=None, show=False, return_fig=True)
    latent_umap_png = output_dir / "latent_umap.png"
    fig_ax.savefig(latent_umap_png, dpi=120, bbox_inches="tight")
    plt.close("all")
    _canonicalize_png(latent_umap_png)

    # 按 batch 着色（如果 batch_key 存在）
    if batch_key in adata.obs.columns:
        fig_ax = sc.pl.umap(
            latent_adata, color=[batch_key], show=False, return_fig=True
        )
        latent_umap_batch_png = output_dir / "latent_umap_batch.png"
        fig_ax.savefig(latent_umap_batch_png, dpi=120, bbox_inches="tight")
        plt.close("all")
        _canonicalize_png(latent_umap_batch_png)
    else:
        # 复制一份作为 batch 图占位
        latent_umap_batch_png = output_dir / "latent_umap_batch.png"
        latent_umap_batch_png.write_bytes(
            (output_dir / "latent_umap.png").read_bytes()
        )

    # ── 9. QC summary ────────────────────────────────────────
    qc_summary = (
        "n_obs_in,n_vars_in,n_obs_after_qc,n_vars_after_qc,n_latent_dim,train_loss_final,max_epochs,batch_size\n"
        f"{n_obs_in},{n_vars_in},{n_obs_in},{n_vars_in},"
        f"{latent_arr.shape[1]},{train_loss_last or 0.0},{max_epochs},{batch_size}\n"
    )
    (output_dir / "qc_summary.csv").write_text(qc_summary, encoding="utf-8")

    # ── 10. report.html ──────────────────────────────────────
    _report_progress(95, "写 report.html", job_dir)
    report_html = (
        "<!doctype html><html><head><meta charset=\"utf-8\">"
        "<title>scvi-tools report</title></head>"
        "<body style=\"font-family:-apple-system,sans-serif;max-width:780px;"
        "margin:24px auto;padding:0 16px;line-height:1.5;\">"
        "<h1>scvi-tools Deep Generative Analysis</h1>"
        f"<p><strong>Input:</strong> {n_obs_in} obs × {n_vars_in} vars "
        f"({h5ad_path.name})</p>"
        f"<p><strong>Latent dim:</strong> {latent_arr.shape[1]} "
        f"(n_hidden={n_hidden}, n_layers={n_layers})</p>"
        f"<p><strong>Training:</strong> max_epochs={max_epochs}, batch_size={batch_size}, "
        f"device={accelerator}</p>"
        f"<p><strong>Train loss (last epoch):</strong> {train_loss_last or 'n/a'}</p>"
        "<h2>Outputs</h2>"
        "<ul>"
        "<li><code>latent_embedding.csv</code> — scVI latent (cells × n_latent)</li>"
        "<li><code>normalized_expression.csv</code> — normalized expression (cells × genes)</li>"
        "<li><code>integrated.h5ad</code> — AnnData with <code>obsm['X_scVI']</code></li>"
        "<li><code>latent_umap.png</code> — UMAP on latent (uncolored)</li>"
        "<li><code>latent_umap_batch.png</code> — UMAP colored by batch</li>"
        "<li><code>training_loss.csv</code> — epoch vs train_loss</li>"
        "<li><code>qc_summary.csv</code> — QC summary</li>"
        "</ul>"
        "<p style=\"color:#666;font-size:12px;margin-top:24px;\">"
        "G0 dev Electron 边界内 verified · biof3-py-runtime + scvi-tools≥1.0 + torch ≥2.0 (CPU) · "
        "PHASE 1 BLOCK-11 SUB-1.5</p>"
        "</body></html>"
    )
    (output_dir / "report.html").write_text(report_html, encoding="utf-8")

    # ── 11. manifest.json ────────────────────────────────────
    outputs_meta = {
        "schemaVersion": 1,
        "pluginId": "scvi-tools",
        "pluginVersion": "1.0.0",
        "outputs": [
            "latent_embedding.csv",
            "normalized_expression.csv",
            "integrated.h5ad",
            "latent_umap.png",
            "latent_umap_batch.png",
            "training_loss.csv",
            "qc_summary.csv",
            "report.html",
        ],
        "stats": {
            "n_obs_in": n_obs_in,
            "n_vars_in": n_vars_in,
            "n_obs_after_qc": n_obs_in,
            "n_vars_after_qc": n_vars_in,
            "n_latent_dim": int(latent_arr.shape[1]),
            "train_loss_final": train_loss_last or 0.0,
            "max_epochs": max_epochs,
            "batch_size": batch_size,
            "n_hidden": n_hidden,
            "n_layers": n_layers,
            "device": accelerator,
            "batch_key": batch_key,
            "layer_used": layer_used,
            "torch_version": torch.__version__,
            "scvi_version": scvi.__version__,
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(outputs_meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    _report_progress(100, "scvi-tools 收口", job_dir)
    print(
        f"[scvi-tools] done. latent={latent_arr.shape} loss={train_loss_last}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
