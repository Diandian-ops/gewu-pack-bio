"""Generate the demo h5ad fixture for celltypist acceptance fixtures.

This script builds a deterministic 270 obs × 100 vars AnnData with 3 simulated
cell types (cell_type A/B/C in obs.cell_type column) — same shape as the
scanpy-advanced demo for cross-plugin consistency. The same file is intentionally
identical to scanpy-advanced/demo-data/pbmc_3k_mini.h5ad (sha256:
e3d6c0f0226479e52c5e8d8d03e2be392466766d25d07c326d0cd0982e2c35bc); we re-build
here with a separate codepath so celltypist fixture stays self-contained.

Run:
    /Applications/anaconda3/envs/biof3-py-runtime/bin/python3 resources/built-in-plugins/celltypist/demo-data/build_demo_h5ad.py

Idempotent: re-running produces byte-identical output (relies on deterministic
RNG seed=42, no time/pickle metadata).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def build() -> Path:
    """Build a deterministic 270 × 100 AnnData with 3 simulated cell types."""
    try:
        import anndata as ad
        import pandas as pd
    except ImportError as exc:
        print(f"[demo-data build] biof3-py-runtime missing dependency: {exc}", file=sys.stderr)
        print("[demo-data build] Run within an environment with anndata installed.", file=sys.stderr)
        sys.exit(1)

    rng = np.random.default_rng(seed=42)
    n_obs = 270
    n_vars = 100
    n_types = 3
    cells_per_type = n_obs // n_types

    # Simulate counts: 3 clusters with shifted mean on a subset of genes.
    X = rng.negative_binomial(5, 0.3, size=(n_obs, n_vars)).astype(np.float32)
    # bump cluster 1 mean on genes 0:30, cluster 2 on genes 30:60, cluster 0 baseline
    X[:cells_per_type, 0:30] += rng.negative_binomial(8, 0.3, size=(cells_per_type, 30)).astype(np.float32)
    X[cells_per_type : 2 * cells_per_type, 30:60] += rng.negative_binomial(8, 0.3, size=(cells_per_type, 30)).astype(
        np.float32
    )

    obs = pd.DataFrame(
        {
            "cell_type": (
                ["A"] * cells_per_type + ["B"] * cells_per_type + ["C"] * (n_obs - 2 * cells_per_type)
            )
        },
        index=[f"cell_{i:04d}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"gene_symbol": [f"gene_{j:03d}" for j in range(n_vars)]},
        index=[f"g{j:03d}" for j in range(n_vars)],
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    out_path = Path(__file__).parent / "pbmc_3k_mini.h5ad"
    adata.write_h5ad(out_path)
    print(f"[demo-data build] wrote {out_path} ({out_path.stat().st_size} bytes)")
    return out_path


if __name__ == "__main__":
    build()
