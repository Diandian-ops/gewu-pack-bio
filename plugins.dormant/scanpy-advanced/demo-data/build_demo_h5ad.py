"""Generate the demo h5ad fixture for scanpy-advanced acceptance fixtures.

This script is intentionally small + deterministic:
  - 270 cells × 100 highly variable genes (PBMC 3k-like compressed)
  - Three simulated cell types via distinct mean expression patterns
  - Seeds scanpy / numpy RNG so the produced h5ad is byte-stable across runs

Output path: resources/built-in-plugins/scanpy-advanced/demo-data/pbmc_3k_mini.h5ad

Run:
    # From BioF3 repo root, with biof3-py-runtime scanpy available
    python resources/built-in-plugins/scanpy-advanced/demo-data/build_demo_h5ad.py

Idempotent: re-running overwrites the same sha256-stable file (modulo Pillow PNG canonicalization
in scripts/scanpy-advanced.py is for figures, not h5ad; h5ad byte-stability relies on this seed).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def build() -> Path:
    """Build a deterministic 270 × 100 AnnData and write pbmc_3k_mini.h5ad."""
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
