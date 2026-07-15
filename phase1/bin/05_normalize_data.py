#!/usr/bin/env python
"""
05_normalize_data.py

Normalize library size, log-transform, and select highly variable genes.
Raw counts and log-normalized data are both kept as named layers so later
steps (scaling, marker-gene testing) can explicitly choose which matrix
they operate on. Python equivalent of 05_normalize_data.R.

Usage:
    python 05_normalize_data.py \
        --in_h5ad kdr_mef2c_wt_e10_04_singlets.h5ad \
        --n_hvgs 2000 \
        --seed 42 \
        --out kdr_mef2c_wt_e10_05_normalized.h5ad
"""

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc

from common import DEFAULT_SEED, log


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Normalize data and select HVGs")
    p.add_argument("--in_h5ad", required=True, type=Path)
    p.add_argument("--n_hvgs", type=int, default=2000,
                    help="Number of highly variable genes to select (default: 2000)")
    p.add_argument("--target_sum", type=float, default=1e4,
                    help="Target total count per cell after normalization (default: 1e4)")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED,
                    help="Random seed for reproducibility (default: 42)")
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)

    log("05", f"reading {args.in_h5ad}")
    adata = sc.read_h5ad(args.in_h5ad)

    # Keep the raw counts as a named layer before any transformation.
    adata.layers["counts"] = adata.X.copy()

    log("05", f"normalizing (target_sum={args.target_sum}) and log1p-transforming")
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)

    # Keep the log-normalized data as a named layer before scaling happens
    # in the next script.
    adata.layers["lognorm"] = adata.X.copy()

    log("05", f"selecting top {args.n_hvgs} highly variable genes "
              "(flavor=seurat_v3, on raw counts)")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=args.n_hvgs, flavor="seurat_v3", layer="counts"
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    log("05", f"saved {args.out}")


if __name__ == "__main__":
    main()
