#!/usr/bin/env python
"""
03_filter_cells.py

Apply gene and cell filters using the same thresholds visualized in
02_visualize_qc.py. Python equivalent of 03_filter_cells.R.

Usage:
    python 03_filter_cells.py \
        --in_h5ad kdr_mef2c_wt_e10_02_with_qc.h5ad \
        --min_features 200 --max_features 6000 \
        --min_counts 1000 --max_counts 50000 \
        --max_percent_mt 10.0 \
        --out kdr_mef2c_wt_e10_03_filtered.h5ad
"""

import argparse
from pathlib import Path

import scanpy as sc

from common import log


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Filter cells and genes on QC thresholds")
    p.add_argument("--in_h5ad", required=True, type=Path)
    p.add_argument("--min_features", type=int, required=True,
                    help="Min genes detected per cell")
    p.add_argument("--max_features", type=int, required=True,
                    help="Max genes detected per cell")
    p.add_argument("--min_counts", type=int, required=True,
                    help="Min total UMI counts per cell")
    p.add_argument("--max_counts", type=int, required=True,
                    help="Max total UMI counts per cell")
    p.add_argument("--max_percent_mt", type=float, required=True,
                    help="Max %% mitochondrial counts per cell")
    p.add_argument("--min_cells", type=int, default=3,
                    help="Min cells a gene must be detected in to be kept (default: 3)")
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    log("03", f"reading {args.in_h5ad}")
    adata = sc.read_h5ad(args.in_h5ad)

    log("03", f"filtering genes (min_cells={args.min_cells})")
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    n_before = adata.n_obs
    log("03", f"filtering cells (min_features={args.min_features}, "
              f"max_features={args.max_features}, min_counts={args.min_counts}, "
              f"max_counts={args.max_counts}, max_percent_mt={args.max_percent_mt})")
    sc.pp.filter_cells(adata, min_genes=args.min_features)
    sc.pp.filter_cells(adata, max_genes=args.max_features)
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.filter_cells(adata, max_counts=args.max_counts)
    adata = adata[adata.obs["pct_counts_mito"] < args.max_percent_mt, :].copy()
    log("03", f"{n_before} -> {adata.n_obs} cells after filtering")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    log("03", f"saved {args.out}")


if __name__ == "__main__":
    main()
