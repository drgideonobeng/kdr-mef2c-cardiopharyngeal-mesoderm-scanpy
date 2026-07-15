#!/usr/bin/env python
"""
02_visualize_qc.py

Compute per-cell QC metrics (genes detected, total counts, %mito) and
plot their distributions with the proposed threshold lines overlaid for
visual reference. This step does NOT filter any cells or genes -- it
only computes, visualizes, and lets you eyeball whether the thresholds
in nextflow.config look reasonable before 03_filter_cells.py applies
them. Python equivalent of 02_visualize_qc.R.

Usage:
    python 02_visualize_qc.py \
        --in_h5ad kdr_mef2c_wt_e10_01_unfiltered.h5ad \
        --sample_id kdr_mef2c_wt_e10 \
        --min_features 200 --max_features 6000 \
        --min_counts 1000 --max_counts 50000 \
        --max_percent_mt 10.0 \
        --out_h5ad kdr_mef2c_wt_e10_02_with_qc.h5ad \
        --out_pdf kdr_mef2c_wt_e10_02_qc_plots.pdf
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

from common import log


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute and plot QC metrics")
    p.add_argument("--in_h5ad", required=True, type=Path)
    p.add_argument("--sample_id", required=True, type=str)
    p.add_argument("--min_features", type=int, required=True)
    p.add_argument("--max_features", type=int, required=True)
    p.add_argument("--min_counts", type=int, required=True)
    p.add_argument("--max_counts", type=int, required=True)
    p.add_argument("--max_percent_mt", type=float, required=True)
    p.add_argument("--mito_prefix", type=str, default="mt-",
                    help="Mitochondrial gene prefix; 'mt-' mouse, 'MT-' human")
    p.add_argument("--out_h5ad", required=True, type=Path)
    p.add_argument("--out_pdf", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    log("02", f"reading {args.in_h5ad}")
    adata = sc.read_h5ad(args.in_h5ad)

    adata.var["mito"] = adata.var_names.str.startswith(args.mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    log("02", "computed n_genes_by_counts, total_counts, pct_counts_mito")

    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mito"]
    thresholds = {
        "n_genes_by_counts": (args.min_features, args.max_features),
        "total_counts": (args.min_counts, args.max_counts),
        "pct_counts_mito": (None, args.max_percent_mt),
    }

    sc.pl.violin(adata, metrics, jitter=0.3, multi_panel=True, show=False)
    fig = plt.gcf()
    for ax, metric in zip(fig.axes, metrics):
        lo, hi = thresholds[metric]
        if lo is not None:
            ax.axhline(lo, color="red", linestyle="--", linewidth=1)
        if hi is not None:
            ax.axhline(hi, color="red", linestyle="--", linewidth=1)
    fig.suptitle(f"{args.sample_id}: proposed QC thresholds (dashed red)")
    fig.savefig(args.out_pdf, bbox_inches="tight")
    plt.close(fig)
    log("02", f"saved {args.out_pdf}")

    args.out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out_h5ad)
    log("02", f"saved {args.out_h5ad}")


if __name__ == "__main__":
    main()
