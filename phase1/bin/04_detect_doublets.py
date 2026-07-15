#!/usr/bin/env python
"""
04_detect_doublets.py

Detect and remove doublets with Scrublet. Must run per-sample, on data
that has already been QC-filtered but NOT yet normalized or pooled with
other samples. Python equivalent of 04_detect_doublets.R (which uses
DoubletFinder; Scrublet is scanpy's native equivalent).

Usage:
    python 04_detect_doublets.py \
        --in_h5ad kdr_mef2c_wt_e10_03_filtered.h5ad \
        --sample_id kdr_mef2c_wt_e10 \
        --seed 42 \
        --out_h5ad kdr_mef2c_wt_e10_04_singlets.h5ad \
        --out_summary kdr_mef2c_wt_e10_04_doublet_summary.csv \
        --out_pdf kdr_mef2c_wt_e10_04_doublet_plots.pdf
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

from common import DEFAULT_SEED, log


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Detect and remove doublets")
    p.add_argument("--in_h5ad", required=True, type=Path)
    p.add_argument("--sample_id", required=True, type=str)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED,
                    help="Random seed for reproducibility (default: 42)")
    p.add_argument("--out_h5ad", required=True, type=Path)
    p.add_argument("--out_summary", required=True, type=Path)
    p.add_argument("--out_pdf", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    log("04", f"reading {args.in_h5ad}")
    adata = sc.read_h5ad(args.in_h5ad)

    log("04", f"running Scrublet (seed={args.seed})")
    sc.pp.scrublet(adata, random_state=args.seed)

    sc.pl.scrublet_score_distribution(adata, show=False)
    plt.gcf().savefig(args.out_pdf, bbox_inches="tight")
    plt.close()
    log("04", f"saved {args.out_pdf} -- inspect for a clear bimodal split "
              "before trusting the threshold")

    n_before = adata.n_obs
    doublet_rate = adata.obs["predicted_doublet"].mean()
    n_doublets = int(adata.obs["predicted_doublet"].sum())
    adata = adata[~adata.obs["predicted_doublet"], :].copy()
    log("04", f"doublet rate: {doublet_rate:.1%} ({n_before} -> {adata.n_obs} cells)")

    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_summary, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample_id", "n_cells_before", "n_doublets",
                          "doublet_rate", "n_cells_after"])
        writer.writerow([args.sample_id, n_before, n_doublets,
                          round(doublet_rate, 4), adata.n_obs])
    log("04", f"saved {args.out_summary}")

    args.out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out_h5ad)
    log("04", f"saved {args.out_h5ad}")


if __name__ == "__main__":
    main()
