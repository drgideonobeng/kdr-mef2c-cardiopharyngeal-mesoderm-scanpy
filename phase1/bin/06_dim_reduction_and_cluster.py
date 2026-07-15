#!/usr/bin/env python
"""
06_dim_reduction_and_cluster.py

Scale the HVGs, run PCA, build the neighbor graph, compute UMAP, cluster
with Leiden, and rank per-cluster marker genes. This is the last script
in Phase 1 -- its output is the frozen per-sample object that Phase 2
(integration) reads in. Python equivalent of 06_dim_reduction_and_cluster.R.

Marker genes are tested on the log-normalized layer, not the scaled
matrix that this script produces in adata.X -- scaled (z-scored) values
would distort test statistics and make fold-changes uninterpretable.

Usage:
    python 06_dim_reduction_and_cluster.py \
        --in_h5ad kdr_mef2c_wt_e10_05_normalized.h5ad \
        --n_pcs 30 --cluster_res 0.3 --seed 42 \
        --out_rds kdr_mef2c_wt_e10_06_phase1.h5ad \
        --out_elbow kdr_mef2c_wt_e10_06_elbow.pdf \
        --out_umap kdr_mef2c_wt_e10_06_umap.pdf
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

from common import DEFAULT_SEED, log


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Scale, reduce dimensions, and cluster")
    p.add_argument("--in_h5ad", required=True, type=Path)
    p.add_argument("--n_pcs", type=int, default=30,
                    help="Number of PCs used to build the neighbor graph (default: 30)")
    p.add_argument("--n_neighbors", type=int, default=15,
                    help="Number of neighbors per cell (default: 15)")
    p.add_argument("--cluster_res", type=float, default=0.3,
                    help="Leiden clustering resolution (default: 0.3)")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED,
                    help="Random seed for reproducibility (default: 42)")
    p.add_argument("--out_rds", required=True, type=Path,
                    help="Output .h5ad path (named --out_rds to mirror the R CLI)")
    p.add_argument("--out_elbow", required=True, type=Path)
    p.add_argument("--out_umap", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    log("06", f"reading {args.in_h5ad}")
    adata = sc.read_h5ad(args.in_h5ad)

    log("06", "scaling highly variable genes (zero_center=True, max_value=10)")
    sc.pp.scale(adata, zero_center=True, max_value=10)

    log("06", "running PCA")
    sc.tl.pca(
        adata, n_comps=50, mask_var="highly_variable",
        svd_solver="arpack", random_state=args.seed,
    )

    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    plt.gcf().savefig(args.out_elbow, bbox_inches="tight")
    plt.close()
    log("06", f"saved {args.out_elbow}")

    log("06", f"computing neighbor graph (n_neighbors={args.n_neighbors}, "
              f"n_pcs={args.n_pcs})")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs,
                     random_state=args.seed)

    log("06", "computing UMAP")
    sc.tl.umap(adata, random_state=args.seed)

    log("06", f"clustering with Leiden (resolution={args.cluster_res})")
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2,
                  resolution=args.cluster_res, random_state=args.seed)
    log("06", f"{adata.obs['leiden'].nunique()} clusters found")

    sc.pl.umap(adata, color=["leiden"], show=False)
    plt.gcf().savefig(args.out_umap, bbox_inches="tight")
    plt.close()
    log("06", f"saved {args.out_umap}")

    log("06", "ranking cluster marker genes (Wilcoxon, on the lognorm layer)")
    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method="wilcoxon",
        layer="lognorm", use_raw=False,
    )

    args.out_rds.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out_rds)
    log("06", f"saved {args.out_rds}")


if __name__ == "__main__":
    main()
