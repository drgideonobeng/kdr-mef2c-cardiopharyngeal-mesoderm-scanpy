#!/usr/bin/env python
"""
01_create_anndata_object.py

Read a 10x filtered_feature_bc_matrix directory into an AnnData object
and tag it with sample-level metadata. Python equivalent of
01_create_seurat_object.R -- CLI shape mirrors it 1:1.

Usage:
    python 01_create_anndata_object.py \
        --tenx_dir raw_data/kdr_Mef2cAHF_wt_e10/filtered_feature_bc_matrix \
        --sample_id kdr_mef2c_wt_e10 \
        --genotype WT \
        --out kdr_mef2c_wt_e10_01_unfiltered.h5ad
"""

import argparse
import sys
from pathlib import Path

import scanpy as sc

from common import log


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create an AnnData object from 10x data")
    p.add_argument("--tenx_dir", required=True, type=Path,
                    help="Path to the 10x filtered_feature_bc_matrix directory")
    p.add_argument("--sample_id", required=True, type=str,
                    help="Sample identifier, e.g. kdr_mef2c_wt_e10")
    p.add_argument("--genotype", required=True, type=str, choices=["WT", "cKO"],
                    help="Genotype label for this sample")
    p.add_argument("--out", required=True, type=Path,
                    help="Path to write the resulting .h5ad file")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.tenx_dir.exists():
        sys.exit(f"[01] ERROR: input directory not found: {args.tenx_dir}")

    log("01", f"reading 10x matrix for {args.sample_id} from {args.tenx_dir}")
    adata = sc.read_10x_mtx(args.tenx_dir, var_names="gene_symbols", make_unique=True)

    adata.obs["sample_id"] = args.sample_id
    adata.obs["genotype"] = args.genotype
    log("01", f"loaded {adata.n_obs} cells x {adata.n_vars} genes")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    log("01", f"saved {args.out}")


if __name__ == "__main__":
    main()
