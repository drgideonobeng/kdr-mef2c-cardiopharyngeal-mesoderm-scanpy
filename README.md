# kdr-mef2c-cardiopharyngeal-mesoderm-scanpy

Single-cell RNA-seq pipeline studying vascular progenitors of the
cardiopharyngeal mesoderm in *Kdr;Mef2cAHF-Cre* WT vs cKO mouse embryos
(E9.5–E10). This is a **Python/scanpy port** of the original R/Seurat
pipeline, built to mirror its architecture and outputs 1:1 while
migrating the implementation to scanpy + Nextflow DSL2 + Docker.

## Biological question

*Kdr* (VEGFR2) marks endothelial and vascular progenitor populations;
*Mef2cAHF-Cre* drives recombination in the anterior heart field lineage.
Conditional loss of function in this compartment is used to identify
vascular progenitors of the pharyngeal arch arteries and the common
cardinal vein, and to characterize how their transcriptional identity
shifts between WT and cKO.

## Data

~4,422 tdTomato+/GFP+ FACS-sorted cells across 2 pooled samples:

| Sample ID            | Genotype | Notes                              |
|-----------------------|----------|-------------------------------------|
| kdr_mef2c_wt_e10      | WT       | tdTomato+/GFP+ FACS-sorted WT pool  |
| kdr_mef2c_cko_e10     | cKO      | tdTomato+/GFP+ FACS-sorted cKO pool |

Raw 10x Genomics `filtered_feature_bc_matrix` output (`barcodes.tsv.gz`,
`features.tsv.gz`, `matrix.mtx.gz`) is expected under `raw_data/`, one
subdirectory per sample. Raw data and pipeline outputs are not tracked
in this repository (see `.gitignore`); see **Reproducing** below.

## Repository structure

```
.
├── Dockerfile              # scanpy + scikit-image + scikit-misc + igraph/leidenalg
├── params/
│   └── samples.yml         # single source of truth: sample IDs, genotypes, paths
├── phase1/                 # per-sample QC through clustering
│   ├── main.nf
│   ├── nextflow.config
│   └── bin/
│       ├── common.py
│       ├── 01_create_anndata_object.py
│       ├── 02_visualize_qc.py
│       ├── 03_filter_cells.py
│       ├── 04_detect_doublets.py
│       ├── 05_normalize_data.py
│       └── 06_dim_reduction_and_cluster.py
├── raw_data/                # not tracked — see Reproducing
└── results/                  # not tracked — pipeline output, regeneratable
```

Later phases (`phase2/` integration, `phase3/` annotation, `phase4/`
perturbation analysis) will follow the same per-phase layout: their own
`main.nf`, `nextflow.config`, and `bin/`, sharing this root-level
`Dockerfile` and `params/samples.yml`.

## Phase 1 — per-sample QC through clustering

One script per concern, chained via explicit `.h5ad` inputs/outputs,
mirroring the R pipeline's `--in_rds`/`--out_rds` discipline:

| Step | Script | What it does |
|---|---|---|
| 01 | `01_create_anndata_object.py` | Read a 10x matrix into an AnnData object, tag `sample_id`/`genotype` |
| 02 | `02_visualize_qc.py` | Compute QC metrics (genes/cell, counts/cell, %mito); plot distributions with proposed thresholds overlaid — no filtering yet |
| 03 | `03_filter_cells.py` | Apply the QC thresholds from step 02 |
| 04 | `04_detect_doublets.py` | Scrublet doublet detection and removal, run per-sample before any pooling |
| 05 | `05_normalize_data.py` | Library-size normalization, log1p, HVG selection (`seurat_v3` on raw counts) |
| 06 | `06_dim_reduction_and_cluster.py` | Scale → PCA → neighbor graph → UMAP → Leiden clustering → marker genes |

Each sample runs through all 6 steps independently. Outputs are
organized **per step** (all samples together, filenames prefixed with
`sample_id`), matching the R pipeline's output convention:

```
results/phase1/06_dim_reduction_and_cluster/
├── kdr_mef2c_cko_e10_06_phase1.h5ad
├── kdr_mef2c_cko_e10_06_umap.pdf
├── kdr_mef2c_wt_e10_06_phase1.h5ad
└── kdr_mef2c_wt_e10_06_umap.pdf
```

QC thresholds are configured in `phase1/nextflow.config` under
`params.defaults`, with optional per-sample overrides under
`params.samples`.

## House style

- **Python**: `argparse`, `pathlib`, explicit `--in_h5ad`/`--out_h5ad`
  I/O, one concern per script, `random_state=42` set on every
  stochastic step, layer discipline (`counts` and `lognorm` kept as
  named `adata.layers` before scaling ever overwrites `adata.X`)
- **Nextflow**: DSL2, Docker-first, per-sample `meta` map
  (`[id, genotype]`) carried through the whole channel chain, named
  `emit:` outputs, per-run provenance (`report`/`timeline`/`trace`/`dag`)
- **Reproducibility**: single seed (`params.seed = 42`) threaded through
  every stochastic step (Scrublet, PCA, UMAP, Leiden)

## Environment

| Tool | Version / details |
|---|---|
| Docker image | `gideondocker98/kdr-mef2c-scanpy:latest` (multi-arch: amd64 + arm64) |
| Base image | `python:3.11-slim` |
| Key Python packages | scanpy, anndata, scikit-image, scikit-misc, python-igraph, leidenalg |
| Nextflow | ≥ 24.04, tested on 26.04 |
| Hardware | macOS M3 Max (profile: `docker,apple_silicon`) |

## Reproducing

```bash
# 1. Build and push the image (multi-arch)
docker buildx build --platform linux/amd64,linux/arm64 \
    -t gideondocker98/kdr-mef2c-scanpy:latest --push .

# 2. Place raw 10x output under raw_data/<sample_id>/filtered_feature_bc_matrix/
#    (see params/samples.yml for exact expected paths)

# 3. Run Phase 1
cd phase1
nextflow run main.nf -profile docker,apple_silicon
```

## Related projects

- **[kdr-mef2c-cardiopharyngeal-mesoderm](https://github.com/drgideonobeng/kdr-mef2c-cardiopharyngeal-mesoderm)**
  — the original R/Seurat implementation of this project. This
  repository is a Python/scanpy port, built to mirror its architecture
  and outputs 1:1.
- **[tbx1-depleted-cardiopharyngeal-mesoderm](https://github.com/drgideonobeng/tbx1-depleted-cardiopharyngeal-mesoderm)**
  — sibling perturbation-arm project (Mesp1-Cre;Tbx1cKO, 4 samples,
  R/Seurat), sharing the same house style and phase structure.

## Status

- [x] Repo scaffolded (Dockerfile, samples.yml, .gitignore)
- [x] Docker image built and pushed
- [x] Phase 1 — per-sample QC through clustering
- [ ] Phase 2 — integration and annotation
- [ ] Phase 3 — cell-type annotation
- [ ] Phase 4 — perturbation analysis