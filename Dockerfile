# Dockerfile — kdr-mef2c-scanpy
# Reproducible Python environment for the Kdr;Mef2cAHF Phase 1 scanpy pipeline.
#
# Build: docker build -t gideondocker98/kdr-mef2c-scanpy:latest .
# Push : docker push gideondocker98/kdr-mef2c-scanpy:latest

FROM python:3.11-slim

LABEL maintainer="Gideon Obeng <obgideon@gmail.com>"
LABEL org.opencontainers.image.source="https://github.com/drgideonobeng/kdr-mef2c-cardiopharyngeal-mesoderm-scanpy"
LABEL description="Kdr;Mef2cAHF Phase 1 scRNA-seq QC/clustering pipeline (scanpy)"

# ── System dependencies ──────────────────────────────────────────────────
# build-essential/libhdf5-dev: compiled deps (igraph, h5py);
# gfortran/git: scikit-misc has no linux/arm64 wheel, so it must compile
#   from source here -- it wraps a Fortran LOESS implementation (gfortran)
#   and its build script probes git metadata during version detection (git);
# procps: provides `ps`, required by Nextflow to collect task metrics.
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential \
      gfortran \
      git \
      libhdf5-dev \
      procps \
    && rm -rf /var/lib/apt/lists/*

# ── Python packages ──────────────────────────────────────────────────────
# scanpy: core analysis; scikit-image: Scrublet automatic thresholding;
# python-igraph + leidenalg: flavor="igraph" Leiden clustering backend;
# scikit-misc: loess regression required by flavor="seurat_v3" HVG selection.
RUN pip install --no-cache-dir \
      scanpy==1.10.* \
      scikit-image \
      scikit-misc \
      python-igraph \
      leidenalg \
      anndata \
      numpy \
      pandas \
      matplotlib

# ── Smoke test — fail the build if key packages are missing ─────────────
RUN python -c "\
import scanpy, igraph, leidenalg, skimage, skmisc; \
print('scanpy', scanpy.__version__)"

WORKDIR /workspace
CMD ["/bin/bash"]