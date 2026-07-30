# Pinned Linux environment where plan("multicore") works — the fork-based
# future strategy unavailable to Windows/macOS natively. Windows users
# (including WSL) run Seurat here to get real multicore parallelism.
FROM rocker/r-ver:4.4.1

# System libraries required to build Seurat and its dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends \
      libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
      libpng-dev libtiff5-dev libjpeg-dev libgeos-dev libgsl-dev \
      libfftw3-dev libglpk-dev libboost-all-dev cmake git \
    && rm -rf /var/lib/apt/lists/*

# Parallelization + timeout deps (also pulled transitively, installed
# explicitly so the harness can run even before Seurat finishes resolving).
RUN R -q -e "install.packages(c('remotes','future','future.apply','R.utils'))"

# Install the Seurat source under test from the build context (repo root),
# so the image exercises the exact code in the PR — not a CRAN release.
WORKDIR /seurat
COPY . /seurat
RUN R -q -e "remotes::install_local('/seurat', dependencies = TRUE, upgrade = 'never')"

# Default: run the smoke harness under multicore.
ENV SEURAT_CI_PLAN=multicore \
    SEURAT_CI_WORKERS=2 \
    SEURAT_CI_TIMEOUT=900
CMD ["Rscript", "/seurat/inst/ci/parallel_smoke.R"]
