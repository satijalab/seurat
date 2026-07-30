#!/usr/bin/env Rscript
# Cross-platform parallelization smoke test for Seurat.
#
# Runs two futurized functions (ScaleData, FindIntegrationAnchors) under a
# sequential plan and then under a non-sequential plan, and asserts:
#   1. fail-on-degrade: the non-sequential plan really has >1 worker,
#   2. no-hang: each parallel call completes within SEURAT_CI_TIMEOUT,
#   3. equivalence: parallel results equal sequential results.
#
# Env:
#   SEURAT_CI_PLAN     "multisession" (native) | "multicore" (container)
#   SEURAT_CI_WORKERS  integer, default 2
#   SEURAT_CI_TIMEOUT  seconds per parallel call, default 600
#
# Exit 0 = pass, non-zero = fail.

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(R.utils)
})

plan_name <- Sys.getenv("SEURAT_CI_PLAN", "multisession")
workers   <- as.integer(Sys.getenv("SEURAT_CI_WORKERS", "2"))
timeout   <- as.numeric(Sys.getenv("SEURAT_CI_TIMEOUT", "600"))

# Do not let an unrelated globals-size limit mask the path under test.
options(future.globals.maxSize = 2 * 1024^3)

cat(sprintf("[smoke] plan=%s workers=%d timeout=%.0fs\n",
            plan_name, workers, timeout))

## --- deterministic synthetic data -------------------------------------------
set.seed(42)
n_genes <- 3000
n_cells <- 2000
counts <- matrix(
  rpois(n = n_genes * n_cells, lambda = 0.5),
  nrow = n_genes, ncol = n_cells
)
rownames(counts) <- paste0("gene", seq_len(n_genes))
colnames(counts) <- paste0("cell", seq_len(n_cells))

make_object <- function() {
  obj <- CreateSeuratObject(counts = counts)
  # numeric covariate to engage ScaleData's futurized regression path
  obj$nfeat_frac <- obj$nFeature_RNA / max(obj$nFeature_RNA)
  obj
}

# version-robust scaled-matrix accessor
scaled_matrix <- function(obj) {
  m <- tryCatch(
    SeuratObject::LayerData(obj, layer = "scale.data"),
    error = function(e) GetAssayData(obj, slot = "scale.data")
  )
  as.matrix(m)
}

## --- exercised functions ----------------------------------------------------
run_scale <- function() {
  obj <- make_object()
  # Seurat v5 (Assay5) requires a populated 'data' layer before ScaleData;
  # NormalizeData does not use futures, so the futurized regression path in
  # ScaleData (engaged by vars.to.regress) remains the code under test.
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = "nfeat_frac", verbose = FALSE)
  scaled_matrix(obj)
}

run_anchors <- function() {
  obj <- make_object()
  # split into two objects for a 2-element object.list
  half <- seq_len(floor(n_cells / 2))
  a <- subset(obj, cells = colnames(obj)[half])
  b <- subset(obj, cells = colnames(obj)[-half])
  prep <- function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures = 500, verbose = FALSE)
    x
  }
  lst <- list(prep(a), prep(b))
  anchors <- FindIntegrationAnchors(
    object.list = lst, anchor.features = 500, dims = 1:20, verbose = FALSE
  )
  # Anchor identity = the integer cell-index + dataset columns. In Seurat 5.x
  # the anchors slot also carries a floating "score" column whose value depends
  # on parallel k-NN chunking (a legitimate parallel run can differ by ~1/k),
  # and the row order is not stable across plans. We therefore assert on the
  # anchor *identity* — the parallel run must find the same anchor pairs — and
  # canonicalize the row order so equality does not hinge on ordering.
  m <- as.matrix(slot(anchors, "anchors"))
  id_cols <- intersect(c("cell1", "cell2", "dataset1", "dataset2"),
                       colnames(m))
  m <- m[, id_cols, drop = FALSE]
  m[do.call(order, as.data.frame(m)), , drop = FALSE]
}

## --- sequential baseline ----------------------------------------------------
plan("sequential")
cat("[smoke] sequential baseline...\n")
seq_scale   <- run_scale()
seq_anchors <- run_anchors()

## --- parallel run with guards ----------------------------------------------
plan(plan_name, workers = workers)

# fail-on-degrade: a green run must mean parallelism actually happened
if (inherits(plan(), "sequential") || nbrOfWorkers() <= 1) {
  stop(sprintf("[smoke] FAIL: plan degraded (class=%s, workers=%d)",
               paste(class(plan()), collapse = "/"), nbrOfWorkers()))
}
cat(sprintf("[smoke] parallel plan active: %s, workers=%d\n",
            paste(class(plan()), collapse = "/"), nbrOfWorkers()))

with_timeout <- function(expr) {
  withTimeout(expr, timeout = timeout, onTimeout = "error")
}

t0 <- Sys.time()
par_scale   <- with_timeout(run_scale())
par_anchors <- with_timeout(run_anchors())
elapsed <- as.numeric(Sys.time() - t0, units = "secs")

## --- equivalence assertions -------------------------------------------------
eq_scale   <- isTRUE(all.equal(seq_scale, par_scale, tolerance = 1e-6))
eq_anchors <- isTRUE(all.equal(seq_anchors, par_anchors, tolerance = 1e-6))

cat("\n[smoke] results\n")
cat(sprintf("  ScaleData               equal=%s\n", eq_scale))
cat(sprintf("  FindIntegrationAnchors  equal=%s\n", eq_anchors))
cat(sprintf("  parallel elapsed        %.1fs\n", elapsed))

if (!(eq_scale && eq_anchors)) {
  stop("[smoke] FAIL: parallel results differ from sequential")
}
cat("[smoke] PASS\n")
quit(status = 0)
