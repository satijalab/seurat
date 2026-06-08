set.seed(42)

# Build a small reference dgCMatrix-backed object and a matching DelayedMatrix
.make_counts <- function(nfeat = 120, ncell = 60, lambda = 1.3) {
  m <- as(matrix(rpois(nfeat * ncell, lambda), nfeat, ncell), "dgCMatrix")
  dimnames(m) <- list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  m
}

context("CheckMatrixSize")

test_that("CheckMatrixSize reports counts and the 2^31 limit", {
  m <- .make_counts()
  res <- CheckMatrixSize(m, warn = FALSE)
  expect_equal(res$limit, .Machine$integer.max)
  expect_equal(res$n, length(slot(m, "x")))
  expect_false(res$exceeds)
  expect_false(res$on.disk)
})

test_that("CheckMatrixSize flags on-disk backends as exempt", {
  skip_if_not_installed("DelayedArray")
  d <- as.DelayedMatrix(.make_counts())
  res <- CheckMatrixSize(d, warn = FALSE)
  expect_true(res$on.disk)
})

context("DelayedMatrix backend")

test_that("as.DelayedMatrix and as.sparse round-trip with dimnames", {
  skip_if_not_installed("DelayedArray")
  m <- .make_counts()
  d <- as.DelayedMatrix(m)
  expect_s4_class(d, "DelayedMatrix")
  expect_identical(dimnames(d), dimnames(m))
  expect_true(all.equal(as.sparse(d), m) == TRUE)
})

test_that("in-memory as.DelayedMatrix uses a 64-bit SVT seed and saves to one rds", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SparseArray")
  m <- .make_counts()
  d <- as.DelayedMatrix(m)
  expect_s4_class(d, "DelayedMatrix")
  expect_true(is(DelayedArray::seed(d), "SVT_SparseMatrix"))
  expect_identical(dimnames(d), dimnames(m))
  obj <- suppressWarnings(CreateSeuratObject(counts = d))
  # plain saveRDS yields a single self-contained file (in-memory, no store)
  f <- tempfile(fileext = ".rds")
  saveRDS(obj, f)
  obj2 <- readRDS(f)
  expect_true(inherits(LayerData(obj2, "counts"), "DelayedMatrix"))
  expect_equal(as.matrix(LayerData(obj2, "counts")), as.matrix(m))
})

test_that("as.DelayedMatrix(ondisk=TRUE) gives a 64-bit on-disk matrix", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("HDF5Array")
  m <- .make_counts()
  d <- as.DelayedMatrix(m, ondisk = TRUE)
  expect_s4_class(d, "DelayedMatrix")
  expect_true(is(d, "TENxMatrix"))
  expect_identical(dimnames(d), dimnames(m))
  expect_equal(as.matrix(d), as.matrix(m))
  # element/index space is 64-bit safe (a sparse > 2^31-element matrix)
  big <- as.DelayedMatrix(
    Matrix::sparseMatrix(i = c(1L, 80000L), j = c(1L, 40000L), x = 1,
                         dims = c(80000L, 40000L)),
    ondisk = TRUE
  )
  expect_gt(length(big), .Machine$integer.max)
})

test_that("DelayedMatrix preprocessing matches dgCMatrix", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("DelayedMatrixStats")
  m <- .make_counts()
  ref <- suppressWarnings(CreateSeuratObject(counts = m))
  dly <- suppressWarnings(CreateSeuratObject(counts = as.DelayedMatrix(m)))
  # nCount / nFeature via .CalcN
  expect_equal(unname(ref$nCount_RNA), unname(dly$nCount_RNA))
  expect_equal(unname(ref$nFeature_RNA), unname(dly$nFeature_RNA))
  # normalization stays on-disk and matches
  ref <- suppressWarnings(NormalizeData(ref, verbose = FALSE))
  dly <- suppressWarnings(NormalizeData(dly, verbose = FALSE))
  expect_s4_class(LayerData(dly, "data"), "DelayedMatrix")
  expect_equal(
    as.matrix(LayerData(ref, "data")),
    as.matrix(LayerData(dly, "data")),
    tolerance = 1e-6
  )
  # variable features: standardized variances match numerically (rank order of
  # near-ties may differ by tiny float differences between the C++ and
  # DelayedMatrixStats paths), and the selected set is the same
  ref <- suppressWarnings(FindVariableFeatures(ref, nfeatures = 30, verbose = FALSE))
  dly <- suppressWarnings(FindVariableFeatures(dly, nfeatures = 30, verbose = FALSE))
  vs.ref <- HVFInfo(ref[["RNA"]])[rownames(ref), "variance.standardized"]
  vs.dly <- HVFInfo(dly[["RNA"]])[rownames(ref), "variance.standardized"]
  expect_equal(vs.ref, vs.dly, tolerance = 1e-3)
  expect_gte(length(intersect(VariableFeatures(ref), VariableFeatures(dly))), 27)
})

test_that("DelayedMatrix-backed object runs PCA", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("DelayedMatrixStats")
  dly <- suppressWarnings(CreateSeuratObject(counts = as.DelayedMatrix(.make_counts())))
  dly <- suppressWarnings(NormalizeData(dly, verbose = FALSE))
  dly <- suppressWarnings(FindVariableFeatures(dly, verbose = FALSE))
  dly <- suppressWarnings(ScaleData(dly, verbose = FALSE))
  dly <- suppressWarnings(RunPCA(dly, npcs = 10, verbose = FALSE))
  expect_true("pca" %in% Reductions(dly))
  expect_equal(ncol(Embeddings(dly, "pca")), 10)
})

test_that("non-block-native ops materialize and match dgCMatrix", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SparseArray")
  m <- .make_counts(ncell = 100)
  ref <- suppressWarnings(NormalizeData(CreateSeuratObject(counts = m), verbose = FALSE))
  dly <- suppressWarnings(NormalizeData(
    CreateSeuratObject(counts = as.DelayedMatrix(m)), verbose = FALSE
  ))
  ref$cond <- rep(c("a", "b"), 50)
  dly$cond <- rep(c("a", "b"), 50)
  # CLR / RC normalization
  for (nm in c("CLR", "RC")) {
    a <- as.matrix(LayerData(suppressWarnings(
      NormalizeData(ref, normalization.method = nm, verbose = FALSE)), "data"))
    b <- as.matrix(LayerData(suppressWarnings(
      NormalizeData(dly, normalization.method = nm, verbose = FALSE)), "data"))
    expect_equal(a, b, tolerance = 1e-6)
  }
  # mean.var.plot variable features
  expect_equal(
    VariableFeatures(suppressWarnings(FindVariableFeatures(
      ref, selection.method = "mvp", verbose = FALSE))),
    VariableFeatures(suppressWarnings(FindVariableFeatures(
      dly, selection.method = "mvp", verbose = FALSE)))
  )
  # pseudobulk summaries
  expect_equal(
    as.matrix(suppressWarnings(AggregateExpression(ref, group.by = "cond", verbose = FALSE))$RNA),
    as.matrix(suppressWarnings(AggregateExpression(dly, group.by = "cond", verbose = FALSE))$RNA)
  )
})

context("Multi-layer / integration on DelayedMatrix")

test_that("split DelayedMatrix layers join and integrate", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("DelayedMatrixStats")
  m1 <- .make_counts(ncell = 80, lambda = 1.2)
  m2 <- .make_counts(ncell = 80, lambda = 1.7)
  colnames(m2) <- paste0("b", seq_len(ncol(m2)))
  obj <- suppressWarnings(CreateSeuratObject(counts = as.DelayedMatrix(cbind(m1, m2))))
  obj$batch <- rep(c("a", "b"), each = 80)
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
  # JoinLayers stitches the DelayedMatrix layers back together
  joined <- JoinLayers(obj)
  expect_s4_class(LayerData(joined, "counts"), "DelayedMatrix")
  expect_equal(ncol(LayerData(joined, "counts")), 160)
  # full integration pipeline
  obj <- suppressWarnings(NormalizeData(obj, verbose = FALSE))
  obj <- suppressWarnings(FindVariableFeatures(obj, verbose = FALSE))
  obj <- suppressWarnings(ScaleData(obj, verbose = FALSE))
  obj <- suppressWarnings(RunPCA(obj, npcs = 10, verbose = FALSE))
  obj <- suppressWarnings(IntegrateLayers(
    obj, method = RPCAIntegration, verbose = FALSE, k.weight = 50
  ))
  expect_true("integrated.dr" %in% Reductions(obj))
})

context("spam coercion")

test_that("as.sparse.spam returns an equivalent dgCMatrix", {
  skip_if_not_installed("spam")
  m <- .make_counts(nfeat = 20, ncell = 10)
  s <- spam::as.spam(as.matrix(m))
  sp <- as.sparse(s)
  expect_s4_class(sp, "dgCMatrix")
  expect_equal(as.matrix(sp), unname(as.matrix(m)))
})

context("AsInMemory and single-file saveRDS")

test_that("AsInMemory materializes on-disk layers for a portable rds", {
  skip_if_not_installed("DelayedArray")
  m <- .make_counts()
  obj <- suppressWarnings(CreateSeuratObject(counts = as.DelayedMatrix(m)))
  expect_s4_class(LayerData(obj, "counts"), "DelayedMatrix")
  obj <- AsInMemory(obj, verbose = FALSE)
  expect_s4_class(LayerData(obj, "counts"), "dgCMatrix")
  expect_equal(LayerData(obj, "counts"), m)
  # single, self-contained file
  f <- tempfile(fileext = ".rds")
  saveRDS(obj, f)
  obj2 <- readRDS(f)
  expect_equal(LayerData(obj2, "counts"), m)
})

context("Portable bundle (SaveSeurat / LoadSeurat)")

test_that("SaveSeurat bundle survives a file move", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("HDF5Array")
  skip_if_not_installed("jsonlite")
  m <- .make_counts()
  obj <- suppressWarnings(CreateSeuratObject(counts = as.DelayedMatrix(m)))
  obj <- suppressWarnings(NormalizeData(obj, verbose = FALSE))
  ref <- as.matrix(LayerData(obj, "counts"))
  f1 <- tempfile(fileext = ".seurat")
  SaveSeurat(obj, f1, verbose = FALSE)
  f2 <- tempfile(fileext = ".seurat")
  file.rename(f1, f2)
  loaded <- LoadSeurat(f2, dir = tempfile(), verbose = FALSE)
  expect_equal(as.matrix(LayerData(loaded, "counts")), ref)
  expect_true("data" %in% Layers(loaded))
})

context("HDF5 container (SaveSeuratH5 / LoadSeuratH5)")

test_that("SaveSeuratH5 round-trips layers, metadata and reductions", {
  skip_if_not_installed("HDF5Array")
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("DelayedMatrixStats")
  m <- .make_counts(ncell = 80)
  obj <- suppressWarnings(CreateSeuratObject(counts = m))
  obj <- suppressWarnings(NormalizeData(obj, verbose = FALSE))
  obj <- suppressWarnings(FindVariableFeatures(obj, verbose = FALSE))
  obj <- suppressWarnings(ScaleData(obj, verbose = FALSE))
  obj <- suppressWarnings(RunPCA(obj, npcs = 10, verbose = FALSE))
  obj$grp <- factor(rep(c("x", "y"), 40))
  obj$flag <- rep(c(TRUE, FALSE), 40)
  ref.counts <- as.matrix(LayerData(obj, "counts"))
  ref.emb <- Embeddings(obj, "pca")
  f <- tempfile(fileext = ".h5")
  SaveSeuratH5(obj, f, verbose = FALSE)
  f2 <- tempfile(fileext = ".h5")
  file.rename(f, f2)
  loaded <- LoadSeuratH5(f2, verbose = FALSE)
  expect_s4_class(LayerData(loaded, "counts"), "DelayedMatrix")
  expect_equal(as.matrix(LayerData(loaded, "counts")), ref.counts)
  expect_equal(as.character(loaded$grp), as.character(obj$grp))
  expect_equal(loaded$flag, obj$flag)
  expect_equal(Embeddings(loaded, "pca"), ref.emb, check.attributes = FALSE)
})
