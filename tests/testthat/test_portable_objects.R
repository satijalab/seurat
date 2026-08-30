set.seed(42)

# Small counts matrix used to back the on-disk layers under test
.portable_counts <- function(nfeat = 60, ncell = 40, lambda = 1.3) {
  m <- as(
    object = matrix(data = rpois(n = nfeat * ncell, lambda = lambda),
                    nrow = nfeat, ncol = ncell),
    Class = "dgCMatrix"
  )
  dimnames(x = m) <- list(
    paste0("g", seq_len(length.out = nfeat)),
    paste0("c", seq_len(length.out = ncell))
  )
  return(m)
}

# Build a Seurat object whose counts layer is a BPCells on-disk matrix stored in
# `dir`; returns the object and the reference in-memory matrix
.portable_bpcells_object <- function(dir, m = .portable_counts()) {
  BPCells::write_matrix_dir(mat = m, dir = dir, overwrite = TRUE)
  mat <- BPCells::open_matrix_dir(dir = dir)
  obj <- suppressWarnings(expr = CreateSeuratObject(counts = mat))
  return(list(object = obj, reference = m))
}

context("On-disk layers and plain saveRDS")

test_that("a plain rds of a BPCells object breaks once the store is gone", {
  skip_if_not_installed("BPCells")
  store <- tempfile(pattern = "bpcells-store")
  built <- .portable_bpcells_object(dir = store)
  f <- tempfile(fileext = ".rds")
  saveRDS(object = built$object, file = f)
  # The rds only records the absolute path of the on-disk store, so sharing the
  # rds on its own -- or moving the store -- leaves the layer unreadable
  unlink(x = store, recursive = TRUE)
  loaded <- readRDS(file = f)
  expect_error(object = as.matrix(x = LayerData(object = loaded, layer = "counts")))
})

context("AsInMemory")

test_that("AsInMemory materializes BPCells layers for a portable rds", {
  skip_if_not_installed("BPCells")
  store <- tempfile(pattern = "bpcells-store")
  built <- .portable_bpcells_object(dir = store)
  obj <- built$object
  expect_true(object = inherits(x = LayerData(object = obj, layer = "counts"),
                                what = "IterableMatrix"))
  obj <- AsInMemory(object = obj, verbose = FALSE)
  expect_s4_class(object = LayerData(object = obj, layer = "counts"),
                  class = "dgCMatrix")
  expect_equal(object = LayerData(object = obj, layer = "counts"),
               expected = built$reference)
  # the object is now self-contained: it survives the store going away
  f <- tempfile(fileext = ".rds")
  saveRDS(object = obj, file = f)
  unlink(x = store, recursive = TRUE)
  loaded <- readRDS(file = f)
  expect_equal(object = LayerData(object = loaded, layer = "counts"),
               expected = built$reference)
})

test_that("AsInMemory leaves in-memory objects untouched", {
  m <- .portable_counts()
  obj <- suppressWarnings(expr = CreateSeuratObject(counts = m))
  obj2 <- AsInMemory(object = obj, verbose = FALSE)
  expect_equal(object = LayerData(object = obj2, layer = "counts"),
               expected = LayerData(object = obj, layer = "counts"))
})

context("Portable bundle (SaveSeurat / LoadSeurat)")

test_that("SaveSeurat bundles BPCells layers and survives a move", {
  skip_if_not_installed("BPCells")
  store <- tempfile(pattern = "bpcells-store")
  built <- .portable_bpcells_object(dir = store)
  obj <- suppressWarnings(expr = NormalizeData(object = built$object, verbose = FALSE))
  f1 <- tempfile(fileext = ".seurat")
  SaveSeurat(object = obj, file = f1, verbose = FALSE)
  # move the bundle and destroy the original store: the bundle must stand alone
  f2 <- tempfile(fileext = ".seurat")
  expect_true(object = file.rename(from = f1, to = f2))
  unlink(x = store, recursive = TRUE)
  loaded <- LoadSeurat(file = f2, dir = tempfile(pattern = "bundle"), verbose = FALSE)
  expect_true(object = inherits(x = LayerData(object = loaded, layer = "counts"),
                                what = "IterableMatrix"))
  expect_equal(
    object = as.matrix(x = LayerData(object = loaded, layer = "counts")),
    expected = as.matrix(x = built$reference)
  )
  expect_true(object = "data" %in% Layers(object = loaded))
})

test_that("SaveSeurat round-trips a fully in-memory object", {
  m <- .portable_counts()
  obj <- suppressWarnings(expr = CreateSeuratObject(counts = m))
  f <- tempfile(fileext = ".seurat")
  SaveSeurat(object = obj, file = f, verbose = FALSE)
  loaded <- LoadSeurat(file = f, dir = tempfile(pattern = "bundle"), verbose = FALSE)
  expect_equal(object = LayerData(object = loaded, layer = "counts"),
               expected = m)
  expect_equal(object = colnames(x = loaded), expected = colnames(x = obj))
})
