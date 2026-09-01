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

test_that("SaveSeurat preserves a full analysed object across a move", {
  skip_if_not_installed("BPCells")
  store <- tempfile(pattern = "bpcells-store")
  built <- .portable_bpcells_object(dir = store, m = .portable_counts(nfeat = 200, ncell = 120))
  obj <- built$object
  obj <- suppressWarnings(expr = NormalizeData(object = obj, verbose = FALSE))
  obj <- suppressWarnings(expr = FindVariableFeatures(object = obj, nfeatures = 50, verbose = FALSE))
  obj <- suppressWarnings(expr = ScaleData(object = obj, verbose = FALSE))
  obj <- suppressWarnings(expr = RunPCA(object = obj, npcs = 10, verbose = FALSE))
  obj <- suppressWarnings(expr = FindNeighbors(object = obj, dims = 1:10, verbose = FALSE))
  obj <- suppressWarnings(expr = FindClusters(object = obj, verbose = FALSE))
  # a second assay, and metadata of several types
  adt <- .portable_counts(nfeat = 12, ncell = 120)
  rownames(x = adt) <- paste0("adt", seq_len(length.out = 12))
  colnames(x = adt) <- colnames(x = obj)
  obj[["ADT"]] <- CreateAssayObject(counts = adt)
  obj$score <- seq_len(length.out = ncol(x = obj)) / 10
  obj$label <- factor(x = rep(c("p", "q"), length.out = ncol(x = obj)))
  obj$flag <- rep(c(TRUE, FALSE), length.out = ncol(x = obj))

  ref.counts <- as.matrix(x = LayerData(object = obj, layer = "counts"))
  ref.emb <- Embeddings(object = obj, reduction = "pca")
  ref.load <- Loadings(object = obj, reduction = "pca")
  ref.vf <- VariableFeatures(object = obj)
  ref.idents <- Idents(object = obj)

  f1 <- tempfile(fileext = ".seurat")
  SaveSeurat(object = obj, file = f1, verbose = FALSE)
  f2 <- tempfile(fileext = ".seurat")
  expect_true(object = file.rename(from = f1, to = f2))
  unlink(x = store, recursive = TRUE)
  loaded <- LoadSeurat(file = f2, dir = tempfile(pattern = "bundle"), verbose = FALSE)

  expect_setequal(object = Assays(object = loaded), expected = c("RNA", "ADT"))
  expect_equal(
    object = as.matrix(x = LayerData(object = loaded, layer = "counts", assay = "RNA")),
    expected = ref.counts
  )
  expect_equal(
    object = as.matrix(x = LayerData(object = loaded, layer = "counts", assay = "ADT")),
    expected = as.matrix(x = adt)
  )
  expect_equal(object = Embeddings(object = loaded, reduction = "pca"), expected = ref.emb)
  expect_equal(object = Loadings(object = loaded, reduction = "pca"), expected = ref.load)
  expect_equal(object = VariableFeatures(object = loaded), expected = ref.vf)
  expect_equal(object = Idents(object = loaded), expected = ref.idents)
  expect_equal(object = loaded$score, expected = obj$score)
  expect_equal(object = as.character(x = loaded$label), expected = as.character(x = obj$label))
  expect_equal(object = loaded$flag, expected = obj$flag)
  # the reloaded object is still usable, not just readable
  expect_silent(object = suppressWarnings(expr = RunPCA(object = loaded, npcs = 5, verbose = FALSE)))
})

test_that("SaveSeurat bundles every on-disk layer of a split assay", {
  skip_if_not_installed("BPCells")
  store <- tempfile(pattern = "bpcells-store")
  built <- .portable_bpcells_object(dir = store, m = .portable_counts(nfeat = 80, ncell = 60))
  obj <- built$object
  obj$batch <- rep(c("a", "b"), each = 30)
  obj[["RNA"]] <- split(x = obj[["RNA"]], f = obj$batch)
  expect_gt(object = length(x = Layers(object = obj)), expected = 1)
  ref <- lapply(
    X = Layers(object = obj),
    FUN = function(l) as.matrix(x = LayerData(object = obj, layer = l))
  )
  names(x = ref) <- Layers(object = obj)

  f <- tempfile(fileext = ".seurat")
  SaveSeurat(object = obj, file = f, verbose = FALSE)
  unlink(x = store, recursive = TRUE)
  loaded <- LoadSeurat(file = f, dir = tempfile(pattern = "bundle"), verbose = FALSE)

  expect_setequal(object = Layers(object = loaded), expected = names(x = ref))
  for (lyr in names(x = ref)) {
    expect_true(object = inherits(
      x = LayerData(object = loaded, layer = lyr), what = "IterableMatrix"
    ))
    expect_equal(
      object = as.matrix(x = LayerData(object = loaded, layer = lyr)),
      expected = ref[[lyr]]
    )
  }
})
