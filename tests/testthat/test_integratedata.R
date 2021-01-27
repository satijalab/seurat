# Tests for integration related fxns
set.seed(42)
pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))

# Setup test objects
ref <- pbmc_small
ref <- FindVariableFeatures(object = ref, verbose = FALSE, nfeatures = 100)
query <- CreateSeuratObject(
  counts = GetAssayData(object = pbmc_small[['RNA']], slot = "counts") + rpois(n = ncol(pbmc_small), lambda = 1)
)
query2 <- CreateSeuratObject(
  counts = GetAssayData(object = pbmc_small[['RNA']], slot = "counts")[, 1:40] + rpois(n = ncol(pbmc_small), lambda = 1)
)
query.list <- list(query, query2)
query.list <- lapply(X = query.list, FUN = NormalizeData, verbose = FALSE)
query.list <- lapply(X = query.list, FUN = FindVariableFeatures, verbose = FALSE, nfeatures = 100)
query.list <- lapply(X = query.list, FUN = ScaleData, verbose = FALSE)
query.list <- suppressWarnings(lapply(X = query.list, FUN = RunPCA, verbose = FALSE, npcs = 20))

anchors2 <- suppressMessages(suppressWarnings(FindIntegrationAnchors(object.list = c(ref, query.list[[1]]), k.filter = NA, verbose = FALSE)))
anchors3 <- suppressMessages(suppressWarnings(FindIntegrationAnchors(object.list = c(ref, query.list), k.filter = NA, verbose = FALSE)))

# Tests for IntegrateEmbeddings
# ------------------------------------------------------------------------------
# context("IntegrateEmbeddings")

# test_that("IntegrateEmbeddings validates properly", {
#   expect_error(IntegrateEmbeddings(anchorset = anchors2))
#   expect_error(IntegrateEmbeddings(anchorset = anchors2, reduction = "pca", k.weight = 100))
#   expect_error(IntegrateEmbeddings(anchorset = anchors2, reduction = c("pca", "pca2"), k.weight = 40))
#   expect_error(IntegrateEmbeddings(anchorset = anchors2, reduction = "pca", k.weight = 40, weight.reduction = c(ref[['pca']])))
#   pca3 <- RenameCells(object = ref[['pca']], new.names = paste0(Cells(ref), "_test"))
#   expect_error(IntegrateEmbeddings(anchorset = anchors2, reduction = "pca", k.weight = 40,
#                                    weight.reduction = c(pca3, ref[['pca']])))
# })
#
# test_that("IntegrateEmbeddings with two objects default works", {
#   skip_on_cran()
#   int2 <- IntegrateEmbeddings(anchorset = anchors2, reduction = "pca", k.weight = 40, verbose = FALSE)
#   expect_equal(Reductions(int2), "integrated_pca")
#   expect_equal(sum(Embeddings(int2[['integrated_pca']])[1,]), -3.13050872287, tolerance = 1e-6)
#   expect_equal(sum(Embeddings(int2[['integrated_pca']])[,1]), -5.78790844887, tolerance = 1e-6)
# })
#
# test_that("IntegrateEmbeddings with three objects default works", {
#   skip_on_cran()
#   int3 <- IntegrateEmbeddings(anchorset = anchors3, reduction = "pca", k.weight = 40, verbose = FALSE)
#   expect_equal(Reductions(int3), "integrated_pca")
#   expect_equal(sum(Embeddings(int3[['integrated_pca']])[1,]), 0.221867815987, tolerance = 1e-6)
#   expect_equal(sum(Embeddings(int3[['integrated_pca']])[,1]), -16.7881409595, tolerance = 1e-6)
# })
#
# test_that("IntegrateEmbeddings works with specified reference objects", {
#   skip_on_cran()
#   anchors4 <- suppressMessages(suppressWarnings(FindIntegrationAnchors(object.list = c(ref, query.list), k.filter = NA, verbose = FALSE, reference = 1)))
#   int4 <- IntegrateEmbeddings(anchorset = anchors4, reduction = "pca", k.weight = 40, verbose = FALSE)
#   expect_equal(Reductions(int4), "integrated_pca")
#   expect_equal(sum(Embeddings(int4[['integrated_pca']])[1,]), -3.13050872287, tolerance = 1e-6)
#   expect_equal(sum(Embeddings(int4[['integrated_pca']])[,1]), 13.1180105492, tolerance = 1e-6)
# })

# Tests for IntegrateData
# ------------------------------------------------------------------------------
context("IntegrateData")

test_that("IntegrateData with two objects default work", {
  expect_error(IntegrateData(anchorset = anchors2))
  int2 <- IntegrateData(anchorset = anchors2, k.weight = 50, verbose = FALSE)
  expect_true(all(Assays(int2) %in% c("integrated", "RNA")))
  expect_equal(Tool(int2), "Integration")
  expect_equal(dim(int2[["integrated"]]), c(133, 160))
  expect_equal(length(VariableFeatures(int2)), 133)
  expect_equal(GetAssayData(int2[["integrated"]], slot = "counts"), new("dgCMatrix"))
  expect_equal(GetAssayData(int2[['integrated']], slot = "scale.data"), matrix())
  expect_equal(sum(GetAssayData(int2[["integrated"]])[1, ]), 44.97355, tolerance = 1e-3)
  expect_equal(sum(GetAssayData(int2[["integrated"]])[, 1]), 78.8965706046, tolerance = 1e-6)
  expect_equal(Tool(object = int2, slot = "Integration")@sample.tree, matrix(c(-1, -2), nrow  = 1))
})

test_that("IntegrateData with three objects default work", {
  expect_error(IntegrateData(anchorset = anchors3, k.weight = 50))
  int3 <- IntegrateData(anchorset = anchors3, k.weight = 25, verbose = FALSE)
  expect_true(all(Assays(int3) %in% c("integrated", "RNA")))
  expect_equal(Tool(int3), "Integration")
  expect_equal(dim(int3[["integrated"]]), c(169, 200))
  expect_equal(length(VariableFeatures(int3)), 169)
  expect_equal(GetAssayData(int3[["integrated"]], slot = "counts"), new("dgCMatrix"))
  expect_equal(GetAssayData(int3[['integrated']], slot = "scale.data"), matrix())
  expect_equal(sum(GetAssayData(int3[["integrated"]])[1, ]), 372.829, tolerance = 1e-6)
  expect_equal(sum(GetAssayData(int3[["integrated"]])[, 1]), 482.5009, tolerance = 1e-6)
  expect_equal(Tool(object = int3, slot = "Integration")@sample.tree, matrix(c(-2, -3, 1, -1), nrow  = 2, byrow = TRUE))
})

test_that("Input validates correctly ", {
  expect_error(anchorset = anchors2, k.weight = 50, features.to.integrate = "BAD")
  expect_error(IntegrateData(anchorset = anchors2, k.weight = 50, normalization.method = "BAD"))
  expect_error(IntegrateData(anchorset = anchors2, k.weight = 50, weight.reduction = "BAD"))
  expect_error(IntegrateData(anchorset = anchors2, reductions.to.integrate = "pca"))
  skip_on_cran()
  #expect_warning(IntegrateData(anchorset = anchors2, k.weight = 50, features = c(rownames(ref), "BAD")))
  #expect_warning(IntegrateData(anchorset = anchors2, k.weight = 50, dims = 1:1000))
})


