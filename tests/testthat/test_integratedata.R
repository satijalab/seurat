# Tests for integration related fxns
set.seed(42)
pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))

# Setup test objects
ref <- pbmc_small
ref <- FindVariableFeatures(object = ref, verbose = FALSE, nfeatures = 100)
query <- CreateSeuratObject(
  counts = as.sparse(
    GetAssayData(
      object = pbmc_small[['RNA']],
      layer = "counts") + rpois(n = ncol(pbmc_small),
      lambda = 1
    )
  )
)
query2 <- CreateSeuratObject(
  counts = as.sparse(
    LayerData(
      object = pbmc_small[['RNA']],
      layer = "counts")[, 1:40] + rpois(n = ncol(pbmc_small),
      lambda = 1
    )
  )
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
  expect_equal(GetAssayData(int2[["integrated"]], layer = "counts"), new("dgCMatrix"))
  expect_equal(GetAssayData(int2[['integrated']], layer = "scale.data"), matrix())
  expect_equal(sum(GetAssayData(int2[["integrated"]], layer = "data")[1, ]), 44.97355, tolerance = 1e-3)
  expect_equal(sum(GetAssayData(int2[["integrated"]], layer = "data")[, 1]), 78.8965706046, tolerance = 1e-6)
  expect_equal(Tool(object = int2, slot = "Integration")@sample.tree, matrix(c(-1, -2), nrow  = 1))
})

test_that("IntegrateData with three objects default work", {
  expect_error(IntegrateData(anchorset = anchors3, k.weight = 50))
  int3 <- IntegrateData(anchorset = anchors3, k.weight = 25, verbose = FALSE)
  expect_true(all(Assays(int3) %in% c("integrated", "RNA")))
  expect_equal(Tool(int3), "Integration")
  expect_equal(dim(int3[["integrated"]]), c(169, 200))
  expect_equal(length(VariableFeatures(int3)), 169)
  expect_equal(GetAssayData(int3[["integrated"]], layer = "counts"), new("dgCMatrix"))
  expect_equal(GetAssayData(int3[['integrated']], layer = "scale.data"), matrix())
  expect_equal(sum(GetAssayData(int3[["integrated"]], layer = "data")[1, ]), 372.829, tolerance = 1e-6)
  expect_equal(sum(GetAssayData(int3[["integrated"]], layer = "data")[, 1]), 482.5009, tolerance = 1e-6)
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

# Tests for IntegrateLayers
# ------------------------------------------------------------------------------
context("IntegrateLayers")
pbmc_small[['RNAv5']] <- CreateAssay5Object(counts = LayerData(pbmc_small[['RNA']], layer = "counts"))

pbmc_small[["RNAv5"]] <- split(pbmc_small[["RNAv5"]], f = pbmc_small$groups)
DefaultAssay(pbmc_small) <- "RNAv5"
pbmc_small <- NormalizeData(pbmc_small)
pbmc_small <- FindVariableFeatures(pbmc_small)
pbmc_small <- ScaleData(pbmc_small)
pbmc_small <- suppressMessages(suppressWarnings(RunPCA(pbmc_small)))


test_that("IntegrateLayers does not work on a v3 assay ", {
  expect_error(IntegrateLayers(object = pbmc_small, method = CCAIntegration,
                               orig.reduction = "pca",
                               assay = "RNA",
                               new.reduction = "integrated.cca"))
})

test_that("IntegrateLayers errors out if incorrect input ", {
  expect_error(IntegrateLayers(object = pbmc_small, method = CCAIntegration,
                               orig.reduction = "pca",
                               assay = "DNA",
                               new.reduction = "integrated.cca"))
  expect_error(IntegrateLayers(object = pbmc_small, method = CCAIntegration,
                               orig.reduction = "lda",
                               new.reduction = "integrated.cca"))
})

#itegration methods
int_cca <- suppressMessages(suppressWarnings(IntegrateLayers(
  object = pbmc_small, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  k.weight=25,
  verbose = FALSE
)))
int_rpca <- suppressMessages(suppressWarnings(IntegrateLayers(
  object = pbmc_small, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  dims = 1:10,
  k.anchor = 10,
  k.weight=10,
  verbose = FALSE
)))

# int_mnn <- suppressMessages(suppressWarnings(IntegrateLayers(
#   object = pbmc_small, method = FastMNNIntegration,
#   new.reduction = "integrated.mnn",
#   k.weight=25,
#   verbose = FALSE
# )))


test_that("IntegrateLayers returns embeddings with correct dimensions ", {
  expect_equal(dim(int_cca[["integrated.cca"]]), c(80, 50))
  expect_equal(dim(int_rpca[["integrated.rpca"]]), c(80, 50))

  int_rpca
  expect_equal(int_cca[["integrated.cca"]]@assay.used, "RNAv5")
  #expect_equal(int_cca[['integrated.cca']]@cell.embeddings, c(3, 4, 5))
})

test_that("IntegrateLayers works with harmony", {
  skip_on_cran()
  skip_if_not_installed("harmony")
  int_harmony <- suppressMessages(suppressWarnings(IntegrateLayers(
    object = pbmc_small, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    k.weight=25,
    verbose = FALSE
  )))
  expect_equal(dim(int_harmony[["harmony"]]), c(80, 50))


})

test_that("group.by ", {
  expect_equal(dim(int_cca[["integrated.cca"]]), c(80, 50))
  expect_equal(int_cca[["integrated.cca"]]@assay.used, "RNAv5")
})


#Harmony integration
# int_2 <- IntegrateLayers(object = pbmc_small, method = CCAIntegration,
#                          group.by = "letter.idents",
#                 orig.reduction = "pca",
#                 assay = "RNAv5",
#                 k.weight = 20,
#                 new.reduction = "integrated.cca")
#
# head(int_2[['integrated.cca']]@cell.embeddings[1:5,1:5])
# head(int_cca[['integrated.cca']]@cell.embeddings[1:5,1:5])

