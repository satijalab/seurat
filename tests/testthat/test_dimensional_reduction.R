context("test-dimensional_reduction")

set.seed(seed = 1)
dummyexpMat <- matrix(
  data = sample(x = c(1:50), size = 1e4, replace = TRUE),
  ncol = 100, nrow = 100
)
colnames(x = dummyexpMat) <- paste0("cell", seq(ncol(x = dummyexpMat)))
row.names(x = dummyexpMat) <- paste0("gene", seq(nrow(x = dummyexpMat)))

# Create Seurat object for testing
obj <- CreateSeuratObject(counts = as.sparse(dummyexpMat))


test_that("different ways of passing distance matrix", {
  # Manually make a distance object to test
  distMat <- dist(t(dummyexpMat))

  expect_equivalent(
    suppressWarnings(expr = RunTSNE(obj, distance.matrix = distMat)),
    suppressWarnings(expr = RunTSNE(obj, distance.matrix = as.matrix(distMat)))
  )
  expect_equivalent(
    suppressWarnings(expr = RunTSNE(obj, distance.matrix = distMat)@reductions$tsne),
    suppressWarnings(expr = RunTSNE(distMat, assay = "RNA"))
  )
  expect_equivalent(
    suppressWarnings(expr = RunTSNE(obj, distance.matrix = distMat)@reductions$tsne),
    suppressWarnings(expr = RunTSNE(as.matrix(distMat), assay = "RNA", is_distance = TRUE))
  )
})

# Normalize, scale, and compute PCA, using RunPCA
obj <- NormalizeData(object = obj, verbose = FALSE)
obj <- ScaleData(object = obj, verbose = FALSE)

pca_result <- suppressWarnings(expr = RunPCA(
  object = obj,
  features = rownames(obj[['RNA']]$counts),
  verbose = FALSE
))

test_that("pca returns total variance (see #982)", {
  # Using stats::prcomp
  scaled_data <- LayerData(object = obj, layer = "scale.data")
  prcomp_result <- stats::prcomp(scaled_data, center = FALSE, scale. = FALSE)

  # Compare
  expect_equivalent(slot(object = pca_result[["pca"]], name = "misc")$total.variance,
                    sum(prcomp_result$sdev^2))

})

test_that("pca is equivalent for BPCells", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  library(BPCells)
  library(Matrix)
  mat_bpcells <- t(x = as(object = t(x = obj[['RNA']]$counts ), Class = "IterableMatrix"))
  obj[['RNAbp']] <- CreateAssay5Object(counts = mat_bpcells)
  DefaultAssay(obj) <- "RNAbp"
  obj <- NormalizeData(object = obj, verbose = FALSE)
  obj <- ScaleData(object = obj, verbose=FALSE)
  pca_result_bp <- suppressWarnings(expr = RunPCA(
    object = obj,
    features = rownames(obj[['RNAbp']]$counts),
    assay = "RNAbp"))
  expect_equivalent(abs(pca_result_bp[['pca']]@cell.embeddings),
                   abs(pca_result[['pca']]@cell.embeddings),
                   tolerance = 1e-5)
})
