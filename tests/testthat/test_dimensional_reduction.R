context("test-dimensional_reduction")

test_that("different ways of passing distance matrix", {
  # Generate dummy data exp matrix
  set.seed(1)
  dummyexpMat <- matrix(data = sample(x = c(1:50), size = 1e4, replace = TRUE),
                        ncol = 100, nrow = 100)
  colnames(dummyexpMat) <- paste0("cell", seq(ncol(dummyexpMat)))
  row.names(dummyexpMat) <- paste0("gene", seq(nrow(dummyexpMat)))

  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)

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

test_that("pca returns total variance (see #982)", {
  # Generate dummy data exp matrix
  set.seed(seed = 1)
  dummyexpMat <- matrix(
    data = sample(x = c(1:50), size = 1e4, replace = TRUE),
    ncol = 100, nrow = 100
  )
  colnames(x = dummyexpMat) <- paste0("cell", seq(ncol(x = dummyexpMat)))
  row.names(x = dummyexpMat) <- paste0("gene", seq(nrow(x = dummyexpMat)))

  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)

  # Scale and compute PCA, using RunPCA
  obj <- ScaleData(object = obj, verbose = FALSE)
  pca_result <- suppressWarnings(expr = RunPCA(
    object = obj,
    features = rownames(x = obj),
    verbose = FALSE
  ))

  # Using stats::prcomp
  scaled_data <- Seurat::GetAssayData(object = obj, slot = "scale.data")
  prcomp_result <- stats::prcomp(scaled_data, center = FALSE, scale. = FALSE)

  # Compare
  expect_equivalent(slot(object = pca_result[["pca"]], name = "misc")$total.variance,
                    sum(prcomp_result$sdev^2))

})
