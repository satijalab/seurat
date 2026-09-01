set.seed(42)

# Reductions failed inside their decompositions, with messages about nu and nv,
# or with "subscript out of bounds", when asked for more than the data holds
data("pbmc_small", package = "SeuratObject", envir = environment())
scaled <- suppressWarnings(ScaleData(pbmc_small, features = rownames(pbmc_small), verbose = FALSE))

test_that("PCA on too few features says so", {
  data <- LayerData(scaled, layer = "scale.data")
  expect_error(
    suppressWarnings(RunPCA(data[1, , drop = FALSE], npcs = 2, verbose = FALSE)),
    "Cannot compute principal components from 1 features"
  )
  expect_error(
    suppressWarnings(RunPCA(data[1:2, , drop = FALSE], npcs = 2, verbose = FALSE)),
    "at least 3 are needed"
  )
  expect_s4_class(
    suppressWarnings(RunPCA(data[1:3, , drop = FALSE], npcs = 2, verbose = FALSE)),
    "DimReduc"
  )
})

test_that("features that were never scaled are reported", {
  unscaled <- setdiff(rownames(pbmc_small), rownames(LayerData(pbmc_small, layer = "scale.data")))
  skip_if(length(unscaled) == 0)
  expect_error(
    suppressWarnings(RunPCA(pbmc_small, features = unscaled[1], npcs = 2, verbose = FALSE)),
    "None of the requested features are in the scale.data layer"
  )
})

test_that("UMAP says when dims go past the reduction", {
  skip_if_not_installed("uwot")
  computed <- ncol(pbmc_small[["pca"]])
  expect_error(
    suppressWarnings(RunUMAP(pbmc_small, dims = 1:(computed + 10), verbose = FALSE)),
    "More dimensions specified in dims than have been computed"
  )
})

test_that("tSNE says when dims go past the reduction", {
  skip_if_not_installed("Rtsne")
  computed <- ncol(pbmc_small[["pca"]])
  expect_error(
    suppressWarnings(RunTSNE(pbmc_small, dims = 1:(computed + 10), verbose = FALSE)),
    "More dimensions specified in dims than have been computed"
  )
})

test_that("dims within the reduction still run", {
  skip_if_not_installed("Rtsne")
  expect_s4_class(
    suppressWarnings(RunTSNE(pbmc_small, dims = 1:5, perplexity = 10, check_duplicates = FALSE, verbose = FALSE))[["tsne"]],
    "DimReduc"
  )
})
