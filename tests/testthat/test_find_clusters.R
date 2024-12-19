path_to_counts <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")


# Builds a `Seurat` instance and annotates it with the requisite data
# structures for running `FindClusters` (i.e. a shared-nearest-neighbor
# (SNN) graph).
get_test_data <- function() {
  raw_counts <- read.table(path_to_counts, sep = "\t", row.names = 1)
  counts <- as.sparse(as.matrix(raw_counts))
  assay <- CreateAssay5Object(counts)
  test_data <- CreateSeuratObject(assay)

  test_data <- NormalizeData(test_data, verbose = FALSE)
  test_data <- FindVariableFeatures(test_data, verbose = FALSE)
  test_data <- ScaleData(test_data, verbose = FALSE)
  # Reduce number of PCs to avoid warning from `irlba` caused by the
  # small size of the dataset being used. Plus, we only want to build our
  # SNN graph using the first 10 PCs to get "interesting" clustering results.
  test_data <- RunPCA(test_data, npcs = 10, verbose = FALSE)
  test_data <- FindNeighbors(
    test_data,
    k.param=10,
    dims = 1:10,
    verbose = FALSE
  )

  return(test_data)
}


context("FindClusters")


test_that("Smoke test for `FindClusters`", {
  test_case <- get_test_data()

  # Spot check cluster assignments with using defaults.
  results <- FindClusters(test_case)$seurat_clusters
  expect_equal(results[[1]], factor(3, levels=0:5))
  expect_equal(results[[15]], factor(4, levels=0:5))
  expect_equal(results[[24]], factor(0, levels=0:5))
  expect_equal(results[[72]], factor(5, levels=0:5))
  expect_equal(results[[length(results)]], factor(2, levels=0:5))

  # Check that every clustering algorithm can be run without errors.
  expect_no_error(FindClusters(test_case, algorithm = 1))
  expect_no_error(FindClusters(test_case, algorithm = 2))
  expect_no_error(FindClusters(test_case, algorithm = 3))
  # The leiden algorithm requires that `random.seed` be greater than 0,
  # which is the default for `FindClusters` so a warning should be raised.
  expect_warning(FindClusters(test_case, algorithm = 4))
  expect_no_warning(FindClusters(test_case, algorithm = 4, random.seed = 1))
})
