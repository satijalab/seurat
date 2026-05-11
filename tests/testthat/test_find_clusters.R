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

  # Validate cluster assignments using default parameters.
  results <- FindClusters(test_case)$seurat_clusters
  # Check that every cell was assigned to a cluster label.
  expect_false(any(is.na(results)))
  # Check that the expected cluster labels were assigned.
  expect_equal(as.numeric(levels(results)), c(0, 1, 2, 3, 4, 5))
  # Check that the cluster sizes match the expected distribution.
  expect_equal(
    as.numeric(sort(table(results))),
    c(9, 10, 10, 11, 20, 20)
  )

  # Check that every clustering algorithm can be run without errors.
  expect_no_error(FindClusters(test_case, algorithm = 1))
  expect_no_error(FindClusters(test_case, algorithm = 2))
  expect_no_error(FindClusters(test_case, algorithm = 3))
  # The leiden algorithm requires that `random.seed` be greater than 0,
  # which is the default for `FindClusters` so a warning should be raised.
  # Test with igraph method
  expect_warning(FindClusters(test_case, algorithm = 4, leiden_method = "igraph"))
  expect_no_warning(FindClusters(test_case, algorithm = 4, leiden_method = "igraph", random.seed = 1))
  
  # Test leidenbase method if available
  skip_if_not_installed("leidenbase")
  expect_warning(FindClusters(test_case, algorithm = 4, leiden_method = "leidenbase"))
  expect_no_warning(FindClusters(test_case, algorithm = 4, leiden_method = "leidenbase", random.seed = 1))
})

test_that("`FindClusters` works if passed a vector of resolutions", {
  test_case <- get_test_data()
  resolutions <- seq(0.4, 0.8, by = 0.1)
  cluster_names <- paste0("resolution_", resolutions)

  # Run FindClusters with multiple resolutions in one call
  clustered <- FindClusters(
    test_case,
    resolution = resolutions,
    cluster.name = cluster_names,
    verbose = FALSE
  )

  # Check that the active identity is set to the last computed clustering
  expect_identical(
    as.character(clustered[[tail(cluster_names, n = 1), drop = TRUE]]),
    as.character(Idents(clustered))
  )

  clustered_loop <- get_test_data()
  for (i in seq_along(resolutions)) {
    clustered_loop <- FindClusters(
      clustered_loop,
      resolution = resolutions[i],
      cluster.name = cluster_names[i],
      verbose = FALSE
    )
  }

  # Check that passing multiple resolutions at once produces identical results to running in a loop
  for (cluster_name in cluster_names) {
    expect_identical(
      as.character(clustered[[cluster_name, drop = TRUE]]),
      as.character(clustered_loop[[cluster_name, drop = TRUE]])
    )
  }
})

test_that("`FindClusters` sorts numeric factor levels correctly", {
  # Helper to verify numeric levels are sorted numerically, not lexicographically
  # We want (1, 2, ..., 9, 10) instead of (1, 10, 2, ... 9)
  check_numeric_levels <- function(factor_col) {
    levels_str <- as.character(levels(factor_col))
    numeric_levels <- levels_str[grepl("^[0-9]+$", levels_str)]
    if (length(numeric_levels) > 0) {
      numeric_values <- as.integer(numeric_levels)
      expect_equal(numeric_levels, as.character(sort(numeric_values)))
    }
  }

  # Test factor levels for default cluster name
  test_case <- get_test_data()
  clustered_default <- FindClusters(test_case, verbose = FALSE)
  check_numeric_levels(clustered_default[["seurat_clusters", drop = TRUE]])

  # and also for the default RNA_snn_res column
  default_res_col <- grep("RNA_snn_res", colnames(clustered_default[[]]), value = TRUE)[1]
  expect_true(!is.na(default_res_col))
  check_numeric_levels(clustered_default[[default_res_col, drop = TRUE]])

  # Test factor levels for custom cluster name
  test_case2 <- get_test_data()
  clustered_custom <- FindClusters(test_case2, cluster.name = "custom_clusters", verbose = FALSE)
  check_numeric_levels(clustered_custom[["custom_clusters", drop = TRUE]])
})
