source(test_path("../testdata/test-objects.R"), local = TRUE)

context("FindClusters")

# Builds a `Seurat` instance and annotates it with the requisite data
# structures for running `FindClusters` (i.e. a shared-nearest-neighbor
# (SNN) graph).
get_findclusters_test_data <- function() {
  counts <- create_pbmc_counts()
  object <- create_seurat_obj(counts = counts)
  create_prep_obj(
    object = object,
    normalize = TRUE,
    variable.features = TRUE,
    scale = TRUE,
    pca = TRUE,
    neighbors = TRUE,
    nfeatures = 20L,
    npcs = 10L,
    dims = 1:10,
    k.param = 10L
  )
}


test_that("`ComputeSNN` handles direct edge cases", {
  nn.ranked <- matrix(
    data = c(
      1L, 1L, 3L, 3L,
      2L, 2L, 4L, 4L
    ),
    nrow = 4,
    ncol = 2
  )

  snn <- ComputeSNN(nn_ranked = nn.ranked, prune = 0, nthreads = 1L)
  expected <- Matrix::bdiag(
    matrix(1, nrow = 2, ncol = 2),
    matrix(1, nrow = 2, ncol = 2)
  )
  expect_equal(as.matrix(snn), as.matrix(expected))

  snn.threaded <- ComputeSNN(nn_ranked = nn.ranked, prune = 0, nthreads = 4L)
  expect_equal(snn, snn.threaded)

  pruned <- ComputeSNN(nn_ranked = nn.ranked, prune = 1.1, nthreads = 2L)
  expect_equal(length(pruned@x), 0)

  invalid <- nn.ranked
  invalid[1, 1] <- 5L
  expect_error(
    ComputeSNN(nn_ranked = invalid, prune = 0, nthreads = 2L),
    regexp = "outside the valid cell range"
  )
})

test_that("Smoke test for `FindClusters`", {
  test_case <- get_findclusters_test_data()

  # Validate cluster assignments using default parameters.
  results <- FindClusters(test_case)$seurat_clusters
  # Check that every cell was assigned to a cluster label.
  expect_false(any(is.na(results)))
  # Check that the expected cluster labels were assigned.
  expect_equal(as.numeric(levels(results)), c(0, 1, 2, 3, 4))
  # Check that the cluster sizes match the expected distribution.
  expect_equal(
    as.numeric(sort(table(results))),
    c(8, 9, 18, 20, 25)
  )

  # Check that every clustering algorithm can be run without errors.
  expect_no_error(FindClusters(test_case, algorithm = 1))
  expect_no_error(FindClusters(test_case, algorithm = 2))
  expect_no_error(FindClusters(test_case, algorithm = 3))
  # The leiden algorithm requires that `random.seed` be greater than 0,
  # so the default `FindClusters` seed warns and is reset.
  # Test with igraph method
  expect_warning(FindClusters(test_case, algorithm = 4, leiden_method = "igraph"))
  expect_no_warning(FindClusters(test_case, algorithm = 4, leiden_method = "igraph", random.seed = 1))
  
  # Test leidenbase method if available
  skip_if_not_installed("leidenbase")
  expect_warning(FindClusters(test_case, algorithm = 4, leiden_method = "leidenbase"))
  expect_no_warning(FindClusters(test_case, algorithm = 4, leiden_method = "leidenbase", random.seed = 1))
})

test_that("`FindClusters` works if passed a vector of resolutions", {
  test_case <- get_findclusters_test_data()
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

  clustered_loop <- get_findclusters_test_data()
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
  test_case <- get_findclusters_test_data()
  clustered_default <- FindClusters(test_case, verbose = FALSE)
  check_numeric_levels(clustered_default[["seurat_clusters", drop = TRUE]])

  # and also for the default RNA_snn_res column
  default_res_col <- grep("RNA_snn_res", colnames(clustered_default[[]]), value = TRUE)[1]
  expect_true(!is.na(default_res_col))
  check_numeric_levels(clustered_default[[default_res_col, drop = TRUE]])

  # Test factor levels for custom cluster name
  test_case2 <- get_findclusters_test_data()
  clustered_custom <- FindClusters(test_case2, cluster.name = "custom_clusters", verbose = FALSE)
  check_numeric_levels(clustered_custom[["custom_clusters", drop = TRUE]])
})

test_that("`FindNeighbors` supports every nn.method", {
  test_case <- get_findclusters_test_data()
  for (method in c("rann", "annoy")) {
    neighbors <- FindNeighbors(
      test_case,
      graph.name = c(paste0("nn_", method), paste0("snn_", method)),
      k.param = 10,
      dims = 1:10,
      nn.method = method,
      verbose = FALSE
    )
    expect_s4_class(neighbors[[paste0("snn_", method)]], "Graph")
  }
  expect_error(
    FindNeighbors(test_case, k.param = 10, dims = 1:10, nn.method = "bogus"),
    "Invalid method"
  )
})

test_that("`FindNeighbors` forwards dots to the neighbor search", {
  # Backend-specific nearest-neighbor arguments are passed through dots.
  embeddings <- Embeddings(get_findclusters_test_data(), reduction = "pca")[, 1:10]

  # Annoy receives search.k from dots.
  expect_error(
    FindNeighbors(
      embeddings, k.param = 10, nn.method = "annoy", search.k = 1,
      return.neighbor = TRUE, verbose = FALSE
    ),
    "Raise search.k"
  )

  expect_no_warning(
    radius <- FindNeighbors(
      embeddings, k.param = 10, nn.method = "rann", searchtype = "radius",
      radius = 0.01, return.neighbor = TRUE, verbose = FALSE
    )
  )
  expect_true(any(Indices(radius) == 0))
  expect_error(
    FindNeighbors(
      embeddings, k.param = 10, nn.method = "rann", searchtype = "radius",
      radius = 0.01, verbose = FALSE
    ),
    "RANN radius search can return fewer than k.param neighbors"
  )

  # Unknown dots are still reported.
  expect_warning(
    FindNeighbors(embeddings, k.param = 10, not_a_real_argument = 1, verbose = FALSE),
    "not used"
  )
})

test_that("`AnnoySearch` matches the RcppAnnoy query loop across metrics and thread counts", {
  old.threads <- getThreads(verbose = FALSE)
  on.exit(setThreads(old.threads, verbose = FALSE), add = TRUE)

  set.seed(123)
  float.data <- matrix(rnorm(160), nrow = 32)
  rownames(float.data) <- paste0("cell", seq_len(nrow(float.data)))
  float.query <- float.data[seq(2, 24, by = 3), , drop = FALSE]
  hamming.data <- matrix(sample(0:1, 32 * 8, replace = TRUE), nrow = 32)
  rownames(hamming.data) <- rownames(float.data)
  hamming.query <- hamming.data[seq(2, 24, by = 3), , drop = FALSE]

  cases <- list(
    euclidean = list(data = float.data, query = float.query),
    cosine = list(data = float.data, query = float.query),
    manhattan = list(data = float.data, query = float.query),
    hamming = list(data = hamming.data, query = hamming.query)
  )

  for (metric in names(x = cases)) {
    index <- AnnoyBuildIndex(
      data = cases[[metric]]$data,
      metric = metric,
      n.trees = 10
    )
    expected <- AnnoySearchR(
      index = index,
      query = cases[[metric]]$query,
      k = 6,
      search.k = -1,
      include.distance = TRUE
    )

    for (nthreads in c(1, 2, 4)) {
      setThreads(nthreads, verbose = FALSE)
      observed <- AnnoySearch(
        index = index,
        query = cases[[metric]]$query,
        k = 6,
        search.k = -1,
        include.distance = TRUE
      )
      expect_identical(
        observed$nn.idx,
        expected$nn.idx,
        info = paste(metric, "indices", nthreads, "threads")
      )
      if (metric == "hamming") { # hamming metric uses exact integer distances
        expect_identical(
          observed$nn.dists,
          expected$nn.dists,
          info = paste(metric, "distances", nthreads, "threads")
        )
      } else { # others use float distances which may not be exactly identical
        expect_equal(
          observed$nn.dists,
          expected$nn.dists,
          tolerance = 1e-6,
          info = paste(metric, "distances", nthreads, "threads")
        )
      }
    }
  }
})

test_that("`AnnoySearchParallel` falls back on file handoff failures", {
  set.seed(456)
  query <- matrix(rnorm(40), nrow = 8)

  save.failures <- list(
    no.permission = function(path) stop("Permission denied"),
    disk.full = function(path) stop("No space left on device"),
    empty.file = function(path) file.create(path),
    missing.file = function(path) { file.create(path); unlink(path) }
  )

  for (failure in names(x = save.failures)) {
    index <- structure(
      list(save = save.failures[[failure]]),
      class = "Rcpp_AnnoyEuclidean"
    )
    capture.output(
      result <- AnnoySearchParallel(index = index, query = query, k = 3, nthreads = 2),
      type = "message"
    )
    expect_null(result, info = failure)
  }
})
