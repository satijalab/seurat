context("NNdist")

# Helper: manual Euclidean distance between two vectors
manual_euclid <- function(a, b) {
  sqrt(sum((a - b)^2))
}

# Helper: old per-cell NNdist implementation for equivalence testing
NNdist_old <- function(nn.idx, embeddings, metric = "euclidean",
                       query.embeddings = NULL, nearest.dist = NULL) {
  if (!is.list(x = nn.idx)) {
    nn.idx <- lapply(X = 1:nrow(x = nn.idx), FUN = function(x) nn.idx[x, ])
  }
  query.embeddings <- query.embeddings %||% embeddings
  nn.dist <- lapply(
    X = 1:nrow(x = query.embeddings),
    FUN = function(x) {
      query_idx <- nn.idx[[x]]
      distances <- Seurat:::fast_dist(
        x = query.embeddings[x, , drop = FALSE],
        y = embeddings,
        n = list(query_idx)
      )[[1]]
      if (!is.null(x = nearest.dist)) {
        r_dist <- distances - nearest.dist[x]
        r_dist[r_dist < 0] <- 0
        return(r_dist)
      } else {
        return(distances)
      }
    }
  )
  return(nn.dist)
}

# --- Test 1: Correctness — NNdist output matches manual reference ---
test_that("NNdist matches manual Euclidean distances", {
  set.seed(42)
  n_cells <- 20
  n_dims <- 5
  k <- 3
  embeddings <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)

  # Build nn.idx as a list: for each cell, pick k random other cells
  nn.idx <- lapply(1:n_cells, function(i) {
    sample(setdiff(1:n_cells, i), k)
  })

  result <- Seurat:::NNdist(nn.idx = nn.idx, embeddings = embeddings)

  expect_length(result, n_cells)
  for (i in 1:n_cells) {
    expect_length(result[[i]], k)
    for (j in seq_along(nn.idx[[i]])) {
      expected <- manual_euclid(embeddings[i, ], embeddings[nn.idx[[i]][j], ])
      expect_equal(result[[i]][j], expected, tolerance = 1e-10)
    }
  }
})

# --- Test 1b: with nearest.dist ---
test_that("NNdist with nearest.dist matches manual calculation", {
  set.seed(42)
  n_cells <- 20
  n_dims <- 5
  k <- 3
  embeddings <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)

  nn.idx <- lapply(1:n_cells, function(i) {
    sample(setdiff(1:n_cells, i), k)
  })
  nearest.dist <- runif(n_cells, min = 0, max = 1)

  result <- Seurat:::NNdist(
    nn.idx = nn.idx, embeddings = embeddings, nearest.dist = nearest.dist
  )

  expect_length(result, n_cells)
  for (i in 1:n_cells) {
    for (j in seq_along(nn.idx[[i]])) {
      raw_dist <- manual_euclid(embeddings[i, ], embeddings[nn.idx[[i]][j], ])
      expected <- max(raw_dist - nearest.dist[i], 0)
      expect_equal(result[[i]][j], expected, tolerance = 1e-10)
    }
  }
})

# --- Test 2: matrix nn.idx input ---
test_that("NNdist works with matrix nn.idx input", {
  set.seed(123)
  n_cells <- 15
  n_dims <- 4
  k <- 3
  embeddings <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)

  # Build nn.idx as a matrix (n_cells x k)
  nn.idx.mat <- t(sapply(1:n_cells, function(i) {
    sample(setdiff(1:n_cells, i), k)
  }))

  result <- Seurat:::NNdist(nn.idx = nn.idx.mat, embeddings = embeddings)

  expect_length(result, n_cells)
  for (i in 1:n_cells) {
    expect_length(result[[i]], k)
    for (j in 1:k) {
      expected <- manual_euclid(embeddings[i, ], embeddings[nn.idx.mat[i, j], ])
      expect_equal(result[[i]][j], expected, tolerance = 1e-10)
    }
  }
})

# --- Test 3: query != reference embeddings ---
test_that("NNdist with separate query and reference embeddings", {
  set.seed(99)
  n_query <- 10
  n_ref <- 25
  n_dims <- 6
  k <- 4
  query.embeddings <- matrix(rnorm(n_query * n_dims), nrow = n_query, ncol = n_dims)
  embeddings <- matrix(rnorm(n_ref * n_dims), nrow = n_ref, ncol = n_dims)

  nn.idx <- lapply(1:n_query, function(i) {
    sample(1:n_ref, k)
  })

  result <- Seurat:::NNdist(
    nn.idx = nn.idx, embeddings = embeddings,
    query.embeddings = query.embeddings
  )

  expect_length(result, n_query)
  for (i in 1:n_query) {
    expect_length(result[[i]], k)
    for (j in seq_along(nn.idx[[i]])) {
      expected <- manual_euclid(query.embeddings[i, ], embeddings[nn.idx[[i]][j], ])
      expect_equal(result[[i]][j], expected, tolerance = 1e-10)
    }
  }
})

# --- Test 4: Performance — batched call is faster than per-cell loop ---
test_that("NNdist is faster than naive per-cell loop", {
  skip_on_cran()
  set.seed(1)
  n_cells <- 2000
  n_dims <- 30
  k <- 199
  embeddings <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)

  nn.idx <- lapply(1:n_cells, function(i) {
    sample(setdiff(1:n_cells, i), k)
  })

  # Time the new (batched) implementation
  t_new <- system.time({
    res_new <- Seurat:::NNdist(nn.idx = nn.idx, embeddings = embeddings)
  })["elapsed"]

  # Time the old (per-cell) implementation
  t_old <- system.time({
    res_old <- NNdist_old(nn.idx = nn.idx, embeddings = embeddings)
  })["elapsed"]

  # New version should be no slower than the old per-cell loop
  # (in practice it is faster due to avoiding per-cell R-to-C++ overhead,
  #  but we use a generous threshold to avoid flaky CI failures)
  expect_lt(t_new, t_old * 1.5)
})

# --- Test 5: Byte-for-byte equivalence with old implementation ---
test_that("NNdist output is identical to old per-cell implementation", {
  skip_on_cran()
  set.seed(1)
  n_cells <- 2000
  n_dims <- 30
  k <- 199
  embeddings <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)

  nn.idx <- lapply(1:n_cells, function(i) {
    sample(setdiff(1:n_cells, i), k)
  })
  nearest.dist <- runif(n_cells, min = 0, max = 2)

  res_new <- Seurat:::NNdist(
    nn.idx = nn.idx, embeddings = embeddings, nearest.dist = nearest.dist
  )
  res_old <- NNdist_old(
    nn.idx = nn.idx, embeddings = embeddings, nearest.dist = nearest.dist
  )

  expect_identical(res_new, res_old)
})
