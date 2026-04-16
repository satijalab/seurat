# Tests for functions in visualization.R

set.seed(42)

# Tests for visualization utilities
# ------------------------------------------------------------------------------
pbmc_small[["tsne_new"]] <- CollapseEmbeddingOutliers(pbmc_small,
     reduction = "tsne", reduction.key = 'tsne_', outlier.sd = 0.5)

test_that("CollapseEmbeddingOutliers works", {
  expect_equal(Embeddings(pbmc_small[["tsne_new"]])[1, 1], -12.59713, tolerance = 1e-6)
  expect_equal(colSums(x = Embeddings(object = pbmc_small[["tsne_new"]])), c(-219.9218, 182.9215), check.attributes = FALSE, tolerance = 1e-5)
})


test_that("DiscretePalette works", {
  isColors <- function(x) {
    all(grepl("#[0-9A-Fa-f]{6}", x))
  }
  expect_true(isColors(DiscretePalette(26)))
  expect_true(isColors(DiscretePalette(32)))
  expect_true(isColors(DiscretePalette(36)))
  expect_warning(DiscretePalette(50), "Not enough colours")
})

test_that("ElbowPlot stdev matches Stdev()", {
  s <- Stdev(pbmc_small, "pca")
  k <- 7L
  ld <- ggplot2::layer_data(ElbowPlot(pbmc_small, ndims = k, reduction = "pca", plot_type = "stdev"), 1L)
  expect_equal(ld$y, s[seq_len(k)], tolerance = 1e-10)
  expect_equal(ld$x, seq_len(k))
})

test_that("ElbowPlot variance and cumulative match PCA variance fractions", {
  s <- Stdev(pbmc_small, "pca")
  den <- sum(s^2)
  pct <- s^2 / den * 100
  k <- 10L
  ldv <- ggplot2::layer_data(ElbowPlot(pbmc_small, ndims = k, plot_type = "variance"), 1L)
  expect_equal(ldv$y, pct[seq_len(k)], tolerance = 1e-10)
  ldc <- ggplot2::layer_data(ElbowPlot(pbmc_small, ndims = k, plot_type = "cumulative_variance"), 1L)
  expect_equal(ldc$y, cumsum(pct)[seq_len(k)], tolerance = 1e-10)
  expect_true(all(diff(ldc$y) >= -1e-12))
  n <- length(s)
  ldall <- ggplot2::layer_data(ElbowPlot(pbmc_small, ndims = n, plot_type = "cumulative_variance"), 1L)
  expect_equal(ldall$y[n], 100, tolerance = 1e-10)
})

test_that("ElbowPlot ndims validation and missing Stdev", {
  expect_error(ElbowPlot(pbmc_small, ndims = 0), "'ndims' must be a finite number")
  expect_error(ElbowPlot(pbmc_small, ndims = -1), "'ndims' must be a finite number")
  expect_error(ElbowPlot(pbmc_small, reduction = "tsne"), "No standard deviation info stored")
  expect_warning(ElbowPlot(pbmc_small, ndims = 50, reduction = "pca"), "only has information for")
})
