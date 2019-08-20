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