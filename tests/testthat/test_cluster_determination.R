# Tests for functions in cluster_determination.R
set.seed(42)

# delete all clustering related slots
pbmc_small@snn <- sparseMatrix(1, 1, x = 1)
pbmc_small@meta.data$res.0.8 <- NULL
pbmc_small@meta.data$res.1 <- NULL
pbmc_small@calc.params$FindClusters.res.0.8 <- NULL
pbmc_small@calc.params$FindClusters.res.1 <- NULL
dist.mat <- as.matrix(dist(t(as.matrix(pbmc_small@data[pbmc_small@var.genes, ])), diag = TRUE, upper = TRUE))
pbmc1 <- FindClusters(object = pbmc_small, reduction.type = "pca", dims.use = 1:4, print.output = 0)
pbmc2 <- FindClusters(object = pbmc_small, reduction.type = "pca", dims.use = 1:4, print.output = 0, save.SNN = TRUE)
pbmc3 <- FindClusters(object = pbmc_small, distance.matrix = dist.mat, print.output = 0, save.SNN = TRUE)


test_that("FindClusters checks calc parameters properly ", {
  expect_warning(FindClusters(object = pbmc1, reduction.type = "pca", dims.use = 1:4, print.output = 0))
  expect_warning(FindClusters(object = pbmc1, reduction.type = "pca", dims.use = 1:4, print.output = 1))
  expect_silent(FindClusters(object = pbmc1, reduction.type = "pca", dims.use = 1:4, print.output = 0, force.recalc = T))
  expect_warning(FindClusters(object = pbmc2, reduction.type = "pca", dims.use = 1:4, print.output = 0))
  expect_silent(FindClusters(object = pbmc2, reduction.type = "pca", dims.use = 1:4, print.output = 0, force.recalc = T))
  expect_silent(FindClusters(object = pbmc2, reduction.type = "pca", dims.use = 1:5, print.output = 0))
  expect_warning(FindClusters(object = pbmc2, reduction.type = "pca", dims.use = 1:4, resolution = 1, print.output = 0))
  expect_silent(FindClusters(object = pbmc3, distance.matrix = dist.mat, print.output = 0, save.SNN = TRUE))
})

pbmc4 <- FindClusters(pbmc_small, reduction.type = "pca", dims.use = 1:4, resolution = 1.5, print.output = 0)

test_that("Singleton grouping done correctly", {
  expect_equal(length(levels(pbmc4@ident)), 12)
  expect_equal(length(which(table(pbmc4@ident) == 1)), 0)
  expect_equal(length(WhichCells(object = pbmc4, ident = 11)), 2)
  expect_equal(length(WhichCells(object = pbmc4, ident = 0)), 18)
})

