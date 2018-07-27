# Tests for functions in snn.R
set.seed(42)

# Tests for building the SNN graph
# --------------------------------------------------------------------------------
context("BuildSNN")

snn.genes <- BuildSNN(pbmc_small, print.output = 0)
dist.mat <- as.matrix(dist(t(as.matrix(pbmc_small@data[pbmc_small@var.genes, ])), diag = TRUE, upper = TRUE))
snn.genes.dist <- BuildSNN(pbmc_small, distance.matrix = dist.mat, print.output = 0)
test_that("SNN builds from genes correctly", {
  expect_equal(diag(snn.genes@snn), rep(1, ncol(snn.genes@snn)))
  expect_equal(snn.genes@snn[2, 1], 0.25)
  expect_equal(snn.genes@snn[40, 39], 1/3)
  expect_equal(snn.genes@snn[79, 80], 0.4285714, tolerance = 1e-6)
  expect_equal(snn.genes@snn, snn.genes.dist@snn)
})

test_that("SNN checks calc parameters properly ", {
  expect_warning(BuildSNN(snn.genes, print.output = 0))
  expect_silent(BuildSNN(snn.genes, genes.use = rownames(snn.genes@data)[1:10], print.output = 0))
  expect_warning(BuildSNN(snn.genes, k.param = 100, print.output = 0))
  # should always rebuild if provided a distance matrix
  expect_silent(BuildSNN(snn.genes.dist, distance.matrix = dist.mat, print.output = 0))
})

snn.pc <- BuildSNN(pbmc_small, reduction.type = "pca", dims.use = 1:10, plot.SNN = TRUE)
test_that("SNN builds from dr (PCs) correctly", {
  expect_equal(diag(snn.pc@snn), rep(1, ncol(snn.pc@snn)))
  expect_equal(snn.pc@snn[2, 1], 1/3)
  expect_equal(snn.pc@snn[40, 39], 0.4285714, tolerance = 1e-6)
  expect_equal(snn.pc@snn[79, 80], 2/3)
})

# Tests for calculating the SNN connectivity
# --------------------------------------------------------------------------------
context("CalcConnectivity")

snn.pc <- FindClusters(snn.pc, reduction.type = "pca", dims.use = 1:10, resolution = 1, print.output = 0)
conn <- CalcConnectivity(snn.pc)
test_that("CalcConnectivity returns expected values", {
  expect_equal(dim(conn), rep(length(unique(snn.pc@ident)), 2))
  expect_equal(unname(diag(conn)), rep(0, length(unique(snn.pc@ident))))
  expect_equal(conn[1, 2], 0.4331914, tolerance = 1e-6)
  expect_equal(conn[3, 4], 0.2730696, tolerance = 1e-6)
})
