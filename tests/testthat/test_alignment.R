# Tests for alignment related functions
set.seed(42)

# Tests for alignment (via AlignSubspace)
# --------------------------------------------------------------------------------
context("AlignSubspace")

pbmc1 <- SubsetData(pbmc_small, cells.use = pbmc_small@cell.names[1:40])
pbmc2 <- SubsetData(pbmc_small, cells.use = pbmc_small@cell.names[41:80])
pbmc1@meta.data$group <- "group1"
pbmc2@meta.data$group <- "group2"
pbmc_cca <- RunCCA(pbmc1, pbmc2)
pbmc_cca <- AlignSubspace(pbmc_cca, reduction.type = "cca", dims.align = 1:5, grouping.var = "group")

test_that("Alignment returns expected values", {
  expect_equal(dim(pbmc_cca@dr$cca.aligned@cell.embeddings), c(80, 5))
  expect_equal(pbmc_cca@dr$cca.aligned@cell.embeddings[1, 1], 0.5337046, tolerance = 1e-6)
  expect_equal(pbmc_cca@dr$cca.aligned@cell.embeddings[5, 3], 1.836985, tolerance = 1e-6)
  expect_equal(pbmc_cca@dr$cca.aligned@cell.embeddings[40, 4], 0.004078839, tolerance = 1e-6)
  expect_equal(pbmc_cca@dr$cca.aligned@cell.embeddings[80, 5], 0.2227416, tolerance = 1e-6)
})

test_that("Alignment score calculated correctly", {
  expect_equal(CalcAlignmentMetric(pbmc_cca, reduction.use = "cca.aligned", dims.use = 1:5, grouping.var = "group", nn = 5), 0.655)
})

pbmc_cca <- CalcVarExpRatio(pbmc_cca, reduction.type = "pca", grouping.var = "group", dims.use = 1:5)
test_that("CalcVarExpRatio performs as expectd", {
  expect_equal(length(pbmc_cca@meta.data$var.ratio.pca), 80)
  expect_equal(pbmc_cca@meta.data$var.ratio.pca[1], 1.818141, tolerance = 1e-6)
  expect_equal(pbmc_cca@meta.data$var.ratio.pca[40], 0.5198426, tolerance = 1e-6)
  expect_equal(pbmc_cca@meta.data$var.ratio.pca[80], 0.9824946, tolerance = 1e-6)
})

test_that("RunMultiCCA works with add.cell.ids", {
  pbmc_multi_cca <- RunMultiCCA(list(pbmc_small, pbmc_small, pbmc_small),
                                add.cell.ids = c("A", "B", "C"))
  expect_s4_class(pbmc_multi_cca, "seurat")
})
