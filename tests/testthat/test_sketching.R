# Tests for sketching related fxns
set.seed(42)
pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))

# Setup test object
pbmc_small <- NormalizeData(pbmc_small, verbose = FALSE)
pbmc_small <- FindVariableFeatures(pbmc_small)

# Tests for SketchData
# ------------------------------------------------------------------------------
context("SketchData")

test_that("SketchData defaults work", {
  pbmc_sketched <- SketchData(pbmc_small, assay = "RNA", ncells = 50, method = "LeverageScore", sketched.assay = "sketch")
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
})

test_that("SketchData with named list works", {
  pbmc_sketched <- SketchData(pbmc_small, assay = "RNA", ncells = c("data" = 50), method = "LeverageScore", sketched.assay = "sketch")
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
})

pbmc_split <- split(pbmc_small, f = pbmc_small$groups)
pbmc_split <- FindVariableFeatures(pbmc_split)

test_that("SketchData with multiple layers works", { # (and one is less than the number of cells in that layer)
  pbmc_sketched <- SketchData(pbmc_split, assay = "RNA", ncells = 10, method = "LeverageScore", sketched.assay = "sketch")
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
})

test_that("SketchData with a different number of cells per layer works", {
  pbmc_sketched <- SketchData(pbmc_split, assay = "RNA", ncells = c(50, 30), method = "LeverageScore", sketched.assay = "sketch")
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  table(pbmc3k_split_sketched@meta.data[colnames(pbmc3k_split_sketched[["sketch"]]),]$random)
})

test_that("SketchData with a different number of cells per layer and a named list works", {
  pbmc_sketched <- SketchData(pbmc_split, assay = "RNA", ncells = c("data.a" = 50, "data.b" = 30), method = "LeverageScore", sketched.assay = "sketch")
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  table(pbmc3k_split_sketched@meta.data[colnames(pbmc3k_split_sketched[["sketch"]]),]$random)
})

test_that("SketchData with specified features works", {
  pbmc_sketched <- SketchData(pbmc_small, assay = "RNA", ncells = 50, method = "LeverageScore", sketched.assay = "sketch", features = VariableFeatures(pbmc_small)[1:100])
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  table(pbmc3k_split_sketched@meta.data[colnames(pbmc3k_split_sketched[["sketch"]]),]$random)
})

test_that("SketchData with specified features and multiple layers works", {
  pbmc_sketched <- SketchData(pbmc_split, assay = "RNA", ncells = c(50, 30), method = "LeverageScore", sketched.assay = "sketch", features = VariableFeatures(pbmc_small)[1:100])
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  table(pbmc3k_split_sketched@meta.data[colnames(pbmc3k_split_sketched[["sketch"]]),]$random)
})

pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))
pbmc_small <- NormalizeData(pbmc_small, verbose = FALSE)
VariableFeatures(pbmc_small) <- rownames(pbmc_small)[1:100]
test_that("SketchData when setting your own variable features and specifying features works", {
  pbmc_sketched <- SketchData(pbmc_small, assay = "RNA", ncells = 50, method = "LeverageScore", sketched.assay = "sketch", features = VariableFeatures(pbmc_small)[1:100])
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  table(pbmc3k_split_sketched@meta.data[colnames(pbmc3k_split_sketched[["sketch"]]),]$random)
})


test_that("SketchData when setting your own variable features and not specifying features errors out", {
  pbmc_sketched <- SketchData(pbmc_small, assay = "RNA", ncells = 50, method = "LeverageScore", sketched.assay = "sketch")
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  table(pbmc3k_split_sketched@meta.data[colnames(pbmc3k_split_sketched[["sketch"]]),]$random)
})

