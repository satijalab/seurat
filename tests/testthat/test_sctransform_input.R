set.seed(42)

# sctransform models log10 of the counts per cell and fits a density over the
# expressed features. Neither is defined for degenerate input, and both failed
# from inside sctransform with messages that name neither the cells nor the
# features involved
data("pbmc_small", package = "SeuratObject", envir = environment())

test_that("cells with no counts are named", {
  skip_if_not_installed("sctransform")
  # subsetting features leaves cells with nothing counted
  few <- pbmc_small[1:8, ]
  expect_gt(sum(Matrix::colSums(LayerData(few, layer = "counts")) == 0), 0)
  expect_error(
    suppressWarnings(SCTransform(few, verbose = FALSE)),
    "cells have no counts for the features given"
  )
  expect_error(
    suppressWarnings(SCTransform(few, verbose = FALSE)),
    "nCount_RNA > 0"
  )
})

test_that("too few detected features are reported", {
  skip_if_not_installed("sctransform")
  counts <- LayerData(pbmc_small, layer = "counts")
  counts[, ] <- 0
  counts[1:2, ] <- 3
  object <- suppressWarnings(CreateSeuratObject(counts = counts))
  expect_error(
    suppressWarnings(SCTransform(object, verbose = FALSE)),
    "features are detected in at least 5 cells"
  )
})

test_that("an ordinary object still transforms", {
  skip_if_not_installed("sctransform")
  transformed <- suppressWarnings(SCTransform(pbmc_small, variable.features.n = 20, verbose = FALSE))
  expect_true("SCT" %in% Assays(transformed))
  expect_identical(ncol(transformed), ncol(pbmc_small))
})
