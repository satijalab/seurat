# Tests for conversions
set.seed(42)

context("Conversions")

test_that("Conversion to anndata", {
  expect_silent(Convert(pbmc_small, "anndata"))
})

test_that("Conversion to anndata, no dr", {
  nodr <- pbmc_small
  nodr@dr <- list()
  expect_silent(Convert(nodr, "anndata"))
})