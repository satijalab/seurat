set.seed(42)

# nv was computed from the object's size and then not used, so LeverageScore()
# always asked irlba for 50 components and failed on anything smaller with
# "max(nu, nv) must be strictly less than min(nrow(A), ncol(A))"
data("pbmc_small", package = "SeuratObject", envir = environment())

test_that("leverage scores can be computed on a small object", {
  scored <- LeverageScore(pbmc_small, verbose = FALSE)
  scores <- scored$leverage.score
  expect_length(scores, ncol(pbmc_small))
  expect_true(all(scores >= 0))
  expect_false(anyNA(scores))
})

test_that("sketching works on a small object", {
  sketched <- suppressWarnings(SketchData(
    pbmc_small,
    ncells = 40,
    method = "LeverageScore",
    sketched.assay = "sketch",
    verbose = FALSE
  ))
  expect_true("sketch" %in% Assays(sketched))
  expect_identical(ncol(sketched[["sketch"]]), 40L)
  expect_true(all(colnames(sketched[["sketch"]]) %in% colnames(pbmc_small)))
})

test_that("asking for more cells than there are keeps them all", {
  sketched <- suppressWarnings(SketchData(
    pbmc_small,
    ncells = ncol(pbmc_small) * 3,
    method = "LeverageScore",
    sketched.assay = "sketch",
    verbose = FALSE
  ))
  expect_identical(ncol(sketched[["sketch"]]), ncol(pbmc_small))
})

test_that("a matrix with too few dimensions is reported", {
  data <- as.sparse(LayerData(pbmc_small, layer = "data")[1, 1, drop = FALSE])
  expect_error(
    LeverageScore(data, verbose = FALSE),
    "at least two of each are needed"
  )
})
