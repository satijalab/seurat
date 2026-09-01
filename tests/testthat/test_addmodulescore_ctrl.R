set.seed(42)

# The control features are drawn from the expression bin each feature sits in,
# so ctrl larger than a bin failed with base R's "cannot take a sample larger
# than the population", which names neither ctrl nor nbin
data("pbmc_small", package = "SeuratObject", envir = environment())

test_that("a ctrl larger than the bins is explained", {
  # 230 features over 24 bins is about 9 per bin, against a default ctrl of 100
  expect_error(
    suppressWarnings(AddModuleScore(pbmc_small, features = list(rownames(pbmc_small)[1:20]))),
    "Cannot sample 100 control features"
  )
  expect_error(
    suppressWarnings(AddModuleScore(pbmc_small, features = list(rownames(pbmc_small)[1:20]))),
    "Lower ctrl to at most"
  )
})

test_that("a ctrl the bins can supply still scores", {
  scored <- suppressWarnings(AddModuleScore(
    pbmc_small,
    features = list(rownames(pbmc_small)[1:20]),
    ctrl = 5,
    name = "module"
  ))
  expect_true("module1" %in% colnames(scored[[]]))
  expect_length(scored$module1, ncol(pbmc_small))
  expect_false(anyNA(scored$module1))
})

test_that("a smaller nbin makes the bins large enough", {
  scored <- suppressWarnings(AddModuleScore(
    pbmc_small,
    features = list(rownames(pbmc_small)[1:20]),
    nbin = 5,
    ctrl = 40,
    name = "module"
  ))
  expect_true("module1" %in% colnames(scored[[]]))
})
