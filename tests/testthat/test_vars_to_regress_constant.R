set.seed(7)
data("pbmc_small")

# A variable holding one value for every cell makes the model matrix
# rank-deficient, and lm() reported "contrasts can be applied only to factors
# with 2 or more levels" without saying which variable it meant
test_that("a constant factor to regress is named", {
  object <- pbmc_small
  object$const <- factor("only")
  expect_error(
    ScaleData(object, vars.to.regress = "const", verbose = FALSE),
    "Cannot regress out .const."
  )
  expect_error(
    suppressWarnings(SCTransform(object, vars.to.regress = "const", verbose = FALSE)),
    "Cannot regress out .const."
  )
})

test_that("every constant variable is named at once", {
  object <- pbmc_small
  object$const <- factor("only")
  object$other <- "same"
  expect_error(
    ScaleData(object, vars.to.regress = c("const", "other"), verbose = FALSE),
    "they hold a single value"
  )
})

# a constant numeric contributes nothing but does not stop the fit, and is
# regularly constant only within one split of a split-by run, so it warns
test_that("a constant numeric warns and carries on", {
  object <- pbmc_small
  object$flat <- 1
  expect_warning(
    scaled <- ScaleData(object, vars.to.regress = "flat", verbose = FALSE),
    "will have no effect"
  )
  expect_identical(ncol(scaled), ncol(object))
})

test_that("variables that vary are untouched", {
  object <- pbmc_small
  object$varies <- rnorm(ncol(object))
  expect_no_error(ScaleData(object, vars.to.regress = "varies", verbose = FALSE))
  expect_no_error(
    suppressWarnings(SCTransform(object, vars.to.regress = "varies", verbose = FALSE))
  )
})
