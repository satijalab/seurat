# GetResidual() on an SCTAssay holding several models could not handle a feature
# that only some of the models cover: the per-model call fell through to a
# variable it never created, and the results were combined with cbind, which
# needs matching row counts.
set.seed(1)

build <- function() {
  data("pbmc_small", package = "SeuratObject")
  first <- suppressWarnings(suppressMessages(
    SCTransform(pbmc_small[, 1:40], verbose = FALSE)
  ))
  second <- suppressWarnings(suppressMessages(
    SCTransform(pbmc_small[, 41:80], verbose = FALSE)
  ))
  merged <- suppressWarnings(merge(first, second))
  partial <- setdiff(VariableFeatures(first), VariableFeatures(second))[1]
  shared <- intersect(VariableFeatures(first), VariableFeatures(second))[1]
  list(object = merged, partial = partial, shared = shared)
}

test_that("the fixture really has a feature only one model covers", {
  b <- build()
  expect_length(levels(b$object[["SCT"]]), 2L)
  expect_true(!is.na(b$partial))
  expect_true(!is.na(b$shared))
})

test_that("na.rm = FALSE returns residuals instead of failing", {
  # previously: "number of rows of matrices must match", or on larger objects
  # "object 'new_residual' not found"
  b <- build()
  got <- expect_no_error(suppressWarnings(GetResidual(
    b$object, features = c(b$partial, b$shared),
    assay = "SCT", na.rm = FALSE, verbose = FALSE
  )))
  scale.data <- LayerData(got[["SCT"]], layer = "scale.data")
  expect_true(b$partial %in% rownames(scale.data))
  expect_true(b$shared %in% rownames(scale.data))
})

test_that("a partially modelled feature gets residuals where it is modelled", {
  b <- build()
  got <- suppressWarnings(GetResidual(
    b$object, features = c(b$partial, b$shared),
    assay = "SCT", na.rm = FALSE, verbose = FALSE
  ))
  values <- LayerData(got[["SCT"]], layer = "scale.data")[b$partial, ]
  models <- slot(b$object[["SCT"]], "SCTModel.list")
  covered <- Cells(models[[1]])
  uncovered <- Cells(models[[2]])
  # finite where the model covers it, NA where it does not
  expect_true(all(is.finite(values[covered])))
  expect_true(all(is.na(values[uncovered])))
})

test_that("na.rm = TRUE keeps only features every model covers", {
  b <- build()
  got <- suppressWarnings(GetResidual(
    b$object, features = c(b$partial, b$shared),
    assay = "SCT", na.rm = TRUE, verbose = FALSE
  ))
  scale.data <- LayerData(got[["SCT"]], layer = "scale.data")
  expect_false(b$partial %in% rownames(scale.data))
  expect_true(b$shared %in% rownames(scale.data))
  expect_false(anyNA(scale.data))
})

test_that("a shared feature is unaffected by the presence of a partial one", {
  b <- build()
  alone <- suppressWarnings(GetResidual(
    b$object, features = b$shared, assay = "SCT", na.rm = FALSE, verbose = FALSE
  ))
  together <- suppressWarnings(GetResidual(
    b$object, features = c(b$partial, b$shared),
    assay = "SCT", na.rm = FALSE, verbose = FALSE
  ))
  expect_equal(
    LayerData(alone[["SCT"]], layer = "scale.data")[b$shared, ],
    LayerData(together[["SCT"]], layer = "scale.data")[b$shared, ]
  )
})

# na.rm = FALSE promises NA for features a model has no fit for, but the
# per-model matrices were combined with cbind, which needs matching rows
test_that("residuals can be fetched with na.rm = FALSE", {
  skip_if_not_installed("sctransform")
  data("pbmc_small", package = "SeuratObject", envir = environment())
  object <- pbmc_small
  object$stim <- rep(c("a", "b"), length.out = ncol(object))
  object[["RNA"]] <- split(object[["RNA"]], f = object$stim)
  object <- suppressWarnings(SCTransform(object, variable.features.n = 30, verbose = FALSE))

  models <- levels(object[["SCT"]])
  modelled <- lapply(
    X = models,
    FUN = function(x) rownames(SCTResults(object[["SCT"]], slot = "feature.attributes", model = x))
  )
  # features one model has and another does not
  partial <- setdiff(modelled[[1]], modelled[[2]])
  skip_if(length(partial) < 2)
  shared <- Reduce(intersect, modelled)
  requested <- c(head(shared, 2), head(partial, 2))

  kept <- suppressWarnings(FetchResiduals(object, features = requested, na.rm = FALSE, verbose = FALSE))
  expect_setequal(rownames(kept), requested)
  # the shared features are complete, the partial ones are NA where no model fits
  expect_false(anyNA(kept[head(shared, 2), ]))
  expect_true(anyNA(kept[head(partial, 2), ]))

  dropped <- suppressWarnings(FetchResiduals(object, features = requested, na.rm = TRUE, verbose = FALSE))
  expect_setequal(rownames(dropped), head(shared, 2))
})
