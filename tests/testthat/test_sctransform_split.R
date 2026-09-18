# Tests for SCTransform's variable features on a split assay
set.seed(42)

# setup shared fixtures
# `pbmc_small` with its assay split in two, so that each layer is modeled
# separately
split.data <- pbmc_small
suppressWarnings(
  split.data[["RNA"]] <- CreateAssay5Object(
    counts = LayerData(split.data, assay = "RNA", layer = "counts")
  )
)
split.data$stim <- rep(c("CTRL", "STIM"), length.out = ncol(split.data))
split.data[["RNA"]] <- split(split.data[["RNA"]], f = split.data$stim)


test_that("a residual variance threshold works on a split assay", {
  # `variable.features.n = NULL` selects by `variable.features.rv.th` instead of
  # by count; on a split assay this used to abort with "'nfeatures' must be a
  # single positive integer"
  object <- expect_no_error(suppressWarnings(SCTransform(
    split.data,
    vst.flavor = "v2",
    variable.features.n = NULL,
    variable.features.rv.th = 1.3,
    verbose = FALSE
  )))
  expect_gt(length(x = VariableFeatures(object)), 0)
  # the layers each selected against the threshold, so the result is the set of
  # features that clear it in at least one of them
  above.threshold <- unique(x = unlist(x = lapply(
    X = levels(x = object[["SCT"]]),
    FUN = function(model) {
      attributes <- SCTResults(
        object[["SCT"]],
        slot = "feature.attributes",
        model = model
      )
      rownames(x = attributes)[attributes$residual_variance >= 1.3]
    }
  )))
  expect_setequal(VariableFeatures(object), above.threshold)
})

test_that("features that could not be scaled are still reported as variable", {
  object <- suppressWarnings(SCTransform(split.data, verbose = FALSE))
  scaled <- rownames(x = LayerData(object[["SCT"]], layer = "scale.data"))
  unscaled <- setdiff(x = VariableFeatures(object), y = scaled)
  # sparsely expressed features are dropped from the layer's model and so have
  # no residuals; they are variable all the same and should not be hidden
  expect_gt(length(x = unscaled), 0)
  expect_true(all(unscaled %in% rownames(x = object[["SCT"]])))
})

test_that("SCTransform warns about variable features it could not scale", {
  warnings <- testthat::capture_warnings(
    object <- SCTransform(split.data, verbose = FALSE)
  )
  unscaled <- setdiff(
    x = VariableFeatures(object),
    y = rownames(x = LayerData(object[["SCT"]], layer = "scale.data"))
  )
  warning <- grep(pattern = "have no residuals", x = warnings, value = TRUE)
  expect_length(warning, 1)
  # the warning names the count, the cause, and the first of the features
  expect_match(warning, paste0("^", length(x = unscaled), " of "))
  expect_match(warning, "`min_cells` (5)", fixed = TRUE)
  expect_match(warning, head(x = unscaled, n = 1L), fixed = TRUE)
})

test_that("the warning reports the `min_cells` that was actually used", {
  warnings <- testthat::capture_warnings(
    SCTransform(split.data, min_cells = 3, verbose = FALSE)
  )
  warning <- grep(pattern = "have no residuals", x = warnings, value = TRUE)
  expect_match(warning, "`min_cells` (3)", fixed = TRUE)
})

test_that("an unsplit assay has no unscalable variable features to warn about", {
  warnings <- testthat::capture_warnings(
    SCTransform(pbmc_small, verbose = FALSE)
  )
  expect_length(grep(pattern = "have no residuals", x = warnings), 0)
})

test_that("VariableFeatures falls back to its default when none are set", {
  object <- suppressWarnings(SCTransform(split.data, verbose = FALSE))
  assay <- object[["SCT"]]
  VariableFeatures(assay) <- character(0)
  # `length()` never returns NULL, so an empty set used to request zero features
  expect_no_error(
    features <- VariableFeatures(assay, use.var.features = FALSE, nfeatures = NULL)
  )
  expect_gt(length(x = features), 0)
})
