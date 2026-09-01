set.seed(42)

# SCTransform on a split assay picks its variable features from the union of
# what the layers chose, but only features present in every layer are scaled.
# Features variable in one layer and absent from another were returned as
# variable while having no residuals, so PrepDR dropped them from PCA.
split_object <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  object <- pbmc_small
  object$stim <- rep(c("CTRL", "STIM"), length.out = ncol(object))
  object[["RNA"]] <- split(object[["RNA"]], f = object$stim)
  object
}

test_that("every variable feature has been scaled", {
  skip_if_not_installed("sctransform")
  object <- suppressWarnings(SCTransform(split_object(), verbose = FALSE))
  scaled <- rownames(LayerData(object[["SCT"]], layer = "scale.data"))
  expect_gt(length(VariableFeatures(object)), 0)
  expect_true(all(VariableFeatures(object) %in% scaled))
  # and PCA runs on all of them, rather than warning that some were not scaled
  warnings <- testthat::capture_warnings(RunPCA(object, npcs = 5, verbose = FALSE))
  expect_false(any(grepl("have not been scaled", warnings)))
})

test_that("a feature missing from one layer is not returned as variable", {
  skip_if_not_installed("sctransform")
  object <- split_object()
  counts <- LayerData(object, layer = "counts.CTRL")
  # silence genes in one batch only, so that layer's model has nothing for them
  silenced <- rownames(counts)[1:20]
  counts[silenced, ] <- 0
  LayerData(object, layer = "counts.CTRL") <- counts

  object <- suppressWarnings(SCTransform(object, verbose = FALSE))
  scaled <- rownames(LayerData(object[["SCT"]], layer = "scale.data"))
  expect_true(all(VariableFeatures(object) %in% scaled))
})

# variable.features.n = NULL selects by residual variance threshold instead of
# by count, and asking for zero features errored
test_that("SCTransform accepts a residual variance threshold on split data", {
  skip_if_not_installed("sctransform")
  object <- expect_no_error(suppressWarnings(SCTransform(
    split_object(),
    vst.flavor = "v2",
    variable.features.n = NULL,
    variable.features.rv.th = 1.3,
    verbose = FALSE
  )))
  scaled <- rownames(LayerData(object[["SCT"]], layer = "scale.data"))
  expect_gt(length(VariableFeatures(object)), 0)
  expect_true(all(VariableFeatures(object) %in% scaled))
})
