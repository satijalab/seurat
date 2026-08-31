set.seed(42)

# FOV-based images (Xenium, Vizgen, VisiumV2) return cell identifiers in a
# 'cell' column with integer row names. FindSpatiallyVariableFeatures used the
# row names to select cells, so it selected none, and RowVar() then aborted the
# session on the resulting empty matrix.
build_fov_object <- function(n = 60, nfeat = 30) {
  cells <- paste0("c", seq_len(n))
  coords <- data.frame(x = runif(n) * 100, y = runif(n) * 100, cell = cells)
  fov <- CreateFOV(CreateCentroids(coords), type = "centroids", assay = "Spatial")
  counts <- as.sparse(matrix(rpois(nfeat * n, lambda = 3), nrow = nfeat))
  dimnames(counts) <- list(paste0("gene", seq_len(nfeat)), cells)
  object <- suppressWarnings(CreateSeuratObject(counts = counts, assay = "Spatial"))
  object[["fov"]] <- fov
  suppressWarnings(NormalizeData(object, verbose = FALSE))
}

test_that("tissue coordinates keep cell names in a column, not row names", {
  object <- build_fov_object()
  coords <- GetTissueCoordinates(object[["fov"]])
  expect_true("cell" %in% colnames(coords))
  expect_identical(rownames(coords), as.character(seq_len(nrow(coords))))
})

test_that("FindSpatiallyVariableFeatures works on an FOV-based object", {
  object <- build_fov_object()
  result <- suppressWarnings(FindSpatiallyVariableFeatures(
    object,
    assay = "Spatial",
    layer = "data",
    features = rownames(object)[1:8],
    selection.method = "markvariogram",
    verbose = FALSE
  ))
  ranked <- SpatiallyVariableFeatures(result, method = "markvariogram")
  expect_length(ranked, 8L)
  expect_true(all(ranked %in% rownames(object)))
})

test_that("coordinates that name no cells give an error, not a crash", {
  object <- build_fov_object()
  coords <- GetTissueCoordinates(object[["fov"]])[, c("x", "y")]
  rownames(coords) <- seq_len(nrow(coords))
  # RowVar() is C++ and aborts on an empty matrix rather than raising
  expect_error(
    suppressWarnings(FindSpatiallyVariableFeatures(
      object[["Spatial"]],
      layer = "data",
      features = rownames(object)[1:5],
      spatial.location = coords,
      selection.method = "markvariogram",
      nfeatures = 5,
      verbose = FALSE
    )),
    "row names must be cell names"
  )
})

test_that("the assay method does not need nfeatures to be supplied", {
  object <- build_fov_object()
  coords <- GetTissueCoordinates(object[["fov"]])
  rownames(coords) <- coords$cell
  coords <- coords[, c("x", "y")]
  # nfeatures defaulted to itself, so omitting it raised
  # "promise already under evaluation"
  expect_no_error(suppressWarnings(FindSpatiallyVariableFeatures(
    object[["Spatial"]],
    layer = "data",
    features = rownames(object)[1:5],
    spatial.location = coords,
    selection.method = "markvariogram",
    verbose = FALSE
  )))
})

test_that("a single feature with variance stays a matrix", {
  object <- build_fov_object()
  data <- LayerData(object, layer = "data")
  # flatten every feature but one so only that one survives the variance filter
  data[2:nrow(data), ] <- 0
  LayerData(object, layer = "data") <- data
  expect_no_error(suppressWarnings(FindSpatiallyVariableFeatures(
    object,
    assay = "Spatial",
    layer = "data",
    features = rownames(object)[1:5],
    selection.method = "markvariogram",
    verbose = FALSE
  )))
})

# selection.method defaults to the whole vector of choices and was never
# resolved, so leaving it out reached switch() with a length-2 value
test_that("selection.method does not have to be supplied", {
  object <- build_fov_object()
  result <- suppressWarnings(FindSpatiallyVariableFeatures(
    object,
    assay = "Spatial",
    layer = "data",
    features = rownames(object)[1:5],
    verbose = FALSE
  ))
  expect_length(SpatiallyVariableFeatures(result, method = "markvariogram"), 5L)
})

test_that("an unknown selection.method is rejected before any work", {
  object <- build_fov_object()
  expect_error(
    suppressWarnings(FindSpatiallyVariableFeatures(
      object,
      assay = "Spatial",
      layer = "data",
      features = rownames(object)[1:5],
      selection.method = "bogus",
      verbose = FALSE
    )),
    "should be one of"
  )
})
