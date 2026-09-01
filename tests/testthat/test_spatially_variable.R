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

# A merged object has one image per sample, and only the default one is used,
# so the analysis silently covered a part of the object
test_that("using one image out of several is reported", {
  first <- build_fov_object(n = 40, nfeat = 20)
  second <- build_fov_object(n = 40, nfeat = 20)
  second <- RenameCells(second, new.names = paste0("second_", Cells(second)))
  names(slot(second, name = "images")) <- "fov2"
  merged <- merge(first, second)

  expect_warning(
    suppressMessages(FindSpatiallyVariableFeatures(
      merged,
      assay = "Spatial",
      layer = "data",
      features = rownames(merged)[1:5],
      selection.method = "moransi",
      verbose = FALSE
    )),
    "covers 40 of the 80 cells"
  )
})

test_that("a single image says nothing about images", {
  object <- build_fov_object(n = 40, nfeat = 20)
  warnings <- testthat::capture_warnings(suppressMessages(FindSpatiallyVariableFeatures(
    object,
    assay = "Spatial",
    layer = "data",
    features = rownames(object)[1:5],
    selection.method = "moransi",
    verbose = FALSE
  )))
  expect_false(any(grepl("images", warnings)))
})
