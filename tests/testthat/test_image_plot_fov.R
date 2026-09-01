set.seed(42)

# "No compatible spatial coordinates present" was raised whether the name given
# was absent, or named an image that is not a field of view, and it named
# neither what was asked for nor what the object has
build_object <- function(n = 20L, nfeat = 10L) {
  cells <- paste0("c", seq_len(n))
  counts <- as.sparse(matrix(rpois(nfeat * n, lambda = 3), nrow = nfeat))
  dimnames(counts) <- list(paste0("gene", seq_len(nfeat)), cells)
  object <- suppressWarnings(CreateSeuratObject(counts = counts, assay = "Spatial"))
  coords <- data.frame(x = runif(n) * 100, y = runif(n) * 100, cell = cells)
  object[["fov"]] <- CreateFOV(CreateCentroids(coords), type = "centroids", assay = "Spatial")
  object
}

test_that("an unknown field of view is named, along with what is present", {
  object <- build_object()
  expect_error(ImageDimPlot(object, fov = "Fov"), "No image named 'Fov'")
  expect_error(ImageDimPlot(object, fov = "Fov"), "Images present: 'fov' \\(FOV\\)")
  expect_error(
    ImageFeaturePlot(object, fov = "Fov", features = "gene1"),
    "No image named 'Fov'"
  )
})

test_that("an image that is not a field of view says so", {
  object <- build_object()
  object[["fov"]] <- NULL
  coordinates <- data.frame(
    x = runif(ncol(object)),
    y = runif(ncol(object)),
    row.names = colnames(object)
  )
  object[["slide"]] <- new(
    Class = "SlideSeq",
    coordinates = coordinates,
    assay = "Spatial",
    key = "slide_"
  )
  expect_error(ImageDimPlot(object, fov = "slide"), "not a field of view")
  expect_error(ImageDimPlot(object, fov = "slide"), "SpatialDimPlot")
})

test_that("an object with no images says that", {
  object <- build_object()
  object[["fov"]] <- NULL
  expect_error(ImageDimPlot(object), "no images")
})

test_that("a field of view that is present still plots", {
  object <- build_object()
  expect_no_error(ImageDimPlot(object, fov = "fov"))
})
