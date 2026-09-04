set.seed(42)

# A VisiumV2 image holds only the pixel position of each spot, so the array row
# and column, and whether the spot is over tissue, were dropped on load. A
# VisiumV1 image kept all five columns in its `coordinates` slot
image.dir <- file.path("../testdata/visium/spatial")

build_object <- function() {
  image <- Read10X_Image(image.dir = image.dir, image.name = "tissue_lowres_image.png")
  cells <- Cells(image)
  counts <- as.sparse(matrix(
    rpois(30 * length(cells), lambda = 3),
    nrow = 30,
    dimnames = list(paste0("gene", seq_len(30)), cells)
  ))
  object <- suppressWarnings(CreateSeuratObject(counts = counts, assay = "Spatial"))
  DefaultAssay(image) <- "Spatial"
  object[["slice1"]] <- image
  object
}

test_that("the array coordinates are kept", {
  object <- build_object()
  added <- Seurat:::.AddArrayCoordinates(object = object, image.dir = image.dir)
  expect_true(all(c("array_row", "array_col", "in_tissue") %in% colnames(added[[]])))

  positions <- Read10X_Coordinates(
    filename = Sys.glob(paths = file.path(image.dir, "*tissue_positions*"))[1],
    filter.matrix = TRUE
  )
  expect_equal(unname(added$array_row), unname(positions[colnames(added), "row"]))
  expect_equal(unname(added$array_col), unname(positions[colnames(added), "col"]))
  expect_true(all(added$in_tissue == 1))
  # and the pixel positions the image holds are unchanged
  expect_identical(
    GetTissueCoordinates(added[["slice1"]]),
    GetTissueCoordinates(object[["slice1"]])
  )
})

test_that("a directory with no tissue positions leaves the object alone", {
  object <- build_object()
  empty <- file.path(tempdir(), "no-positions")
  dir.create(path = empty, showWarnings = FALSE)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  expect_identical(
    colnames(Seurat:::.AddArrayCoordinates(object = object, image.dir = empty)[[]]),
    colnames(object[[]])
  )
})

test_that("Load10X_Spatial keeps the array coordinates", {
  skip_if_not_installed("hdf5r")
  object <- suppressWarnings(Load10X_Spatial(
    data.dir = file.path("../testdata/visium"),
    filename = "filtered_feature_bc_matrix.h5"
  ))
  expect_true(all(c("array_row", "array_col", "in_tissue") %in% colnames(object[[]])))
  expect_false(anyNA(object$array_row))
})
