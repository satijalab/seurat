context("Read10X")
# These tests were added to ensure Seurat was forwards and backwards compatible for 3.0 data

dname = "../testdata/cr3.0"
test.data <- Read10X(dname)
test.data2 <- Read10X(c(dname, dname))

test_that("Cell Ranger 3.0 Data Parsing", {
  expect_is(test.data, "list")
  expect_equal(ncol(test.data$`Gene Expression`), .5 * ncol(test.data2$`Gene Expression`))
  expect_equal(ncol(test.data$`Antibody Capture`), .5 * ncol(test.data2$`Antibody Capture`))
  expect_equal(colnames(test.data2[[1]])[6], "2_AAAGTAGCACAGTCGC-1")
  expect_equal(test.data$`Gene Expression`[2,2], 1000)
})

# Tests of Pre-3.0 Data
test.data3 <- Read10X("../testdata/")
test_that("Read10X creates sparse matrix", {
  expect_is(test.data3, "dgCMatrix")
  expect_equal(colnames(test.data3)[1], "ATGCCAGAACGACT-1")
  expect_equal(rownames(test.data3)[1], "MS4A1")
})

test_that("Read10X handles missing files properly", {
  expect_error(Read10X("."))
  expect_error(Read10X("./notadir/"))
  expect_error(Read10X(dname, gene.column = 10))
})


context("Load10X_Spatial")

# setup test fixtures
path.to.data <- file.path("../testdata")
path.to.visium <- file.path(path.to.data, "visium")
path.to.visium.hd <- file.path(path.to.data, "visium_hd")

test_that("Read10X_h5 works as expected", {
  skip_if_not_installed("hdf5r")

  path.to.counts <- file.path(
    path.to.visium,
    "filtered_feature_bc_matrix.h5"
  )
  counts <- Read10X_h5(path.to.counts)

  # check that the shape of the returned matrix is correct
  expect_equal(ncol(counts), 2695)
  expect_equal(nrow(counts), 100)
  # spot-check a few of the values
  expect_equal(colnames(counts)[[151]], "AATGCAACCGGGTACC-1")
  expect_equal(counts[1, 151], 1)
  expect_equal(colnames(counts)[[9]], "AAACCGTTCGTCCAGG-1")
  expect_equal(counts[5, 9], 1)
  expect_equal(colnames(counts)[[2328]], "TCTCGAGGAGGTTCGC-1")
  expect_equal(counts[99, 2328], 1)
})

test_that("Read10X_Image works as expected", {
  path.to.images <- file.path(path.to.visium, "spatial")
  coordinate.filenames <-  "tissue_positions_list.csv"
  # only test HD the dataset containing a parquet file if `arrow` is installed
  if (requireNamespace("arrow", quietly = TRUE)) {
    path.to.images <- c(
      path.to.images,
      file.path(
        path.to.visium.hd,
        "binned_outputs/square_008um/spatial"
      )
    )
    coordinate.filenames <- c(
      coordinate.filenames,
      "tissue_positions.parquet"
    )
  }

  for (i in seq_along(path.to.images)) {
    path.to.image <- path.to.images[[i]]
    coordinate.filename <- coordinate.filenames[[i]]

    # read in the coordinates as a data.frame - only keep the relevant columns
    # and rename them to match the processed output
    coordinates.expected <- Read10X_Coordinates(
      file.path(path.to.image, coordinate.filename),
      filter.matrix = TRUE
    )[, c("imagerow", "imagecol")]
    colnames(coordinates.expected) <- c("x", "y")
    # read in the scale factors as an S3 object
    scale.factors.expected <- Read10X_ScaleFactors(
      file.path(path.to.image, "scalefactors_json.json")
    )
    # default/lowres scaling
    image.lowres <- Read10X_Image(
      path.to.image, 
      image.name = "tissue_lowres_image.png"
    )
    coordinates <- GetTissueCoordinates(image.lowres, scale = "lowres")
    spot.radius <- Radius(image.lowres, scale = "lowres")
    scale.factors <- ScaleFactors(image.lowres)
    # check that the scale factors were read in as expected
    expect_true(identical(scale.factors, scale.factors.expected))
    # check that `coordinates` contains values scaled for the low resolution PNG
    expect_equal(
      coordinates[, c("x", "y")] / scale.factors[["lowres"]], 
      coordinates.expected
    )
    # check that the spot size is similarly scaled
    expect_equal(
      (spot.radius / scale.factors[["lowres"]] * max(dim(image.lowres))),
      scale.factors.expected[["spot"]],
    )

    # hires scaling
    image.hires <- Read10X_Image(
      path.to.image, 
      image.name = "tissue_hires_image.png"
    )
    coordinates <- GetTissueCoordinates(image.hires, scale = "hires")
    spot.radius <- Radius(image.hires, scale = "hires")
    scale.factors <- ScaleFactors(image.hires)
    # check that the scale factors were read in as expected
    expect_true(identical(scale.factors, scale.factors.expected))
    # check that `coordinates` contains values scaled for the high resolution PNG
    expect_equal(
      coordinates[, c("x", "y")] / scale.factors[["hires"]], 
      coordinates.expected
    )
    # check that the spot size is similarly scaled
    expect_equal(
      (spot.radius / scale.factors[["hires"]] * max(dim(image.hires))),
      scale.factors.expected[["spot"]]
    )
    # the size of the two images should be different
    expect_false(all(dim(image.hires) == dim(image.lowres)))

    # `VisiumV1` image
    image.v1 <- Read10X_Image(
      path.to.image,
      image.name = "tissue_lowres_image.png",
      image.type = "VisiumV1"
    )
    coordinates <- GetTissueCoordinates(image.v1, scale = "lowres")
    spot.radius <- Radius(image.v1, scale = "lowres")
    scale.factors <- ScaleFactors(image.v1)
    # check that the scale factors were read in as expected
    expect_true(identical(scale.factors, scale.factors.expected))
    # check that `coordinates` contains values scaled for the low resolution PNG
    # also make sure that it has the expected column names
    coordinates.expected.v1 <- coordinates.expected
    colnames(coordinates.expected.v1) <- c("imagerow", "imagecol")
    expect_equal(
      coordinates[, c("imagerow", "imagecol")] / scale.factors[["lowres"]],
      coordinates.expected.v1
    )
    # check that the spot size is similarly scaled
    expect_equal(
      (spot.radius / scale.factors[["lowres"]] * max(dim(image.lowres))),
      scale.factors.expected[["spot"]],
    )
  }
})

test_that("Load10X_Spatial works with SD data", {
  skip_if_not_installed("hdf5r")

  path.to.counts <- file.path(
    path.to.visium,
    "filtered_feature_bc_matrix.h5"
  )
  path.to.image <- file.path(path.to.visium, "spatial")

  # load the expected counts matrix
  counts.expected <- Read10X_h5(path.to.counts)

  # load the expected image
  image.expected <- Read10X_Image(path.to.image)
  # set the image's key
  Key(image.expected) <- "slice1_"
  # update the expected image's assay to match the default
  DefaultAssay(image.expected) <- "Spatial"
  # align the expected image's identifiers with the expected count matrix
  image.expected <- image.expected[colnames(counts.expected)]

  spatial <- Load10X_Spatial(path.to.visium)

  # check that `spatial` contains the expected counts matrix
  expect_true(
    identical(
      LayerData(spatial, assay = "Spatial", layer = "counts"),
      counts.expected
    )
  )
  # check that `spatial` contains the expected image
  expect_true(
    identical(
      spatial[["slice1"]],
      image.expected
    )
  )
})

test_that("Load10X_Spatial works with HD data", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("arrow")

  # since the "HD" test data is actually just the standard definition test
  # data re-structured to look like it's been binned at multiple resolutions
  # the `merge` call inside `Load10X_Spatial` throws a warning
  spatial <- suppressWarnings(Load10X_Spatial(path.to.visium.hd))

  bin.size <- c(16, 8)
  for (i in seq_along(bin.size)) {
    bin.size.pretty <- paste0(sprintf("%03d", bin.size[[i]]), "um")
    assay.name <- paste0("Spatial.", bin.size.pretty)
    image.key <- paste0("slice1", bin.size.pretty, "_")
    image.key.pretty <- paste0("slice1.", bin.size.pretty)

    path.to.bin <- file.path(
      path.to.visium.hd,
      paste0("binned_outputs/square_", bin.size.pretty)
    )
    path.to.counts <- file.path(
      path.to.bin,
      "filtered_feature_bc_matrix.h5"
    )
    path.to.image <- file.path(path.to.bin, "spatial")

    # load the expected counts matrix
    counts.expected <- Read10X_h5(path.to.counts)
    # accomodate for the way cell identifiers will be reset during
    # the call to `merge` inside `Load10X_Spatial`
    colnames(counts.expected) <- paste0(colnames(counts.expected), "_", i)
    # load the expected image
    image.expected <- Read10X_Image(path.to.image, assay = assay.name)
    # set the image's key
    Key(image.expected) <- image.key
    # again, accomodate for the way cell identifiers will be reset during
    # the call to `merge` inside `Load10X_Spatial`
    image.expected <- RenameCells(image.expected, paste0(Cells(image.expected), "_", i))
    # align the expected image's identifiers with the expected count matrix
    image.expected <- image.expected[colnames(counts.expected)]

    # check that `spatial` contains the expected counts matrix
    expect_true(
      identical(
        LayerData(spatial, assay = assay.name, layer = "counts"),
        counts.expected
      )
    )
    # check that `spatial` contains the expected image
    expect_true(
      identical(
        spatial[[image.key.pretty]],
        image.expected
      )
    )
  }
})

test_that("Read10X_Spatial handles missing files properly", {
  expect_error(Load10X_Spatial(data.dir = "."))
  expect_error(Load10X_Spatial(data.dir = "./notadir/"))
})
