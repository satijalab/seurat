# setup shared test fixtures
path_to_counts <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")
path_to_image = file.path("../testdata/visium/spatial")


build_spatial_data <- function(assay_name, image_name, id_prefix) {
  raw_counts <- read.table(path_to_counts, sep = "\t", row.names = 1)

  image <- Read10X_Image(
    path_to_image,
    assay = assay_name,
    slice = image_name
  )
  cell_names <- Cells(image)

  counts <- do.call(cbind, replicate(34, raw_counts, simplify = FALSE))
  counts <- counts[1:length(cell_names)]
  counts <- as.sparse(as.matrix(counts))
  colnames(counts) <- cell_names

  test_data <- CreateSeuratObject(counts, assay = assay_name)
  test_data[[image_name]] <- image
  test_data <- RenameCells(
    test_data,
    add.cell.id = id_prefix
  )

  return (test_data)
}

test_render <- function(plot) {
  grDevices::pdf(NULL)
  print(plot)
  dev.off()
}

equivalent_plots <- function(plot1, plot2) {
  if (length(plot1$layers) != length(plot2$layers)) {
    return(FALSE)
  }

  for (i in seq_along(plot1$layers)) {
    layer1 <- plot1$layers[[i]]
    layer2 <- plot2$layers[[i]]

    if (class(layer1$geom)[1] != class(layer2$geom)[1]) {
      return(FALSE)
    }

    if (!identical(layer1$data, layer2$data)) {
      return(FALSE)
    }
  }
  return(TRUE)
}

test.data.1 <- build_spatial_data(
  assay_name = "Spatial.A",
  image_name = "slice1.A",
  id_prefix = "test-data-1"
)

test.data.2 <- build_spatial_data(
  assay_name = "Spatial.A",
  image_name = "slice2.A",
  id_prefix = "test-data-2"
)

test.data.3 <- build_spatial_data(
  assay_name = "Spatial.B",
  image_name = "slice2.B",
  id_prefix = "test-data-3"
)

context("SpatialFeaturePlot")

test_that("SpatialFeaturePlot works with a single assay/image", {
  test.case <- test.data.1

  plot.1 <- SpatialFeaturePlot(test.case, features = "nCount_Spatial.A")
  plot.2 <- SpatialFeaturePlot(
    test.case,
    images = "slice1.A",
    features = "nCount_Spatial.A"
  )

  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plot.1, plot.2))
})

test_that("SpatialFeaturePlot works with multiple layers & images", {
  test.case <- merge(test.data.1, test.data.2)

  plots <- SpatialFeaturePlot(test.case, features = "nCount_Spatial.A")
  plot.1 <- SpatialFeaturePlot(
    test.case,
    images = "slice1.A",
    features = "nCount_Spatial.A"
  )
  plot.2 <- SpatialFeaturePlot(
    test.case,
    images = "slice2.A",
    features = "nCount_Spatial.A"
  )
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.1))
  expect_true(equivalent_plots(plots[[2]], plot.2))
})

test_that("SpatialFeaturePlot works with multiple overlapping images", {
  skip_on_cran()
  skip_if_not_installed("sf")

  test.case <- test.data.1
  suppressWarnings(
    test.case[["slice1.crop"]] <- Crop(
      test.case[["slice1.A"]],
      x = c(0, 5000),
      y = c(0, 5000)
    )
  )

  plots <- SpatialFeaturePlot(test.case, features = "nCount_Spatial.A")
  plot.1 <- SpatialFeaturePlot(
    test.case,
    images = "slice1.A",
    features = "nCount_Spatial.A"
  )
  plot.2 <- SpatialFeaturePlot(
    test.case,
    images = "slice1.crop",
    features = "nCount_Spatial.A"
  )
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.1))
  expect_true(equivalent_plots(plots[[2]], plot.2))
})

test_that("SpatialFeaturePlot works with multiple assays & images", {
  test.case <- merge(test.data.1, test.data.3)

  DefaultAssay(test.case) <- "Spatial.A"
  plot.1 <- SpatialFeaturePlot(test.case, features = "nCount_Spatial.A")
  plot.2 <- SpatialFeaturePlot(
    test.case,
    features = "nCount_Spatial.A",
    images = "slice2.B"
  )
  plots <- SpatialFeaturePlot(
    test.case,
    features = "nCount_Spatial.A",
    images = c("slice1.A", "slice2.B")
  )
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.1))
  expect_true(equivalent_plots(plots[[2]], plot.2))

  DefaultAssay(test.case) <- "Spatial.B"
  plot.1 <- SpatialFeaturePlot(test.case, features = "nCount_Spatial.B")
  plot.2 <- SpatialFeaturePlot(
    test.case,
    features = "nCount_Spatial.B",
    images = "slice1.A"
  )
  plots <- SpatialFeaturePlot(
    test.case,
    features = "nCount_Spatial.B",
    images = c("slice1.A", "slice2.B")
  )
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[2]], plot.1))
  expect_true(expect_true(equivalent_plots(plots[[1]], plot.2)))
})

test_that("SpatialFeaturePlot works with multiple assays, layers, & images", {
  test.case <- merge(test.data.1, c(test.data.2, test.data.3))

  DefaultAssay(test.case) <- "Spatial.A"
  plots <- SpatialFeaturePlot(test.case, features = "nCount_Spatial.A")
  plot.1 <- SpatialFeaturePlot(
    test.case,
    images = "slice1.A",
    features = "nCount_Spatial.A"
  )
  plot.2 <- SpatialFeaturePlot(
    test.case,
    images = "slice2.A",
    features = "nCount_Spatial.A"
  )
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.1))
  expect_true(equivalent_plots(plots[[2]], plot.2))

  DefaultAssay(test.case) <- "Spatial.B"
  plot.1 <- SpatialFeaturePlot(test.case, features = "nCount_Spatial.B")
  plot.2 <- SpatialFeaturePlot(
    test.case,
    features = "nCount_Spatial.B",
    images = "slice2.A"
  )
  plots <- SpatialFeaturePlot(
    test.case,
    features = "nCount_Spatial.B",
    images = c("slice2.A", "slice2.B")
  )
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.2))
  expect_true(equivalent_plots(plots[[2]], plot.1))
})


context("SpatialDimPlot")

test_that("SpatialDimPlot works with a single assay/image", {
  test.case <- test.data.1

  plot.1 <- SpatialDimPlot(test.case)
  plot.2 <- SpatialDimPlot(
    test.case,
    images = "slice1.A"
  )
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plot.1, plot.2))
})

test_that("SpatialDimPlot works with multiple layers/images", {
  test.case <- merge(test.data.1, test.data.2)

  plots <- SpatialDimPlot(test.case)
  plot.1 <- SpatialDimPlot(
    test.case,
    images = "slice1.A"
  )
  plot.2 <- SpatialDimPlot(
    test.case,
    images = "slice2.A"
  )
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.1))
  expect_true(equivalent_plots(plots[[2]], plot.2))
})

test_that("SpatialDimPlot works with multiple overlapping images", {
  skip_on_cran()
  skip_if_not_installed("sf")
  
  test.case <- test.data.1
  suppressWarnings(
    test.case[["slice1.crop"]] <- Crop(
      test.case[["slice1.A"]],
      x = c(0, 5000),
      y = c(0, 5000)
    )
  )

  plots <- SpatialDimPlot(test.case)
  plot.1 <- SpatialDimPlot(
    test.case,
    images = "slice1.A"
  )
  plot.2 <- SpatialDimPlot(
    test.case,
    images = "slice1.crop"
  )
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.1))
  expect_true(equivalent_plots(plots[[2]], plot.2))
})

test_that("SpatialDimPlot works with multiple assays/images", {
  test.case <- merge(test.data.1, test.data.3)

  DefaultAssay(test.case) <- "Spatial.A"
  plot.1 <- SpatialDimPlot(test.case)
  plot.2 <- SpatialDimPlot(
    test.case,
    images = "slice2.B"
  )
  plots <- SpatialDimPlot(
    test.case,
    images = c("slice1.A", "slice2.B")
  )

  expect_equal(length(plots), 2)
  expect_true(equivalent_plots(plots[[1]], plot.1))
  expect_true(equivalent_plots(plots[[2]], plot.2))

  DefaultAssay(test.case) <- "Spatial.B"
  plot.1 <- SpatialDimPlot(test.case)
  plot.2 <- SpatialDimPlot(
    test.case,
    images = "slice1.A"
  )
  plots <- SpatialDimPlot(test.case, images = c("slice1.A", "slice2.B"))

  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[2]], plot.1))
  expect_true(expect_true(equivalent_plots(plots[[1]], plot.2)))
})

test_that("SpatialDimPlot works with multiple assays, layers, & images", {
  test.case <- merge(test.data.1, c(test.data.2, test.data.3))

  DefaultAssay(test.case) <- "Spatial.A"
  plots <- SpatialDimPlot(test.case)
  plot.1 <- SpatialDimPlot(test.case, images = "slice1.A")
  plot.2 <- SpatialDimPlot(test.case, images = "slice2.A")
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.1))
  expect_true(equivalent_plots(plots[[2]], plot.2))

  DefaultAssay(test.case) <- "Spatial.B"
  plot.1 <- SpatialDimPlot(test.case)
  plot.2 <- SpatialDimPlot(test.case, images = "slice2.A")
  plots <- SpatialDimPlot(test.case, images = c("slice2.A", "slice2.B"))
  expect_equal(length(plots), 2)
  expect_no_error(test_render(plots))
  expect_no_error(test_render(plot.1))
  expect_no_error(test_render(plot.2))
  expect_true(equivalent_plots(plots[[1]], plot.2))
  expect_true(equivalent_plots(plots[[2]], plot.1))
})

test_that("FindSpatiallyVariableFeatures uses FOV cell names for the requested assay", {
  skip_on_cran()
  skip_if_not_installed("spatstat.explore")
  skip_if_not_installed("spatstat.geom")

  set.seed(42)
  test.case <- merge(test.data.1, test.data.3)
  segment.cells <- colnames(x = test.case[["Spatial.A"]])
  bin.cells <- colnames(x = test.case[["Spatial.B"]])
  features <- rownames(x = test.case[["Spatial.B"]])[1:6]

  test.case[["slice1.A"]] <- NULL
  test.case[["slice2.B"]] <- NULL
  test.case[["segmentation"]] <- CreateFOV(
    CreateSegmentation(data.frame(
      x = rep(c(0, 1, 1, 0), 2L) + rep(1:2, each = 4L),
      y = rep(c(0, 0, 1, 1), 2L),
      cell = rep(segment.cells[1:2], each = 4L)
    )),
    assay = "Spatial.A"
  )
  test.case[["bins"]] <- CreateFOV(
    CreateCentroids(data.frame(
      x = runif(length(x = bin.cells)) * 100,
      y = runif(length(x = bin.cells)) * 100,
      cell = bin.cells
    )),
    type = "centroids",
    assay = "Spatial.B"
  )
  DefaultAssay(test.case) <- "Spatial.A"
  test.case <- suppressWarnings(NormalizeData(test.case, assay = "Spatial.B", verbose = FALSE))

  result <- suppressWarnings(FindSpatiallyVariableFeatures(
    test.case,
    assay = "Spatial.B",
    layer = "data",
    features = features,
    selection.method = "markvariogram",
    verbose = FALSE
  ))
  expect_length(SpatiallyVariableFeatures(result[["Spatial.B"]], method = "markvariogram"), 6L)

  coords <- GetTissueCoordinates(test.case[["bins"]])[, c("x", "y")]
  rownames(coords) <- seq_len(nrow(coords))
  expect_error(
    FindSpatiallyVariableFeatures(
      test.case[["Spatial.B"]],
      layer = "data",
      features = features,
      spatial.location = coords,
      selection.method = "markvariogram",
      nfeatures = length(x = features),
      verbose = FALSE
    ),
    "row names in 'spatial.location' do not match cells"
  )

  expect_error(
    suppressWarnings(FindSpatiallyVariableFeatures(
      test.case,
      assay = "Spatial.A",
      layer = "counts",
      image = "segmentation",
      features = features,
      selection.method = "markvariogram",
      verbose = FALSE
    )),
    "centroid-based boundary"
  )
})

#' Returns a small assay plus matching coordinates for the tests below.
build_svf_case <- function(n_features_vary = 0L) {
  cells <- paste0("cell", seq_len(12L))
  counts <- matrix(3L, nrow = 4L, ncol = length(x = cells))
  dimnames(x = counts) <- list(paste0("gene", 1:4), cells)
  for (i in seq_len(n_features_vary)) {
    counts[i, ] <- sample(x = seq_along(cells))
  }
  object <- suppressWarnings(CreateSeuratObject(counts = as.sparse(counts)))
  object <- suppressWarnings(NormalizeData(object, verbose = FALSE))
  list(
    assay = object[["RNA"]],
    coordinates = data.frame(
      x = seq_len(length(x = cells)),
      y = rev(seq_len(length(x = cells))),
      row.names = cells
    )
  )
}

test_that("FindSpatiallyVariableFeatures handles bad inputs", {
  skip_on_cran()
  skip_if_not_installed("spatstat.explore")
  skip_if_not_installed("spatstat.geom")
  skip_if(!requireNamespace("Rfast2", quietly = TRUE) && !requireNamespace("ape", quietly = TRUE))

  set.seed(42)
  case <- build_svf_case(n_features_vary = 4L)
  # error when there are too few cells passed to spatial.location
  expect_error(
    FindSpatiallyVariableFeatures(
      case$assay,
      layer = "counts",
      features = rownames(x = case$assay),
      spatial.location = case$coordinates[1, , drop = FALSE],
      selection.method = "moransi",
      verbose = FALSE
    ),
    "at least two"
  )

  # test that no errors are thrown when identifying a single varying feature
  for (method in c("moransi", "markvariogram")) {
    set.seed(42)
    case <- build_svf_case(n_features_vary = 1L)
    result <- suppressWarnings(FindSpatiallyVariableFeatures(
      case$assay,
      layer = "counts",
      features = rownames(x = case$assay),
      spatial.location = case$coordinates,
      selection.method = method,
      nfeatures = 4L,
      verbose = FALSE
    ))
    expect_equal(SpatiallyVariableFeatures(result, method = method), "gene1")
  }

  # warn when no features vary
  case <- build_svf_case(n_features_vary = 0L)
  expect_warning(
    result <- FindSpatiallyVariableFeatures(
      case$assay,
      layer = "counts",
      features = rownames(x = case$assay),
      spatial.location = case$coordinates,
      selection.method = "markvariogram",
      verbose = FALSE
    ),
    "None of the requested features vary"
  )
  expect_identical(result, case$assay)
})

test_that("RunMarkVario returns one named entry per feature", {
  skip_on_cran()
  skip_if_not_installed("spatstat.explore")
  skip_if_not_installed("spatstat.geom")

  set.seed(42)
  n.cells <- 12L
  positions <- data.frame(x = seq_len(n.cells), y = rev(seq_len(n.cells)))
  # markvario() hands back a bare 'fv' for a single mark and a named list of
  # them for several; the caller relies on always getting the list.
  for (n.features in c(1L, 3L)) {
    marks <- matrix(
      data = rnorm(n = n.features * n.cells),
      nrow = n.features,
      dimnames = list(
        paste0("gene", seq_len(n.features)),
        paste0("cell", seq_len(n.cells))
      )
    )
    mv <- RunMarkVario(spatial.location = positions, data = marks)
    expect_named(mv, rownames(x = marks))
    expect_true(all(vapply(mv, inherits, logical(1L), what = "fv")))
  }
})

context("Spatial image orientation")

image_orientation_fixture <- function() {
  cells <- c("cell-a", "cell-b", "cell-c")
  counts <- as.sparse(matrix(1:6, nrow = 2, dimnames = list(c("GeneA", "GeneB"), cells)))
  list(
    matrix = list(`Gene Expression` = counts),
    centroids = data.frame(x = c(2, 12, 22), y = c(102, 202, 302), cell = cells),
    segmentations = data.frame(
      cell = rep(cells, each = 4),
      x = rep(c(0, 0, 4, 4), 3) + rep(c(0, 10, 20), each = 4),
      y = rep(c(100, 104, 104, 100), 3) + rep(c(0, 100, 200), each = 4)
    ),
    microns = data.frame(x = c(1, 11, 21), y = c(101, 201, 301), gene = "GeneA")
  )
}

image_orientation_object <- function(fixture, fov = "fov", assay = "Xenium") {
  object <- CreateSeuratObject(fixture$matrix[["Gene Expression"]], assay = assay)
  image <- CreateFOV(list(centroids = CreateCentroids(fixture$centroids),
    segmentations = CreateSegmentation(fixture$segmentations)), assay = assay,
    molecules = CreateMolecules(fixture$microns))
  object[[fov]] <- image
  object
}

test_that("Image plots orient FOVs without changing stored boundaries or molecules", {
  # check native and legacy flipped plot axes without changing stored boundaries or molecules
  fixture <- image_orientation_fixture()
  object <- NormalizeData(image_orientation_object(fixture,
    fov = "sample_native_fov"), verbose = FALSE)
  # boundary names and FOV names may both contain underscores
  object[["sample_native_fov"]][["cell_boundary"]] <- object[["sample_native_fov"]][["segmentations"]]
  for (orientation in list("vertical", NULL)) {
    attr(object@images[["sample_native_fov"]], "coords_x_orientation") <- orientation
    before <- serialize(object@images, NULL)
    for (flip in list(NULL, FALSE, TRUE)) {
      swapped <- flip %||% identical(orientation, "vertical")
      for (plotter in list(ImageDimPlot, ImageFeaturePlot)) {
        args <- list(object = object, boundaries = "cell_boundary", molecules = "GeneA",
          mols.cols = c(GeneA = "red"), flip_xy = flip, combine = FALSE)
        if (identical(plotter, ImageFeaturePlot)) args$features <- "GeneA"
        plot <- do.call(plotter, args)[[1]]
        layers <- ggplot2::ggplot_build(plot)$data
        expect_equal(sort(layers[[1]]$x), sort(if (swapped) plot$data$y else plot$data$x))
        expect_equal(sort(layers[[1]]$y), sort(if (swapped) plot$data$x else -plot$data$y))
        expect_equal(sort(layers[[2]]$x), sort(if (swapped) fixture$microns$y else fixture$microns$x))
        expect_equal(sort(layers[[2]]$y), sort(if (swapped) fixture$microns$x else -fixture$microns$y))
      }
    }
    expect_identical(serialize(object@images, NULL), before)
  }
})

test_that("split ImageFeaturePlot panels retain their orientation", {
  # check that split labels preserve the selected axis direction without changing stored data
  object <- NormalizeData(image_orientation_object(image_orientation_fixture()), verbose = FALSE)
  object$split <- c("one", "two", "one")
  before <- serialize(object@images, NULL)
  for (flip in c(FALSE, TRUE)) {
    plots <- suppressWarnings(ImageFeaturePlot(object, features = c("GeneA", "GeneB"),
      split.by = "split", flip_xy = flip, combine = FALSE))
    for (plot in plots) {
      data <- ggplot2::ggplot_build(plot)$data[[1]]
      expect_equal(data$x, if (flip) plot$data$y else plot$data$x)
      expect_equal(data$y, if (flip) plot$data$x else -plot$data$y)
    }
  }
  expect_identical(serialize(object@images, NULL), before)
})

context("Visium plotting")

visium_plot_fixture <- function(type = "VisiumV2", polygons = FALSE, rgba = FALSE) {
  cells <- letters[1:4]
  coords <- data.frame(x = c(2, 4, 6, 8), y = c(1, 3, 2, 4), cell = cells)
  pixels <- array(seq_len(8 * 10 * 3) / (8 * 10 * 3), dim = c(8, 10, 3))
  if (rgba) {
    pixels <- array(c(pixels, rep(.6, 80)), dim = c(8, 10, 4))
  }
  factors <- scalefactors(spot = 2, fiducial = 3, hires = 1, lowres = .5)
  if (type == "VisiumV1") {
    image <- new("VisiumV1", image = pixels, scale.factors = factors,
      coordinates = data.frame(imagecol = coords$x, imagerow = coords$y, row.names = cells),
      spot.radius = .1, assay = "Spatial", key = "slice_")
  } else {
    boundaries <- list(centroids = CreateCentroids(coords))
    if (polygons) {
      boundaries$segmentations <- CreateSegmentation(data.frame(
        x = rep(c(1, 1, 3, 3), 4) + rep(c(0, 2, 4, 6), each = 4),
        y = rep(c(.5, 1.5, 1.5, .5), 4) + rep(c(0, 2, 1, 3), each = 4),
        cell = rep(cells, each = 4)
      ))
    }
    image <- new("VisiumV2", image = pixels, scale.factors = factors,
      boundaries = boundaries, assay = "Spatial", key = "slice_",
      coords_x_orientation = "horizontal")
  }
  counts <- as.sparse(matrix(c(1, 0, 2, 10, 4, 3, 8, 2), nrow = 2,
    dimnames = list(c("GeneA", "GeneB"), cells)))
  object <- CreateSeuratObject(counts, assay = "Spatial")
  object$group1 <- factor(c("a", "b", "b", "a"))
  Idents(object) <- "group1"
  object[["slice"]] <- image
  object
}

test_that("Visium plotting requires horizontal coordinate orientation", {
  # (5.4.0) check that images with missing or legacy orientation metadata request an object update
  object <- visium_plot_fixture()
  for (orientation in list(character(0), "vertical")) {
    slot(object[["slice"]], "coords_x_orientation") <- orientation
    expect_error(SpatialFeaturePlot(object, features = "GeneA", slot = "counts"),
      "Please run `UpdateSeuratObject`", fixed = TRUE)
  }
  attr(object@images[["slice"]], "coords_x_orientation") <- NULL
  expect_error(SpatialFeaturePlot(object, features = "GeneA", slot = "counts"),
    "Please run `UpdateSeuratObject`", fixed = TRUE)
})

test_that("Visium radius matches the raster-derived value", {
  # compare spot radii with the original raster calculation for both image classes
  for (type in c("VisiumV1", "VisiumV2")) {
    for (rgba in c(FALSE, TRUE)) {
      image <- visium_plot_fixture(type, rgba = rgba)[["slice"]]
      raster.size <- max(dim(GetImage(image)$raster))
      for (scale in c("lowres", "hires")) {
        expect_equal(Radius(image, scale = scale), ScaleFactors(image)$spot * ScaleFactors(image)[[scale]] / raster.size)
      }
    }
  }
})


test_that("Visium panels render without changing stored images or coordinates", {
  # check that feature panels share coordinates and leave stored image data unchanged
  for (type in c("VisiumV1", "VisiumV2")) {
    object <- visium_plot_fixture(type)
    before <- serialize(object@images, NULL)
    plots <- SpatialFeaturePlot(object, features = c("GeneA", "GeneB"),
      slot = "counts", image.alpha = .4, combine = FALSE)
    expect_no_error(lapply(plots, ggplot2::ggplotGrob))
    layers <- lapply(plots, function(plot) ggplot2::ggplot_build(plot)$data[[1]])
    expect_equal(layers[[1]][, c("x", "y")], layers[[2]][, c("x", "y")])
    expect_identical(serialize(object@images, NULL), before)
  }
})

test_that("Visium polygons select boundaries explicitly and retain feature limits", {
  # check polygon selection and feature scales without changing stored image boundaries
  object <- visium_plot_fixture(polygons = TRUE)
  before <- serialize(object@images, NULL)
  for (scale in c("feature", "all")) {
    plots <- SpatialFeaturePlot(object, features = c("GeneA", "GeneB"), slot = "counts",
      plot_segmentations = TRUE, keep.scale = scale, image.alpha = .4, combine = FALSE)
    # Standard polygons include a closing vertex for each of the four cells.
    expect_equal(nrow(plots[[1]]$layers[[2]]$data), 20L)
    expect_equal(plots[[1]]$scales$get_scales("fill")$limits, c(NA, if (scale == "feature") 8 else 10))
    expect_equal(plots[[2]]$scales$get_scales("fill")$limits, c(NA, 10))
    expect_no_error(lapply(plots, ggplot2::ggplotGrob))
  }
  expect_identical(serialize(object@images, NULL), before)
})


context("Spatial file loading")

test_that("STARmap, Curio Seeker, and Slide-seq agree with base file readers", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # compare accelerated spatial loading with base readers for counts and stored coordinates
  writeLines(c("1,2", "3,4"), file.path(bundle, "cell_barcode_count.csv"))
  writeLines(c("GeneA", "GeneB"), file.path(bundle, "genes.csv"))
  writeLines(c("1\t10\t1", "1\t12\t2", "2\t20\t3", "2\t22\t4"),
    file.path(bundle, "qhulls.tsv"))
  writeLines(c("11\t1.5", "21\t3.5"), file.path(bundle, "centroids.tsv"))
  counts <- Matrix::Matrix(matrix(c(1, 2, 3, 4), 2L), sparse = TRUE)
  Matrix::writeMM(counts, file.path(bundle, "MoleculesPerMatchedBead.mtx"))
  writeLines(c("0001", "0002"), file.path(bundle, "barcodes.tsv"))
  writeLines(c("GeneA", "GeneB"), file.path(bundle, "genes.tsv"))
  coordinates <- file.path(bundle, "MatchedBeadLocation.csv")
  writeLines(c("barcode,x,y", "0001,2,12", "0002,4,24"), coordinates)
  read <- function() {
    starmap <- expect_warning(LoadSTARmap(bundle), "Coercing to dgCMatrix")
    list(starmap, LoadCurioSeeker(bundle), ReadSlideSeq(coordinates))
  }
  actual <- read()
  expected <- local({
    testthat::local_mocked_bindings(readMM_fast = Matrix::readMM,
      .ReadLoaderTable = function(file, header = TRUE, sep = ",", check.names = FALSE, ...) {
        read.table(file, header = header, sep = sep, stringsAsFactors = FALSE,
          check.names = check.names, row.names = NULL, comment.char = "", ...)
      }, .package = "Seurat")
    read()
  })
  for (i in 1:2) {
    expect_equal(LayerData(actual[[i]], layer = "counts"), LayerData(expected[[i]], layer = "counts"))
    expect_equal(actual[[i]]@images, expected[[i]]@images)
  }
  expect_equal(actual[[3]], expected[[3]])
  expect_identical(Cells(actual[[2]]), c("0001", "0002"))
  expect_equal(GetTissueCoordinates(actual[[2]][["Slice"]])$y, c(-12, -24))
})

write_xenium_fixture <- function(bundle) {
  dir.create(bundle, recursive = TRUE, showWarnings = FALSE)
  matrix.dir <- file.path(bundle, "cell_feature_matrix")
  dir.create(matrix.dir)
  counts <- Matrix::Matrix(matrix(c(1, 0, 2, 3), nrow = 2,
    dimnames = list(c("A", "B"), c("cell1", "cell2"))), sparse = TRUE)
  Matrix::writeMM(counts, file.path(matrix.dir, "matrix.mtx"))
  write.table(c("cell1", "cell2"), file.path(matrix.dir, "barcodes.tsv"),
    quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(data.frame(c("idA", "idB"), c("A", "B"), "Gene Expression"),
    file.path(matrix.dir, "features.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
  for (path in list.files(matrix.dir, full.names = TRUE)) {
    con <- gzfile(paste0(path, ".gz"), "wb")
    writeBin(readBin(path, "raw", n = file.info(path)$size), con)
    close(con)
  }
  cells <- data.frame(cell_id = c("cell1", "cell2"), x_centroid = c(1, 3),
    y_centroid = c(4, 6), segmentation_method = "cell")
  molecules <- data.frame(x_location = c(1, 2, 3, 4, 5),
    y_location = c(4, 5, 6, 7, 8), feature_name = c("A", "B", "A", "B", "A"),
    qv = c(10, 20, 30, NA, 40), unused = "unused")
  for (name in c("cells", "transcripts")) {
    con <- gzfile(file.path(bundle, paste0(name, ".csv.gz")), "wt")
    write.csv(if (name == "cells") cells else molecules, con, row.names = FALSE)
    close(con)
  }
  bundle
}

test_that("Xenium CSV molecules filter genes, QV, and flipped regions", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check transcript selection by gene, quality score, and region in either coordinate orientation
  bundle <- write_xenium_fixture(bundle)
  read <- function(...) ReadXenium(bundle, outs = "microns", type = character(), ...)$microns
  expect_equal(read(genes = "A"), data.frame(x = c(3, 5), y = c(6, 8), gene = "A"))
  expect_equal(read(genes = "B"), data.frame(x = 2, y = 5, gene = "B"))
  expect_equal(read(region = list(x = c(2, 3), y = c(5, 6))),
    data.frame(x = c(2, 3), y = c(5, 6), gene = c("B", "A")))
  expect_equal(read(flip.xy = TRUE, region = list(x = c(5, 6), y = c(2, 3))),
    data.frame(x = c(5, 6), y = c(2, 3), gene = c("B", "A")))
  expect_equal(nrow(read(genes = character())), 0L)
  expect_equal(nrow(read(genes = "unknown")), 0L)
  expect_error(read(region = list(x = c(3, 1), y = c(1, 2))), "region")
  expect_error(read(mols.qv.threshold = NA_real_), "mols.qv.threshold")
  expect_error(read(genes = NA_character_), "genes")
})

test_that("Xenium Parquet filters agree with CSV including binary feature names", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check that columnar and text transcript files produce the same filtered molecules
  skip_if_not_installed("arrow")
  bundle <- write_xenium_fixture(bundle)
  expected <- ReadXenium(bundle, outs = "microns", type = character(), genes = "A")$microns
  transcripts <- read.csv(gzfile(file.path(bundle, "transcripts.csv.gz")))
  arrow::write_parquet(transcripts, file.path(bundle, "transcripts.parquet"), chunk_size = 2L)
  expect_equal(ReadXenium(bundle, outs = "microns", type = character(), genes = "A")$microns, expected)
  transcripts$feature_name <- lapply(transcripts$feature_name, charToRaw)
  tab <- arrow::Table$create(transcripts, schema = arrow::schema(
    x_location = arrow::float64(), y_location = arrow::float64(), feature_name = arrow::binary(),
    qv = arrow::float64(), unused = arrow::utf8()))
  arrow::write_parquet(tab, file.path(bundle, "transcripts.parquet"))
  expect_equal(ReadXenium(bundle, outs = "microns", type = character(), genes = "A")$microns, expected)
})

test_that("Visium positions filter barcodes and tissue consistently across formats", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check that supported position file formats retain the same cells and coordinates
  positions <- data.frame(barcode = c("0001", "0002", "0003"), in_tissue = c(1, 0, 1),
    array_row = 1:3, array_col = 4:6, pxl_row_in_fullres = 7:9, pxl_col_in_fullres = 10:12)
  parquet <- file.path(bundle, "tissue_positions.parquet")
  csv <- file.path(bundle, "tissue_positions.csv")
  legacy <- file.path(bundle, "tissue_positions_list.csv")
  if (requireNamespace("arrow", quietly = TRUE)) arrow::write_parquet(positions, parquet)
  write.csv(positions, csv, row.names = FALSE)
  write.table(positions, legacy, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  expected <- data.frame(tissue = 1, row = 3, col = 6, imagerow = 9, imagecol = 12, row.names = "0003")
  for (path in c(csv, legacy, if (file.exists(parquet)) parquet)) {
    expect_equal(Read10X_Coordinates(path, TRUE, cells = c("0002", "0003")), expected)
    expect_equal(nrow(Read10X_Coordinates(path, FALSE)), 3L)
    expect_equal(nrow(Read10X_Coordinates(path, TRUE, cells = character())), 0L)
  }
})

write_xenium_boundaries <- function(bundle) {
  vertices <- data.frame(cell_id = rep(c("cell1", "cell2"), each = 5L),
    vertex_x = c(0, 2, 2, 0, 0, 2, 4, 4, 2, 2),
    vertex_y = c(3, 3, 5, 5, 3, 5, 5, 7, 7, 5), unused = "drop")
  for (name in c("cell", "nucleus")) {
    con <- gzfile(file.path(bundle, paste0(name, "_boundaries.csv.gz")), "wt")
    write.csv(vertices, con, row.names = FALSE)
    close(con)
  }
  vertices
}

test_that("Xenium boundary readers filter cells and preserve coordinate orientation", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check cell selection and coordinate ordering for cell and nucleus boundary files
  skip_if_not_installed("sf")
  bundle <- write_xenium_fixture(bundle)
  vertices <- write_xenium_boundaries(bundle)
  for (kind in c("segmentations", "nucleus_segmentations")) {
    selected <- ReadXenium(bundle, outs = character(), type = kind, cells = "cell2",
      flip.xy = TRUE)[[kind]]
    expect_equal(selected$cell, rep("cell2", 5L))
    expect_equal(range(selected$x), c(5, 7))
    expect_equal(range(selected$y), c(2, 4))
  }
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(vertices, file.path(bundle, "cell_boundaries.parquet"))
    selected <- ReadXenium(bundle, outs = character(), type = "segmentations", cells = "cell1")$segmentations
    expect_equal(selected$cell, rep("cell1", 5L))
    expect_named(selected, c("cell", "x", "y"))
  }
})

test_that("Xenium retries CSV after Parquet failures and caches cell metadata", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check schema and read failures across xenium tables while reading cell metadata only once
  skip_if_not_installed("arrow")
  bundle <- write_xenium_fixture(bundle)
  write_xenium_boundaries(bundle)
  read <- function() ReadXenium(bundle, outs = c("segmentation_method", "microns"),
    type = c("centroids", "segmentations", "nucleus_segmentations"))
  expected <- read()
  read.csv <- Seurat:::.ReadLoaderTable
  cell.reads <- 0L
  testthat::local_mocked_bindings(.ReadLoaderTable = function(file, ...) {
    if (basename(file) == "cells.csv.gz") cell.reads <<- cell.reads + 1L
    read.csv(file, ...)
  }, .package = "Seurat")
  testthat::local_mocked_bindings(read_parquet = function(...) stop("unsupported parquet"),
    .package = "arrow")
  for (name in c("cells", "transcripts", "cell_boundaries", "nucleus_boundaries")) {
    path <- file.path(bundle, paste0(name, ".parquet"))
    if (name %in% c("cells", "transcripts")) {
      writeLines("invalid parquet", path)
    } else {
      csv <- utils::read.csv(gzfile(file.path(bundle, paste0(name, ".csv.gz"))))
      arrow::write_parquet(csv, path)
    }
  }
  expect_equal(read(), expected)
  expect_equal(cell.reads, 1L)
})

write_visium_segmentations <- function(bundle) {
  dir.create(file.path(bundle, "segmented_outputs"), recursive = TRUE)
  polygons <- sf::st_sf(cell_id = 1:2, geometry = sf::st_sfc(lapply(0:1, function(i) {
    sf::st_polygon(list(cbind(c(0, 2, 2, 0, 0) + i * 4, c(0, 0, 2, 2, 0))))
  })))
  sf::write_sf(polygons, file.path(bundle, "segmented_outputs", "cell_segmentations.geojson"))
  bundle
}

test_that("Visium GeoJSON filters cell identifiers", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check that polygon selection uses formatted cell identifiers and handles empty selections
  skip_if_not_installed("sf")
  bundle <- write_visium_segmentations(bundle)
  ids <- Format10X_GeoJson_CellID(1:2)
  all <- Read10X_HD_GeoJson(bundle)
  expect_equal(all$barcodes, ids)
  selected <- Read10X_HD_GeoJson(bundle, cells = ids[2L])
  expect_equal(selected$barcodes, ids[2L])
  expect_equal(nrow(Read10X_HD_GeoJson(bundle, cells = character())), 0L)
  expect_equal(nrow(Read10X_HD_GeoJson(bundle, cells = "invalid")), 0L)

})


write_xenium_custom_probes <- function(bundle, format = c("csv", "parquet")) {
  format <- match.arg(format)
  bundle <- write_xenium_fixture(bundle)
  probes <- c("Dash-probe", "Under_score", "Pipe|probe", "Control_probe", "Protein_probe")
  counts <- Matrix::Matrix(matrix(seq_len(10L), 5L,
    dimnames = list(probes, c("cell1", "cell2"))), sparse = TRUE)
  types <- c(rep("Gene Expression", 3L), "Negative Control Probe", "Protein Expression")
  matrix.dir <- file.path(bundle, "cell_feature_matrix")
  Matrix::writeMM(counts, file.path(matrix.dir, "matrix.mtx"))
  write.table(data.frame(probes, probes, types), file.path(matrix.dir, "features.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  for (name in c("matrix.mtx", "features.tsv")) {
    path <- file.path(matrix.dir, name)
    con <- gzfile(paste0(path, ".gz"), "wb")
    writeBin(readBin(path, "raw", n = file.info(path)$size), con)
    close(con)
  }
  transcripts <- data.frame(x_location = seq_along(probes), y_location = seq_along(probes),
    feature_name = probes, qv = 30)
  if (format == "parquet") {
    # Exercise the native HDF5 counts reader and binary Parquet gene names.
    BPCells::write_matrix_10x_hdf5(methods::as(counts, "IterableMatrix"),
      file.path(bundle, "cell_feature_matrix.h5"), feature_types = types, type = "double")
    transcripts$feature_name <- lapply(transcripts$feature_name, charToRaw)
    arrow::write_parquet(arrow::Table$create(transcripts, schema = arrow::schema(
      x_location = arrow::int32(), y_location = arrow::int32(),
      feature_name = arrow::binary(), qv = arrow::float64())), file.path(bundle, "transcripts.parquet"))
  } else {
    con <- gzfile(file.path(bundle, "transcripts.csv.gz"), "wt")
    write.csv(transcripts, con, row.names = FALSE)
    close(con)
  }
  bundle
}

test_that("Xenium custom-probe molecule names match assays for initial reads", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check that molecule names match sanitized assay features, including singleton control assays
  for (format in c("csv", "parquet")) {
    if (format == "parquet" && !all(vapply(c("arrow", "BPCells", "hdf5r"),
      requireNamespace, logical(1), quietly = TRUE))) next
    fixture <- write_xenium_custom_probes(file.path(bundle, format), format)
    load <- function(...) suppressWarnings(suppressMessages(LoadXenium(fixture, ...)))
    expected <- c("Dash-probe", "Under-score", "Pipe-probe", "Control-probe", "Protein-probe")
    eager <- load()
    expect_s4_class(eager[["ControlProbe"]], "Assay")
    expect_s4_class(eager[["ProteinExpression"]], "Assay")
    features <- unlist(lapply(c("Xenium", "ControlProbe", "ProteinExpression"),
      function(assay) rownames(eager[[assay]])), use.names = FALSE)
    expect_identical(features, expected)
    molecules <- GetTissueCoordinates(eager, which = "molecules")
    expect_setequal(molecules$molecule, features)
    expect_identical(GetTissueCoordinates(eager, which = "molecules", features = "Under-score")$molecule, "Under-score")
    raw <- ReadXenium(fixture, outs = "microns", type = character())$microns
    expect_identical(raw$gene, c("Dash-probe", "Under_score", "Pipe|probe", "Control_probe", "Protein_probe"))

    selected <- load(genes = c("Under-score", "Pipe|probe"))
    expect_setequal(GetTissueCoordinates(selected[["fov"]], which = "molecules")$molecule,
      c("Under-score", "Pipe-probe"))
  }
})

test_that("Visium segmented loads honor assay, slice, and uppercase settings", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check custom names and uppercase features while preserving counts and segmentation storage
  skip_if_not_installed("sf")
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("BPCells")
  bundle <- write_visium_segmentations(bundle)
  spatial <- file.path(bundle, "segmented_outputs", "spatial")
  dir.create(spatial)
  png::writePNG(array(1, dim = c(10, 10, 3)), file.path(spatial, "tissue_lowres_image.png"))
  jsonlite::write_json(list(tissue_lowres_scalef = 1, tissue_hires_scalef = 1,
    fiducial_diameter_fullres = 1, spot_diameter_fullres = 1),
    file.path(spatial, "scalefactors_json.json"), auto_unbox = TRUE)
  counts <- Matrix::Matrix(matrix(1:4, 2, dimnames = list(c("geneA", "geneB"),
    Format10X_GeoJson_CellID(1:2))), sparse = TRUE)
  BPCells::write_matrix_10x_hdf5(methods::as(counts, "IterableMatrix"),
    file.path(bundle, "segmented_outputs", "filtered_feature_cell_matrix.h5"), type = "double")
  for (compact in c(FALSE, TRUE)) {
    object <- Load10X_Spatial(bundle, bin.size = "polygons", assay = "Custom", slice = "section1",
      to.upper = TRUE, compact = compact)
    image <- object[["section1.polygons"]]
    expect_identical(DefaultAssay(image), "Custom.Polygons")
    expect_identical(Key(image), Key("section1.polygons", quiet = TRUE))
    expect_identical(rownames(object), c("GENEA", "GENEB"))
    expect_equal(unname(as.matrix(LayerData(object, layer = "counts"))), unname(as.matrix(counts)))
    expect_true("segmentations" %in% Boundaries(image))
    expect_identical(length(slot(image[["segmentations"]], "polygons")) == 0L, compact)
    expect_s4_class(LayerData(object, layer = "counts"), "dgCMatrix")
  }
})


test_that("Xenium count failures report corrupt inputs and allow MEX fallback", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check that failed count reads identify the input and fall back to a valid matrix directory
  read <- function() ReadXenium(bundle, outs = "matrix", type = character())$matrix
  expect_error(read(), "missing cell_feature_matrix")
  writeLines("invalid HDF5", file.path(bundle, "cell_feature_matrix.h5"))
  dir.create(file.path(bundle, "cell_feature_matrix"))
  expect_error(read(), "cell_feature_matrix.h5:")
  expect_error(read(), "cell_feature_matrix:")
  valid <- write_xenium_fixture(file.path(bundle, "valid"))
  file.copy(list.files(file.path(valid, "cell_feature_matrix"), full.names = TRUE),
    file.path(bundle, "cell_feature_matrix"))
  expected <- Read10X(file.path(valid, "cell_feature_matrix"))
  expect_equal(read(), list(`Gene Expression` = expected))
})

test_that("LoadXenium preserves native coordinates and deprecated calls still swap them", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check native loading and plotting while retaining the deprecated coordinate swap behavior
  bundle <- write_xenium_fixture(bundle)
  expect_no_warning(native <- suppressMessages(LoadXenium(bundle)))
  lifecycle::expect_deprecated(
    flipped <- suppressMessages(LoadXenium(bundle, flip.xy = TRUE)), regexp = "flip.xy")
  centroids <- GetTissueCoordinates(native[["fov"]], which = "centroids")
  swapped <- GetTissueCoordinates(flipped[["fov"]], which = "centroids")
  expect_equal(centroids$x, c(1, 3))
  expect_equal(centroids$y, c(4, 6))
  before <- serialize(native@images, NULL)
  plot <- ImageDimPlot(native, combine = FALSE)[[1]]
  data <- ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(data$x, centroids$x)
  expect_equal(data$y, -centroids$y)
  expect_identical(serialize(native@images, NULL), before)
  expect_equal(swapped$x, centroids$y)
  expect_equal(swapped$y, centroids$x)
  molecules <- GetTissueCoordinates(native[["fov"]], which = "molecules")
  swapped <- GetTissueCoordinates(flipped[["fov"]], which = "molecules")
  expect_equal(swapped$x, molecules$y)
  expect_equal(swapped$y, molecules$x)
  plot <- ImageDimPlot(flipped, combine = FALSE)[[1]]
  data <- ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(data$x, centroids$x)
  expect_equal(data$y, centroids$y)
  create.fov <- SeuratObject::CreateFOV
  testthat::local_mocked_bindings(CreateFOV = function(...) {
    object <- create.fov(...)
    attr(object, "coords_x_orientation") <- NULL
    object
  }, .package = "Seurat")
  lifecycle::expect_deprecated(
    expect_error(suppressMessages(LoadXenium(bundle, flip.xy = TRUE)),
      "SeuratObject >= 5.3.0. Please update SeuratObject and restart R", fixed = TRUE),
    regexp = "flip.xy")
})

test_that("initial Xenium loading filters molecules by gene and region", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check that initial molecule loading applies gene and region selection in either orientation
  bundle <- write_xenium_fixture(bundle)
  for (flip in c(FALSE, TRUE)) {
    region <- if (flip) list(x = c(5, 6), y = c(2, 3)) else list(x = c(2, 3), y = c(5, 6))
    lifecycle::expect_deprecated(
      object <- suppressMessages(LoadXenium(bundle, genes = "A", region = region,
        flip.xy = flip)), regexp = "flip.xy")
    molecules <- GetTissueCoordinates(object[["fov"]], which = "molecules")
    expect_identical(molecules$molecule, "A")
    expect_equal(molecules$x, if (flip) 6 else 3)
    expect_equal(molecules$y, if (flip) 3 else 6)
    expect_length(object@misc, 0L)
  }
})

test_that("Xenium can skip transcript reads without retaining molecule state", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check that count and cell loading succeeds without reading or retaining transcripts
  bundle <- write_xenium_fixture(bundle)
  file.remove(file.path(bundle, "transcripts.csv.gz"))
  object <- suppressMessages(LoadXenium(bundle, molecule.coordinates = FALSE))
  expect_length(Molecules(object[["fov"]]), 0L)
  expect_equal(Cells(object), c("cell1", "cell2"))
  expect_length(object@misc, 0L)
})

test_that("transcript projection tolerates missing QV in CSV and Parquet", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check that transcript projection accepts files without quality scores in either format
  transcripts <- data.frame(x_location = c(1, 2), y_location = c(3, 4),
    feature_name = c("A", "B"), unused = c("large", "column"))
  con <- gzfile(file.path(bundle, "transcripts.csv.gz"), "wt")
  write.csv(transcripts, con, row.names = FALSE)
  close(con)
  expected <- data.frame(x = 1, y = 3, gene = "A")
  expect_equal(Seurat:::.XeniumReadMolecules(bundle, genes = "A"), expected)
  skip_if_not_installed("arrow")
  arrow::write_parquet(transcripts, file.path(bundle, "transcripts.parquet"))
  expect_equal(Seurat:::.XeniumReadMolecules(bundle, genes = "A"), expected)
})

context("Atera file loading")

write_atera_fixture <- function(bundle, multiplier = 1) {
  skip_if_not_installed("blosc")
  skip_if_not_installed("data.table")
  skip_if(Sys.which("zip") == "", "zip executable is required")
  dir.create(bundle, recursive = TRUE, showWarnings = FALSE)
  archive <- function(name, populate) {
    root <- tempfile(tmpdir = bundle)
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    array <- function(path, values) {
      folder <- file.path(root, path)
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      shape <- dim(values) %||% length(values)
      dtype <- if (is.character(values)) "|O" else if (is.integer(values)) "<i4" else "<f8"
      bytes <- if (is.character(values)) {
        c(writeBin(as.integer(length(values)), raw(), endian = "little"),
          do.call(c, lapply(values, function(value) {
            value <- charToRaw(enc2utf8(value))
            c(writeBin(as.integer(length(value)), raw(), endian = "little"), value)
          })))
      } else {
        if (is.matrix(values)) values <- as.vector(t(values))
        writeBin(values, raw(), endian = "little")
      }
      jsonlite::write_json(list(zarr_format = 2L, shape = as.list(shape),
        chunks = as.list(shape), dtype = dtype, order = "C", fill_value = 0,
        compressor = list(id = "blosc")), file.path(folder, ".zarray"), auto_unbox = TRUE)
      key <- paste(rep(0, length(shape)), collapse = ".")
      writeBin(blosc::blosc_compress(bytes), file.path(folder, key))
    }
    populate(array, root)
    files <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    previous <- setwd(root)
    on.exit(setwd(previous), add = TRUE, after = FALSE)
    utils::zip(file.path(bundle, name), files, flags = "-0 -q")
  }
  genes <- c("Under_score", "Pipe|probe", "Empty", "Ctrl_one", "Ctrl|two")
  archive("csc_cell_feature_matrix.zarr.zip", function(array, root) {
    array("var/feature_name", genes)
    array("var/feature_type/categories", c("Gene Expression", "Negative Control Probe"))
    array("var/feature_type/codes", c(0L, 0L, 0L, 1L, 1L))
    array("obs/cell_id", 1:2)
    array("X/indptr", seq.int(0L, 10L, 2L))
    array("X/indices", rep(0:1, 5L))
    array("X/data", as.double(1:10) * multiplier)
  })
  archive("cells.zarr.zip", function(array, root) {
    array("cell_id", 1:2)
    array("cell_summary", matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE))
    for (kind in 0:1) {
      path <- paste0("polygon_sets/", kind, "/")
      array(paste0(path, "cell_index"), 0:1)
      array(paste0(path, "num_vertices"), c(5L, 5L))
      array(paste0(path, "vertices"), rbind(c(0, 1, 0, 3, 2, 3, 2, 1, 0, 1),
        c(2, 3, 2, 5, 4, 5, 4, 3, 2, 3)))
    }
  })
  archive("transcripts.zarr.zip", function(array, root) {
    jsonlite::write_json(list(gene_names = genes), file.path(root, ".zattrs"))
    array("grid/0/gene_offset", matrix(c(0L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 4L), ncol = 2L, byrow = TRUE))
    array("grid/0/location", matrix(rep(c(1, 2), 4L), ncol = 2L, byrow = TRUE))
    array("grid/0/quality_score", rep(30, 4L))
    array("grid/0/codeword_identity", rep(1L, 4L))
  })
  bundle
}

test_that("Atera reads counts, centroids, controls, and cropped boundaries from archives", {
  skip_on_cran()
  bundle <- tempfile(pattern = "seurat-test-")
  dir.create(bundle)
  # check archive loading, native plot axes, unchanged stored data, and boundary selection
  bundle <- write_atera_fixture(bundle)
  object <- suppressWarnings(LoadAtera(bundle, bpcells.dir = FALSE))
  expect_identical(Cells(object), c("1", "2"))
  expect_identical(rownames(object), c("Under-score", "Pipe-probe", "Empty"))
  expect_equal(unname(as.matrix(LayerData(object, layer = "counts"))), matrix(1:6, nrow = 3L, byrow = TRUE))
  before <- serialize(object@images, NULL)
  coordinates <- GetTissueCoordinates(object[["fov"]], which = "centroids")
  plot <- ImageDimPlot(object, combine = FALSE)[[1]]
  data <- ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(data$x, coordinates$x)
  expect_equal(data$y, -coordinates$y)
  expect_identical(serialize(object@images, NULL), before)
  object <- suppressWarnings(subset(object, cells = "2"))
  for (kind in c("cell", "nucleus")) {
    result <- LoadAteraSegmentations(object, bundle, "fov", segmentations = kind)
    expect_setequal(GetTissueCoordinates(result[["fov"]][["segmentations"]])$cell, "2")
  }
})
