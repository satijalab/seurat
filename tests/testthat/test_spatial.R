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
