# setup shared test fixtures
path.to.data = file.path("../testdata/visium")

test.data.1 <- Load10X_Spatial(
  path.to.data,
  assay = "Spatial.A",
  slice = "slice1.A"
)
test.data.1 <- RenameCells(
  test.data.1,
  add.cell.id = "test.data.1"
)

test.data.2 <- Load10X_Spatial(
  path.to.data,
  assay = "Spatial.A",
  slice = "slice2.A"
)
test.data.2 <- RenameCells(
  test.data.2,
  add.cell.id = "test.data.2"
)

test.data.3 <- Load10X_Spatial(
  path.to.data,
  assay = "Spatial.B",
  slice = "slice2.B"
)
test.data.3 <- RenameCells(
  test.data.3,
  add.cell.id = "test.data.3"
)

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
