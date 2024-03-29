set.seed(42)

pbmc.file <- system.file('extdata', 'pbmc_raw.txt', package = 'Seurat')
pbmc.test <- as.sparse(x = as.matrix(read.table(pbmc.file, sep = "\t", row.names = 1)))

meta.data <- data.frame(
  a = rep(as.factor(c('a', 'b', 'c')), length.out = ncol(pbmc.test)),
  row.names = colnames(pbmc.test)
)

object <- CreateSeuratObject(
  counts = pbmc.test,
  min.cells = 10,
  min.features = 30,
  meta.data = meta.data
)
object <- NormalizeData(object)
object <- SetIdent(object, value = 'a')

group.by = "a"
data <- FetchData(object = object, vars = rev(x = group.by))
data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = F]
category.matrix.avg <- CreateCategoryMatrix(labels = data, method = 'average')
category.matrix.sum <- CreateCategoryMatrix(labels = data, method = 'aggregate')


test_that("CreateCategoryMatrix works for average and aggregate", {
  expect_equal(unname(colSums(category.matrix.avg)), c(1, 1, 1))
  expect_equal(unname(colSums(category.matrix.sum)), c(27, 26, 24))
})

test_that("AverageExpression works for different layers", {
  #average expression on data layer is equal to log of average exponentiated data
  suppressWarnings(average.expression <- AverageExpression(object, layer = 'data')$RNA)
  counts.from.data.avg <- expm1(object[['RNA']]$data) %*% category.matrix.avg
  expect_equivalent(
    log1p(counts.from.data.avg),
    average.expression,
    tolerance = 1e-6
  )
  #average expression on counts layer is equal to average of counts
  suppressWarnings(average.counts <- AverageExpression(object, layer = 'counts')$RNA)
  avg.counts <- object[['RNA']]$data %*% category.matrix.avg
  expect_equivalent(
    avg.counts,
    average.counts,
    tolerance = 1e-6
  )
  #average expression on scale.data layer is equal to average of scale.data
  object <- ScaleData(object, features = rownames(object[['RNA']]$data))
  suppressWarnings(average.scale.data <- AverageExpression(object, layer = 'scale.data')$RNA)
  avg.scale <- object[['RNA']]$scale.data %*% category.matrix.avg
  expect_equivalent(
    average.scale.data,
    avg.scale,
    tolerance = 1e-6
  )
  })

test_that("AverageExpression handles features properly", {
  features <- rownames(x = object)[1:10]

  # check that the average expression is calcualted for the specifed features
  average.expression <- AverageExpression(
    object, 
    layer = "data", 
    features = features
  )$RNA
  expect_equal(rownames(x = average.expression), features)
  
  # check that an error is raised if none of the specified features are present
  expect_warning(AverageExpression(object, layer = 'data', features = "BAD"))
  # check that an error is raised if any of the features are missing
  expect_warning(
    AverageExpression(
      object, 
      layer = "data", 
      features = c(features, "BAD")
    )
  )

  # check that features can be specified as a simple vector even when
  # `layer="scale.data"` and `return.seurat=TRUE`
  object <- ScaleData(object = object, verbose = FALSE)
  avg.scale <- AverageExpression(
    object, 
    layer = "scale.data", 
    return.seurat = TRUE, 
    features = features,
    verbose = FALSE
  )$RNA
  expect_equal(rownames(avg.scale), features)
})

test_that("AverageExpression with return.seurat", {
  # counts
  avg.counts <- AverageExpression(object, layer = "counts", return.seurat = TRUE, verbose = FALSE)
  avg.counts.calc <- object[['RNA']]$counts %*% category.matrix.avg
  #test that counts are indeed equal to average counts
  expect_equivalent(
    as.matrix(avg.counts[['RNA']]$counts),
    as.matrix(avg.counts.calc),
    tolerance = 1e-6
  )
  expect_s4_class(object = avg.counts, "Seurat")
  avg.counts.mat <- AverageExpression(object, layer = 'counts')$RNA
  expect_equal(unname(as.matrix(LayerData(avg.counts[["RNA"]], layer = "counts"))),
               unname(as.matrix(avg.counts.mat)))
  avg.data <- LayerData(avg.counts[["RNA"]], layer = "data")
  #test that data returned is log1p of average counts
  expect_equivalent(
    as.matrix(log1p(avg.counts.mat)),
    as.matrix(avg.data),
    tolerance = 1e-6
  )
  #test that scale.data returned is scaled data
  avg.scale <- LayerData(avg.counts[["RNA"]], layer = "scale.data")
  expect_equal(
    avg.scale,
    ScaleData(avg.counts)[['RNA']]$scale.data,
    tolerance = 1e-6
  )
  # data
  avg.data <- AverageExpression(object, layer = "data", return.seurat = TRUE, verbose = FALSE)
  expect_s4_class(object = avg.data, "Seurat")
  avg.data.mat <- AverageExpression(object, layer = 'data')$RNA
  expect_equal(unname(as.matrix(LayerData(avg.data[["RNA"]], layer = "counts"))),
               unname(as.matrix(avg.data.mat)))
  expect_equal(unname(as.matrix(LayerData(avg.data[["RNA"]], layer = "data"))),
               as.matrix(unname(log1p(x = avg.data.mat))))
  avg.scale <- LayerData(avg.data[["RNA"]], layer = "scale.data")
  expect_equal(
    avg.scale['MS4A1', ],
    c(a = -0.07823997, b = 1.0368218, c = -0.9585818),
    tolerance = 1e-6
  )
  expect_equal(
    avg.scale['SPON2', ],
    c(a = 0.1213127, b = 0.9338096, c = -1.0551222),
    tolerance = 1e-6
  )
  # scale.data
  object <- ScaleData(object = object, verbose = FALSE)
  avg.scale <- AverageExpression(object, layer = "scale.data", return.seurat = TRUE, verbose = FALSE)
  expect_s4_class(object = avg.scale, "Seurat")
  avg.scale.mat <- AverageExpression(object, layer = 'scale.data')$RNA
  expect_equal(unname(as.matrix(LayerData(avg.scale[["RNA"]], layer = "scale.data"))), unname(as.matrix(avg.scale.mat)))
})

test.dat <- LayerData(object = object, layer = "data")
rownames(x = test.dat) <- paste0("test-", rownames(x = test.dat))
object[["TEST"]] <- CreateAssayObject(data = test.dat)

test_that("AverageExpression with multiple assays", {
  avg.test <- AverageExpression(object = object, assays = "TEST", layer = "data")
  expect_equal(names(x = avg.test), "TEST")
  expect_equal(length(x = avg.test), 1)
  expect_equivalent(
    avg.test[[1]]['test-KHDRBS1', 1:3],
    c(a = 10.329153, b = 92.287109, c = 5.620942),
    tolerance = 1e-6
  )
  expect_equivalent(
    avg.test[[1]]['test-DNAJB1', 1:3] ,
    c(a = 42.32240, b = 15.94807, c = 15.96319),
    tolerance = 1e-6
  )
  avg.all <- AverageExpression(object = object, layer = "data")
  expect_equal(names(x = avg.all), c("RNA", "TEST"))
  expect_equal(length(x = avg.all), 2)
})


meta.data.2 <- data.frame(
  b = rep(as.factor(c('c', 'd', 'e')), length.out = ncol(pbmc.test)),
  row.names = colnames(pbmc.test)
)
object <- AddMetaData(object, meta.data.2)
if(class(object[['RNA']]) == "Assay5")  {
  test_that("AggregateExpression works with multiple layers", {
    object.split <- split(object, f = object$b)
    aggregate.split <- AggregateExpression(object.split, assay = "RNA")
    aggregate <- AggregateExpression(object, assay = "RNA")
    expect_equivalent(
      aggregate.split$RNA,
      aggregate$RNA,
      tolerance = 1e-6
    )
    avg.split <- AverageExpression(object.split, assay = "RNA")
    avg <- AverageExpression(object, assay = "RNA")
    expect_equivalent(
      avg.split$RNA,
      avg$RNA,
      tolerance = 1e-6
    )
  })
}

test_that("PercentAbove works as expected", {
  vals <- c(1, 1, 2, 2, NA)
  expect_equal(PercentAbove(vals, threshold = 1), 0.4)
})


context("BuildNicheAssay")

test_that("BuildNicheAssay works as expected", {
  test.data <- pbmc_small

  # generate fake coordinates arranging the cells from pbmc_small into a grid
  test.coordinates <- data.frame(
    cell = Cells(test.data), 
    x = rep(1:4, times = 20), 
    y = rep(1:4, each = 20)
  )
  # associate the coordinates and counts with a FOV
  fov <- CreateFOV(
    test.coordinates,
    type = "centroids",
    assay = "RNA"
  )
  test.data[["fov"]] <- fov

  # dividing the grid into 4 along each axis creates 16 regions - label
  # each cell with the region containing it's x, y position
  x.regions <- cut(test.coordinates[["x"]], breaks = 4, labels = FALSE)
  y.regions <- cut(test.coordinates[["y"]], breaks = 4, labels = FALSE)
  test.data[["test_labels"]] <- ((y.regions - 1) * 4 + x.regions)

  results <- BuildNicheAssay(
    test.data,
    fov = "fov",
    group.by = "test_labels",
    assay = "niche"
  )
  # the new niche assay should contain the same number of cells as the input
  # and a feature for each unique label from the specified `group.by` variable
  expect_equal(
    c(16, ncol(test.data[["RNA"]])),
    dim(results[["niche"]])
  )
  # exaclty a quarter of the cells should be assigned to each niche 
  for (niche in 1:4) {
    expect_equal(
      20,
      length(results[["niches"]][results[["niches"]]["niches"] == niche])
    )
  }
})
