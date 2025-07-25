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

test_that("AverageExpression with return.seurat behaves as expected", {
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

test_that("Aggregate Expression correctly sums counts", {
  test.obj <- CreateSeuratObject(
    counts = pbmc.test,
    min.cells = 10,
    min.features = 30,
    meta.data = meta.data
  )
  # test normalization method
  aggregate.obj <- AggregateExpression(test.obj, assay = "RNA", group.by = "a", return.seurat = TRUE,
                                       normalization.method = "CLR")
  aggregate.counts.calc <- object[['RNA']]$counts %*% category.matrix.sum
  #test that counts are indeed equal to average counts
  expect_equivalent(
    as.matrix(aggregate.obj[['RNA']]$counts),
    as.matrix(aggregate.counts.calc),
    tolerance = 1e-6
  )
})

test_that("Aggregate Expression return.seurat=TRUE returns normalized data", {
  test.obj <- CreateSeuratObject(
    counts = pbmc.test,
    min.cells = 10,
    min.features = 30,
    meta.data = meta.data
  )
  # test normalization method
  aggregate.obj <- AggregateExpression(test.obj, assay = "RNA", group.by = "a", return.seurat = TRUE,
                                       normalization.method = "CLR")
  aggregate.obj.2 <- NormalizeData(aggregate.obj, normalization.method = "CLR")
  expect_equal(as.matrix(LayerData(aggregate.obj, assay = "RNA", layer = "data")),
               as.matrix(LayerData(aggregate.obj.2, assay = "RNA", layer = "data"))
  )
  # test margin
  aggregate.obj <- AggregateExpression(test.obj, assay = "RNA", group.by = "a", return.seurat = TRUE,
                                       margin = 2, normalization.method = "CLR" )
  aggregate.obj.2 <- NormalizeData(aggregate.obj, margin = 2, normalization.method = "CLR")
  expect_equal(as.matrix(LayerData(aggregate.obj, assay = "RNA", layer = "data")),
               as.matrix(LayerData(aggregate.obj.2, assay = "RNA", layer = "data"))
  )
  # test scale.factor
  aggregate.obj <- AggregateExpression(test.obj, assay = "RNA", group.by = "a", return.seurat = TRUE)
  aggregate.obj.2 <- NormalizeData(aggregate.obj)
  expect_equal(as.matrix(LayerData(aggregate.obj, assay = "RNA", layer = "data")),
               as.matrix(LayerData(aggregate.obj.2, assay = "RNA", layer = "data")))
  # test scale.factor
  aggregate.obj <- AggregateExpression(test.obj, assay = "RNA", group.by = "a", return.seurat = TRUE,
                                       scale.factor = 1e6)
  aggregate.obj.2 <- NormalizeData(aggregate.obj, scale.factor = 1e6)
  expect_equal(as.matrix(LayerData(aggregate.obj, assay = "RNA", layer = "data")),
               as.matrix(LayerData(aggregate.obj.2, assay = "RNA", layer = "data"))
  )
})

test_that("Aggregate Expression with return.seurat=TRUE returns normalized data in multi-assay case", {
  test.obj <- CreateSeuratObject(
    counts = pbmc.test,
    min.cells = 10,
    min.features = 30,
    meta.data = meta.data
  )
  test.obj[['RNA2']] <- CreateAssay5Object(counts = pbmc.test,  min.cells = 10,
                                          min.features = 30)
  # test normalization method
  aggregate.obj <- AggregateExpression(test.obj,  group.by = "a", return.seurat = TRUE,
                                       normalization.method = "CLR")
  aggregate.obj.2 <- NormalizeData(aggregate.obj, assay = "RNA2", normalization.method = "CLR")
  expect_equal(as.matrix(LayerData(aggregate.obj, assay = "RNA2", layer = "data")),
               as.matrix(LayerData(aggregate.obj.2, assay = "RNA2", layer = "data"))
  )
  # test margin
  aggregate.obj <- AggregateExpression(test.obj, group.by = "a", return.seurat = TRUE,
                                       margin = 2, normalization.method = "CLR" )
  aggregate.obj.2 <- NormalizeData(aggregate.obj, assay = "RNA2", margin = 2, normalization.method = "CLR")
  expect_equal(as.matrix(LayerData(aggregate.obj, assay = "RNA2", layer = "data")),
               as.matrix(LayerData(aggregate.obj.2, assay = "RNA2", layer = "data"))
  )
  # test scale.factor
  aggregate.obj <- AggregateExpression(test.obj,group.by = "a", return.seurat = TRUE,
                                       scale.factor = 1e6)
  aggregate.obj.2 <- NormalizeData(aggregate.obj, assay = "RNA2", scale.factor = 1e6)
  expect_equal(as.matrix(LayerData(aggregate.obj, assay = "RNA2", layer = "data")),
               as.matrix(LayerData(aggregate.obj.2, assay = "RNA2", layer = "data"))
  )
})

test_that("Aggregate Expression return.seurat=TRUE keeps group.by variables as meta data", {
  test.obj <- CreateSeuratObject(
    counts = pbmc.test,
    min.cells = 10,
    min.features = 30,
    meta.data = meta.data
  )
  test.obj$b <- rep(as.factor(c('c', 'd', 'e')), length.out = ncol(test.obj))
  aggregate.obj <- AggregateExpression(test.obj,  group.by = c("a", "b"), return.seurat = TRUE)
  expect_true(all(c("a", "b", "orig.ident") %in% colnames(aggregate.obj[[]])))
})

test_that("Aggregate Expression return.seurat=TRUE passes arguments to CreateSeuratObject", {
  test.obj <- CreateSeuratObject(
    counts = pbmc.test,
    min.cells = 10,
    min.features = 30,
    meta.data = meta.data
  )
  test.obj$b <- rep(as.factor(c('c', 'd', 'e')), length.out = ncol(test.obj))
  # check that min.cells works
  aggregate.obj <- AggregateExpression(test.obj,  group.by = c("a", "b"), return.seurat = TRUE,
                                       min.cells = 5)
  aggregate.obj.2 <- CreateSeuratObject(counts = aggregate.obj[['RNA']]$counts, min.cells=5)
  expect_equal(nrow(LayerData(aggregate.obj, layer = "counts")),
               nrow(LayerData(aggregate.obj.2, layer = "counts"))
               )
  # check that min.features works
  aggregate.obj <- AggregateExpression(test.obj,  group.by = c("a", "b"), return.seurat = TRUE,
                                       min.features = 50)
  aggregate.obj.2 <- CreateSeuratObject(counts = aggregate.obj[['RNA']]$counts, min.features=50)
  expect_equal(ncol(LayerData(aggregate.obj, layer = "counts")),
               ncol(LayerData(aggregate.obj.2, layer = "counts"))
  )
})

test.dat <- LayerData(object = object, layer = "data")[, 1:10]
rownames(x = test.dat) <- paste0("test-", rownames(x = test.dat))
suppressWarnings(object[["TEST"]] <- CreateAssayObject(data = test.dat))
test_that("AverageExpression with multiple assays", {
  avg.test <- AverageExpression(object = object, assays = "TEST", layer = "data")
  expect_equal(names(x = avg.test), "TEST")
  expect_equal(length(x = avg.test), 1)
  expect_equivalent(
    avg.test[[1]]['test-KHDRBS1', 1:3],
    c(a = 20.66116, b = 751.54445, c = 38.31418),
    tolerance = 1e-6
  )
  expect_equivalent(
    avg.test[[1]]['test-DNAJB1', 1:3] ,
    c(a = 230.50887, b = 50.50505, c = 65.35948),
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
    assay = "niche",
    cluster.name = "niches"
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

test_that("BuildNicheAssay works with FOV and VisiumV2 instances", {
  skip_if_not_installed("hdf5r")

  path.to.data = file.path("../testdata/visium")

  test.case <- Load10X_Spatial(
    path.to.data,
    assay = "Spatial",
    slice = "slice"
  )

  # populate a meta.data column with random labels
  random_labels <- sample(1:3, size = ncol(test.case), replace = TRUE)
  test.case[["random_labels"]] <- random_labels

  fov <- CreateFOV(
    GetTissueCoordinates(test.case[["slice"]]),
    type = "centroids",
    radius = Radius(test.case[["slice"]]),
    assay = "Spatial",
    key = "fov"
  )

  test.case[["fov"]] <- fov

  left <- BuildNicheAssay(
    test.case,
    fov = "fov",
    group.by = "random_labels"
  )

  right <- BuildNicheAssay(
    test.case,
    fov = "fov",
    group.by = "random_labels"
  )

  expect_equal(
    LayerData(left, layer = "scale.data"),
    LayerData(left, layer = "scale.data")
  )
})

test_that("AddModuleScore works in the multi-layer case", {
  cd_features <- list(c('CD79B','CD79A','CD3D','CD2','CD3E','CD7',
                        'CD14','CD68','CD247'))
  object1 <- AddModuleScore(object = object,
                            features = cd_features,
                            ctrl = 5,
                            name = 'CD_Features')
  
  split_object <- SplitObject(object, split.by = 'a')
  split_object1 <- AddModuleScore(object = split_object[[1]],
                                  features = cd_features,
                                  ctrl = 5,
                                  name = 'CD_Features')
  
  object2 <- split(object, f = object$a)
  object2 <- AddModuleScore(object = object2,
                            features = cd_features,
                            ctrl = 5,
                            name = 'CD_Features')
  
  expect_equal(
    object2$CD_Features1[object2$a=='a'],
    split_object1$CD_Features1
  )
  
  expect_false(
    all(object1$CD_Features1 == object2$CD_Features1)
  )
})
