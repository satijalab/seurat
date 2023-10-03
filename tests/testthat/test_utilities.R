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
  average.expression <- AverageExpression(object, layer = 'data', features = features)$RNA
  expect_equal(rownames(x = average.expression), features)
  expect_warning(AverageExpression(object, layer = 'data', features = "BAD"))
  expect_warning(AverageExpression(object, layer = "data", features = c(features, "BAD")))
})

test_that("AverageExpression with return.seurat", {
  # counts
  avg.counts <- AverageExpression(object, layer = "counts", return.seurat = TRUE, verbose = FALSE)
  expect_s4_class(object = avg.counts, "Seurat")
  avg.counts.mat <- AverageExpression(object, layer = 'counts')$RNA
  expect_equal(unname(as.matrix(LayerData(avg.counts[["RNA"]], layer = "counts"))),
               unname(as.matrix(avg.counts.mat)))
  avg.data <- LayerData(avg.counts[["RNA"]], layer = "data")

  expect_equivalent(
    as.matrix(NormalizeData(avg.counts.mat)),
    as.matrix(avg.data),
    tolerance = 1e-6
  )

  avg.scale <- LayerData(avg.counts[["RNA"]], layer = "scale.data")
  expect_equal(
    avg.scale['MS4A1', ],
    c(a = -0.8141426, b = 1.1162108, c = -0.3020683),
    tolerance = 1e-6
  )
  expect_equal(
    avg.scale['SPON2', ],
    c(a = 0.3387626, b = 0.7866155, c = -1.1253781),
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

