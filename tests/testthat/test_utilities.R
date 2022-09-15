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
object <- SetIdent(object, value = 'a')

test_that("AverageExpression works for different slots", {
  average.expression <- AverageExpression(object, slot = 'data')$RNA
  expect_equivalent(
    average.expression['KHDRBS1', 1:3],
    c(a = 7.278237e-01, b = 1.658166e+14, c = 1.431902e-01),
    tolerance = 1e-6
  )
  expect_equivalent(
    average.expression['DNAJB1', 1:3] ,
    c(a = 1.374079e+00, b = 5.100840e-01, c = 5.011655e-01),
    tolerance = 1e-6
  )
  avg.counts <- AverageExpression(object, slot = 'counts')$RNA
  expect_equal(
    avg.counts['MS4A1', ],
    c(a = 0.37037037, b = 0.3461538, c = 0.3333333),
    tolerance = 1e-6
  )
  expect_equal(
    avg.counts['SPON2', ],
    c(a = 0.5185185, b = 0.6153846, c = 0.08333333),
    tolerance = 1e-6
  )
  expect_warning(AverageExpression(object, slot = 'scale.data'))
  object <- ScaleData(object = object, verbose = FALSE)
  avg.scale <- AverageExpression(object, slot = "scale.data")$RNA
  expect_equal(
    avg.scale['MS4A1', ],
    c(a = 0.02092088, b = -0.004769018, c = -0.018369549),
    tolerance = 1e-6
  )
  expect_equal(
    avg.scale['SPON2', ],
    c(a = 0.1052434, b = 0.2042827, c = -0.3397051),
    tolerance = 1e-6
  )
})

test_that("AverageExpression handles features properly", {
  features <- rownames(x = object)[1:10]
  average.expression <- AverageExpression(object, slot = 'data', features = features)$RNA
  expect_equal(rownames(x = average.expression), features)
  expect_warning(AverageExpression(object, slot = 'data', features = "BAD"))
  expect_warning(AverageExpression(object, slot = "data", features = c(features, "BAD")))
})

test_that("AverageExpression with return.seurat", {
  # counts
  avg.counts <- AverageExpression(object, slot = "counts", return.seurat = TRUE, verbose = FALSE)
  expect_s4_class(object = avg.counts, "Seurat")
  avg.counts.mat <- AverageExpression(object, slot = 'counts')$RNA
  expect_equal(as.matrix(GetAssayData(avg.counts[["RNA"]], slot = "counts")), avg.counts.mat)
  avg.data <- GetAssayData(avg.counts[["RNA"]], slot = "data")
  expect_equal(
    avg.data['MS4A1', ],
    c(a = 0.31508105, b = 0.2972515, c = 0.2876821),
    tolerance = 1e-6
  )
  expect_equal(
    avg.data['SPON2', ],
    c(a = 0.4177352, b = 0.4795731, c = 0.08004271),
    tolerance = 1e-6
  )
  avg.scale <- GetAssayData(avg.counts[["RNA"]], slot = "scale.data")
  expect_equal(
    avg.scale['MS4A1', ],
    c(a = 1.0841908, b = -0.1980056, c = -0.8861852),
    tolerance = 1e-6
  )
  expect_equal(
    avg.scale['SPON2', ],
    c(a = 0.4275778, b = 0.7151260, c = -1.1427038),
    tolerance = 1e-6
  )
  # data
  avg.data <- AverageExpression(object, slot = "data", return.seurat = TRUE, verbose = FALSE)
  expect_s4_class(object = avg.data, "Seurat")
  avg.data.mat <- AverageExpression(object, slot = 'data')$RNA
  expect_equal(as.matrix(GetAssayData(avg.data[["RNA"]], slot = "counts")), avg.data.mat)
  expect_equal(unname(as.matrix(GetAssayData(avg.data[["RNA"]], slot = "data"))), unname(log1p(x = avg.data.mat)))
  avg.scale <- GetAssayData(avg.data[["RNA"]], slot = "scale.data")
  expect_equal(
    avg.scale['MS4A1', ],
    c(a = 0.721145238, b = -1.1415734, c = 0.4204281),
    tolerance = 1e-6
  )
  expect_equal(
    avg.scale['SPON2', ],
    c(a = 0.08226771, b = 0.9563249, c = -1.0385926),
    tolerance = 1e-6
  )
  # scale.data
  object <- ScaleData(object = object, verbose = FALSE)
  avg.scale <- AverageExpression(object, slot = "scale.data", return.seurat = TRUE, verbose = FALSE)
  expect_s4_class(object = avg.scale, "Seurat")
  avg.scale.mat <- AverageExpression(object, slot = 'scale.data')$RNA
  expect_equal(unname(as.matrix(GetAssayData(avg.scale[["RNA"]], slot = "scale.data"))), unname(avg.scale.mat))
  expect_true(all(is.na(GetAssayData(avg.scale[["RNA"]], slot = "data"))))
  expect_equal(GetAssayData(avg.scale[["RNA"]], slot = "counts"), matrix())
})

test.dat <- GetAssayData(object = object, slot = "data")
rownames(x = test.dat) <- paste0("test-", rownames(x = test.dat))
object[["TEST"]] <- CreateAssayObject(data = test.dat)

test_that("AverageExpression with multiple assays", {
  avg.test <- AverageExpression(object = object, assays = "TEST")
  expect_equal(names(x = avg.test), "TEST")
  expect_equal(length(x = avg.test), 1)
  expect_equivalent(
    avg.test[[1]]['test-KHDRBS1', 1:3],
    c(a = 7.278237e-01, b = 1.658166e+14, c = 1.431902e-01),
    tolerance = 1e-6
  )
  expect_equivalent(
    avg.test[[1]]['test-DNAJB1', 1:3] ,
    c(a = 1.374079e+00, b = 5.100840e-01, c = 5.011655e-01),
    tolerance = 1e-6
  )
  avg.all <- AverageExpression(object = object)
  expect_equal(names(x = avg.all), c("RNA", "TEST"))
  expect_equal(length(x = avg.all), 2)
})
