set.seed(42)

pbmc.file <- system.file('extdata', 'pbmc_raw.txt', package = 'Seurat')
pbmc.test <- as(as.matrix(read.table(pbmc.file, sep = "\t", row.names = 1)), "dgCMatrix")

meta.data <- data.frame(
  a = rep(as.factor(c('a', 'b', 'c')), length.out = ncol(pbmc.test)),
  row.names = colnames(pbmc.test)
)

object.filtered <- CreateSeuratObject(
  counts = pbmc.test,
  min.cells = 10,
  min.features = 30,
  meta.data = meta.data
)

test_that("AverageExpression", {
  object <- SetIdent(object.filtered, value = 'a')
  average.expression <- AverageExpression(object, slot = 'data')$RNA

  expect_equivalent(average.expression['KHDRBS1', 1:3],
                    c(a = 7.278237e-01, b = 1.658166e+14, c = 1.431902e-01),
                    tolerance = 1e-6
                    )
  expect_equivalent(average.expression['DNAJB1', 1:3] ,
                    c(a = 1.374079e+00, b = 5.100840e-01, c = 5.011655e-01),
                    tolerance = 1e-6)
})
