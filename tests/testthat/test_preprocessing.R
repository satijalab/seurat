# Tests for functions dependent on a seurat object
set.seed(42)

pbmc.file <- system.file('extdata', 'pbmc_raw.txt', package = 'Seurat')
pbmc.test <- as(as.matrix(read.table(pbmc.file, sep = "\t", row.names = 1)), "dgCMatrix")

# Tests for object creation (via CreateSeuratObject)
# --------------------------------------------------------------------------------
context("Object creation")

fake.meta.data <- data.frame(rep(1, ncol(pbmc.test)))
rownames(fake.meta.data) <- colnames(pbmc.test)
colnames(fake.meta.data) <- "FMD"
object <- CreateSeuratObject(raw.data = pbmc.test,
                             normalization.method = "LogNormalize",
                             do.scale = T,
                             meta.data = fake.meta.data,
                             # save.raw = F,
                             display.progress = F)
test_that("object initialization actually creates seurat object", {
  expect_is(object, "seurat")
})

# test_that("save.raw option handled properly", {
#   expect_equal(dim(object@raw.data), c(1, 1))
#   expect_equal(object@raw.data[1, 1], NA)
# })

test_that("meta.data slot generated correctly", {
  expect_equal(dim(object@meta.data), c(80, 4))
  expect_equal(colnames(object@meta.data), c("nGene", "nUMI", "FMD", "orig.ident"))
  expect_equal(rownames(object@meta.data), colnames(object@data))
  expect_equal(object@meta.data$nGene[1:5], c(47, 52, 50, 56, 53))
  expect_equal(object@meta.data$nUMI[75:80], c(228, 527, 202, 157, 150, 233))
})

test_that("normalization and scaling run during object creation process", {
  expect_equal(object@data[2, 1], 4.968821, tolerance = 1e-6)
  expect_equal(object@data[174, 80], 5.554937, tolerance = 1e-6)
  expect_equal(object@scale.data[2, 1], 1.917418, tolerance = 1e-6)
  expect_equal(object@scale.data[174, 80], 1.998957, tolerance = 1e-6)
})

object <- CreateSeuratObject(raw.data = pbmc.test,
                             is.expr = 2,
                             min.cells = 3,
                             min.genes = 10)

test_that("Expression threshold zeros out proper entries in expression matrix", {
  expect_equal(nnzero(object@data), 1939)
  expect_equal(nnzero(object@raw.data), 3197)
})

test_that("Genes are filtered out based on min.cells", {
  expect_equal(nrow(object@raw.data), 162)
  expect_equal(nrow(object@data), 162)
})

# Tests for Read10X
# --------------------------------------------------------------------------------
context("Read10X")

test_that("Read10X handles missing files properly", {
  expect_error(Read10X("."))
  expect_error(Read10X("./notadir/"))
})

test.data <- Read10X("../testdata/")
test_that("Read10X creates sparse matrix", {
  expect_is(test.data, "dgTMatrix")
})


# Tests for NormalizeData
# --------------------------------------------------------------------------------
context("NormalizeData")

test.object <- object
test.object@raw.data <- NULL

test_that("NormalizeData error handling", {
  expect_error(NormalizeData(object, assay.type = "FAKE"))
  expect_equal(object@data, NormalizeData(object, normalization.method = NULL)@data)
  expect_error(NormalizeData(test.object))
})

object <- NormalizeData(object, scale.factor = 1e6, display.progress = F)
test_that("NormalizeData scales properly", {
  expect_equal(object@data[2, 2], 9.304742, tolerance = 1e-6)
  expect_equal(object@data[161, 55], 7.659003, tolerance = 1e-6)
  expect_equal(object@calc.params$NormalizeData$scale.factor, 1e6)
  expect_equal(object@calc.params$NormalizeData$normalization.method, "LogNormalize")
})

normalized.data <- LogNormalize(data = object@raw.data)
test_that("LogNormalize normalizes properly", {
  expect_equal(LogNormalize(data = object@raw.data, display.progress = F), LogNormalize(data = as.data.frame(as.matrix(object@raw.data)), display.progress = F))
})

# Tests for ScaleData
# --------------------------------------------------------------------------------
context("ScaleData")
object <- ScaleData(object, do.cpp = F, display.progress = F)
test_that("Old R implementation (ScaleDataR) works properly", {
  expect_equal(object@scale.data[1, 1], -0.2995232, tolerance = 1e-6)
  expect_equal(object@scale.data[75, 25], 1.993555, tolerance = 1e-6)
  expect_equal(object@scale.data[162, 59], -0.5480965, tolerance = 1e-6)
})

object <- ScaleData(object, display.progress = F)
test_that("ScaleData returns expected values when input is a sparse matrix", {
  expect_equal(object@scale.data[1, 1], -0.2995232, tolerance = 1e-6)
  expect_equal(object@scale.data[75, 25], 1.993555, tolerance = 1e-6)
  expect_equal(object@scale.data[162, 59], -0.5480965, tolerance = 1e-6)
})

object@data <- as.matrix(object@data)
object <- ScaleData(object)
test_that("ScaleData returns expected values when input is not sparse", {
  expect_equal(object@scale.data[1, 1], -0.2995232, tolerance = 1e-6)
  expect_equal(object@scale.data[75, 25], 1.993555, tolerance = 1e-6)
  expect_equal(object@scale.data[162, 59], -0.5480965, tolerance = 1e-6)
})

# Tests for various regression techniques
context("Regression")

object <- ScaleData(object,
                    vars.to.regress = "nUMI",
                    genes.use = rownames(object@data)[1:10],
                    display.progress = F,
                    model.use = "linear")

test_that("Linear regression works as expected", {
  expect_equal(dim(object@scale.data), c(10, 59))
  expect_equal(object@scale.data[1, 1], -0.4039399, tolerance = 1e-6)
  expect_equal(object@scale.data[5, 25], -0.9216946, tolerance = 1e-6)
  expect_equal(object@scale.data[10, 59], -0.5475258, tolerance = 1e-6)
})

object <- ScaleData(object,
                    vars.to.regress = "nUMI",
                    genes.use = rownames(object@data)[1:10],
                    display.progress = F,
                    model.use = "negbinom")

test_that("Negative binomial regression works as expected", {
  expect_equal(dim(object@scale.data), c(10, 59))
  expect_equal(object@scale.data[1, 1], -0.4150756, tolerance = 1e-6)
  expect_equal(object@scale.data[5, 25], -0.6586565, tolerance = 1e-6)
  expect_equal(object@scale.data[10, 59], -0.4537495, tolerance = 1e-6)
})

object <- suppressWarnings(RegressOutNB(object = object,
                                        latent.vars = "nUMI",
                                        genes.regress = rownames(object@data)[1:10]))

test_that("Other negative binomial regression works as expected", {
  expect_equal(dim(object@scale.data), c(10, 59))
  expect_equal(object@scale.data[1, 1], -0.274358, tolerance = 1e-6)
  expect_equal(object@scale.data[5, 25], -0.5623909, tolerance = 1e-6)
  expect_equal(object@scale.data[10, 59], -0.3456492, tolerance = 1e-6)
})

test_that("Regression error handling checks out", {
  expect_error(ScaleData(object, vars.to.regress = "nUMI", model.use = "not.a.model"))
})

object <- ScaleData(object,
                    vars.to.regress = "nUMI",
                    genes.use = rownames(object@data)[1:10],
                    display.progress = F,
                    model.use = "poisson")

test_that("Poisson regression works as expected", {
  expect_equal(dim(object@scale.data), c(10, 59))
  expect_equal(object@scale.data[1, 1], -0.6115097, tolerance = 1e-6)
  expect_equal(object@scale.data[5, 25], -0.5971585, tolerance = 1e-6)
  expect_equal(object@scale.data[10, 59], -0.4533085, tolerance = 1e-6)
})

if (detectCores() > 1) {
  object <- ScaleData(object,
  vars.to.regress = "nUMI",
  genes.use = rownames(object@data)[1:10],
  display.progress = F,
  model.use = "linear",
  do.par = TRUE,
  num.cores = 2)

  test_that("Parallelization works", {
    expect_equal(dim(object@scale.data), c(10, 59))
    expect_equal(object@scale.data[1, 1], -0.4039399, tolerance = 1e-6)
    expect_equal(object@scale.data[5, 25], -0.9216946, tolerance = 1e-6)
    expect_equal(object@scale.data[10, 59], -0.5475258, tolerance = 1e-6)
  })
}


# Tests for SampleUMI
# --------------------------------------------------------------------------------
context("SampleUMI")

downsampled.umis <- SampleUMI(object@raw.data, max.umi = 100, progress.bar = F)
downsampled.umis.p.cell <- SampleUMI(object@raw.data, max.umi = seq(50, 630, 10), progress.bar = F, upsample = T)
test_that("SampleUMI gives reasonable downsampled/upsampled UMI counts", {
  expect_true(!any(colSums(downsampled.umis) < 80, colSums(downsampled.umis) > 120))
  expect_error(SampleUMI(object@raw.data, max.umi = rep(1, 5)))
  expect_true(!is.unsorted(colSums(downsampled.umis.p.cell)))
})

# Tests for FindVariableGenes
# --------------------------------------------------------------------------------
context("FindVariableGenes")

object <- FindVariableGenes(object, display.progress = F, do.cpp = F, do.plot = F)
test_that("R implementation of FindVariableGenes returns expected values", {
  expect_equal(object@var.genes[1:2], c("MS4A1", "CD2"))
  expect_equal(length(object@var.genes), 10)
  expect_equal(object@hvg.info$gene.mean[1:2], c(8.856202, 10.472897), tolerance = 1e6)
  expect_equal(object@hvg.info$gene.dispersion[1:2], c(12.41696, 12.23218), tolerance = 1e6)
  expect_equal(as.numeric(object@hvg.info$gene.dispersion.scaled[1:2]), c(1.7506589, 1.1963021), tolerance = 1e6)
  expect_true(!is.unsorted(rev(object@hvg.info$gene.dispersion)))
})

object <- FindVariableGenes(object, display.progress = F, do.plot = F)
test_that("C++ implementation of FindVariableGenes returns expected values", {
  expect_equal(object@var.genes[1:2], c("MS4A1", "CD2"))
  expect_equal(length(object@var.genes), 10)
  expect_equal(object@hvg.info$gene.mean[1:2], c(8.856202, 10.472897), tolerance = 1e6)
  expect_equal(object@hvg.info$gene.dispersion[1:2], c(12.41696, 12.23218), tolerance = 1e6)
  expect_equal(as.numeric(object@hvg.info$gene.dispersion.scaled[1:2]), c(1.7506589, 1.1963021), tolerance = 1e6)
  expect_true(!is.unsorted(rev(object@hvg.info$gene.dispersion)))
  expect_warning(FindVariableGenes(object, display.progress = F, do.plot = F, mean.function = mean))
  expect_warning(FindVariableGenes(object, display.progress = F, do.plot = F, dispersion.function = ExpSD))
})

var.genes <- FindVariableGenes(object, set.var.genes = F)
test_that("Option to only return vector of genes works", {
  expect_equal(length(var.genes), 10)
  expect_equal(length(var.genes), length(object@var.genes))
  expect_equal(var.genes[1:2], c("MS4A1", "CD2"))
})

# Tests for FilterCells
# --------------------------------------------------------------------------------
context("FilterCells")

object.filtered <- FilterCells(object, subset.names = c("nGene", "nUMI"), low.thresholds = c(20, 100))

test_that("FilterCells low thresholds work properly", {
  expect_equal(length(object.filtered@cell.names), 31)
  expect_true(!any(object.filtered@meta.data$nGene < 20))
  expect_true(!any(object.filtered@meta.data$nUMI < 100))
})

object.filtered <- FilterCells(object, subset.names = c("nGene", "nUMI"), high.thresholds = c(30, 300))

test_that("FilterCells high thresholds work properly", {
  expect_equal(length(object.filtered@cell.names), 38)
  expect_true(!any(object.filtered@meta.data$nGene > 30))
  expect_true(!any(object.filtered@meta.data$nUMI > 300))
})

test_that("FilterCells handles input correctly", {
  expect_error(FilterCells(object, subset.names = c("nGene", "nUMI"), high.thresholds = 30))
  expect_error(FilterCells(object, subset.names = c("nGene", "nUMI"), low.thresholds = 20))
  expect_error(FilterCells(object, subset.names = c("nGene"), high.thresholds = c(30, 300)))
})
