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
                             meta.data = fake.meta.data)
test_that("object initialization actually creates seurat object", {
  expect_is(object, "Seurat")
})

test_that("meta.data slot generated correctly", {
  expect_equal(dim(object[]), c(80, 4))
  expect_equal(colnames(object[]), c("orig.ident", "nUMI", "nFeature_RNA", "FMD"))
  expect_equal(rownames(object[]), colnames(object))
  expect_equal(object["nFeature_RNA"][1:5, ], c(47, 52, 50, 56, 53))
  expect_equal(object["nUMI"][75:80, ], c(228, 527, 202, 157, 150, 233))
})

object.filtered <- CreateSeuratObject(
  raw.data = pbmc.test,
  min.cells = 10,
  min.features = 30
)

test_that("Filtering handled properly", {
  expect_equal(nrow(x = GetAssayData(object = object.filtered, slot = "raw.data")), 162)
  expect_equal(ncol(x = GetAssayData(object = object.filtered, slot = "raw.data")), 75)
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

test_that("NormalizeData error handling", {
  expect_error(NormalizeData(object = object, assay.use = "FAKE"))
  expect_equal(GetAssayData(object = object, slot = "raw.data"),
               GetAssayData(object = NormalizeData(
                  object = object,
                  normalization.method = NULL),
                slot = "data"))
})

object <- NormalizeData(object = object, verbose = FALSE, scale.factor = 1e6)
test_that("NormalizeData scales properly", {
  expect_equal(GetAssayData(object = object, slot = "data")[2, 1], 9.567085, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object, slot = "data")[161, 55], 8.415309, tolerance = 1e-6)
  expect_equal(Command(object = object, command = "NormalizeData.RNA", value = "scale.factor"), 1e6)
  expect_equal(Command(object = object, command = "NormalizeData.RNA", value = "normalization.method"), "LogNormalize")
})

normalized.data <- LogNormalize(data = GetAssayData(object = object[["RNA"]], slot = "raw.data"))
test_that("LogNormalize normalizes properly", {
  expect_equal(
    LogNormalize(data = GetAssayData(object = object[["RNA"]], slot = "raw.data"), verbose = FALSE),
    LogNormalize(data = as.data.frame(as.matrix(GetAssayData(object = object[["RNA"]], slot = "raw.data"))), verbose = FALSE)
  )
})

# Tests for ScaleData
# --------------------------------------------------------------------------------
context("ScaleData")
object <- ScaleData(object, verbose = FALSE)
test_that("ScaleData returns expected values when input is a sparse matrix", {
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[1, 1], -0.4148587, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[75, 25], -0.2562305, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[162, 59], -0.4363939, tolerance = 1e-6)
})

object[["RNA"]] <- SetAssayData(
  object = object[["RNA"]],
  slot = "data",
  new.data = as.matrix(GetAssayData(object = object[["RNA"]], slot = "data"))
)
object <- ScaleData(object = object)
test_that("ScaleData returns expected values when input is not sparse", {
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[1, 1], -0.4148587, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[75, 25], -0.2562305, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[162, 59], -0.4363939, tolerance = 1e-6)
})

# Tests for various regression techniques
context("Regression")

object <- ScaleData(
  object = object,
  vars.to.regress = "nUMI",
  features.use = rownames(x = object)[1:10],
  verbose = FALSE,
  model.use = "linear")

test_that("Linear regression works as expected", {
  expect_equal(dim(x = GetAssayData(object = object[["RNA"]], slot = "scale.data")), c(10, 80))
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[1, 1], -0.6436435, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[5, 25], -0.09035383, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[10, 80], -0.2723782, tolerance = 1e-6)
})

object <- ScaleData(
  object,
  vars.to.regress = "nUMI",
  features.use = rownames(x = object)[1:10],
  verbose = FALSE,
  model.use = "negbinom")

test_that("Negative binomial regression works as expected", {
  expect_equal(dim(x = GetAssayData(object = object[["RNA"]], slot = "scale.data")), c(10, 80))
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[1, 1], -0.5888811, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[5, 25], -0.2553394, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[10, 80], -0.1921429, tolerance = 1e-6)
})

test_that("Regression error handling checks out", {
  expect_error(ScaleData(object, vars.to.regress = "nUMI", model.use = "not.a.model"))
})

object <- ScaleData(
  object,
  vars.to.regress = "nUMI",
  features.use = rownames(x = object)[1:10],
  verbose = FALSE,
  model.use = "poisson")

test_that("Poisson regression works as expected", {
  expect_equal(dim(x = GetAssayData(object = object[["RNA"]], slot = "scale.data")), c(10, 80))
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[1, 1], -1.011717, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[5, 25], 0.05575307, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[10, 80], -0.1662119, tolerance = 1e-6)
})


# Tests for SampleUMI
# --------------------------------------------------------------------------------
# context("SampleUMI")
#
# downsampled.umis <- SampleUMI(object@raw.data, max.umi = 100, progress.bar = F)
# downsampled.umis.p.cell <- SampleUMI(object@raw.data, max.umi = seq(50, 630, 10), progress.bar = F, upsample = T)
# test_that("SampleUMI gives reasonable downsampled/upsampled UMI counts", {
#   expect_true(!any(colSums(downsampled.umis) < 80, colSums(downsampled.umis) > 120))
#   expect_error(SampleUMI(object@raw.data, max.umi = rep(1, 5)))
#   expect_true(!is.unsorted(colSums(downsampled.umis.p.cell)))
# })

# Tests for FindVariableGenes
# --------------------------------------------------------------------------------
# context("FindVariableGenes")
#
# object <- FindVariableGenes(object, display.progress = F, do.cpp = F, do.plot = F)
# test_that("R implementation of FindVariableGenes returns expected values", {
#   expect_equal(object@var.genes[1:2], c("MS4A1", "CD2"))
#   expect_equal(length(object@var.genes), 10)
#   expect_equal(object@hvg.info$gene.mean[1:2], c(8.856202, 10.472897), tolerance = 1e6)
#   expect_equal(object@hvg.info$gene.dispersion[1:2], c(12.41696, 12.23218), tolerance = 1e6)
#   expect_equal(as.numeric(object@hvg.info$gene.dispersion.scaled[1:2]), c(1.7506589, 1.1963021), tolerance = 1e6)
#   expect_true(!is.unsorted(rev(object@hvg.info$gene.dispersion)))
# })
#
# object <- FindVariableGenes(object, display.progress = F, do.plot = F)
# test_that("C++ implementation of FindVariableGenes returns expected values", {
#   expect_equal(object@var.genes[1:2], c("MS4A1", "CD2"))
#   expect_equal(length(object@var.genes), 10)
#   expect_equal(object@hvg.info$gene.mean[1:2], c(8.856202, 10.472897), tolerance = 1e6)
#   expect_equal(object@hvg.info$gene.dispersion[1:2], c(12.41696, 12.23218), tolerance = 1e6)
#   expect_equal(as.numeric(object@hvg.info$gene.dispersion.scaled[1:2]), c(1.7506589, 1.1963021), tolerance = 1e6)
#   expect_true(!is.unsorted(rev(object@hvg.info$gene.dispersion)))
#   expect_warning(FindVariableGenes(object, display.progress = F, do.plot = F, mean.function = mean))
#   expect_warning(FindVariableGenes(object, display.progress = F, do.plot = F, dispersion.function = ExpSD))
# })
#
# var.genes <- FindVariableGenes(object, set.var.genes = F)
# test_that("Option to only return vector of genes works", {
#   expect_equal(length(var.genes), 10)
#   expect_equal(length(var.genes), length(object@var.genes))
#   expect_equal(var.genes[1:2], c("MS4A1", "CD2"))
# })
#
# object2 <- FindVariableGenes(object, display.progress = F, do.plot = F, do.recalc = FALSE)
# test_that("do.recalc doesn't change vector of variable genes", {
#   expect_equal(intersect(object@var.genes, object2@var.genes), object@var.genes)
# })

# Tests for FilterCells
# --------------------------------------------------------------------------------
context("FilterCells")

object.filtered <- FilterCells(
  object = object,
  subset.names = c("nFeature_RNA", "nUMI"),
  low.thresholds = c(20, 100)
)

test_that("FilterCells low thresholds work properly", {
  expect_equal(ncol(x = object.filtered), 62)
  expect_true(!any(object.filtered["nFeature_RNA"] < 20))
  expect_true(!any(object.filtered["nUMI"] < 100))
})

object.filtered <- FilterCells(
  object = object,
  subset.names = c("nFeature_RNA", "nUMI"),
  high.thresholds = c(50, 300)
)

test_that("FilterCells high thresholds work properly", {
  expect_equal(ncol(x = object.filtered), 35)
  expect_true(!any(object.filtered["nFeature_RNA"] > 50))
  expect_true(!any(object.filtered["nUMI"] > 300))
})

test_that("FilterCells handles input correctly", {
  expect_error(FilterCells(object, subset.names = c("nGene", "nUMI"), high.thresholds = 30))
  expect_error(FilterCells(object, subset.names = c("nGene", "nUMI"), low.thresholds = 20))
  expect_error(FilterCells(object, subset.names = c("nGene"), high.thresholds = c(30, 300)))
})
