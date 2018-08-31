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
object <- CreateSeuratObject(counts = pbmc.test,
                             meta.data = fake.meta.data)
test_that("object initialization actually creates seurat object", {
  expect_is(object, "Seurat")
})

test_that("meta.data slot generated correctly", {
  expect_equal(dim(object[[]]), c(80, 4))
  expect_equal(colnames(object[[]]), c("orig.ident", "nCount_RNA", "nFeature_RNA", "FMD"))
  expect_equal(rownames(object[[]]), colnames(object))
  expect_equal(object[["nFeature_RNA"]][1:5, ], c(47, 52, 50, 56, 53))
  expect_equal(object[["nCount_RNA"]][75:80, ], c(228, 527, 202, 157, 150, 233))
})

object.filtered <- CreateSeuratObject(
  counts = pbmc.test,
  min.cells = 10,
  min.features = 30
)

test_that("Filtering handled properly", {
  expect_equal(nrow(x = GetAssayData(object = object.filtered, slot = "counts")), 162)
  expect_equal(ncol(x = GetAssayData(object = object.filtered, slot = "counts")), 75)
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
  expect_equal(GetAssayData(object = object, slot = "counts"),
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

normalized.data <- LogNormalize(data = GetAssayData(object = object[["RNA"]], slot = "counts"))
test_that("LogNormalize normalizes properly", {
  expect_equal(
    LogNormalize(data = GetAssayData(object = object[["RNA"]], slot = "counts"), verbose = FALSE),
    LogNormalize(data = as.data.frame(as.matrix(GetAssayData(object = object[["RNA"]], slot = "counts"))), verbose = FALSE)
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
  vars.to.regress = "nCount_RNA",
  features = rownames(x = object)[1:10],
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
  vars.to.regress = "nCount_RNA",
  features = rownames(x = object)[1:10],
  verbose = FALSE,
  model.use = "negbinom")

test_that("Negative binomial regression works as expected", {
  expect_equal(dim(x = GetAssayData(object = object[["RNA"]], slot = "scale.data")), c(10, 80))
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[1, 1], -0.5888811, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[5, 25], -0.2553394, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[10, 80], -0.1921429, tolerance = 1e-6)
})

test_that("Regression error handling checks out", {
  expect_error(ScaleData(object, vars.to.regress = "nCount_RNA", model.use = "not.a.model"))
})

object <- ScaleData(
  object,
  vars.to.regress = "nCount_RNA",
  features = rownames(x = object)[1:10],
  verbose = FALSE,
  model.use = "poisson")

test_that("Poisson regression works as expected", {
  expect_equal(dim(x = GetAssayData(object = object[["RNA"]], slot = "scale.data")), c(10, 80))
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[1, 1], -1.011717, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[5, 25], 0.05575307, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[10, 80], -0.1662119, tolerance = 1e-6)
})


#Tests for SampleUMI
#--------------------------------------------------------------------------------
context("SampleUMI")

downsampled.umis <- SampleUMI(
  data = GetAssayData(object = object, slot = "counts"),
  max.umi = 100,
  verbose = FALSE
)
downsampled.umis.p.cell <- SampleUMI(
  data = GetAssayData(object = object, slot = "counts"),
  max.umi = seq(50, 840, 10),
  verbose = FALSE,
  upsample = TRUE
)
test_that("SampleUMI gives reasonable downsampled/upsampled UMI counts", {
  expect_true(!any(colSums(x = downsampled.umis) < 30, colSums(x = downsampled.umis) > 120))
  expect_error(SampleUMI(data = GetAssayData(object = object, slot = "raw.data"), max.umi = rep(1, 5)))
  expect_true(!is.unsorted(x = colSums(x = downsampled.umis.p.cell)))
  expect_error(SampleUMI(
    data = GetAssayData(object = object, slot = "counts"),
    max.umi = seq(50, 900, 10),
    verbose = FALSE,
    upsample = TRUE
  ))
})

# Tests for FindVariableFeautres
# --------------------------------------------------------------------------------
context("FindVariableFeatures")

object <- FindVariableFeatures(object = object, selection.method = "mean.var.plot", verbose = FALSE)
test_that("mean.var.plot selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PTGDR", "SATB1", "ZNF330", "S100B"))
  expect_equal(length(x = VariableFeatures(object = object)), 20)
  expect_equal(HVFInfo(object = object[["RNA"]])$mean[1:2], c(8.328927, 8.444462), tolerance = 1e6)
  expect_equal(HVFInfo(object = object[["RNA"]])$dispersion[1:2], c(10.552507, 10.088223), tolerance = 1e6)
  expect_equal(as.numeric(HVFInfo(object = object[["RNA"]])$dispersion.scaled[1:2]), c(0.1113214, -0.1113214), tolerance = 1e6)
})

object <- FindVariableFeatures(object, selection.method = "dispersion", verbose = FALSE)
test_that("dispersion selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PCMT1", "PPBP", "LYAR", "VDAC3"))
  expect_equal(length(x = VariableFeatures(object = object)), 230)
  expect_equal(HVFInfo(object = object[["RNA"]])$mean[1:2], c(8.328927, 8.444462), tolerance = 1e6)
  expect_equal(HVFInfo(object = object[["RNA"]])$dispersion[1:2], c(10.552507, 10.088223), tolerance = 1e6)
  expect_equal(as.numeric(HVFInfo(object = object[["RNA"]])$dispersion.scaled[1:2]), c(0.1113214, -0.1113214), tolerance = 1e6)
  expect_true(!is.unsorted(rev(HVFInfo(object = object[["RNA"]])[VariableFeatures(object = object), "dispersion"])))
})

object <- FindVariableFeatures(object, selection.method = "vst", verbose = FALSE)
test_that("vst selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PPBP", "IGLL5", "VDAC3", "CD1C"))
  expect_equal(length(x = VariableFeatures(object = object)), 230)
  expect_equal(unname(object[["RNA"]][["variance", drop = TRUE]][1:2]), c(1.0251582, 1.2810127), tolerance = 1e6)
  expect_equal(unname(object[["RNA"]][["variance.expected", drop = TRUE]][1:2]), c(1.1411616, 2.7076228), tolerance = 1e6)
  expect_equal(unname(object[["RNA"]][["variance.standardized", drop = TRUE]][1:2]), c(0.8983463, 0.4731134), tolerance = 1e6)
  expect_true(!is.unsorted(rev(object[["RNA"]][["variance.standardized", drop = TRUE]][VariableFeatures(object = object)])))
})

# # Tests for FilterCells
# # --------------------------------------------------------------------------------
# context("FilterCells")
#
# object.filtered <- FilterCells(
#   object = object,
#   subset.names = c("nFeature_RNA", "nCount_RNA"),
#   low.thresholds = c(20, 100)
# )
#
# test_that("FilterCells low thresholds work properly", {
#   expect_equal(ncol(x = object.filtered), 62)
#   expect_true(!any(object.filtered["nFeature_RNA"] < 20))
#   expect_true(!any(object.filtered["nCount_RNA"] < 100))
# })
#
# object.filtered <- FilterCells(
#   object = object,
#   subset.names = c("nFeature_RNA", "nCount_RNA"),
#   high.thresholds = c(50, 300)
# )
#
# test_that("FilterCells high thresholds work properly", {
#   expect_equal(ncol(x = object.filtered), 35)
#   expect_true(!any(object.filtered["nFeature_RNA"] > 50))
#   expect_true(!any(object.filtered["nCount_RNA"] > 300))
# })
#
# test_that("FilterCells handles input correctly", {
#   expect_error(FilterCells(object, subset.names = c("nGene", "nCount_RNA"), high.thresholds = 30))
#   expect_error(FilterCells(object, subset.names = c("nGene", "nCount_RNA"), low.thresholds = 20))
#   expect_error(FilterCells(object, subset.names = c("nGene"), high.thresholds = c(30, 300)))
# })
