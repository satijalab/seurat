# Tests for functions dependent on a seurat object
set.seed(42)

pbmc.file <- system.file('extdata', 'pbmc_raw.txt', package = 'Seurat')
pbmc.test <- as.sparse(x = as.matrix(read.table(pbmc.file, sep = "\t", row.names = 1)))

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

#this should be moved to seurat object
# test_that("meta.data slot generated correctly", {
#   expect_equal(dim(object[[]]), c(80, 4))
#   expect_equal(colnames(object[[]]), c("orig.ident", "nCount_RNA", "nFeature_RNA", "FMD"))
#   expect_equal(rownames(object[[]]), colnames(object))
#   expect_equal(object[["nFeature_RNA"]][1:5, ], c(47, 52, 50, 56, 53))
#   expect_equal(object[["nCount_RNA"]][75:80, ], c(228, 527, 202, 157, 150, 233))
# })

object.filtered <- CreateSeuratObject(
  counts = pbmc.test,
  min.cells = 10,
  min.features = 30
)

test_that("Filtering handled properly", {
  expect_equal(nrow(x = LayerData(object = object.filtered, layer = "counts")), 163)
  expect_equal(ncol(x = LayerData(object = object.filtered, layer = "counts")), 77)
})

#this should be moved to seurat object
# test_that("Metadata check errors correctly", {
#   pbmc.md <- pbmc_small[[]]
#   pbmc.md.norownames <- as.matrix(pbmc.md)
#   rownames(pbmc.md.norownames) <- NULL
#   expect_error(CreateSeuratObject(counts = pbmc.test, meta.data = pbmc.md.norownames),
#                "Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
# })

# Tests for NormalizeData
# --------------------------------------------------------------------------------
context("NormalizeData")
test_that("NormalizeData error handling", {
  expect_error(NormalizeData(object = object, assay = "FAKE"))
  expect_equal(
    object = LayerData(
      object = NormalizeData(
        object = object,
        normalization.method = NULL,
        verbose = FALSE
      ),
      layer = "data"
    ),
    expected = LayerData(object = object, layer = "counts")
  )
})

object <- NormalizeData(object = object, verbose = FALSE, scale.factor = 1e6)
test_that("NormalizeData scales properly", {
  expect_equal(LayerData(object = object, layer = "data")[2, 1], 9.567085, tolerance = 1e-6)
  expect_equal(LayerData(object = object, layer = "data")[161, 55], 8.415309, tolerance = 1e-6)
  expect_equal(Command(object = object, command = "NormalizeData.RNA", value = "scale.factor"), 1e6)
  expect_equal(Command(object = object, command = "NormalizeData.RNA", value = "normalization.method"), "LogNormalize")
})

normalized.data <- LogNormalize(data = GetAssayData(object = object[["RNA"]], layer = "counts"), verbose = FALSE)
test_that("LogNormalize normalizes properly", {
  expect_equal(
    as.matrix(LogNormalize(data = GetAssayData(object = object[["RNA"]], layer = "counts"), verbose = FALSE)),
    as.matrix(LogNormalize(data = as.data.frame(as.matrix(GetAssayData(object = object[["RNA"]], layer = "counts"))), verbose = FALSE))
  )
})

clr.counts <- NormalizeData(object = pbmc.test, normalization.method = "CLR", verbose = FALSE)
test_that("CLR normalization returns expected values", {
  expect_equal(dim(clr.counts), c(dim(pbmc.test)))
  expect_equal(clr.counts[2, 1], 0.5517828, tolerance = 1e-6)
  expect_equal(clr.counts[228, 76], 0.5971381, tolerance = 1e-6)
  expect_equal(clr.counts[230, 80], 0)
})

rc.counts <- NormalizeData(object = pbmc.test, normalization.method = "RC", verbose = FALSE)
test_that("Relative count normalization returns expected values", {
  expect_equal(rc.counts[2, 1], 142.8571, tolerance = 1e-6)
  expect_equal(rc.counts[228, 76], 18.97533, tolerance = 1e-6)
  expect_equal(rc.counts[230, 80], 0)
  rc.counts <- NormalizeData(object = pbmc.test, normalization.method = "RC", verbose = FALSE, scale.factor = 1e6)
  expect_equal(rc.counts[2, 1], 14285.71, tolerance = 1e-6)
})

denseMatrix <- as.matrix(pbmc.test)  # Matrix to test LogNormalize.V3Matrix and RelativeCounts methods
test_that("LogNormalize.V3Matrix computes median scale factor correctly", {
  expectedMedian <- median(colSums(denseMatrix))
  resultFromExpectedMedian <- LogNormalize.V3Matrix(data = denseMatrix, scale.factor = expectedMedian, margin = 2L, verbose = FALSE)
  resultFromScaleFactorSetToMedian <- LogNormalize.V3Matrix(data = denseMatrix, scale.factor = "median", margin = 2L, verbose = FALSE)
  expect_equal(as.matrix(resultFromExpectedMedian), as.matrix(resultFromScaleFactorSetToMedian), tolerance = 1e-6)
})

test_that("RelativeCounts computes median scale factor correctly", {
  expectedMedian <- median(colSums(denseMatrix))
  resultFromExpectedMedian <- RelativeCounts(data = denseMatrix, scale.factor = expectedMedian, verbose = FALSE)
  resultFromScaleFactorSetToMedian <- RelativeCounts(data = denseMatrix, scale.factor = "median", verbose = FALSE)
  expect_equal(as.matrix(resultFromExpectedMedian), as.matrix(resultFromScaleFactorSetToMedian), tolerance = 1e-6)
})

# Tests for v5 NormalizeData
# --------------------------------------------------------------------------------
context("v5 NormalizeData")

if(class(object[['RNA']]) == "Assay5")  {
  fake.groups <- c(rep(1, floor(ncol(pbmc.test)/2)),
                   rep(2, ncol(pbmc.test) - (floor(ncol(pbmc.test)/2))) )
  object$groups <- fake.groups
  object.split <- CreateSeuratObject(split(object[["RNA"]], f = object$groups))
  object.split <-  NormalizeData(object = object.split)

  group1 <- subset(object, groups==1)
  group1 <- NormalizeData(group1)

  test_that("Normalization is performed for each layer", {
    expect_equal(Layers(object.split),c("counts.1", "counts.2", "data.1", "data.2"))
    expect_equal(group1[['RNA']]$data, LayerData(object.split, layer="data.1"))
  })

  object.split <- NormalizeData(object = object.split, normalization.method = "CLR", verbose = FALSE)
  group1 <- NormalizeData(object = group1, normalization.method = "CLR", verbose = FALSE)
  test_that("CLR normalization works with multiple layers", {
    expect_equal(Layers(object.split),c("counts.1", "counts.2", "data.1", "data.2"))
    expect_equal(group1[['RNA']]$data, LayerData(object.split, layer="data.1"))
  })

  object.split <- NormalizeData(object = object.split, normalization.method = "RC", verbose = FALSE)
  group1 <- NormalizeData(object = group1, normalization.method = "RC", verbose = FALSE)
  test_that("RC normalization works with multiple layers", {
    expect_equal(Layers(object.split),c("counts.1", "counts.2", "data.1", "data.2"))
    expect_equal(group1[['RNA']]$data, LayerData(object.split, layer="data.1"))
  })
}



test_that("NormalizeData scales properly for BPcells", {
  # Tests for BPCells NormalizeData
  # --------------------------------------------------------------------------------

  skip_on_cran()
  library(Matrix)
  skip_if_not_installed("BPCells")
  library(BPCells)
  mat_bpcells <- t(as(t(object[['RNA']]$counts ), "IterableMatrix"))
  object[['RNAbp']] <- CreateAssay5Object(counts = mat_bpcells)

  object <- NormalizeData(object = object, verbose = FALSE, scale.factor = 1e6, assay = "RNAbp")
  object <- NormalizeData(object = object, verbose = FALSE, scale.factor = 1e6, assay = "RNA")

  expect_equal(as.matrix(object[['RNAbp']]$data), as.matrix(object[['RNA']]$data), tolerance = 1e-6)
  expect_equal(Command(object = object, command = "NormalizeData.RNAbp", value = "scale.factor"), 1e6)
  expect_equal(Command(object = object, command = "NormalizeData.RNAbp", value = "normalization.method"), "LogNormalize")
})



test_that("LogNormalize normalizes properly for BPCells", {
  skip_on_cran()
  library(Matrix)
  skip_if_not_installed("BPCells")
  library(BPCells)
  mat_bpcells <- t(as(t(object[['RNA']]$counts ), "IterableMatrix"))
  object[['RNAbp']] <- CreateAssay5Object(counts = mat_bpcells)

  object <- NormalizeData(object = object, verbose = FALSE, scale.factor = 1e6, assay = "RNAbp")
  object <- NormalizeData(object = object, verbose = FALSE, scale.factor = 1e6, assay = "RNA")

  normalized.data.bp <- LogNormalize(data = GetAssayData(object = object[["RNAbp"]], layer = "counts"), verbose = FALSE)
  normalized.data <- LogNormalize(data = GetAssayData(object = object[["RNA"]], layer = "counts"), verbose = FALSE)

  expect_equal(
    as.matrix(normalized.data.bp),
    as.matrix(normalized.data),
    tolerance = 1e-6
  )
})

test_that("LogNormalize.IterableMatrix computes median scale factor correctly", {
  skip_on_cran()
  library(Matrix)
  skip_if_not_installed("BPCells")
  library(BPCells)
  mat_bpcells <- t(as(t(object[['RNA']]$counts ), "IterableMatrix"))
  expectedMedian <- median(colSums(mat_bpcells))
  resultFromExpectedMedian <- LogNormalize.IterableMatrix(data = mat_bpcells, scale.factor = expectedMedian, margin = 2L, verbose = FALSE)
  resultFromScaleFactorSetToMedian <- LogNormalize.IterableMatrix(data = mat_bpcells, scale.factor = "median", margin = 2L, verbose = FALSE)
  expect_equal(as.matrix(resultFromExpectedMedian), as.matrix(resultFromScaleFactorSetToMedian), tolerance = 1e-6)
})

denseMatrix <- as.matrix(pbmc.test)  # Matrix to test LogNormalize.default when scale.factor is set to "median"
test_that("LogNormalize.default computes median scale factor correctly for both margin values", {
  expectedMedianForMargin1L <- median(rowSums(denseMatrix))
  expectedMedianForMargin2L <- median(colSums(denseMatrix))
  
  resultFromExpectedMedianForMargin1L <- LogNormalize.default(data = denseMatrix, scale.factor = expectedMedianForMargin1L, margin = 1L, verbose = FALSE)
  resultFromExpectedMedianForMargin2L <- LogNormalize.default(data = denseMatrix, scale.factor = expectedMedianForMargin2L, margin = 2L, verbose = FALSE)
  
  resultsFromScaleFactorSetToMedianForMargin1L <- LogNormalize.default(data = denseMatrix, scale.factor = "median", margin = 1L, verbose = FALSE)#if the normalization is across rows (genes)
  resultsFromScaleFactorSetToMedianForMargin2L <- LogNormalize.default(data = denseMatrix, scale.factor = "median", margin = 2L, verbose = FALSE)#if the normalization is across columns (cells)
  
  expect_equal(as.matrix(resultFromExpectedMedianForMargin1L), as.matrix(resultsFromScaleFactorSetToMedianForMargin1L), tolerance = 1e-6)
  expect_equal(as.matrix(resultFromExpectedMedianForMargin2L), as.matrix(resultsFromScaleFactorSetToMedianForMargin2L), tolerance = 1e-6)
})

theSparseMatrix <- as.sparse(denseMatrix) # Sparse Matrix to test .SparseNormalize computes median scale factor correctly
test_that("LogNormalize.default computes median scale factor correctly for both margin values", {
  expectedMedian <- median(colSums(theSparseMatrix))
  
  resultFromExpectedMedian <- .SparseNormalize(data = theSparseMatrix, scale.factor = expectedMedian, verbose = FALSE)
  resultsFromScaleFactorSetToMedian <- .SparseNormalize(data = theSparseMatrix, scale.factor = "median", verbose = FALSE)
  
  expect_equal(resultFromExpectedMedian, resultsFromScaleFactorSetToMedian, tolerance = 1e-6)
})


# Tests for ScaleData
# --------------------------------------------------------------------------------
context("ScaleData")
object <- ScaleData(object, verbose = FALSE)
test_that("ScaleData returns expected values when input is a sparse matrix", {
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[1, 1], -0.4148587, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[75, 25], -0.2562305, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[162, 59], -0.4363939, tolerance = 1e-6)
})

new.data <- as.matrix(GetAssayData(object = object[["RNA"]], layer = "data"))
new.data[1, ] <- rep(x = 0, times = ncol(x = new.data))
object2 <- object

object2 <- SetAssayData(
  object = object,
  assay = "RNA",
  slot = "data",
  new.data = new.data
)
object2 <- ScaleData(object = object2, verbose = FALSE)

object <- ScaleData(object = object, verbose = FALSE)
test_that("ScaleData returns expected values when input is not sparse", {
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[75, 25], -0.2562305, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[162, 59], -0.4363939, tolerance = 1e-6)
})

test_that("ScaleData handles zero variance features properly", {
  expect_equal(GetAssayData(object = object2[["RNA"]], layer = "scale.data")[1, 1], 0)
  expect_equal(GetAssayData(object = object2[["RNA"]], layer = "scale.data")[1, 80], 0)
})

ng1 <- rep(x = "g1", times = round(x = ncol(x = object) / 2))
object$group <- c(ng1, rep(x = "g2", times = ncol(x = object) - length(x = ng1)))
g1 <- subset(x = object, group == "g1")
g1 <- ScaleData(object = g1, features = rownames(x = g1), verbose = FALSE)
g2 <- subset(x = object, group == "g2")
g2 <- ScaleData(object = g2, features = rownames(x = g2), verbose = FALSE)
object <- ScaleData(object = object, features = rownames(x = object), verbose = FALSE, split.by = "group")

#move to SeuratObject
# test_that("split.by option works", {
#   expect_equal(GetAssayData(object = object, layer = "scale.data")[, Cells(x = g1)],
#                GetAssayData(object = g1, layer = "scale.data"))
#   expect_equal(GetAssayData(object = object, layer = "scale.data")[, Cells(x = g2)],
#                GetAssayData(object = g2, layer = "scale.data"))
# })

g1 <- ScaleData(object = g1, features = rownames(x = g1), vars.to.regress = "nCount_RNA", verbose = FALSE)
g2 <- ScaleData(object = g2, features = rownames(x = g2), vars.to.regress = "nCount_RNA", verbose = FALSE)
object <- ScaleData(object = object, features = rownames(x = object), verbose = FALSE, split.by = "group", vars.to.regress = "nCount_RNA")
test_that("split.by option works with regression", {
  expect_equal(LayerData(object = object, layer = "scale.data")[, Cells(x = g1)],
               LayerData(object = g1, layer = "scale.data"))
  expect_equal(LayerData(object = object, layer = "scale.data")[, Cells(x = g2)],
               LayerData(object = g2, layer = "scale.data"))
})


# Tests for various regression techniques
context("Regression")

suppressWarnings({
  object <- ScaleData(
  object = object,
  vars.to.regress = "nCount_RNA",
  features = rownames(x = object)[1:10],
  verbose = FALSE,
  model.use = "linear")
  })

test_that("Linear regression works as expected", {
  expect_equal(dim(x = GetAssayData(object = object[["RNA"]], layer = "scale.data")), c(10, 80))
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[1, 1], -0.6436435, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[5, 25], -0.09035383, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[10, 80], -0.2723782, tolerance = 1e-6)
})

object <- ScaleData(
  object,
  vars.to.regress = "nCount_RNA",
  features = rownames(x = object)[1:10],
  verbose = FALSE,
  model.use = "negbinom")

test_that("Negative binomial regression works as expected", {
  expect_equal(dim(x = GetAssayData(object = object[["RNA"]], layer = "scale.data")), c(10, 80))
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[1, 1], -0.5888811, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[5, 25], -0.2553394, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[10, 80], -0.1921429, tolerance = 1e-6)
})

test_that("Regression error handling checks out", {
  expect_error(ScaleData(object, vars.to.regress = "nCount_RNA", model.use = "not.a.model", verbose = FALSE))
})

object <- ScaleData(
  object,
  vars.to.regress = "nCount_RNA",
  features = rownames(x = object)[1:10],
  verbose = FALSE,
  model.use = "poisson")

test_that("Poisson regression works as expected", {
  expect_equal(dim(x = GetAssayData(object = object[["RNA"]], layer = "scale.data")), c(10, 80))
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[1, 1], -1.011717, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[5, 25], 0.05575307, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], layer = "scale.data")[10, 80], -0.1662119, tolerance = 1e-6)
})


#Tests for SampleUMI
#--------------------------------------------------------------------------------
context("SampleUMI")

downsampled.umis <- SampleUMI(
  data = LayerData(object = object, layer = "counts"),
  max.umi = 100,
  verbose = FALSE
)
downsampled.umis.p.cell <- SampleUMI(
  data = LayerData(object = object, layer = "counts"),
  max.umi = seq(50, 1640, 20),
  verbose = FALSE,
  upsample = TRUE
)
test_that("SampleUMI gives reasonable downsampled/upsampled UMI counts", {
  expect_true(!any(colSums(x = downsampled.umis) < 30, colSums(x = downsampled.umis) > 120))
  expect_error(SampleUMI(data = LayerData(object = object, layer = "counts"), max.umi = rep(1, 5)))
  expect_true(!is.unsorted(x = colSums(x = downsampled.umis.p.cell)))
  expect_error(SampleUMI(
    data = LayerData(object = object, layer = "counts"),
    max.umi = seq(50, 900, 10),
    verbose = FALSE,
    upsample = TRUE
  ))
})

# Tests for FindVariableFeatures
# --------------------------------------------------------------------------------
context("FindVariableFeatures")

object <- FindVariableFeatures(object = object, selection.method = "mean.var.plot", verbose = FALSE)
test_that("mean.var.plot selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PTGDR", "SATB1", "ZNF330", "S100B"))
  expect_equal(length(x = VariableFeatures(object = object)), 20)
  hvf_info <- HVFInfo(object = object[["RNA"]], method = 'mvp')
  expect_equal(hvf_info[[grep("mean$", colnames(hvf_info), value = TRUE)]][1:2], c(8.328927, 8.444462), tolerance = 1e-6)
  expect_equal(hvf_info[[grep("dispersion$", colnames(hvf_info), value = TRUE)]][1:2], c(10.552507, 10.088223), tolerance = 1e-6)
  expect_equal(as.numeric(hvf_info[[grep("dispersion.scaled$", colnames(hvf_info), value = TRUE)]][1:2]), c(0.1113214, -0.1332181523), tolerance = 1e-6)
})

object <- FindVariableFeatures(object, selection.method = "dispersion", verbose = FALSE)
test_that("dispersion selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PCMT1", "PPBP", "LYAR", "VDAC3"))
  expect_equal(length(x = VariableFeatures(object = object)), 230)
  hvf_info <- HVFInfo(object = object[["RNA"]], method = 'mvp')
  expect_equal(hvf_info[[grep("mean$", colnames(hvf_info), value = TRUE)]][1:2], c(8.328927, 8.444462), tolerance = 1e-6)
  expect_equal(hvf_info[[grep("dispersion$", colnames(hvf_info), value = TRUE)]][1:2], c(10.552507, 10.088223), tolerance = 1e-6)
  expect_equal(as.numeric(hvf_info[[grep("dispersion.scaled$", colnames(hvf_info), value = TRUE)]][1:2]), c(0.1113214, -0.1332181523), tolerance = 1e-6)
  expect_true(!is.unsorted(rev(hvf_info[VariableFeatures(object = object), "dispersion"])))
})

object <- FindVariableFeatures(object, selection.method = "vst", verbose = FALSE)
test_that("vst selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PPBP", "IGLL5", "VDAC3", "CD1C"))
  expect_equal(length(x = VariableFeatures(object = object)), 230)
  hvf_info <- HVFInfo(object = object[["RNA"]], method = 'vst')
  expect_equal(hvf_info[[grep("variance$", colnames(hvf_info), value = TRUE)]][1:2], c(1.0251582, 1.2810127), tolerance = 1e-6)
  expect_equal(hvf_info[[grep("variance.standardized$", colnames(hvf_info), value = TRUE)]][1:2], c(0.8983463, 0.4731134), tolerance = 1e-6)
  expect_true(!is.unsorted(rev(hvf_info[VariableFeatures(object = object), grep("variance.standardized$", colnames(hvf_info))])))
})

#object <- FindVariableFeatures(object, assay = "RNAbp")
#this breaks currently

# Tests for internal functions
# ------------------------------------------------------------------------------
norm.fxn <- function(x) {x / mean(x)}
test_that("CustomNormalize works as expected", {
  expect_equal(
    CustomNormalize(data = pbmc.test, custom_function = norm.fxn, margin = 2),
    apply(X = pbmc.test, MARGIN = 2, FUN = norm.fxn)
  )
  expect_equal(
    CustomNormalize(data = as.matrix(pbmc.test), custom_function = norm.fxn, margin = 2),
    apply(X = pbmc.test, MARGIN = 2, FUN = norm.fxn)
  )
  expect_equal(
    CustomNormalize(data = as.data.frame(as.matrix(pbmc.test)), custom_function = norm.fxn, margin = 2),
    apply(X = pbmc.test, MARGIN = 2, FUN = norm.fxn)
  )
  expect_equal(
    CustomNormalize(data = pbmc.test, custom_function = norm.fxn, margin = 1),
    t(apply(X = pbmc.test, MARGIN = 1, FUN = norm.fxn))
  )
  expect_error(CustomNormalize(data = pbmc.test, custom_function = norm.fxn, margin = 10))
})

# Tests for SCTransform
# --------------------------------------------------------------------------------
context("SCTransform")
object <- suppressWarnings(SCTransform(object = object, verbose = FALSE, vst.flavor = "v1",  seed.use = 1448145))

test_that("SCTransform v1 works as expected", {
  expect_true("SCT" %in% names(object))
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "scale.data"))[1]), 11.40288448)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "scale.data"))[5]), 0)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "data"))[1]), 57.7295742, tolerance = 1e-6)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "data"))[5]), 11.74403719, tolerance = 1e-6)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "counts"))[1]), 129)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "counts"))[5]), 28)
  expect_equal(length(VariableFeatures(object[["SCT"]])), 220)
  fa <- SCTResults(object = object, assay = "SCT", slot = "feature.attributes")
  expect_equal(fa["MS4A1", "detection_rate"], 0.15)
  expect_equal(fa["MS4A1", "gmean"], 0.2027364, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "variance"], 1.025158, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "residual_mean"], 0.2362887, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "residual_variance"], 2.875761, tolerance = 1e-6)
})

suppressWarnings(RNGversion(vstr = "3.5.0"))
object <- suppressWarnings(SCTransform(object = object, vst.flavor = "v1", ncells = 80, verbose = FALSE, seed.use =  42))
test_that("SCTransform ncells param works", {
  expect_true("SCT" %in% names(object))
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "scale.data"))[1]), 11.40288, tolerance = 1e-6)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "scale.data"))[5]), 0)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "data"))[1]), 57.72957, tolerance = 1e-6)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "data"))[5]), 11.74404, tolerance = 1e-6)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "counts"))[1]), 129)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "counts"))[5]), 28)
  expect_equal(length(VariableFeatures(object[["SCT"]])), 220)
  fa <- SCTResults(object = object, assay = "SCT", slot = "feature.attributes")
  expect_equal(fa["MS4A1", "detection_rate"], 0.15)
  expect_equal(fa["MS4A1", "gmean"], 0.2027364, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "variance"], 1.025158, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "residual_mean"], 0.2362887, tolerance = 1e-3)
  expect_equal(fa["MS4A1", "residual_variance"], 2.875761, tolerance = 1e-3)
})

suppressWarnings(object[["SCT_SAVE"]] <- object[["SCT"]])
object[["SCT"]] <- suppressWarnings({SetAssayData(object = object[["SCT"]], slot = "scale.data", new.data = GetAssayData(object = object[["SCT"]], layer = "scale.data")[1:100, ])})
object <- GetResidual(object = object, features = rownames(x = object), verbose = FALSE)
test_that("GetResidual works", {
  expect_equal(dim(GetAssayData(object = object[["SCT"]], layer = "scale.data")), c(220, 80))
  expect_equal(
    GetAssayData(object = object[["SCT"]], layer = "scale.data"),
    GetAssayData(object = object[["SCT_SAVE"]], layer = "scale.data")
  )
  expect_warning(GetResidual(object, features = "asd"))
})


test_that("SCTransform v2 works as expected", {
  skip_on_cran()
  skip_if_not_installed("glmGamPoi")
  object <- suppressWarnings(SCTransform(object = object, verbose = FALSE, vst.flavor = "v2",  seed.use = 1448145))

  expect_true("SCT" %in% names(object))
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "scale.data"))[1]), 24.5813, tolerance = 1e-4)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "scale.data"))[5]), 0)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "data"))[1]), 58.65829, tolerance = 1e-6)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "data"))[5]), 13.75449, tolerance = 1e-6)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], layer = "counts"))[1]), 141)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], layer = "counts"))[5]), 40)
  expect_equal(length(VariableFeatures(object[["SCT"]])), 220)
  fa <- SCTResults(object = object, assay = "SCT", slot = "feature.attributes")
  expect_equal(fa["MS4A1", "detection_rate"], 0.15)
  expect_equal(fa["MS4A1", "gmean"], 0.2027364, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "variance"], 1.025158, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "residual_mean"], 0.2763993, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "residual_variance"], 3.023062, tolerance = 1e-6)
  expect_equal(fa["FCER2", "theta"], Inf)
})

test_that("SCTransform `clip.range` param works as expected", {
  # make a copy of the testing data
  test.data <- object
  # override defaults for ease of testing
  clip.min <- -0.1
  clip.max <- 0.1

  # for some reason, the clipping seems to be a little fuzzy at the upper end,
  # since this is expected behaviour we'll need to accomodate the difference
  clip.max.tolerance <- 0.1

  test.result <- suppressWarnings(
      SCTransform(
      test.data,
      clip.range = c(clip.min, clip.max),
    )
  )
  scale.data <- LayerData(test.result[["SCT"]], layer = "scale.data")
  expect_true(min(scale.data) >= clip.min)
  expect_true(max(scale.data) <= (clip.max + clip.max.tolerance))

  # when `ncells` is less than the size of the dataset the residuals will get 
  # re-clipped in batches, make sure this clipping is done correctly as well
  test.result <- suppressWarnings(
    SCTransform(
      test.data,
      clip.range = c(clip.min, clip.max),
      ncells = 40
    )
  )
  scale.data <- LayerData(test.result[["SCT"]], layer = "scale.data")
  expect_true(min(scale.data) >= clip.min)
  expect_true(max(scale.data) <= (clip.max + clip.max.tolerance))
})

test_that("SCTransform `vars.to.regress` param works as expected", {
  # make a copy of the testing data
  test.data <- object
  # add a fake mitochondrial gene to the counts matrix
  counts <- LayerData(test.data, assay = "RNA", layer = "counts")
  counts <- rbind(counts, 5)
  rownames(counts)[nrow(counts)] <- "MT-TEST"
  # use the fake feature to populate a new meta.data column
  test.data[[ "percent.mt" ]] <- PercentageFeatureSet(
    test.data,
    pattern="^MT-"
  )

  # make sure that `ncells` is smaller than the datset being transformed 
  # so tha the regression model is trained on a subset of the data - make sure 
  # the regression is applied to the entire dataset
  left <- suppressWarnings(
      SCTransform(
      test.data,
      vars.to.regress = NULL,
      ncells = ncol(test.data) / 2,
      verbose = FALSE
    )
  )
  right <- suppressWarnings(
      SCTransform(
      test.data,
      vars.to.regress = "percent.mt",
      ncells = ncol(test.data) / 2,
      verbose = FALSE
    )
  )
  expect_false(identical(left[["SCT"]]$scale.data, right[["SCT"]]$scale.data))

  # if the `assay` points to an `Assay5` instance the regression is handled
  # using separate logic
  test.data[["RNAv5"]] <- CreateAssay5Object(
    counts = LayerData(
      test.data,
      assay = "RNA",
      layer = "counts"
    )
  )
  left <- suppressWarnings(
    SCTransform(
      test.data,
      assay = "RNAv5",
      vars.to.regress = NULL,
      ncells = ncol(test.data) / 2,
      verbose = FALSE
    )
  )
  right <- suppressWarnings(
    SCTransform(
      test.data,
      assay = "RNAv5",
      vars.to.regress = "percent.mt",
      ncells = ncol(test.data) / 2,
      verbose = FALSE
    )
  )
  expect_false(identical(left[["SCT"]]$scale.data, right[["SCT"]]$scale.data))
})

test_that("SCTransform is equivalent for BPcells ", {
  skip_on_cran()
  skip_on_cran()
  skip_if_not_installed("glmGamPoi")

  library(Matrix)
  skip_if_not_installed("BPCells")
  library(BPCells)
  mat_bpcells <- t(as(t(object[['RNA']]$counts ), "IterableMatrix"))
  object[['RNAbp']] <- CreateAssay5Object(counts = mat_bpcells)
  object <- suppressWarnings(SCTransform(object = object, assay = "RNA", new.assay.name = "SCT",
                                         verbose = FALSE, vst.flavor = "v2",  seed.use = 1448145))

  object <- suppressWarnings(SCTransform(object = object, assay = "RNAbp", new.assay.name = "SCTbp",
                                         verbose = FALSE, vst.flavor = "v2",  seed.use = 1448145))

  expect_equal(as.matrix(LayerData(object = object[["SCT"]], layer = "data")),
               as.matrix(LayerData(object = object[["SCTbp"]], layer = "data")),
               tolerance = 1e-6)
})
