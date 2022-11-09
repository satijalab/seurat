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
  expect_equal(nrow(x = GetAssayData(object = object.filtered, slot = "counts")), 163)
  expect_equal(ncol(x = GetAssayData(object = object.filtered, slot = "counts")), 77)
})

test_that("Metadata check errors correctly", {
  pbmc.md <- pbmc_small[[]]
  pbmc.md.norownames <- as.matrix(pbmc.md)
  rownames(pbmc.md.norownames) <- NULL
  expect_error(CreateSeuratObject(counts = pbmc.test, meta.data = pbmc.md.norownames),
               "Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
})

# Tests for NormalizeData
# --------------------------------------------------------------------------------
context("NormalizeData")
test_that("NormalizeData error handling", {
  expect_error(NormalizeData(object = object, assay = "FAKE"))
  expect_equal(
    object = GetAssayData(
      object = NormalizeData(
        object = object,
        normalization.method = NULL,
        verbose = FALSE
      ),
      slot = "data"
    ),
    expected = GetAssayData(object = object, slot = "counts")
  )
})

object <- NormalizeData(object = object, verbose = FALSE, scale.factor = 1e6)
test_that("NormalizeData scales properly", {
  expect_equal(GetAssayData(object = object, slot = "data")[2, 1], 9.567085, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object, slot = "data")[161, 55], 8.415309, tolerance = 1e-6)
  expect_equal(Command(object = object, command = "NormalizeData.RNA", value = "scale.factor"), 1e6)
  expect_equal(Command(object = object, command = "NormalizeData.RNA", value = "normalization.method"), "LogNormalize")
})

normalized.data <- LogNormalize(data = GetAssayData(object = object[["RNA"]], slot = "counts"), verbose = FALSE)
test_that("LogNormalize normalizes properly", {
  expect_equal(
    LogNormalize(data = GetAssayData(object = object[["RNA"]], slot = "counts"), verbose = FALSE),
    LogNormalize(data = as.data.frame(as.matrix(GetAssayData(object = object[["RNA"]], slot = "counts"))), verbose = FALSE)
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

# Tests for ScaleData
# --------------------------------------------------------------------------------
context("ScaleData")
object <- ScaleData(object, verbose = FALSE)
test_that("ScaleData returns expected values when input is a sparse matrix", {
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[1, 1], -0.4148587, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[75, 25], -0.2562305, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[162, 59], -0.4363939, tolerance = 1e-6)
})

new.data <- as.matrix(GetAssayData(object = object[["RNA"]], slot = "data"))
new.data[1, ] <- rep(x = 0, times = ncol(x = new.data))
object2 <- object

object2[["RNA"]] <- SetAssayData(
  object = object[["RNA"]],
  slot = "data",
  new.data = new.data
)
object2 <- ScaleData(object = object2, verbose = FALSE)

object <- ScaleData(object = object, verbose = FALSE)
test_that("ScaleData returns expected values when input is not sparse", {
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[75, 25], -0.2562305, tolerance = 1e-6)
  expect_equal(GetAssayData(object = object[["RNA"]], slot = "scale.data")[162, 59], -0.4363939, tolerance = 1e-6)
})

test_that("ScaleData handles zero variance features properly", {
  expect_equal(GetAssayData(object = object2[["RNA"]], slot = "scale.data")[1, 1], 0)
  expect_equal(GetAssayData(object = object2[["RNA"]], slot = "scale.data")[1, 80], 0)
})

ng1 <- rep(x = "g1", times = round(x = ncol(x = object) / 2))
object$group <- c(ng1, rep(x = "g2", times = ncol(x = object) - length(x = ng1)))
g1 <- subset(x = object, group == "g1")
g1 <- ScaleData(object = g1, features = rownames(x = g1), verbose = FALSE)
g2 <- subset(x = object, group == "g2")
g2 <- ScaleData(object = g2, features = rownames(x = g2), verbose = FALSE)
object <- ScaleData(object = object, features = rownames(x = object), verbose = FALSE, split.by = "group")

test_that("split.by option works", {
  expect_equal(GetAssayData(object = object, slot = "scale.data")[, Cells(x = g1)],
               GetAssayData(object = g1, slot = "scale.data"))
  expect_equal(GetAssayData(object = object, slot = "scale.data")[, Cells(x = g2)],
               GetAssayData(object = g2, slot = "scale.data"))
})

g1 <- ScaleData(object = g1, features = rownames(x = g1), vars.to.regress = "nCount_RNA", verbose = FALSE)
g2 <- ScaleData(object = g2, features = rownames(x = g2), vars.to.regress = "nCount_RNA", verbose = FALSE)
object <- ScaleData(object = object, features = rownames(x = object), verbose = FALSE, split.by = "group", vars.to.regress = "nCount_RNA")
test_that("split.by option works with regression", {
  expect_equal(GetAssayData(object = object, slot = "scale.data")[, Cells(x = g1)],
               GetAssayData(object = g1, slot = "scale.data"))
  expect_equal(GetAssayData(object = object, slot = "scale.data")[, Cells(x = g2)],
               GetAssayData(object = g2, slot = "scale.data"))
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
  expect_error(ScaleData(object, vars.to.regress = "nCount_RNA", model.use = "not.a.model", verbose = FALSE))
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
  max.umi = seq(50, 1640, 20),
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

# Tests for FindVariableFeatures
# --------------------------------------------------------------------------------
context("FindVariableFeatures")

object <- FindVariableFeatures(object = object, selection.method = "mean.var.plot", verbose = FALSE)
test_that("mean.var.plot selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PTGDR", "SATB1", "ZNF330", "S100B"))
  expect_equal(length(x = VariableFeatures(object = object)), 20)
  expect_equal(HVFInfo(object = object[["RNA"]], selection.method = 'mvp')$mean[1:2], c(8.328927, 8.444462), tolerance = 1e-6)
  expect_equal(HVFInfo(object = object[["RNA"]], selection.method = 'mvp')$dispersion[1:2], c(10.552507, 10.088223), tolerance = 1e-6)
  expect_equal(as.numeric(HVFInfo(object = object[["RNA"]], selection.method = 'mvp')$dispersion.scaled[1:2]), c(0.1113214, -0.1332181523), tolerance = 1e-6)
})

object <- FindVariableFeatures(object, selection.method = "dispersion", verbose = FALSE)
test_that("dispersion selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PCMT1", "PPBP", "LYAR", "VDAC3"))
  expect_equal(length(x = VariableFeatures(object = object)), 230)
  expect_equal(HVFInfo(object = object[["RNA"]], selection.method = 'mvp')$mean[1:2], c(8.328927, 8.444462), tolerance = 1e-6)
  expect_equal(HVFInfo(object = object[["RNA"]], selection.method = 'mvp')$dispersion[1:2], c(10.552507, 10.088223), tolerance = 1e-6)
  expect_equal(as.numeric(HVFInfo(object = object[["RNA"]], selection.method = 'mvp')$dispersion.scaled[1:2]), c(0.1113214, -0.1332181523), tolerance = 1e-6)
  expect_true(!is.unsorted(rev(HVFInfo(object = object[["RNA"]], selection.method = 'mvp')[VariableFeatures(object = object), "dispersion"])))
})

object <- FindVariableFeatures(object, selection.method = "vst", verbose = FALSE)
test_that("vst selection option returns expected values", {
  expect_equal(VariableFeatures(object = object)[1:4], c("PPBP", "IGLL5", "VDAC3", "CD1C"))
  expect_equal(length(x = VariableFeatures(object = object)), 230)
  expect_equal(unname(object[["RNA"]][["vst.variance", drop = TRUE]][1:2]), c(1.0251582, 1.2810127), tolerance = 1e-6)
  expect_equal(unname(object[["RNA"]][["vst.variance.expected", drop = TRUE]][1:2]), c(1.1411616, 2.7076228), tolerance = 1e-6)
  expect_equal(unname(object[["RNA"]][["vst.variance.standardized", drop = TRUE]][1:2]), c(0.8983463, 0.4731134), tolerance = 1e-6)
  expect_true(!is.unsorted(rev(object[["RNA"]][["vst.variance.standardized", drop = TRUE]][VariableFeatures(object = object)])))
})

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
# ------------------------------------------------------------------------------
context("SCTransform")
object <- suppressWarnings(SCTransform(object = object, verbose = FALSE))

test_that("SCTransform wrapper works as expected", {
  expect_true("SCT" %in% names(object))
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], slot = "scale.data"))[1]), 11.40288448)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], slot = "scale.data"))[5]), 0)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], slot = "data"))[1]), 57.7295742, tolerance = 1e-6)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], slot = "data"))[5]), 11.74403719, tolerance = 1e-6)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], slot = "counts"))[1]), 129)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], slot = "counts"))[5]), 28)
  expect_equal(length(VariableFeatures(object[["SCT"]])), 220)
  fa <- SCTResults(object = object, assay = "SCT", slot = "feature.attributes")
  expect_equal(fa["MS4A1", "detection_rate"], 0.15)
  expect_equal(fa["MS4A1", "gmean"], 0.2027364, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "variance"], 1.025158, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "residual_mean"], 0.2362887, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "residual_variance"], 2.875761, tolerance = 1e-6)
})

suppressWarnings(RNGversion(vstr = "3.5.0"))
object <- suppressWarnings(SCTransform(object = object, ncells = 40, verbose = FALSE, seed.use =  42))
test_that("SCTransform ncells param works", {
  expect_true("SCT" %in% names(object))
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], slot = "scale.data"))[1]), 12.02126, tolerance = 1e6)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], slot = "scale.data"))[5]), 0)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], slot = "data"))[1]), 60.65299, tolerance = 1e-6)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], slot = "data"))[5]), 11.74404, tolerance = 1e-6)
  expect_equal(as.numeric(colSums(GetAssayData(object = object[["SCT"]], slot = "counts"))[1]), 136)
  expect_equal(as.numeric(rowSums(GetAssayData(object = object[["SCT"]], slot = "counts"))[5]), 28)
  expect_equal(length(VariableFeatures(object[["SCT"]])), 220)
  fa <- SCTResults(object = object, assay = "SCT", slot = "feature.attributes")
  expect_equal(fa["MS4A1", "detection_rate"], 0.15)
  expect_equal(fa["MS4A1", "gmean"], 0.2027364, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "variance"], 1.025158, tolerance = 1e-6)
  expect_equal(fa["MS4A1", "residual_mean"], 0.2829672, tolerance = 1e-3)
  expect_equal(fa["MS4A1", "residual_variance"], 3.674079, tolerance = 1e-3)
})

suppressWarnings(object[["SCT_SAVE"]] <- object[["SCT"]])
object[["SCT"]] <- SetAssayData(object = object[["SCT"]], slot = "scale.data", new.data = GetAssayData(object = object[["SCT"]], slot = "scale.data")[1:100, ])
object <- GetResidual(object = object, features = rownames(x = object), verbose = FALSE)
test_that("GetResidual works", {
  expect_equal(dim(GetAssayData(object = object[["SCT"]], slot = "scale.data")), c(220, 80))
  expect_equal(
    GetAssayData(object = object[["SCT"]], slot = "scale.data"),
    GetAssayData(object = object[["SCT_SAVE"]], slot = "scale.data")
  )
  expect_warning(GetResidual(object, features = "asd"))
})

