# Tests for integration/transfer related fxns
set.seed(42)
pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))

# Setup test objects
ref <- pbmc_small
query <- CreateSeuratObject(
  counts = as.sparse(
    GetAssayData(
      object = pbmc_small[['RNA']],
      layer = "counts") + rpois(n = ncol(pbmc_small),
      lambda = 1
    )
  )
)
query <- NormalizeData(object = query, verbose = FALSE)
query <- FindVariableFeatures(object = query, verbose = FALSE, nfeatures = 100)
ref <- FindVariableFeatures(object = ref, verbose = FALSE, nfeatures = 100)

# Tests for FindTransferAnchors
# ------------------------------------------------------------------------------
context("FindTransferAnchors")

test_that("FindTransferAnchors defaults work", {
  anchors <- FindTransferAnchors(reference = ref, query = query, k.filter = 50)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]], layer ="data")[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")[1, 3], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 0.4840944592, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], 0.2103563963, tolerance = 1e-7)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], 0.05175486778, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], 0.2103563963, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(128, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(5, 5, 0.08361970218), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

test_that("FindTransferAnchors catches bad input", {
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.assay = "BAD", k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, query.assay = "BAD", k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, normalization.method = "BAD", k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, reduction = "BAD", k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, npcs = NULL, k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, npcs = NULL, reference.reduction = "BAD", k.filter = 50))
  expect_error(suppressWarngings(FindTransferAnchors(reference = ref, query = query, dims = 1:100, k.filter = 50)))
  expect_error(suppressWarnings(FindTransferAnchors(reference = ref, query = query, dims = 1:100, project.query = TRUE, k.filter = 50)))
  expect_error(FindTransferAnchors(reference = ref, query = query, k.anchor = 80, k.filter = 50))
  expect_warning(FindTransferAnchors(reference = ref, query = query, k.filter = 81))
  expect_error(FindTransferAnchors(reference = ref, query = query, k.filter = 50, k.score = 80))
  expect_error(suppressWarnings(FindTransferAnchors(reference = ref, query = query, k.filter = 50, features = "BAD")))
  expect_error(FindTransferAnchors(reference = ref, query = query, k.filter = 50, reduction = "cca", project.query = TRUE))
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.reduction = "BAD", k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.reduction = "BAD", project.query = TRUE, k.filter = 50))
})

ref <- ScaleData(ref, verbose = FALSE)
ref <- suppressWarnings(RunPCA(ref, npcs = 30, verbose = FALSE))

test_that("FindTransferAnchors allows reference.reduction to be precomputed", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, k.filter = 50, reference.reduction = "pca")
  expect_error(FindTransferAnchors(reference = ref, query = query, k.filter = 50, reference.reduction = "pca", reduction = "cca"))
  expect_error(FindTransferAnchors(reference = ref, query = query, k.filter = 50, reference.reduction = "pca", project.query = TRUE))
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]])[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")[1, 3], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 0.4840944592, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], 0.2103563963, tolerance = 1e-7)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], 0.05175486778, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], 0.2103563963, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(128, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(5, 5, 0.08361970218), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

test_that("FindTransferAnchors with cca defaults work", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, reduction = "cca", k.filter = 50)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("cca", "cca.l2"))
  expect_equal(GetAssayData(co[["RNA"]])["PPBP", 3], 0)
  expect_equal(GetAssayData(co[["RNA"]])["PPBP", 1], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")["PPBP", 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")["PPBP", 1], 0)
  expect_equal(dim(co[['cca']]), c(160, 30))
  expect_equal(Embeddings(co[['cca']])[1, 1], 0.04611130861, tolerance = 1e-7)
  expect_equal(Loadings(co[['cca']], projected = T)["PPBP", 1], 12.32379661, tolerance = 1e-7)
  expect_equal(dim(co[['cca.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['cca.l2']])[1, 1], 0.06244169641, tolerance = 1e-7)
  expect_equal(Loadings(co[['cca.l2']], projected = T)["PPBP", 1], 12.32379661, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(324, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.8211091234), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

test_that("FindTransferAnchors with project.query defaults work", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, project.query = TRUE, k.filter = 50)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]], layer = "data")["PPBP", 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "data")["PPBP", 1], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")["PPBP", 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")["PPBP", 1], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 1.577959404, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)["PPBP", 1], 0.1145472305, tolerance = 1e-7)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], 0.1358602536, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)["PPBP", 1], 0.1145472305, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(208, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 10, 0.4984040128), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "GZMA")
  expect_equal(anchors@neighbors, list())
})

query <- ScaleData(query, verbose = FALSE)
query <- suppressWarnings(RunPCA(query, npcs = 30, verbose = FALSE))

test_that("FindTransferAnchors with project.query and reference.reduction works", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, k.filter = 50, reference.reduction = "pca", project.query = TRUE)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]])["PPBP", 3], 0)
  expect_equal(GetAssayData(co[["RNA"]])["PPBP", 1], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")["PPBP", 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")["PPBP", 1], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 1.577959404, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)["PPBP", 1], 0.1145472305, tolerance = 1e-7)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], 0.1358602536, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)["PPBP", 1], 0.1145472305, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(208, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 10, 0.4984040128), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "GZMA")
  expect_equal(anchors@neighbors, list())
})

ref <- FindNeighbors(object = ref, reduction = "pca", dims = 1:30, return.neighbor = TRUE, k.param = 31, verbose = FALSE, l2.norm = TRUE, nn.method = "annoy")
test_that("FindTransferAnchors with reference.neighbors precomputed works", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, reference.neighbors = "RNA.nn", k.filter = 50)
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.neighbors = "BAD", k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.neighbors = "RNA.nn", k.filter = 50, k.score = 31))
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.neighbors = "RNA.nn", k.filter = 50, k.anchor = 31))
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]])[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")[1, 3], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 0.4840944592, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], 0.2103563963, tolerance = 1e-7)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], 0.05175486778, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], 0.2103563963, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(128, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(5, 5, 0.08361970218), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

test_that("FindTransferAnchors with no l2 works", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, l2.norm = FALSE, k.filter = 50)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("pcaproject"))
  expect_equal(GetAssayData(co[["RNA"]])[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")[1, 3], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 0.4840944592, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], 0.2103563963, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(115, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(5, 5, 0.2950654582), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

# SCTransform tests V1
query <- suppressWarnings(SCTransform(object = query, verbose = FALSE,vst.flavor = 'v1'))
ref <- suppressWarnings(SCTransform(object = ref, verbose = FALSE,vst.flavor = 'v1'))

test_that("FindTransferAnchors with default SCT works", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", k.filter = 50)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(220, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(DefaultAssay(co), "SCT")
  expect_equal(GetAssayData(co[["SCT"]], layer = "scale.data"), new(Class = "matrix"))
  expect_equal(GetAssayData(co[["SCT"]])[1, 1], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], -1.852491719, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], -0.1829401539, tolerance = 1e-7)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], -0.1971047407, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], -0.1829401539, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(256, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.688195991), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 220)
  expect_equal(anchors@anchor.features[1], "NKG7")
  expect_equal(anchors@neighbors, list())
})

test_that("Mixing SCT and non-SCT assays fails", {
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.assay = "SCT", query.assay = "RNA", k.filter = 50))
  ref.0 <- ref
  ref.2 <- ref
  ref.0[["SCT"]]@SCTModel.list <- list()
  ref.2[["SCT"]]@SCTModel.list$model2 <-  ref.2[["SCT"]]@SCTModel.list$model1
  expect_error(FindTransferAnchors(reference = ref.0, query = query, reference.assay = "SCT", query.assay = "RNA", k.filter = 50, normalization.method = "SCT"))
  #expect_error(FindTransferAnchors(reference = ref.2, query = query, reference.assay = "SCT", query.assay = "RNA", k.filter = 50, normalization.method = "SCT")) # references with multiple sct models are now supported
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.assay = "RNA", query.assay = "SCT", k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.assay = "RNA", query.assay = "SCT", k.filter = 50, normalization.method = "SCT"))
})

test_that("FindTransferAnchors with default SCT works", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", reduction = "cca", k.filter = 50)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(220, 160))
  expect_equal(Reductions(co), c("cca", "cca.l2"))
  expect_equal(DefaultAssay(co), "SCT")
  expect_equal(GetAssayData(co[["SCT"]])[1, 1], 0)
  expect_equal(dim(co[['cca']]), c(160, 30))
  expect_equal(Embeddings(co[['cca']])[1, 1], 0.0459135444, tolerance = 1e-7)
  expect_equal(Loadings(co[['cca']], projected = T)["NKG7", 1], 8.51477973, tolerance = 1e-7)
  expect_equal(dim(co[['cca.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['cca.l2']])[1, 1], 0.0625989664, tolerance = 1e-7)
  expect_equal(Loadings(co[['cca.l2']], projected = T)["NKG7", 1], 8.51477973, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(313, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.616858238), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 220)
  expect_equal(anchors@anchor.features[1], "NKG7")
  expect_equal(anchors@neighbors, list())
})

test_that("FindTransferAnchors with SCT and project.query work", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", project.query = TRUE, k.filter = 50, recompute.residuals = FALSE)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(220, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(DefaultAssay(co), "SCT")
  expect_equal(GetAssayData(co[["SCT"]])[1, 1], 0)
  expect_equal(GetAssayData(co[["SCT"]], slot = "scale.data"), new("matrix"))
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 0.3049308, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], 0.05788217444, tolerance = 1e-7)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], 0.04334884, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], 0.05788217444, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(290, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.6315789), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 220)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

test_that("FindTransferAnchors with SCT and l2.norm FALSE work", {
  skip_on_cran()
  anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", l2.norm = FALSE, k.filter = 50)
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(220, 160))
  expect_equal(Reductions(co), c("pcaproject"))
  expect_equal(DefaultAssay(co), "SCT")
  expect_equal(GetAssayData(co[["SCT"]])[1, 1], 0)
  expect_equal(GetAssayData(co[["SCT"]], layer = "scale.data"), new("matrix"))
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], -1.852491719, tolerance = 1e-7)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], -0.1829401539, tolerance = 1e-7)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(249, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.760589319), tolerance = 1e-7)
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 220)
  expect_equal(anchors@anchor.features[1], "NKG7")
  expect_equal(anchors@neighbors, list())
})

# Added test for multi-layer query 
pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))
reference <- pbmc_small
reference <- NormalizeData(reference, verbose = FALSE)
reference <- FindVariableFeatures(reference, verbose = FALSE)
reference <- ScaleData(reference, verbose = FALSE)
reference <- suppressWarnings(RunPCA(reference, npcs = 30, verbose = FALSE))

query <- CreateSeuratObject(
  counts = as.sparse(
    GetAssayData(
      object = pbmc_small[['RNA']],
      layer = "counts"
    ) + rpois(n = ncol(pbmc_small), lambda = 1)
  )
)
multilayer_query <- query
multilayer_query$dummy_group <- ifelse(
  test = as.numeric(factor(Cells(multilayer_query))) %% 2 == 0,
  yes = "Even",
  no  = "Odd"
)
multilayer_query[["RNA"]] <- split(x = multilayer_query[["RNA"]], f = multilayer_query$dummy_group)
multilayer_query <- NormalizeData(multilayer_query, verbose = FALSE)

test_that("FindTransferAnchors handles multi-layer queries when mapping.score,k is set", {
  anchors <- FindTransferAnchors(
    reference = reference,
    query = multilayer_query,
    reference.reduction = "pca",
    mapping.score.k = 30
  )
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(230, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]], layer ="data")[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")[1, 3], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(multilayer_query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_true("query.neighbors" %in% names(anchors@neighbors))
  query_layers <- multilayer_query@assays$RNA@layers
  expect_equal(length(anchors@neighbors$query.neighbors), length(query_layers)/2)
})

test_that("FindTransferAnchors handles multi-layer queries when mapping.score,k is NOT set", {
  anchors <- FindTransferAnchors(
    reference = reference,
    query = multilayer_query,
    reference.reduction = "pca",
  )
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(230, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]], layer ="data")[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], layer = "counts")[1, 3], 0)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(multilayer_query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_false("query.neighbors" %in% names(anchors@neighbors))
})

# Tests for MappingScore
context("MappingScore")

test_that("MappingScore works as expected", {
  # Calculate leverage score with and without precomputed query neighbors.
  anchors_vanilla <- FindTransferAnchors(
    reference           = reference,
    query               = multilayer_query,
    reference.reduction = "pca"
  )
  anchors_with_neighbors <- FindTransferAnchors(
    reference = reference,
    query = multilayer_query,
    reference.reduction = "pca",
    mapping.score.k = 10
  )
  for (anchors in c(anchors_vanilla, anchors_with_neighbors)) {
    scores <- MappingScore(
      anchors,
      ndim = 30,
      ksmooth = 10,
      ksnn = 5,
      kanchor = 10
    )
    # Scores should be numeric values between 0 and 1.
    expect_setequal(names(scores), Cells(multilayer_query))
    expect_true(is.numeric(scores))
    expect_false(any(is.na(scores)))
    expect_true(all(scores >= 0 & scores <= 1))
  }
})
