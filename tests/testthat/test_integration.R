# Tests for integration/transfer related fxns
set.seed(42)
pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))

# Setup test objects
ref <- pbmc_small
query <- CreateSeuratObject(
  counts = GetAssayData(object = pbmc_small[['RNA']], slot = "counts") + rpois(n = ncol(pbmc_small), lambda = 1)
)
query <- NormalizeData(object = query, verbose = FALSE)
query <- FindVariableFeatures(object = query, verbose = FALSE, nfeatures = 100)
ref <- FindVariableFeatures(object = ref, verbose = FALSE, nfeatures = 100)

# Tests for FindTransferAnchors
# ------------------------------------------------------------------------------
context("FindTransferAnchors")

anchors <- FindTransferAnchors(reference = ref, query = query, k.filter = 50)
test_that("FindTransferAnchors defaults work", {
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]])[1, 3], 4.753094626)
  expect_equal(GetAssayData(co[["RNA"]], slot = "counts")[1, 3], 1)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 0.4840944592)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], 0.2103563963)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], 0.05175486778)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], 0.2103563963)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(128, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(5, 5, 0.08361970218))
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

test_that("FindTransferAnchors catches bad input", {
  expect_error(FindTransferAnchors(reference = ref, query = query, reference.assay = "BAD", k.filter = 50)) # could have better msgs
  expect_error(FindTransferAnchors(reference = ref, query = query, query.assay = "BAD", k.filter = 50)) # could have better msgs
  expect_error(FindTransferAnchors(reference = ref, query = query, reduction = "BAD", k.filter = 50))
  expect_error(FindTransferAnchors(reference = ref, query = query, npcs = NULL, k.filter = 50)) # needs better error message
  expect_error(FindTransferAnchors(reference = ref, query = query, dims = 1:100,k.filter = 50)) # needs better error message
  expect_error(FindTransferAnchors(reference = ref, query = query, k.anchor = 80, k.filter = 50)) # needs better error message
  expect_error(FindTransferAnchors(reference = ref, query = query, k.filter = 81)) # needs better error message
  expect_error(FindTransferAnchors(reference = ref, query = query, k.filter = 50, k.score = 80)) # needs better error message
  expect_error(FindTransferAnchors(reference = ref, query = query, k.filter = 50, npcs = NULL)) # needs better error message
  # should probably error
  # FindTransferAnchors(reference = ref, query = query, k.filter = 50, reduction = "cca", project.query = TRUE)
})

anchors <- FindTransferAnchors(reference = ref, query = query, reduction = "cca", k.filter = 50)
test_that("FindTransferAnchors with cca defaults work", {
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(230, 160))
  expect_equal(Reductions(co), c("cca", "cca.l2"))
  expect_equal(GetAssayData(co[["RNA"]])[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]])[2, 1], 4.968820744)
  expect_equal(GetAssayData(co[["RNA"]], slot = "counts")[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], slot = "counts")[2, 1], 1)
  expect_equal(dim(co[['cca']]), c(160, 30))
  expect_equal(Embeddings(co[['cca']])[1, 1], 0.04611130861)
  expect_equal(Loadings(co[['cca']], projected = T)[1, 1], -1.712007073)
  expect_equal(dim(co[['cca.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['cca.l2']])[1, 1], 0.06244169641)
  expect_equal(Loadings(co[['cca.l2']], projected = T)[1, 1], -1.712007073)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(324, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.8211091234))
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

anchors <- FindTransferAnchors(reference = ref, query = query, project.query = TRUE, k.filter = 50)
test_that("FindTransferAnchors with project.query defaults work", {
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(230, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(GetAssayData(co[["RNA"]])[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]])[2, 1], 4.968820744)
  expect_equal(GetAssayData(co[["RNA"]], slot = "counts")[1, 3], 0)
  expect_equal(GetAssayData(co[["RNA"]], slot = "counts")[2, 1], 1)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], -3.403044358)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], -0.04359905081)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], -0.2365379481)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], -0.04359905081)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(117, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(4, 52, 0.06493506494))
  expect_equal(max(anchor.mat[, 2]), 79)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

anchors <- FindTransferAnchors(reference = ref, query = query, l2.norm = FALSE, k.filter = 50)
test_that("FindTransferAnchors with no l2 works", {
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(100, 160))
  expect_equal(Reductions(co), c("pcaproject"))
  expect_equal(GetAssayData(co[["RNA"]])[1, 3], 4.753094626)
  expect_equal(GetAssayData(co[["RNA"]], slot = "counts")[1, 3], 1)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], 0.4840944592)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], 0.2103563963)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(115, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(5, 5, 0.2950654582))
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 100)
  expect_equal(anchors@anchor.features[1], "PPBP")
  expect_equal(anchors@neighbors, list())
})

# SCTransform tests
query <- suppressWarnings(SCTransform(object = query, verbose = FALSE))
ref <- suppressWarnings(SCTransform(object = ref, verbose = FALSE))

anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", k.filter = 50)
test_that("FindTransferAnchors with default SCT works", {
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(220, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(DefaultAssay(co), "SCT")
  expect_equal(GetAssayData(co[["SCT"]], slot = "scale.data"), new(Class = "matrix"))
  expect_equal(GetAssayData(co[["SCT"]])[1, 1], -0.542716846)
  expect_equal(GetAssayData(co[["SCT"]], slot = "counts")[1, 1], -0.542716846)
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], -1.941335809)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], -0.1778101795)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], -0.1980063984)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], -0.1778101795)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(302, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.4497248624))
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 220)
  expect_equal(anchors@anchor.features[1], "NKG7")
  expect_equal(anchors@neighbors, list())
})

anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", reduction = "cca", k.filter = 50)
test_that("FindTransferAnchors with default SCT works", {
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(220, 160))
  expect_equal(Reductions(co), c("cca", "cca.l2"))
  expect_equal(DefaultAssay(co), "SCT")
  expect_equal(GetAssayData(co[["SCT"]])[1, 1], -0.5427168460)
  expect_equal(as.matrix(GetAssayData(co[["SCT"]])), GetAssayData(co[["SCT"]], slot = "scale.data"))
  expect_equal(GetAssayData(co[["SCT"]]), GetAssayData(co[["SCT"]], slot = "counts"))
  expect_equal(dim(co[['cca']]), c(160, 30))
  expect_equal(Embeddings(co[['cca']])[1, 1], 0.04945926847)
  expect_equal(Loadings(co[['cca']], projected = T)[1, 1], 12.32750457)
  expect_equal(dim(co[['cca.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['cca.l2']])[1, 1], 0.07510796779)
  expect_equal(Loadings(co[['cca.l2']], projected = T)[1, 1], 12.32750457)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(319, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.9565217391))
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 220)
  expect_equal(anchors@anchor.features[1], "NKG7")
  expect_equal(anchors@neighbors, list())
})

anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", project.query = TRUE, k.filter = 50)
test_that("FindTransferAnchors with SCT and project.query work", {
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(220, 160))
  expect_equal(Reductions(co), c("pcaproject", "pcaproject.l2"))
  expect_equal(DefaultAssay(co), "SCT")
  expect_equal(GetAssayData(co[["SCT"]])[1, 1], -0.5427168460)
  expect_equal(GetAssayData(co[["SCT"]], slot = "scale.data"), new("matrix"))
  expect_equal(GetAssayData(co[["SCT"]]), GetAssayData(co[["SCT"]], slot = "counts"))
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], -0.206189831)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], 0.04942085688)
  expect_equal(dim(co[['pcaproject.l2']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject.l2']])[1, 1], -0.02461344272)
  expect_equal(Loadings(co[['pcaproject.l2']], projected = T)[1, 1], 0.04942085688)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(287, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.6141124587))
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 220)
  expect_equal(anchors@anchor.features[1], "NKG7")
  expect_equal(anchors@neighbors, list())
})

anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", l2.norm = FALSE, k.filter = 50)
test_that("FindTransferAnchors with SCT and project.query work", {
  co <- anchors@object.list[[1]]
  expect_equal(dim(co), c(220, 160))
  expect_equal(Reductions(co), c("pcaproject"))
  expect_equal(DefaultAssay(co), "SCT")
  expect_equal(GetAssayData(co[["SCT"]])[1, 1], -0.5427168460)
  expect_equal(GetAssayData(co[["SCT"]], slot = "scale.data"), new("matrix"))
  expect_equal(GetAssayData(co[["SCT"]]), GetAssayData(co[["SCT"]], slot = "counts"))
  expect_equal(dim(co[['pcaproject']]), c(160, 30))
  expect_equal(Embeddings(co[['pcaproject']])[1, 1], -1.941335809)
  expect_equal(Loadings(co[['pcaproject']], projected = T)[1, 1], -0.1778101795)
  ref.cells <- paste0(Cells(ref), "_reference")
  query.cells <- paste0(Cells(query), "_query")
  expect_equal(anchors@reference.cells, ref.cells)
  expect_equal(anchors@query.cells, query.cells)
  expect_equal(anchors@reference.objects, logical())
  anchor.mat <- anchors@anchors
  expect_equal(dim(anchor.mat), c(279, 3))
  expect_equal(as.vector(anchor.mat[1, ]), c(1, 1, 0.6428571429))
  expect_equal(max(anchor.mat[, 2]), 80)
  expect_null(anchors@offsets)
  expect_equal(length(anchors@anchor.features), 220)
  expect_equal(anchors@anchor.features[1], "NKG7")
  expect_equal(anchors@neighbors, list())
})

# input cases to improve
# - npcs too large
# - npcs too small
# remove pcaqueryproject option - handled through project.query parameter
# running log norm method with SCT assays should warn/error
