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

anchors <- FindTransferAnchors(reference = ref, query = query, k.filter = 50)

# Tests for TransferData
# ------------------------------------------------------------------------------
context("TransferData")

preds.standard <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, verbose = FALSE)
test_that("TransferData default work", {
  # categorical metadata
  expect_equal(dim(preds.standard), c(80, 5))
  expect_equal(colnames(preds.standard)[c(1, 5)], c("predicted.id", "prediction.score.max"))
  expect_equal(rownames(preds.standard), Cells(query))
  expect_equal(preds.standard[1, 1], "1")
  expect_equal(preds.standard[1, 5], 0.4280746, tolerance = 1e-6)
  expect_equal(as.vector(rowSums(as.matrix(preds.standard[, 2:4]))), rep(1, times = ncol(query)))
  expect_true(inherits(preds.standard, "data.frame"))
  # continuous assay data
  pred.assay <- TransferData(anchorset = anchors, refdata = GetAssayData(ref[["RNA"]]), verbose = FALSE)
  expect_equal(dim(pred.assay), c(230, 80))
  expect_equal(GetAssayData(pred.assay, slot = "counts"), new("matrix"))
  expect_equal(GetAssayData(pred.assay, slot = "scale.data"), new("matrix"))
  expect_equal(colnames(pred.assay), Cells(query))
  expect_equal(rownames(pred.assay), rownames(ref[["RNA"]]))
  expect_equal(sum(GetAssayData(pred.assay)[1, ]), 64.46388, tolerance = 1e-6)
  expect_equal(sum(GetAssayData(pred.assay)[, 1]), 281.0306, tolerance = 1e-6)
  expect_true(inherits(pred.assay, "Assay"))
  expect_equal(pred.assay@var.features, logical(0))
  expect_equal(ncol(pred.assay@meta.features), 0)
})

test_that("TransferData can return predictions assay, ", {
  pred.assay <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, prediction.assay = TRUE, verbose = FALSE)
  expect_true(inherits(pred.assay, "Assay"))
  expect_equal(dim(pred.assay), c(4, 80))
  expect_equal(GetAssayData(pred.assay, slot = "counts"), new("matrix"))
  expect_equal(GetAssayData(pred.assay, slot = "scale.data"), new("matrix"))
  expect_equal(colnames(pred.assay), Cells(query))
  expect_equal(pred.assay@var.features, logical(0))
  expect_equal(ncol(pred.assay@meta.features), 0)
  expect_equal(sum(GetAssayData(pred.assay)[1, ]), 26.59365, tolerance = 1e-6)
  expect_equal(sum(GetAssayData(pred.assay)[, 1]), 1.428075, tolerance = 1e-6)
  expect_equal(as.vector(colSums(GetAssayData(pred.assay)[1:3, ])), rep(1, ncol(query)))
})

test_that("TransferData handles weight.reduction properly, ", {
  skip_on_cran()
  # test for custom dimreduc
  custom.dr <- anchors@object.list[[1]][["pcaproject"]]
  custom.dr <- subset(x = custom.dr, cells = anchors@query.cells)
  custom.dr <- RenameCells(object = custom.dr, new.names = sapply(X = Cells(custom.dr), FUN = function(x){
    x <- gsub(pattern = "_query", replacement = "", x = x)
  }))
  expect_error(TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, weight.reduction = custom.dr, dims = 1:100))
  preds <-TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, verbose = FALSE)
  cdr.preds <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, weight.reduction = custom.dr, verbose = FALSE, dims = 1:30)
  expect_equal(preds, cdr.preds)
  # weight.reduction = "pca
  pca.preds <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, query = query, weight.reduction = "pca", verbose = FALSE)
  expect_true(inherits(pca.preds, "Seurat"))
  expect_equal(sum(GetAssayData(pca.preds[['prediction.score.id']])[1, ]), 27.83330252, tolerance = 1e-6)
  # weight.reduction = "cca"
  anchors.cca <- FindTransferAnchors(reference = ref, query = query, k.filter = 50, reduction = "cca")
  cca.preds <- TransferData(anchorset = anchors.cca, refdata = ref$RNA_snn_res.1, weight.reduction = "cca", verbose = FALSE)
  expect_true(inherits(cca.preds, "data.frame"))
  expect_equal(sum(cca.preds[, 2]), 43.61738383, tolerance = 1e-6)
})

test_that("TransferData with l2.norm works", {
  skip_on_cran()
  preds <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, l2.norm = TRUE, verbose = FALSE)
  expect_equal(dim(preds), c(80, 5))
  expect_equal(colnames(preds)[c(1, 5)], c("predicted.id", "prediction.score.max"))
  expect_equal(rownames(preds), Cells(query))
  expect_equal(preds[1, 1], "0")
  expect_equal(preds[1, 5], 0.3973124793, tolerance = 1e-6)
  expect_equal(as.vector(rowSums(as.matrix(preds[, 2:4]))), rep(1, times = ncol(query)))
  expect_true(inherits(preds, "data.frame"))
})

test_that("TransferData with other k.weight works", {
  skip_on_cran()
  preds <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, k.weight = 10, verbose = FALSE)
  expect_equal(dim(preds), c(80, 5))
  expect_equal(colnames(preds)[c(1, 5)], c("predicted.id", "prediction.score.max"))
  expect_equal(rownames(preds), Cells(query))
  expect_equal(preds[1, 1], "2")
  expect_equal(preds[1, 5], 0.6145459065, tolerance = 1e-6)
  expect_equal(as.vector(rowSums(as.matrix(preds[, 2:4]))), rep(1, times = ncol(query)))
  expect_true(inherits(preds, "data.frame"))
})

test_that("TransferData with reference specified works", {
  skip_on_cran()
  pred2 <- TransferData(anchorset = anchors, refdata = "RNA_snn_res.1", reference = ref, verbose = FALSE)
  expect_equal(preds.standard, pred2)
})

test_that("TransferData throws expected errors ", {
  expect_error(TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, weight.reduction = "BAD")) # better message
  expect_error(TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, weight.reduction = "cca")) # better message
  expect_error(TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, dims = 1:100))
  expect_error(ransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, k.weight = 1000))
  expect_error(suppressWarnings(TransferData(anchorset = anchors, refdata = "RNA_snn_res.1")))
  expect_error(suppressWarnings(TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1[1:10])))
  expect_error(TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, query = subset(x = query, cells = Cells(query)[1:10])))
})

test_that("TransferData with multiple items to transfer works ", {
  skip_on_cran()
  preds <- TransferData(anchorset = anchors, refdata = list(
    ids = ref$RNA_snn_res.1, groups = ref$groups, dat = GetAssayData(ref[["RNA"]])),
    verbose = FALSE)
  expect_equal(length(preds), 3)
  expect_equal(preds[[1]], preds.standard)
})

test_that("TransferData can return a modified query object ", {
  query <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, query = query, verbose = FALSE)
  expect_true("prediction.score.id" %in% Assays(query))
  expect_true("predicted.id" %in% colnames(query[[]]))
  expect_true("predicted.id.score" %in% colnames(query[[]]))
  query <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, query = query, store.weights = TRUE, verbose = FALSE)
  expect_equal(dim(Tool(query, slot = "TransferData")$weights.matrix), c(128, 80))
})

