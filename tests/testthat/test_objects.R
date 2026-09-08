# Tests for functions in objects.R

# Tests for SCE conversion
# ------------------------------------------------------------------------------
test_that("as.SingleCellExperiment works", {
  skip_on_cran()
  if (requireNamespace('SingleCellExperiment', quietly = TRUE)) {
    mat <- pbmc_small[["RNA"]]$counts
    seuratObj <- Seurat::CreateSeuratObject(mat)
    sce <- suppressWarnings(as.SingleCellExperiment(seuratObj))

    expect_equal(ncol(sce), 80)
    expect_equal(nrow(sce), 230)
    # expect_equal(length(SingleCellExperiment::altExps(sce)), 0)
    # expect_equal(SingleCellExperiment::mainExpName(sce), 'RNA')

    seuratObj <- Seurat::CreateSeuratObject(mat)
    seuratObj[['ADT']] <- CreateAssayObject(mat)
    sce <- suppressWarnings(as.SingleCellExperiment(seuratObj))
    expect_equal(ncol(sce), 80)
    expect_equal(nrow(sce), 230)
    # expect_equal(names(SingleCellExperiment::altExps(sce)), 'ADT')
    # expect_equal(SingleCellExperiment::mainExpName(sce), 'RNA')
  }
})

test_that("as.Seurat handles a SingleCellExperiment with duplicated names", {
  skip_on_cran()
  skip_if_not_installed('SingleCellExperiment')
  mat <- as.matrix(x = pbmc_small[["RNA"]]$counts[1:20, 1:6])
  # a multi-sample aggregate repeats barcodes across samples
  colnames(x = mat) <- rep(x = c('AAA', 'BBB', 'CCC'), times = 2)
  rownames(x = mat) <- c(paste0('g', 1:19), 'g1')
  mat <- as.sparse(x = mat)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat, logcounts = mat)
  )

  warns <- character(length = 0L)
  obj <- withCallingHandlers(
    expr = as.Seurat(x = sce),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(c = w))
      invokeRestart(r = 'muffleWarning')
    }
  )
  expect_true(object = any(grepl(pattern = 'Non-unique features', x = warns)))
  expect_true(object = any(grepl(pattern = 'Non-unique cell names', x = warns)))
  expect_s4_class(object = obj, class = 'Seurat')
  expect_equal(object = ncol(x = obj), expected = 6)
  expect_equal(object = nrow(x = obj), expected = 20)
  expect_equal(
    object = colnames(x = obj),
    expected = c('AAA', 'BBB', 'CCC', 'AAA.1', 'BBB.1', 'CCC.1')
  )
  expect_equal(
    object = rownames(x = obj),
    expected = c(paste0('g', 1:19), 'g1.1')
  )
  # the second matrix must survive the round trip, not just `counts`
  expect_equal(
    object = as.matrix(x = LayerData(object = obj, layer = 'data')),
    expected = as.matrix(x = LayerData(object = obj, layer = 'counts'))
  )
  expect_equal(object = rownames(x = obj[[]]), expected = colnames(x = obj))
})
