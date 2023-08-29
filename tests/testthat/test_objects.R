# Tests for functions in objects.R

# Tests for SCE conversion
# ------------------------------------------------------------------------------
test_that("as.SingleCellExperiment works", {
  skip_on_cran()
  if (requireNamespace('SingleCellExperiment', quietly = TRUE)) {
    mat <- matrix(1:100, ncol = 10)
    colnames(mat) <- LETTERS[1:10]
    rownames(mat) <- LETTERS[1:10]
    seuratObj <- Seurat::CreateSeuratObject(mat)
    sce <- as.SingleCellExperiment(seuratObj)

    expect_equal(ncol(sce), 10)
    expect_equal(nrow(sce), 10)
    # expect_equal(length(SingleCellExperiment::altExps(sce)), 0)
    # expect_equal(SingleCellExperiment::mainExpName(sce), 'RNA')

    seuratObj <- Seurat::CreateSeuratObject(mat)
    seuratObj[['ADT']] <- CreateAssayObject(mat)
    sce <- as.SingleCellExperiment(seuratObj)
    expect_equal(ncol(sce), 10)
    expect_equal(nrow(sce), 10)
    # expect_equal(names(SingleCellExperiment::altExps(sce)), 'ADT')
    # expect_equal(SingleCellExperiment::mainExpName(sce), 'RNA')
  }
})
