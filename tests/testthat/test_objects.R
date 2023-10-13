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
