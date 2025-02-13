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


# Tests for merge SCT
# ------------------------------------------------------------------------------

pbmc_small_sct <- pbmc_small
pbmc_small_sct[['RNA']] <- CreateAssay5Object(counts = LayerData(pbmc_small_sct[['RNA']], layer = "counts"))
pbmc_small_sct <- SCTransform(pbmc_small)
test_that("merge.SCTAssay throws as a warning if attempting to merge with an non-SCT assay",  {
  expect_warning(merge(pbmc_small_sct[['SCT']], pbmc_small[['RNA']]))
})

test_that("merge.SCTAssay works for multi-layer objects",  {
  g1 <- LayerData(subset(pbmc_small, groups=="g1"), layer="counts")
  colnames(g1) <- paste0(colnames(g1), "-g1")
  g1 <- CreateSeuratObject(counts = g1)
  g2 <- LayerData(subset(pbmc_small, groups=="g2"), layer="counts")
  colnames(g2) <- paste0(colnames(g2), "-g2")
  g2 <- CreateSeuratObject(counts = g2)
  m1 <- merge(g1, g2)
  m1 <- SCTransform(m1)
  m2 <- merge(pbmc_small_sct, m1)
  expect_equal(length(m2[['SCT']]@SCTModel.list), 3)
})
