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

test_that("as.Seurat converts alternative experiments that lack an assay", {
  skip_on_cran()
  skip_if_not_installed("SingleCellExperiment")
  counts <- pbmc_small[["RNA"]]$counts
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts, logcounts = log1p(counts))
  )
  # spike-ins with counts only, and ADTs with neither counts nor logcounts
  SingleCellExperiment::altExp(sce, "Spikes") <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts[1:10, ])
  )
  SingleCellExperiment::altExp(sce, "ADT") <- SingleCellExperiment::SingleCellExperiment(
    assays = list(clr = log1p(counts[11:15, ]))
  )
  warnings <- capture_warnings(object <- as.Seurat(sce))
  expect_match(warnings, "'Spikes' has no 'logcounts' assay; converting the rest", all = FALSE)
  expect_match(warnings, "'ADT' has no 'counts' or 'logcounts' assay; skipping it", all = FALSE)
  expect_setequal(Assays(object), c("originalexp", "Spikes"))
  expect_equal(object[["originalexp"]]$data, log1p(counts))
  expect_equal(object[["Spikes"]]$counts, counts[1:10, ])
  # a missing assay in the main experiment is still an error
  sce.counts <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  expect_error(as.Seurat(sce.counts), "No data in provided assay - logcounts")
})
