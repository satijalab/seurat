# PCAPlot/TSNEPlot/UMAPPlot used to infer their reduction from the name of the
# calling function in the call stack, so any call that was not a bare
# `PCAPlot(obj)` looked up the wrong reduction.
obj <- suppressWarnings(pbmc_small)

# A DimPlot's data columns are named from the reduction's key (PC_1, tSNE_1),
# so they identify which reduction was actually used
used_reduction <- function(plot) {
  colnames(plot$data)[1]
}

test_that("PCAPlot uses the pca reduction however it is called", {
  expected <- used_reduction(suppressWarnings(PCAPlot(obj)))
  expect_true(nzchar(expected))
  expect_identical(used_reduction(suppressWarnings(Seurat::PCAPlot(obj))), expected)
  aliased <- PCAPlot
  expect_identical(used_reduction(suppressWarnings(aliased(obj))), expected)
  expect_identical(
    used_reduction(suppressWarnings(do.call(PCAPlot, list(obj)))),
    expected
  )
  expect_identical(
    used_reduction(suppressWarnings(lapply(list(obj), PCAPlot)[[1]])),
    expected
  )
})

test_that("TSNEPlot uses the tsne reduction however it is called", {
  expected <- used_reduction(suppressWarnings(TSNEPlot(obj)))
  expect_true(nzchar(expected))
  expect_identical(used_reduction(suppressWarnings(Seurat::TSNEPlot(obj))), expected)
  expect_identical(
    used_reduction(suppressWarnings(do.call(TSNEPlot, list(obj)))),
    expected
  )
})

test_that("PCAPlot and TSNEPlot pick different reductions", {
  expect_false(identical(
    used_reduction(suppressWarnings(PCAPlot(obj))),
    used_reduction(suppressWarnings(TSNEPlot(obj)))
  ))
})

test_that("extra arguments still reach DimPlot", {
  plot <- suppressWarnings(PCAPlot(obj, group.by = "letter.idents"))
  expect_s3_class(plot, "ggplot")
  expect_true("letter.idents" %in% colnames(plot$data))
  expect_setequal(as.character(unique(plot$data$letter.idents)), c("A", "B"))
  # and the reduction is still the pca one
  expect_identical(used_reduction(plot), used_reduction(suppressWarnings(PCAPlot(obj))))

  subset <- suppressWarnings(PCAPlot(obj, cells = colnames(obj)[1:10]))
  expect_equal(nrow(subset$data), 10L)
})
