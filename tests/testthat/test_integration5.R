# Tests for IntegrateLayers
set.seed(42)


# checks that the absolute value of `x` and `y` are within `tolerance`
expect_abs_equal <- function(x, y, tolerance = 1.0e-04) {
  expect_equal(abs(x), abs(y), tolerance = tolerance)
}


# setup shared fixtures
# update `pbmc_small` to use `Assay5` instances
test.data <- pbmc_small
suppressWarnings(
  test.data[["RNA"]] <- CreateAssay5Object(
    counts = LayerData(
      test.data,
      assay = "RNA",
      layer = "counts"
    )
  )
)
# split the assay into multiple layers
test.data[["RNA"]] <- split(test.data[["RNA"]], f = test.data$groups)


context("IntegrateLayers")

# setup fixtures for standard integration workflow
test.data.std <- NormalizeData(test.data, verbose = FALSE)
test.data.std <- FindVariableFeatures(test.data.std, verbose = FALSE)
test.data.std <- ScaleData(test.data.std, verbose = FALSE)
test.data.std <- suppressWarnings(
  RunPCA(test.data.std, verbose = FALSE)
)

test_that("IntegrateLayers works with HarmonyIntegration", {
  skip_on_cran()
  skip_if_not_installed("harmony")

  integrated <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE
    )
  )
  # the integrated reduction should have the same dimensions as the original 
  expect_equal(dim(integrated[["integrated"]]), dim(integrated[["pca"]]))
  # spot-check a few of the integrated values - since the integrated 
  # reductions sporadically flip sign only compare absolute values
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[5, 5],
    0.3912
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[40, 25],
    0.6668
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[75, 45],
    0.7248
  )
})

test_that("IntegrateLayers works with CCAIntegration", {
  integrated <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10
    )
  )
  # the integrated reduction should have the same dimensions as the original 
  expect_equal(dim(integrated[["integrated"]]), dim(integrated[["pca"]]))
  # spot-check a few of the integrated values - since the integrated 
  # reductions sporadically flip sign only compare absolute values
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[5, 5],
    0.9174
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[40, 25],
    1.4885
    )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[75, 45],
    0.5442
  )

  integrated_sub <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      dims.to.integrate = 1:10
    )
  )
  # check that the integrated reduction has the specified number of
  # `dims.to.integrate`
  expect_equal(ncol(integrated_sub[["integrated"]]), 10)

  integrated_overflow <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      dims.to.integrate = 1:100
    )
  )
  # check that the integrated reduction is the same as you'd get if you
  # didn't specify `dims.to.integrate` (i.e. the same size as the initial
  # reduction)
  expect_equal(Embeddings(integrated_overflow), Embeddings(integrated))
})

test_that("StandardizeBPCells matches C++ Standardize", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  set.seed(42)
  m <- matrix(rnorm(200), nrow = 20, ncol = 10)
  colnames(m) <- paste0("cell", 1:10)
  rownames(m) <- paste0("gene", 1:20)
  std_cpp <- Standardize(mat = m, display_progress = FALSE)
  bm <- as(as(m, "dgCMatrix"), "IterableMatrix")
  std_bp <- as.matrix(StandardizeBPCells(mat = bm))
  expect_equal(unname(std_bp), unname(std_cpp), tolerance = 1e-10)
})

test_that("RunCCA.IterableMatrix matches RunCCA.default", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  set.seed(42)
  m1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  m2 <- matrix(rnorm(600), nrow = 50, ncol = 12)
  colnames(m1) <- paste0("a_", 1:10)
  colnames(m2) <- paste0("b_", 1:12)
  rownames(m1) <- rownames(m2) <- paste0("gene", 1:50)
  res_dense <- RunCCA(object1 = m1, object2 = m2, num.cc = 5, seed.use = 42)
  bm1 <- as(as(m1, "dgCMatrix"), "IterableMatrix")
  bm2 <- as(as(m2, "dgCMatrix"), "IterableMatrix")
  res_bp <- RunCCA(object1 = bm1, object2 = bm2, num.cc = 5, seed.use = 42)
  expect_equal(abs(res_bp$ccv), abs(res_dense$ccv), tolerance = 1e-6)
  expect_equal(res_bp$d, res_dense$d, tolerance = 1e-6)
})

test_that("CCAIntegration works with BPCells on-disk data", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  test.data.bp <- test.data.std
  for (lyr in Layers(test.data.bp[["RNA"]])) {
    mat <- LayerData(test.data.bp, assay = "RNA", layer = lyr)
    bp_mat <- as(as(as.matrix(mat), "dgCMatrix"), "IterableMatrix")
    LayerData(test.data.bp, assay = "RNA", layer = lyr) <- bp_mat
  }
  # BPCells integration should run without converting to dgCMatrix
  integrated_bp <- suppressWarnings(IntegrateLayers(
    test.data.bp, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated",
    verbose = FALSE, k.weight = 10
  ))
  integrated_mem <- suppressWarnings(IntegrateLayers(
    test.data.std, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated",
    verbose = FALSE, k.weight = 10
  ))
  # the integrated reduction should have the same dimensions
  expect_equal(dim(integrated_bp[["integrated"]]), dim(integrated_mem[["integrated"]]))
  expect_equal(
    rownames(Embeddings(integrated_bp[["integrated"]])),
    rownames(Embeddings(integrated_mem[["integrated"]]))
  )
  expect_equal(
    colnames(Embeddings(integrated_bp[["integrated"]])),
    colnames(Embeddings(integrated_mem[["integrated"]]))
  )
})

test_that("IntegrateLayers works with RPCAIntegration", {
  integrated <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10
    )
  )
  # the integrated reduction should have the same dimensions as the original 
  expect_equal(dim(integrated[["integrated"]]), dim(integrated[["pca"]]))
  # spot-check a few of the integrated values - since the integrated 
  # reductions sporadically flip sign only compare absolute values
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[5, 5],
    0.1785
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[40, 25],
    0.5832
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[75, 45],
    0.5442
  )

  # check that the integrated reduction has the specified number of
  # `dims.to.integrate`
  integrated_sub <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      dims.to.integrate = 1:10
    )
  )
  # check that dims.to.integrate is not being overwritten
  expect_equal(ncol(integrated_sub[["integrated"]]), 10)

  integrated_overflow <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      dims.to.integrate = 1:100
    )
  )
  # check that the integrated reduction is the same as you'd get if you
  # didn't specify `dims.to.integrate` (i.e. the same size as the initial
  # reduction)
  expect_equal(Embeddings(integrated_overflow), Embeddings(integrated))
})

test_that("IntegrateLayers works with JointPCAIntegration", {
  integrated <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = JointPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10
    )
  )
  # the integrated reduction should have the same dimensions as the original 
  expect_equal(dim(integrated[["integrated"]]), dim(integrated[["pca"]]))
  # spot-check a few of the integrated values - since the integrated 
  # reductions sporadically flip sign only compare absolute values
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[5, 5],
    0.4092
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[40, 25],
    0.3246
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[75, 45],
    0.5442
  )
  # check that the integrated reduction has the specified number of
  # `dims.to.integrate`
  integrated_sub <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = JointPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      dims.to.integrate = 1:10
    )
  )
  # check that dims.to.integrate is not being overwritten
  expect_equal(ncol(integrated_sub[["integrated"]]), 10)

  integrated_overflow <- suppressWarnings(
    IntegrateLayers(
      test.data.std,
      method = JointPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      dims.to.integrate = 1:100
    )
  )
  # check that the integrated reduction is the same as you'd get if you
  # didn't specify `dims.to.integrate` (i.e. the same size as the initial
  # reduction)
  expect_equal(Embeddings(integrated_overflow), Embeddings(integrated))
})

test_that("IntegrateLayers fails when expected", {
  # an error should be raised if a v3 assay is passed in
  expect_error(
    IntegrateLayers(
      pbmc_small,
      method = CCAIntegration,
      orig.reduction = "pca",
      assay = "RNA",
      new.reduction = "integrated"
    )
  )

  # an error should be raised if a nonexistent `assay` is specified
  expect_error(
    IntegrateLayers(
      test.data.std,
      method = CCAIntegration,
      orig.reduction = "pca",
      assay = "DNA",
      new.reduction = "integrated"
    )
  )

  # an error should be raised if a nonexistent `orig.reduction` is specified
  expect_error(
    IntegrateLayers(
      test.data.std,
      method = CCAIntegration,
      orig.reduction = "lda",
      new.reduction = "integrated"
    )
  )
})


context("IntegrateData with SCTransform")

# setup fixtures for SCTransform workflow
test.data.sct <- suppressWarnings(
  SCTransform(
    test.data, 
    # use v1 to avoid potentially different
    # return values depending on if `glmGamPoi`
    # is installed or not
    vst.flavor="v1", 
    # set seed for reproducibility
    seed.use = 12345, 
    verbose = FALSE
  )
)
test.data.sct <- suppressWarnings(
  RunPCA(test.data.sct, verbose = FALSE)
)

test_that("IntegrateLayers works with HarmonyIntegration & SCTransform", {
  skip_on_cran()
  skip_if_not_installed("harmony")

  integrated <- suppressWarnings(
    IntegrateLayers(
      test.data.sct,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE
    )
  )
  # the integrated reduction should have the same dimensions as the original 
  expect_equal(dim(integrated[["integrated"]]), dim(integrated[["pca"]]))
  # spot-check a few of the integrated values - since the integrated 
  # reductions sporadically flip sign only compare absolute values
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[5, 5],
    1.1520
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[40, 25],
    1.0302
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[75, 45],
    0.1886
  )
})

test_that("IntegrateLayers works with CCAIntegration & SCTransform", {
  integrated <- suppressWarnings(
    IntegrateLayers(
      test.data.sct,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      # `CCAIntegration` needs to know how the data was normalized
      normalization.method = "SCT"
    )
  )
  # the integrated reduction should have the same dimensions as the original 
  expect_equal(dim(integrated[["integrated"]]), dim(integrated[["pca"]]))
  # spot-check a few of the integrated values - since the integrated 
  # reductions sporadically flip sign only compare absolute values
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[5, 5],
    1.6113
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[40, 25],
    0.6927
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[75, 45],
    0.0855
  )
})

test_that("IntegrateLayers works with RPCAIntegration & SCTransform", {
  integrated <- suppressWarnings(
    IntegrateLayers(
      test.data.sct,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      # `RPCAIntegration` needs to know how the data was normalized
      normalization.method = "SCT"
    )
  )
  # the integrated reduction should have the same dimensions as the original 
  expect_equal(dim(integrated[["integrated"]]), dim(integrated[["pca"]]))
  # spot-check a few of the integrated values - since the integrated 
  # reductions sporadically flip sign only compare absolute values
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[5, 5],
    1.6492
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[40, 25],
    0.7343
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[75, 45],
    0.0855
  )
})

test_that("IntegrateLayers works with JointPCAIntegration & SCTransform", {
  integrated <- suppressWarnings(
    IntegrateLayers(
      test.data.sct,
      method = JointPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated",
      verbose = FALSE,
      # since `k.weight` must be less than the number of samples in the
      # smallest layer being integrated, it must be set to accommodate the
      # small dataset used for testing
      k.weight = 10,
      # `JointPCAIntegration` needs to know how the data was normalized
      normalization.method = "SCT"
    )
  )
  # the integrated reduction should have the same dimensions as the original 
  expect_equal(dim(integrated[["integrated"]]), dim(integrated[["pca"]]))
  # spot-check a few of the integrated values - since the integrated 
  # reductions sporadically flip sign only compare absolute values
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[5, 5],
    0.3427
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[40, 25],
    0.1015
  )
  expect_abs_equal(
    Embeddings(integrated[["integrated"]])[75, 45],
    0.0855
  )
})
