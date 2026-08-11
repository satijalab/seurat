# Tests for IntegrateLayers
set.seed(42)
is_not_cran_submission <- isTRUE(as.logical(Sys.getenv("NOT_CRAN")))


# checks that the absolute value of `x` and `y` are within `tolerance`
expect_abs_equal <- function(x, y, tolerance = 1.0e-04) {
  expect_equal(abs(x), abs(y), tolerance = tolerance)
}

expected_harmony_embeddings <- function(object, assay = NULL, orig.reduction = "pca", key = "harmony_") {
  assay <- assay %||% DefaultAssay(object = object)
  assay.object <- object[[assay]]
  is.sct <- inherits(x = assay.object, what = "SCTAssay")

  # use the exact options & default parameters used in HarmonyIntegration
  advanced_options <- harmony::harmony_options(
    tau = 0,
    block.size = 0.05,
    max.iter.cluster = 20L,
    epsilon.cluster = 1e-05,
    epsilon.harmony = 0.01
  )
  embeddings <- harmony::RunHarmony(
    data_mat = Embeddings(object = object[[orig.reduction]]),
    meta_data = CreateIntegrationGroups(
      object = assay.object,
      layers = if (is.sct) "data" else Layers(object = object, assay = assay, search = "data"),
      scale.layer = if (is.sct) "scale.data" else Layers(object = object, search = "scale.data")
    ),
    vars_use = 'group',
    theta = NULL,
    lambda = NULL,
    sigma = 0.1,
    nclust = NULL,
    max_iter = 10L,
    return_object = FALSE,
    verbose = FALSE,
    .options = advanced_options
  )
  rownames(x = embeddings) <- Cells(x = object[[orig.reduction]])
  colnames(x = embeddings) <- paste0(key, seq_len(length.out = ncol(x = embeddings)))
  return(embeddings)
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

if (is_not_cran_submission) {
  test_that("IntegrateLayers works with HarmonyIntegration", {
    skip_if_not_installed("harmony")

    set.seed(seed = 42)
    expected <- expected_harmony_embeddings(test.data.std)

    set.seed(seed = 42)
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
    # but should actually perform integration and produce different values than the original reduction
    expect_false(isTRUE(all.equal(
      Embeddings(object = integrated[["integrated"]]),
      Embeddings(object = integrated[["pca"]])
    )))
    # check that the integrated reduction is the same as calling harmony::RunHarmony directly
    expect_equal(
      Embeddings(object = integrated[["integrated"]]),
      expected,
      tolerance = 1e-06
    )
  })
}

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


if (is_not_cran_submission) {
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
    skip_if_not_installed("harmony")

    set.seed(seed = 42)
    expected <- expected_harmony_embeddings(test.data.sct)

    set.seed(seed = 42)
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
    # but should actually perform integration and produce different values than the original reduction
    expect_false(isTRUE(all.equal(
      Embeddings(object = integrated[["integrated"]]),
      Embeddings(object = integrated[["pca"]])
    )))
    # check that the integrated reduction is the same as calling harmony::RunHarmony directly
    expect_equal(
      Embeddings(object = integrated[["integrated"]]),
      expected,
      tolerance = 1e-06
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
}


