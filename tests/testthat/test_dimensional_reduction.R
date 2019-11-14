context("test-dimensional_reduction")

test_that("different ways of passing distance matrix", {
  # Generate dummy data exp matrix
  set.seed(1)
  dummyexpMat <- matrix(data = sample(x = c(1:50), size = 1e4, replace = TRUE),
                        ncol = 100, nrow = 100)
  colnames(dummyexpMat) <- paste0("cell", seq(ncol(dummyexpMat)))
  row.names(dummyexpMat) <- paste0("gene", seq(nrow(dummyexpMat)))

  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)

  # Manually make a distance object to test
  distMat <- dist(t(dummyexpMat))

  expect_equivalent(
    suppressWarnings(expr = RunTSNE(obj, distance.matrix = distMat)),
    suppressWarnings(expr = RunTSNE(obj, distance.matrix = as.matrix(distMat)))
  )
  expect_equivalent(
    suppressWarnings(expr = RunTSNE(obj, distance.matrix = distMat)@reductions$tsne),
    suppressWarnings(expr = RunTSNE(distMat, assay = "RNA"))
  )
  expect_equivalent(
    suppressWarnings(expr = RunTSNE(obj, distance.matrix = distMat)@reductions$tsne),
    suppressWarnings(expr = RunTSNE(as.matrix(distMat), assay = "RNA", is_distance = TRUE))
  )
})

test_that("pca returns total variance (see #982)", {
  # Generate dummy data exp matrix
  set.seed(seed = 1)
  dummyexpMat <- matrix(
    data = sample(x = c(1:50), size = 1e4, replace = TRUE),
    ncol = 100, nrow = 100
  )
  colnames(x = dummyexpMat) <- paste0("cell", seq(ncol(x = dummyexpMat)))
  row.names(x = dummyexpMat) <- paste0("gene", seq(nrow(x = dummyexpMat)))

  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)

  # Scale and compute PCA, using RunPCA
  obj <- ScaleData(object = obj, verbose = FALSE)
  for (scale in c(TRUE, FALSE)) {
    pca_result <- suppressWarnings(expr = RunPCA(
      object = obj,
      features = rownames(x = obj),
      verbose = FALSE,
      scale = scale
    ))

    # Using stats::prcomp
    scaled_data <- Seurat::GetAssayData(object = obj, slot = "scale.data")
    prcomp_result <- stats::prcomp(scaled_data, center = FALSE, scale. = scale)

    # Compare
    expect_equivalent(
      slot(object = pca_result[["pca"]], name = "misc")$total.variance,
      sum(prcomp_result$sdev^2),
      label = paste0("RunPCA with `scale=", scale, "`"))
  }

})

expect_point_symmetric <- function(object, expected, ..., info = NULL, label = NULL, expected.label = NULL) {
  act <- testthat::quasi_label(rlang::enquo(object), label, arg = "object")
  exp <- testthat::quasi_label(rlang::enquo(expected), expected.label, arg = "expected")
  object_mirrored_if_needed <- t(t(act$val) * sign(act$val[1,]) * sign(exp$val[1,]))
  comp <- testthat::compare(object_mirrored_if_needed, exp$val, ...)
  expect(comp$equal, sprintf("%s not point symmetric to %s.\n%s",
                             act$lab, exp$lab, comp$message), info = info)
  invisible(act$val)
}

test_that("PCA of scaled data with and withou scaling is identical", {
  #test new behavior
  Seurat.RunPCA.use.correct.scaling.bak <- options(Seurat.RunPCA.use.correct.scaling = TRUE)
  on.exit(options(Seurat.RunPCA.use.correct.scaling.bak))

  # Generate dummy data exp matrix
  set.seed(seed = 1)
  dummyexpMat <- matrix(
    data = stats::rexp(n = 1e4, rate = 1),
    ncol = 100, nrow = 100
  )
  npcs <- 10
  colnames(x = dummyexpMat) <- paste0("cell", seq(ncol(x = dummyexpMat)))
  row.names(x = dummyexpMat) <- paste0("gene", seq(nrow(x = dummyexpMat)))

  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)

  # Normalize and compute PCA
  obj <- NormalizeData(object = obj, verbose = FALSE)

  # Scale
  obj <- ScaleData(object = obj, verbose = FALSE)

  # compute PCA with different values of related parameters
  for (approx in list(list(approx=TRUE), list(approx=FALSE))) {
    pars <- c(approx)
    pars.str <- paste0(deparse(pars), collapse="")
    pca.res.rescale <- do.call(
      what = RunPCA,
      args = c(
        list(object = obj, features = rownames(x = obj), npcs = npcs, verbose = FALSE),
        pars,
        list(scale=TRUE, center=TRUE, slot="scale.data")
      )
    )[["pca"]]
    pca.res.onscale <- do.call(
      what = RunPCA,
      args = c(
        list(object = obj, features = rownames(x = obj), npcs = npcs, verbose = FALSE),
        pars,
        list(scale=FALSE, center=FALSE, slot="scale.data")
      )
    )[["pca"]]

    expect_point_symmetric(
      object = slot(object = pca.res.rescale, name = "cell.embeddings"),
      expected = slot(object = pca.res.onscale, name = "cell.embeddings"),
      label = paste0("RunPCA rescale with ", pars.str),
      tolerance = 1E-7
    )
    expect_point_symmetric(
      object = slot(object = pca.res.rescale, name = "feature.loadings"),
      expected = slot(object = pca.res.onscale, name = "feature.loadings"),
      label = paste0("RunPCA rescale with ", pars.str),
      tolerance = 1E-7
    )
    expect_equal(
      object = slot(object = pca.res.rescale, name = "stdev"),
      expected = slot(object = pca.res.onscale, name = "stdev"),
      label = paste0("RunPCA rescale with ", pars.str),
      tolerance = 1E-7
    )
  }
})

test_that("pca with implicit scaling yields identical result", {
  #test new behavior
  Seurat.RunPCA.use.correct.scaling.bak <- options(Seurat.RunPCA.use.correct.scaling = TRUE)
  on.exit(options(Seurat.RunPCA.use.correct.scaling.bak))
  # Generate dummy data exp matrix
  set.seed(seed = 1)
  dummyexpMat <- matrix(
    data = stats::rexp(n = 1e4, rate = 1),
    ncol = 100, nrow = 100
  )
  npcs <- 10
  colnames(x = dummyexpMat) <- paste0("cell", seq(ncol(x = dummyexpMat)))
  row.names(x = dummyexpMat) <- paste0("gene", seq(nrow(x = dummyexpMat)))

  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)

  # Normalize and compute PCA
  obj <- NormalizeData(object = obj, verbose = FALSE)


  # scale and compute PCA with different values of related parameters
  for (approx in c(TRUE, FALSE)) {
    for (scale in c(TRUE, FALSE)) {
      for (center in c(TRUE, FALSE)) {
        pars <- list(approx = approx, scale = scale, center = center)
        pars.str <- paste0(deparse(pars), collapse="")

        # Scale
        obj <- ScaleData(object = obj, do.scale = scale, do.center = center, verbose = FALSE)

        pca.res.onnorm <-  RunPCA(object = obj, features = rownames(x = obj), npcs = npcs,
                                  verbose = FALSE, scale = scale, center = center, slot="data")[["pca"]]

        pca.res.onscale <- RunPCA(object = obj, features = rownames(x = obj), npcs = npcs,
                                 verbose = FALSE, scale = FALSE, center = FALSE, slot="scale.data")[["pca"]]

        expect_point_symmetric(
          object = slot(object = pca.res.onnorm, name = "cell.embeddings"),
          expected = slot(object = pca.res.onscale, name = "cell.embeddings"),
          label = paste0("RunPCA with ", pars.str),
          tolerance = 1E-7
        )
        expect_point_symmetric(
          object = slot(object = pca.res.onnorm, name = "feature.loadings"),
          expected = slot(object = pca.res.onscale, name = "feature.loadings"),
          label = paste0("RunPCA with ", pars.str),
          tolerance = 1E-7
        )
        expect_equal(
          object = slot(object = pca.res.onnorm, name = "stdev"),
          expected = slot(object = pca.res.onscale, name = "stdev"),
          label = paste0("RunPCA with ", pars.str),
          tolerance = 1E-7
        )
        expect_equal(
          object = slot(object = pca.res.onnorm, name = "misc"),
          expected = slot(object = pca.res.onscale, name = "misc"),
          label = paste0("RunPCA with ", pars.str),
          tolerance = 1E-7
        )
      }
    }
  }
})

test_that("pca embedding weighting works", {
  #test new behavior
  Seurat.RunPCA.use.correct.scaling.bak <- options(Seurat.RunPCA.use.correct.scaling = TRUE)
  on.exit(options(Seurat.RunPCA.use.correct.scaling.bak))
  # Generate dummy data exp matrix
  set.seed(seed = 1)
  npcs <- 50
  dummyexpMat <- matrix(
    data = stats::rexp(n = 2e4, rate = 1),
    ncol = 200, nrow = 100
  )
  colnames(x = dummyexpMat) <- paste0("cell", seq(ncol(x = dummyexpMat)))
  row.names(x = dummyexpMat) <- paste0("gene", seq(nrow(x = dummyexpMat)))

  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)

  # Normalize
  obj <- NormalizeData(object = obj, verbose = FALSE)
  # Scale
  obj <- ScaleData(object = obj, verbose = FALSE)

  # un(weighted/scaled)
  # compute PCA
  obj <- suppressWarnings(expr = RunPCA(
    object = obj,
    features = rownames(x = obj),
    verbose = FALSE,
    reduction.name = "pca.prcomp.unscaled",
    npcs = npcs,
    weight.by.var = FALSE,
    approx = FALSE
  ))

  obj <- suppressWarnings(expr = RunPCA(
    object = obj,
    features = rownames(x = obj),
    verbose = FALSE,
    reduction.name = "pca.irlba.unscaled",
    npcs = npcs,
    weight.by.var = FALSE,
    approx = TRUE
  ))

  # Compare
  expect_equivalent(
    diag(x = cov(x = slot(object = obj[["pca.prcomp.unscaled"]], name = "cell.embeddings"))),
    rep(x = 1/(ncol(x = obj)-1), times=npcs)
  )
  expect_equivalent(
    diag(x = cov(x = slot(object = obj[["pca.irlba.unscaled"]], name = "cell.embeddings"))),
    rep(x = 1/(ncol(x = obj)-1), times=npcs)
  )

  # weighted/scaled
  # compute PCA
  obj <- suppressWarnings(expr = RunPCA(
    object = obj,
    features = rownames(x = obj),
    verbose = FALSE,
    reduction.name = "pca.prcomp.var_scaled",
    npcs = npcs,
    weight.by.var = TRUE,
    approx = FALSE
  ))

  obj <- suppressWarnings(expr = RunPCA(
    object = obj,
    features = rownames(x = obj),
    verbose = FALSE,
    reduction.name = "pca.irlba.var_scaled",
    npcs = npcs,
    weight.by.var = TRUE,
    approx = TRUE
  ))

  # Compare
  expect_equivalent(
    diag(x = cov(x = slot(object = obj[["pca.prcomp.var_scaled"]], name = "cell.embeddings"))),
    slot(object = obj[["pca.prcomp.var_scaled"]], name = "stdev")[1:npcs]^2
  )
  expect_equivalent(
    diag(x = cov(x = slot(object = obj[["pca.irlba.var_scaled"]], name = "cell.embeddings"))),
    slot(object = obj[["pca.prcomp.var_scaled"]], name = "stdev")[1:npcs]^2
  )
})

test_that("pca reduction behaves as previously", {
  #test default behavior
  Seurat.RunPCA.use.correct.scaling.bak <- options(Seurat.RunPCA.use.correct.scaling = NULL) # change to FALSE, when changing default to TRUE
  on.exit(options(Seurat.RunPCA.use.correct.scaling.bak))
  # Generate dummy data exp matrix
  set.seed(seed = 1)
  npcs <- 3
  dummyexpMat <- matrix(
    data = stats::rexp(n = 2e4, rate = 1),
    ncol = 200, nrow = 100
  )
  colnames(x = dummyexpMat) <- paste0("cell", seq(ncol(x = dummyexpMat)))
  row.names(x = dummyexpMat) <- paste0("gene", seq(nrow(x = dummyexpMat)))

  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)

  # Normalize
  obj <- NormalizeData(object = obj, verbose = FALSE)
  # Scale
  obj <- ScaleData(object = obj, verbose = FALSE)

  # compute PCA with different values of related parameters
  for (approx in list(list(approx=TRUE), list(approx=FALSE))) {
    for (weight.by.var in list(list(weight.by.var=TRUE), list(weight.by.var=FALSE), list())) {
      pars <- c(approx, weight.by.var)
      pars.str <- paste0(deparse(pars), collapse="")
      expect_known_value(
        object = suppressWarnings(
          expr = do.call(
            what = RunPCA,
            args = c(
              list(object = obj, features = rownames(x = obj), npcs = npcs),
              pars
            )
          )
        )[["pca"]],
        file = paste0("pca", make.names(pars.str), ".rds"),
        label = paste0("RunPCA with ", pars.str)
      )
    }
  }
})

