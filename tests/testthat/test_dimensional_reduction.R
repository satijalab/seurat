source(test_path("../testdata/test-objects.R"), local = TRUE)

#' Checks that the specified dimensional reduction `method` returns equivalent
#' results for each test case in `inputs`.
test_dimensional_reduction <- function(inputs, method, ...) {
  # Use explicit reduction names.
  reduction_name = "test_reduction"

  # Run `method` on each test case in `inputs`.
  outputs <- lapply(
    inputs,
    # Use all features from the input for each dimensional reduction.
    \(input, ...) method(input, features = rownames(input), ...),
    reduction.name = reduction_name,
    verbose = FALSE,
    ...
  )

  # Fetch the embeddings for each dimensional reduction in `outputs`.
  embeddings_all <- lapply(
    outputs,
    Embeddings,
    reduction = reduction_name
  )
  embeddings_1 <- embeddings_all[[1]]
  # Check that the first set of embeddings has the expected row names.
  expect_true(all.equal(colnames(inputs[[1]]), rownames(embeddings_1)))
  # Check that all of the embeddings are equivalent.
  for (embeddings_i in embeddings_all[-1]) {
    expect_equivalent(
      abs(embeddings_1), 
      abs(embeddings_i), 
      tolerance = 1e-4
    )
  }

  # Fetch the feature loadings for each dimensional reduction in `outputs`.
  loadings_all <- lapply(
    outputs,
    Loadings,
    reduction = reduction_name
  )
  loadings_1 <- loadings_all[[1]]
  # Check that the first set of feature loadings has the expected row names.
  expect_true(all.equal(rownames(inputs[[1]]), rownames(loadings_1)))
  # Check that all of the feature loadings are equivalent.
  for (loadings_i in loadings_all[-1]) {
    expect_equivalent(
      abs(loadings_1), 
      abs(loadings_i), 
      tolerance = 1e-4
    )
  }
}

context("RunCCA")

test_that("CCA C++ implicit multiply matches crossprod", {
  set.seed(123)
  left <- matrix(rnorm(40 * 25), nrow = 40)
  right <- matrix(rnorm(40 * 30), nrow = 40)
  x <- rnorm(ncol(right))
  expected <- as.numeric(crossprod(x = left, y = right %*% x))

  expect_equal(CcaCrossprodMultiply(left, right, x), expected, tolerance = 1e-12)
})

test_that("CCA C++ implicit multiply validates dimensions", {
  expect_error(
    CcaCrossprodMultiply(matrix(1, 3, 2), matrix(1, 4, 2), rep(1, 2)),
    "same number of rows"
  )
  expect_error(
    CcaCrossprodMultiply(matrix(1, 3, 2), matrix(1, 3, 2), rep(1, 3)),
    "must match the number of columns"
  )
})

test_that("RunCCA RSpectra backend matches explicit CCA within tolerance", {
  set.seed(123)
  object1 <- matrix(rnorm(80 * 35), nrow = 80)
  object2 <- matrix(rnorm(80 * 32), nrow = 80)
  colnames(object1) <- paste0("cell_a_", seq_len(ncol(object1)))
  colnames(object2) <- paste0("cell_b_", seq_len(ncol(object2)))

  cca.irlba <- RunCCA(
    object1 = object1,
    object2 = object2,
    num.cc = 8,
    svd.method = "irlba"
  )
  cca.rspectra <- RunCCA(
    object1 = object1,
    object2 = object2,
    num.cc = 8,
    svd.method = "rspectra"
  )

  expect_equal(cca.rspectra$d, cca.irlba$d, tolerance = 1e-5)
  expect_equal(abs(cca.rspectra$ccv), abs(cca.irlba$ccv), tolerance = 1e-5)
  expect_equal(dim(cca.rspectra$ccv), c(ncol(object1) + ncol(object2), 8))
  expect_equal(rownames(cca.rspectra$ccv), c(colnames(object1), colnames(object2)))
  expect_equal(colnames(cca.rspectra$ccv), paste0("CC", seq_len(8)))
  expect_true(all(is.finite(cca.rspectra$d)))
})

test_that("RunCCA RSpectra backend matches explicit CCA without standardization", {
  set.seed(123)
  object1 <- matrix(rnorm(80 * 35), nrow = 80)
  object2 <- matrix(rnorm(80 * 32), nrow = 80)
  colnames(object1) <- paste0("cell_a_", seq_len(ncol(object1)))
  colnames(object2) <- paste0("cell_b_", seq_len(ncol(object2)))

  cca.irlba <- RunCCA(
    object1 = object1,
    object2 = object2,
    standardize = FALSE,
    num.cc = 8,
    svd.method = "irlba"
  )
  cca.rspectra <- RunCCA(
    object1 = object1,
    object2 = object2,
    standardize = FALSE,
    num.cc = 8,
    svd.method = "rspectra"
  )

  expect_equal(cca.rspectra$d, cca.irlba$d, tolerance = 1e-5)
  expect_equal(abs(cca.rspectra$ccv), abs(cca.irlba$ccv), tolerance = 1e-5)
  expect_equal(dim(cca.rspectra$ccv), c(ncol(object1) + ncol(object2), 8))
  expect_equal(rownames(cca.rspectra$ccv), c(colnames(object1), colnames(object2)))
  expect_equal(colnames(cca.rspectra$ccv), paste0("CC", seq_len(8)))
  expect_true(all(is.finite(cca.rspectra$d)))
})

test_that("RunCCA RSpectra backend supports sparse inputs without standardization", {
  set.seed(123)
  object1 <- Matrix::rsparsematrix(80, 35, density = 0.2)
  object2 <- Matrix::rsparsematrix(80, 32, density = 0.2)
  colnames(object1) <- paste0("cell_a_", seq_len(ncol(object1)))
  colnames(object2) <- paste0("cell_b_", seq_len(ncol(object2)))

  cca.irlba <- RunCCA(
    object1 = object1,
    object2 = object2,
    standardize = FALSE,
    num.cc = 8,
    svd.method = "irlba"
  )
  cca.rspectra <- RunCCA(
    object1 = object1,
    object2 = object2,
    standardize = FALSE,
    num.cc = 8,
    svd.method = "rspectra"
  )

  expect_equal(cca.rspectra$d, cca.irlba$d, tolerance = 1e-5)
  expect_equal(abs(cca.rspectra$ccv), abs(cca.irlba$ccv), tolerance = 1e-5)
})

test_that("RunCCA default backend switches only for large dense CCA inputs", {
  set.seed(123)
  small1 <- matrix(rnorm(20 * 35), nrow = 20)
  small2 <- matrix(rnorm(20 * 32), nrow = 20)
  large1 <- matrix(rnorm(20 * 1000), nrow = 20)
  large2 <- matrix(rnorm(20 * 1000), nrow = 20)
  run_cca <- function(x, y, svd.method = NULL) {
    RunCCA(x, y, standardize = FALSE, num.cc = 3, svd.method = svd.method)
  }
  expect_equal(
    expect_message(run_cca(small1, small2), NA),
    run_cca(small1, small2, "irlba")
  )
  expect_equal(
    expect_message(run_cca(large1, large2), regexp = "Using RSpectra for CCA"),
    run_cca(large1, large2, "rspectra")
  )
})

context("RunPCA")

test_that("`RunPCA` returns total variance", {
  # For the motivation behind this test see https://github.com/satijalab/seurat/issues/982.
  counts <- create_random_counts()
  object <- create_seurat_obj(counts = counts)
  test_case <- create_prep_obj(
    object = object,
    normalize = TRUE,
    scale = TRUE
  )
  counts_scaled <- LayerData(test_case, layer = "scale.data")

  # Calculate the expected total variance using `prcomp`
  prcomp_result <- stats::prcomp(
    counts_scaled, 
    center = FALSE, 
    scale. = FALSE
  )
  expected_total_variance <- sum(prcomp_result$sdev^2)

  pca_result <- suppressWarnings(
    RunPCA(
      test_case,
      features = rownames(counts_scaled),
      verbose = FALSE
    )
  )

  expect_equivalent(
    expected_total_variance,
    slot(object = pca_result[["pca"]], name = "misc")$total.variance,
  )
})

test_that("`RunPCA` does not drop features when there are zero-variance features outside the supplied set", {
  set.seed(123)
  
  mat <- create_random_counts()
  # Set some features outside the supplied set to have zero variance
  mat[1:50, ] <- 0

  object <- create_seurat_obj(counts = mat, assay_version = "v5")
  test_case <- create_prep_obj(
    object = object,
    normalize = TRUE,
    scale = TRUE
  )
  supplied_features <- rownames(mat)[51:80]

  result <- suppressWarnings(RunPCA(test_case, features = supplied_features, npcs = 10, verbose = FALSE))

  loadings_features <- rownames(Loadings(result[["pca"]]))
  expect_equal(length(loadings_features), length(supplied_features))
  expect_equal(loadings_features, supplied_features)
})

test_that("`RunPCA` drops zero-variance features only within the supplied feature set", {
  set.seed(123)

  mat <- create_random_counts()
  # Set some features outside the supplied set to have zero variance
  mat[1:50, ] <- 0
  # Set some features within the supplied set to have zero variance
  mat[51:60, ] <- 0

  object <- create_seurat_obj(counts = mat, assay_version = "v5")
  test_case <- create_prep_obj(
    object = object,
    normalize = TRUE,
    scale = TRUE
  )

  supplied_features <- rownames(mat)[51:100]

  result <- suppressWarnings(RunPCA(test_case, features = supplied_features, npcs = 10, verbose = FALSE))

  loadings_features <- rownames(Loadings(result[["pca"]]))
  expected_features <- setdiff(supplied_features, rownames(mat)[51:60])

  expect_equal(loadings_features, expected_features)
  expect_false(any(rownames(mat)[51:60] %in% loadings_features))
})

test_that("RunPCA preserves RNG state and is reproducible with a seed", {
  set.seed(1)
  counts <- create_random_counts()
  object <- create_seurat_obj(counts = counts, assay_version = "v5")
  test_case <- create_prep_obj(object = object, normalize = TRUE, scale = TRUE)
  features <- rownames(x = LayerData(object = test_case, layer = "scale.data"))[1:30]

  before <- .Random.seed
  pca1 <- {result <- suppressWarnings(RunPCA(test_case, features = features, npcs = 10, seed.use = 42, verbose = FALSE)); expect_identical(.Random.seed, before); result}
  pca2 <- {result <- suppressWarnings(RunPCA(test_case, features = features, npcs = 10, seed.use = 42, verbose = FALSE)); expect_identical(.Random.seed, before); result}
  expect_identical(.Random.seed, before)
  expect_equal(abs(Embeddings(pca1[["pca"]])), abs(Embeddings(pca2[["pca"]])))
})

test_that("RunUMAP preserves RNG state and is reproducible with a seed", {
  set.seed(1)
  object <- create_seurat_obj(counts = create_random_counts(), assay_version = "v5")
  object <- create_prep_obj(object = object, normalize = TRUE, scale = TRUE)
  object <- suppressWarnings(RunPCA(object, npcs = 10, verbose = FALSE))

  before <- .Random.seed
  umap1 <- {result <- RunUMAP(object = object, dims = 1:10, seed.use = 42, n.epochs = 0L, verbose = FALSE); expect_identical(.Random.seed, before); result}
  umap2 <- {result <- RunUMAP(object = object, dims = 1:10, seed.use = 42, n.epochs = 0L, verbose = FALSE); expect_identical(.Random.seed, before); result}
  expect_identical(.Random.seed, before)
  expect_equal(Embeddings(umap1[["umap"]]), Embeddings(umap2[["umap"]]))
})

test_that("RunTSNE preserves RNG state and is reproducible with a seed", {
  set.seed(1)
  object <- create_seurat_obj(counts = create_random_counts(), assay_version = "v5")
  object <- create_prep_obj(object = object, normalize = TRUE, scale = TRUE)
  object <- suppressWarnings(RunPCA(object, npcs = 10, verbose = FALSE))
  pca <- object[["pca"]]

  before <- .Random.seed
  tsne1 <- {result <- RunTSNE(pca, dims = 1:5, seed.use = 42, dim.embed = 2, verbose = FALSE); expect_identical(.Random.seed, before); result}
  tsne2 <- {result <- RunTSNE(pca, dims = 1:5, seed.use = 42, dim.embed = 2, verbose = FALSE); expect_identical(.Random.seed, before); result}
  expect_identical(.Random.seed, before)
  expect_equal(Embeddings(tsne1), Embeddings(tsne2))
})

context("RunICA")

test_that("`RunPCA` works as expected", {
  counts <- create_random_counts()
  object_v3 <- create_seurat_obj(counts, assay_version = "v3")
  input_v3 <- create_prep_obj(
    object = object_v3,
    normalize = TRUE,
    scale = TRUE
  )
  object_v5 <- create_seurat_obj(counts, assay_version = "v5")
  input_v5 <- create_prep_obj(
    object = object_v5,
    normalize = TRUE,
    scale = TRUE
  )
  inputs <- c(input_v3, input_v5)

  # If `BPCells` is installed, add `IterableMatrix` inputs to the set of 
  # equivalent inputs.
  if (requireNamespace("BPCells", quietly = TRUE)) {
    counts_bpcells <- t(as(t(counts), Class = "IterableMatrix"))
    object_bpcells <- create_seurat_obj(counts_bpcells, assay_version = "v5")
    input_bpcells <- create_prep_obj(
      object = object_bpcells,
      normalize = TRUE,
      scale = TRUE
    )
    inputs <- c(inputs, input_bpcells)
  }

  # Check that `RunPCA` returns equivalent results for every input in `inputs`.
  test_dimensional_reduction(
    inputs = inputs,
    method = RunPCA, 
    # Use fewer PCs for the small fixture.
    npcs = 10
  )
})
