#' Returns a random counts matrix.
get_random_counts <- function() {
  # Populate a 100 by 100 matrix with random integers from 1 to 50.
  counts <- matrix(
    data = sample(c(1:50), size = 1e4, replace = TRUE),
    ncol = 100,
    nrow = 100
  )

  # Assign column and row names to the matrix to label cells and genes.
  colnames(counts) <- paste0("cell", seq(ncol(counts)))
  row.names(counts) <- paste0("gene", seq(nrow(counts)))

  # Convert `counts` to a `dgCMatrix`.
  counts_sparse <- as.sparse(counts)

  return(counts_sparse)
}

#' Returns a `Seurat` instance containing the specified `assay_version` and
#' populated with `counts` which is also preprocessed (normalized + scaled).
get_test_data <- function(
  counts = get_random_counts(),
  assay_version = getOption("Seurat.object.assay.version")
) {
  # Use the `assay_version` param to choose the correct assay builder.
  create_assay <- switch(assay_version,
    v3 = CreateAssayObject,
    v5 = CreateAssay5Object,
    stop("`assay_version` should be one of 'v3', 'v5'")
  )
  # And then instantiate the specified assay type.
  assay <- create_assay(counts)

  # Instantiate a `Seurat` instance using the default assay name.
  test_data <- CreateSeuratObject(assay)

  # Normalize, and then scale the input data.
  test_data <- NormalizeData(test_data, verbose = FALSE)
  test_data <- ScaleData(test_data, verbose = FALSE)

  return(test_data)
}

#' Checks that the specified dimensional reduction `method` returns equivalent
#' results for each test case in `inputs`.
test_dimensional_reduction <- function(inputs, method, ...) {
  # Avoid replying on default reduction names.
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

context("RunPCA")

test_that("`RunPCA` returns total variance", {
  # For the motivation behind this test see https://github.com/satijalab/seurat/issues/982.
  test_case <- get_test_data()
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
  
  mat <- get_random_counts()
  # Set some features outside the supplied set to have zero variance
  mat[1:50, ] <- 0

  test_case <- get_test_data(counts = mat, assay_version = "v5")
  supplied_features <- rownames(mat)[51:80]

  result <- suppressWarnings(RunPCA(test_case, features = supplied_features, npcs = 10, verbose = FALSE))

  loadings_features <- rownames(Loadings(result[["pca"]]))
  expect_equal(length(loadings_features), length(supplied_features))
  expect_equal(loadings_features, supplied_features)
})

test_that("`RunPCA` drops zero-variance features only within the supplied feature set", {
  set.seed(123)

  mat <- get_random_counts()
  # Set some features outside the supplied set to have zero variance
  mat[1:50, ] <- 0
  # Set some features within the supplied set to have zero variance
  mat[51:60, ] <- 0

  test_case <- get_test_data(counts = mat, assay_version = "v5")

  supplied_features <- rownames(mat)[51:100]

  result <- suppressWarnings(RunPCA(test_case, features = supplied_features, npcs = 10, verbose = FALSE))

  loadings_features <- rownames(Loadings(result[["pca"]]))
  expected_features <- setdiff(supplied_features, rownames(mat)[51:60])

  expect_equal(loadings_features, expected_features)
  expect_false(any(rownames(mat)[51:60] %in% loadings_features))
})

context("RunICA")

test_that("`RunPCA` works as expected", {
  counts <- get_random_counts()
  input_v3 <- get_test_data(counts, assay_version = "v3")
  input_v5 <- get_test_data(counts, assay_version = "v5")
  inputs <- c(input_v3, input_v5)

  # If `BPCells` is installed, add `IterableMatrix` inputs to the set of 
  # equivalent inputs.
  if (requireNamespace("BPCells", quietly = TRUE)) {
    counts_bpcells <- t(as(t(counts), Class = "IterableMatrix"))
    input_bpcells <- get_test_data(counts_bpcells, assay_version = "v5")
    inputs <- c(inputs, input_bpcells)
  }

  # Check that `RunPCA` returns equivalent results for every input in `inputs`.
  test_dimensional_reduction(
    inputs = inputs,
    method = RunPCA, 
    # Reduce number of PCs from the default of 20 to avoid warning from 
    # `irlba` caused by the small size of dataset being used.
    npcs = 10
  )
})

# Determinism ------------------------------------------------------------------

set.seed(42)
determinism_object <- suppressWarnings(CreateSeuratObject(counts = get_random_counts()))
determinism_object <- suppressWarnings(NormalizeData(determinism_object, verbose = FALSE))
determinism_object <- suppressWarnings(
  FindVariableFeatures(determinism_object, verbose = FALSE)
)
determinism_object <- suppressWarnings(ScaleData(determinism_object, verbose = FALSE))
determinism_object <- suppressWarnings(
  RunPCA(determinism_object, npcs = 10, verbose = FALSE)
)

test_that("RunPCA uses a reproducible sign convention", {
  loadings <- Loadings(determinism_object[["pca"]])
  pivots <- apply(abs(loadings), 2, which.max)
  pivot.values <- loadings[cbind(pivots, seq_len(ncol(loadings)))]
  expect_true(all(pivot.values > 0))
})

test_that("RunPCA sign convention only flips signs", {
  legacy <- suppressWarnings(RunPCA(
    determinism_object,
    npcs = 10,
    verbose = FALSE,
    deterministic = FALSE
  ))
  signs <- sign(diag(crossprod(
    Loadings(legacy[["pca"]]),
    Loadings(determinism_object[["pca"]])
  )))
  expect_equal(
    sweep(Embeddings(determinism_object[["pca"]]), 2, signs, `*`),
    Embeddings(legacy[["pca"]])
  )
  expect_equal(Stdev(determinism_object[["pca"]]), Stdev(legacy[["pca"]]))
})

test_that("PC sign conventions do not change downstream results", {
  legacy <- suppressWarnings(RunPCA(
    determinism_object,
    npcs = 10,
    verbose = FALSE,
    deterministic = FALSE
  ))
  expect_equal(
    as.matrix(FindNeighbors(determinism_object, dims = 1:10, verbose = FALSE)[["RNA_snn"]]),
    as.matrix(FindNeighbors(legacy, dims = 1:10, verbose = FALSE)[["RNA_snn"]])
  )
  expect_equal(
    Embeddings(suppressWarnings(RunUMAP(
      determinism_object, dims = 1:10, n.neighbors = 10L, verbose = FALSE
    ))[["umap"]]),
    Embeddings(suppressWarnings(RunUMAP(
      legacy, dims = 1:10, n.neighbors = 10L, verbose = FALSE
    ))[["umap"]])
  )
})

test_that("RunUMAP does not depend on the ambient random stream", {
  set.seed(1)
  first <- suppressWarnings(RunUMAP(
    determinism_object, dims = 1:10, n.neighbors = 10L, verbose = FALSE
  ))
  set.seed(9999)
  invisible(runif(n = 37))
  second <- suppressWarnings(RunUMAP(
    determinism_object, dims = 1:10, n.neighbors = 10L, verbose = FALSE
  ))
  expect_identical(
    Embeddings(first[["umap"]]),
    Embeddings(second[["umap"]])
  )
})

test_that("RunUMAP does not depend on the number of workers", {
  skip_if_not_installed("future")
  sequential <- suppressWarnings(RunUMAP(
    determinism_object, dims = 1:10, n.neighbors = 10L, n.threads = 1L, verbose = FALSE
  ))
  parallel <- suppressWarnings(RunUMAP(
    determinism_object, dims = 1:10, n.neighbors = 10L, n.threads = 4L, verbose = FALSE
  ))
  expect_identical(
    Embeddings(sequential[["umap"]]),
    Embeddings(parallel[["umap"]])
  )
})

test_that("RunUMAP warns when asked for nondeterministic SGD", {
  expect_warning(
    RunUMAP(
      determinism_object,
      dims = 1:10,
      n.neighbors = 10L,
      uwot.sgd = TRUE,
      n.sgd.threads = 4L,
      verbose = FALSE
    ),
    "nondeterministic order"
  )
})

test_that("stochastic reductions leave the caller's random stream intact", {
  expect_stream_preserved <- function(expr) {
    set.seed(123)
    expected <- runif(n = 3)
    set.seed(123)
    force(expr)
    expect_identical(runif(n = 3), expected)
  }
  expect_stream_preserved(suppressWarnings(
    RunPCA(determinism_object, npcs = 10, verbose = FALSE)
  ))
  expect_stream_preserved(suppressWarnings(
    RunUMAP(determinism_object, dims = 1:10, n.neighbors = 10L, verbose = FALSE)
  ))
  expect_stream_preserved(suppressWarnings(RunTSNE(
    determinism_object,
    dims = 1:10,
    perplexity = 10,
    check_duplicates = FALSE
  )))
})

test_that("FindNeighbors is reproducible", {
  first <- FindNeighbors(determinism_object, dims = 1:10, verbose = FALSE)
  set.seed(7)
  invisible(runif(n = 11))
  second <- FindNeighbors(determinism_object, dims = 1:10, verbose = FALSE)
  expect_identical(
    as.matrix(first[["RNA_snn"]]),
    as.matrix(second[["RNA_snn"]])
  )
})
