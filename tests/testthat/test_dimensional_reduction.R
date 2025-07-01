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
