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
