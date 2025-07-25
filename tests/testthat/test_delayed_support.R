# Test DelayedArray support in Seurat
# 
# This file contains tests for DelayedArray integration

test_that("DelayedArray support detection works", {
  # Test CheckDelayedArraySupport function
  result <- CheckDelayedArraySupport(verbose = FALSE)
  expect_is(result, "logical")
  expect_length(result, 1)
})

test_that("Matrix format optimization works", {
  # Skip if DelayedArray not available
  skip_if_not_installed("DelayedArray")
  
  # Create a test matrix
  test_matrix <- matrix(rnorm(100), 10, 10)
  rownames(test_matrix) <- paste0("Gene_", 1:10)
  colnames(test_matrix) <- paste0("Cell_", 1:10)
  
  # Test auto optimization
  optimized <- OptimizeMatrixFormat(test_matrix, prefer_format = "sparse")
  expect_s4_class(optimized, "dgCMatrix")
  
  # Test DelayedMatrix conversion
  if (requireNamespace("DelayedArray", quietly = TRUE)) {
    delayed <- OptimizeMatrixFormat(test_matrix, prefer_format = "delayed")
    expect_s4_class(delayed, "DelayedMatrix")
  }
})

test_that("DelayedMatrix conversion works", {
  skip_if_not_installed("DelayedArray")
  
  # Create test matrix
  test_matrix <- matrix(rnorm(100), 10, 10)
  
  # Test as.DelayedMatrix
  delayed <- as.DelayedMatrix(test_matrix)
  expect_s4_class(delayed, "DelayedMatrix")
  expect_equal(dim(delayed), dim(test_matrix))
})

test_that("VST works with DelayedMatrix", {
  skip_if_not_installed("DelayedArray")
  
  # Create test DelayedMatrix
  test_matrix <- matrix(rpois(200, 5), 20, 10)
  rownames(test_matrix) <- paste0("Gene_", 1:20)
  delayed_matrix <- as.DelayedMatrix(test_matrix)
  
  # Test VST.DelayedMatrix
  hvf_result <- VST.DelayedMatrix(delayed_matrix, nselect = 5)
  
  expect_is(hvf_result, "data.frame")
  expect_equal(nrow(hvf_result), 20)
  expect_true(sum(hvf_result$variable) <= 5)
  expect_true(all(c("mean", "variance", "variable", "rank") %in% colnames(hvf_result)))
})

test_that("as.sparse.DelayedMatrix works", {
  skip_if_not_installed("DelayedArray")
  
  # Create test DelayedMatrix
  test_matrix <- Matrix::sparseMatrix(
    i = c(1, 3, 5),
    j = c(1, 2, 3),
    x = c(1, 2, 3),
    dims = c(10, 10)
  )
  delayed_matrix <- as.DelayedMatrix(test_matrix)
  
  # Test conversion back to sparse
  sparse_result <- as.sparse.DelayedMatrix(delayed_matrix)
  expect_s4_class(sparse_result, "dgCMatrix")
  expect_equal(dim(sparse_result), dim(delayed_matrix))
})

test_that("MemoryEfficientMatrix class union works", {
  # Test that different matrix types are recognized
  dense_matrix <- matrix(1:12, 3, 4)
  sparse_matrix <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1:3)
  
  expect_true(inherits(dense_matrix, "MemoryEfficientMatrix"))
  expect_true(inherits(sparse_matrix, "MemoryEfficientMatrix"))
  
  if (requireNamespace("DelayedArray", quietly = TRUE)) {
    delayed_matrix <- as.DelayedMatrix(dense_matrix)
    expect_true(inherits(delayed_matrix, "MemoryEfficientMatrix"))
  }
  
  if (requireNamespace("spam", quietly = TRUE)) {
    spam_matrix <- spam::as.spam(dense_matrix)
    expect_true(inherits(spam_matrix, "MemoryEfficientMatrix"))
  }
})

test_that("Large matrix warnings work", {
  # Test memory threshold warnings in CheckMatrixSize
  expect_message(
    CheckMatrixSize(nrows = 10000, ncols = 20000),  # 200M elements
    "Large matrix detected"
  )
  
  # Test that smaller matrices don't trigger warnings
  expect_silent(CheckMatrixSize(nrows = 100, ncols = 100))
})

test_that("DelayedArray utilities work with row operations", {
  skip_if_not_installed("DelayedArray")
  
  # Create test matrix
  test_matrix <- matrix(rnorm(100), 10, 10)
  delayed_matrix <- as.DelayedMatrix(test_matrix)
  
  # Test row operations
  row_means <- RowMeanSparse(delayed_matrix)
  row_sums <- RowSumSparse(delayed_matrix)
  row_vars <- RowVarSparse(delayed_matrix)
  
  expect_length(row_means, 10)
  expect_length(row_sums, 10)
  expect_length(row_vars, 10)
  
  # Compare with base R results
  expected_means <- rowMeans(test_matrix)
  expected_sums <- rowSums(test_matrix)
  
  expect_equal(row_means, expected_means, tolerance = 1e-10)
  expect_equal(row_sums, expected_sums, tolerance = 1e-10)
})
