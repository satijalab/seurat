# Tests for spam matrix support
# Tests for functions with spam matrix support

# Skip tests if spam is not available
if (!requireNamespace("spam", quietly = TRUE)) {
  skip("spam package not available")
}

context("Spam matrix support")

# Create a simple spam matrix for testing
test_spam_matrix <- function() {
  spam::spam(0, nrow = 100, ncol = 50)
}

test_that("CheckMatrixSize correctly identifies large matrices", {
  # Test with normal sized matrix
  expect_false(CheckMatrixSize(1000, 1000))
  
  # Test with large matrix that exceeds limits
  expect_true(CheckMatrixSize(50000, 50000))
  
  # Test with large nnz
  expect_true(CheckMatrixSize(1000, 1000, nnz = 2^31))
})

test_that("RowMeanSparse works with spam matrices", {
  if (requireNamespace("spam", quietly = TRUE)) {
    mat <- test_spam_matrix()
    # Add some values
    mat[1:10, 1:10] <- 1:100
    
    means <- RowMeanSparse(mat)
    expect_equal(length(means), nrow(mat))
    expect_true(all(means >= 0))
  }
})

test_that("RowSumSparse works with spam matrices", {
  if (requireNamespace("spam", quietly = TRUE)) {
    mat <- test_spam_matrix()
    # Add some values
    mat[1:10, 1:10] <- 1:100
    
    sums <- RowSumSparse(mat)
    expect_equal(length(sums), nrow(mat))
    expect_true(all(sums >= 0))
  }
})

test_that("as.sparse works with spam matrices", {
  if (requireNamespace("spam", quietly = TRUE)) {
    mat <- test_spam_matrix()
    mat[1:5, 1:5] <- 1:25
    
    sparse_mat <- as.sparse(mat)
    expect_s4_class(sparse_mat, "dgCMatrix")
    expect_equal(dim(sparse_mat), dim(mat))
  }
})

test_that("VST works with spam matrices", {
  if (requireNamespace("spam", quietly = TRUE)) {
    mat <- test_spam_matrix()
    # Add some realistic gene expression values
    for (i in 1:50) {
      for (j in 1:25) {
        if (runif(1) < 0.1) {  # 10% sparsity
          mat[i, j] <- rpois(1, lambda = 2)
        }
      }
    }
    
    hvf_info <- VST(mat, nselect = 10)
    expect_s3_class(hvf_info, "data.frame")
    expect_equal(nrow(hvf_info), nrow(mat))
    expect_true("variable" %in% colnames(hvf_info))
    expect_true(sum(hvf_info$variable) <= 10)
  }
})

test_that("CreateSeuratObject works with spam matrices", {
  if (requireNamespace("spam", quietly = TRUE)) {
    mat <- test_spam_matrix()
    # Add some gene expression values
    for (i in 1:20) {
      for (j in 1:20) {
        if (runif(1) < 0.2) {  # 20% sparsity
          mat[i, j] <- rpois(1, lambda = 1)
        }
      }
    }
    
    # Set row and column names
    rownames(mat) <- paste0("Gene_", 1:nrow(mat))
    colnames(mat) <- paste0("Cell_", 1:ncol(mat))
    
    # This should work without errors
    expect_no_error({
      obj <- CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
    })
    
    obj <- CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
    expect_s4_class(obj, "Seurat")
    expect_equal(ncol(obj), ncol(mat))
  }
})
