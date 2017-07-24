# Tests for functions in data_manipulation.cpp
set.seed(42)
library(Matrix)

# Tests for row merging 
# --------------------------------------------------------------------------------
context("Row Merging")

m1 <- rsparsematrix(10, 10, 0.1)
m2 <- rsparsematrix(10, 10, 0.1)
m1.names <- paste0("row", sample(1:10, size = 10))
m2.names <- paste0("row", sample(1:20, size = 10))
all.names <- union(m1.names, m2.names)
rownames(m1) <- m1.names
rownames(m2) <- m2.names
m1 <- as(m1, "RsparseMatrix")
m2 <- as(m2, "RsparseMatrix")

test_that("Row merging done correctly", {
  m3 <- RowMergeMatrices(mat1 = m1, mat2 = m2, mat1_rownames = m1.names, mat2_rownames = m2.names, 
                  all_rownames = all.names)
  expect_equal(m3[1, 14], -0.17)
  expect_equal(m3[3, 2], -1.4)
  expect_equal(m3[14, 18], -0.43)
  expect_equal(length(m3), 280)
})

# Tests for log normalization
# --------------------------------------------------------------------------------
context("Log Normalization")

mat <- as(matrix(1:16, ncol = 4, nrow = 4), "sparseMatrix")

test_that("Log Normalization returns expected values", {
  mat.norm.r <- log1p(sweep(mat, 2, Matrix::colSums(mat), FUN = "/") * 1e4)
  mat.norm <- LogNorm(mat, 1e4, display_progress = F)
  expect_equal(mat.norm[1, ], mat.norm.r[1, ])
  expect_equal(mat.norm[4, 4], mat.norm.r[4, 4])
})

# Tests for matrix multiply
# --------------------------------------------------------------------------------
context("Matrix Multiply")

mat <- as.matrix(mat)

test_that("Fast implementation of matrix multiply returns as expected", {
  expect_equal(mat %*% mat, FastMatMult(mat, mat))
  mat[1, 1] <- NA
  expect_equal(mat %*% mat, FastMatMult(mat, mat))
  mat[1, 1] <- NaN
  expect_equal(mat %*% mat, FastMatMult(mat, mat))
})

# Tests for scaling data
# --------------------------------------------------------------------------------
context("Fast Scale Data Functions")

mat <- matrix(seq(0.001, 0.1, 0.001), nrow = 10, ncol = 10)

# should be the equivalent of t(scale(t(mat)))
test_that("Fast implementation of row scaling returns expected values", {
  expect_equal(t(scale(t(mat))[1:10, 1:10]), FastRowScale(mat, display_progress = FALSE))
  expect_equal(t(scale(t(mat), center = FALSE))[1:10, 1:10], 
               FastRowScale(mat, center = FALSE, display_progress = FALSE))
  expect_equal(t(scale(t(mat), scale = FALSE))[1:10, 1:10], 
               FastRowScale(mat, scale = FALSE, display_progress = FALSE))
  expect_equal(t(scale(t(mat), scale = FALSE, center = F))[1:10, 1:10], 
               FastRowScale(mat, scale = FALSE, center = F, display_progress = FALSE))
})

# should be the equivalent of scale(mat, TRUE, apply(mat, 2, sd))
test_that("Standardize returns expected values", {
  expect_equal(Standardize(mat, display_progress = FALSE), scale(mat, TRUE, apply(mat, 2, sd)), 
               check.attributes = FALSE)
})

# should be the equivalent of t(scale(t(mat)))
mat <- rsparsematrix(10, 15, 0.1)
test_that("Fast implementation of row scaling returns expected values", {
  expect_equal(t(scale(t(as.matrix(mat))))[1:10, 1:15], FastSparseRowScale(mat, display_progress = FALSE), 
               check.attributes = FALSE)
  expect_equal(t(scale(t(as.matrix(mat)), center = FALSE))[1:10, 1:15], 
               FastSparseRowScale(mat, center = FALSE, display_progress = FALSE), 
               check.attributes = FALSE)
  expect_equal(t(scale(t(as.matrix(mat)), scale = FALSE))[1:10, 1:15], 
               FastSparseRowScale(mat, scale = FALSE, display_progress = FALSE),
               check.attributes = FALSE)
  expect_equal(t(scale(t(as.matrix(mat)), scale = FALSE, center = F))[1:10, 1:15], 
               FastSparseRowScale(mat, scale = FALSE, center = F, display_progress = FALSE),
               check.attributes = FALSE)
})

# Tests for fast basic stats functions
# --------------------------------------------------------------------------------
set.seed(42)
mat <- replicate(10, rchisq(10, 4))
fcv <- FastCov(mat)
cv <- cov(mat)
test_that("Fast implementation of covariance returns expected values", {
  expect_equal(fcv[1,1], 9.451051142)
  expect_equal(fcv[10,10], 5.6650068)
  expect_equal(fcv, cv)
})

merged.mat <- FastRBind(mat, fcv)
test_that("Fast implementation of rbind returns expected values", {
  expect_equal(merged.mat, rbind(mat, fcv))
  expect_equal(mat[1,1], merged.mat[1,1])
  expect_equal(fcv[10,10], merged.mat[20,10])
})
