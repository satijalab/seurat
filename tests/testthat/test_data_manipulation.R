# Tests for functions in data_manipulation.cpp
# change in random number generation in R3.6, this ensures tests will pass under older and newer Rs
suppressWarnings(RNGversion(vstr = "3.5.3"))
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

#test_that("Row merging with a list done correctly", {
#  m3 <- RowMergeMatricesList(mat_list = list(m1, m2), mat_rownames = list(m1.names, m2.names), all_rownames = all.names)
#  expect_equal(m3[1, 14], -0.17)
#  expect_equal(m3[3, 2], -1.4)
#  expect_equal(m3[14, 18], -0.43)
#  expect_equal(length(m3), 280)
#})

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

# Tests for scaling data
# --------------------------------------------------------------------------------
context("Fast Scale Data Functions")

mat <- matrix(rnorm(n = 10*15), nrow = 10, ncol = 15)

# should be the equivalent of t(scale(t(mat)))
test_that("Fast implementation of row scaling returns expected values", {
  expect_equal(t(scale(t(mat)))[1:10, 1:15], FastRowScale(mat))
  expect_equal(t(scale(t(mat), center = FALSE))[1:10, 1:15],
               FastRowScale(mat, center = FALSE))
  expect_equal(t(scale(t(mat), scale = FALSE))[1:10, 1:15],
               FastRowScale(mat, scale = FALSE))
  expect_equal(t(scale(t(mat), scale = FALSE, center = F))[1:10, 1:15],
               FastRowScale(mat, scale = FALSE, center = F))
  mat.clipped <- FastRowScale(mat, scale_max = 0.2)
  expect_true(max(mat.clipped, na.rm = T) >= 0.2)
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
  mat.clipped <- FastSparseRowScale(mat, scale_max = 0.2, display_progress = F)
  expect_true(max(mat.clipped, na.rm = T) >= 0.2)
})

mat <- as.sparse(x = matrix(rnorm(100), nrow = 10, ncol = 10))

test_that("Row scaling with known stats works", {
  mat.rowmeans <- rowMeans(x = mat)
  mat.sd <- apply(X = mat, MARGIN = 1, FUN = sd)
  expect_equal(
    t(scale(t(as.matrix(mat)), center = mat.rowmeans, scale = mat.sd)),
    FastSparseRowScaleWithKnownStats(mat = mat, mu = mat.rowmeans, sigma = mat.sd, scale = TRUE, center = TRUE, scale_max = 10, display_progress = FALSE),
    check.attributes = FALSE
  )
  expect_equal(
    t(scale(t(as.matrix(mat)), center = FALSE, scale = mat.sd)),
    FastSparseRowScaleWithKnownStats(mat = mat, mu = mat.rowmeans, sigma = mat.sd, scale = TRUE, center = FALSE, scale_max = 10, display_progress = FALSE),
    check.attributes = FALSE
  )
  expect_equal(
    t(scale(t(as.matrix(mat)), center = mat.rowmeans, scale = FALSE)),
    FastSparseRowScaleWithKnownStats(mat = mat, mu = mat.rowmeans, sigma = mat.sd, scale = FALSE, center = TRUE, scale_max = 10, display_progress = FALSE),
    check.attributes = FALSE
  )
  mat.clipped <- FastSparseRowScaleWithKnownStats(mat = mat, mu = mat.rowmeans, sigma = mat.sd, scale = FALSE, center = TRUE, scale_max = 0.2, display_progress = FALSE)
  expect_true(max(mat.clipped, na.rm = T) >= 0.2)
})


# Tests for fast basic stats functions
# --------------------------------------------------------------------------------
context("Fast Basic Stats Functions")

set.seed(42)
mat <- replicate(10, rchisq(10, 4))
fcv <- FastCov(mat)
cv <- cov(mat)
test_that("Fast implementation of covariance returns expected values", {
  expect_equal(fcv[1,1], 9.451051142)
  expect_equal(fcv[10,10], 5.6650068)
  expect_equal(fcv, cv)
})

mat2 <- replicate(10, rchisq(10, 4))
fcv <- FastCovMats(mat1 = mat, mat2 = mat2)
cv <- cov(mat, mat2)
test_that("Fast implementation of covariance returns expected values for matrices", {
  expect_equal(fcv[1,1], 1.523417, tolerance = 1e-6)
  expect_equal(fcv[10,10], -0.6031694, tolerance = 1e-6)
  expect_equal(fcv, cv)
})


merged.mat <- FastRBind(mat, fcv)
test_that("Fast implementation of rbind returns expected values", {
  expect_equal(merged.mat, rbind(mat, fcv))
  expect_equal(mat[1,1], merged.mat[1,1])
  expect_equal(fcv[10,10], merged.mat[20,10])
})

mat <- as.sparse(mat)
test_that("Fast implementation of ExpMean returns expected values",{
  expect_equal(ExpMean(mat[1,]), FastExpMean(mat, display_progress = F)[1])
  expect_equal(ExpMean(mat[5,]), FastExpMean(mat, display_progress = F)[5])
  expect_equal(ExpMean(mat[10,]), FastExpMean(mat, display_progress = F)[10])
  expect_equal(length(FastExpMean(mat, display_progress = F)), nrow(mat))
  expect_error(FastExpMean(mat[1, ], display_progress = F))
  expect_equal(FastExpMean(mat[1, ,drop = F], display_progress = F), ExpMean(mat[1,]))
  expect_equal(FastExpMean(mat, display_progress = F)[1], 6.493418, tolerance = 1e-6)
  expect_equal(FastExpMean(mat, display_progress = F)[5], 6.255206, tolerance = 1e-6)
  expect_equal(FastExpMean(mat, display_progress = F)[10], 7.84965, tolerance = 1e-6)
})
test_that("Fast implementation of LogVMR returns expected values", {
  expect_equal(LogVMR(mat[1,]), FastLogVMR(mat, display_progress = F)[1])
  expect_equal(LogVMR(mat[5,]), FastLogVMR(mat, display_progress = F)[5])
  expect_equal(LogVMR(mat[10,]), FastLogVMR(mat, display_progress = F)[10])
  expect_equal(length(FastExpMean(mat, display_progress = F)), nrow(mat))
  expect_error(FastLogVMR(mat[1, ], display_progress = F))
  expect_equal(FastLogVMR(mat[1, ,drop = F], display_progress = F), LogVMR(mat[1,]))
  expect_equal(FastLogVMR(mat, display_progress = F)[1], 7.615384, tolerance = 1e-6)
  expect_equal(FastLogVMR(mat, display_progress = F)[5], 7.546768, tolerance = 1e-6)
  expect_equal(FastLogVMR(mat, display_progress = F)[10], 10.11755, tolerance = 1e-6)
})

test_that("Row variance calculations for sparse matrices work", {
  expect_equal(apply(X = mat, MARGIN = 1, FUN = var), SparseRowVar(mat = mat, display_progress = FALSE), tolerance = 1e-6)
  expect_equal(apply(X = mat2, MARGIN = 1, FUN = var), SparseRowVar(mat = as.sparse(x = mat2), display_progress = FALSE), tolerance = 1e-6)
})

# Tests for data structure manipulations
# --------------------------------------------------------------------------------
context("Data structure manipulations")

mat <- rsparsematrix(nrow = 10, ncol = 100, density = 0.1)
mat2 <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
cols.to.replace1 <- 1:10
cols.to.replace2 <- 10:1
cols.to.replace3 <- 91:100
cols.to.replace4 <- c(10, 15, 33, 2, 6, 99, 55, 30, 25, 42)

ReplaceCols <- function(mat, cols, replace){
  mat[, cols] <- replace
  return(mat)
}

test_that("Replacing columns works", {
  expect_equal(ReplaceColsC(mat = mat, col_idx = cols.to.replace1 - 1, replacement = mat2),
               ReplaceCols(mat = mat, cols = cols.to.replace1, replace = mat2))
  expect_equal(ReplaceColsC(mat = mat, col_idx = cols.to.replace2 - 1, replacement = mat2),
               ReplaceCols(mat = mat, cols = cols.to.replace2, replace = mat2))
  expect_equal(ReplaceColsC(mat = mat, col_idx = cols.to.replace3 - 1, replacement = mat2),
               ReplaceCols(mat = mat, cols = cols.to.replace3, replace = mat2))
  expect_equal(ReplaceColsC(mat = mat, col_idx = cols.to.replace4 - 1, replacement = mat2),
               ReplaceCols(mat = mat, cols = cols.to.replace4, replace = mat2))
})

test_that("Cpp implementation of row variance is correct", {
  expect_equal(apply(X = mat, MARGIN = 1, FUN = var), RowVar(as.matrix(mat)))
  expect_equal(apply(X = merged.mat, MARGIN = 1, FUN = var), RowVar(as.matrix(merged.mat)))
})
