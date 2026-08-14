# Ported from scSLIDE. Uses testthat edition 3 (local_mocked_bindings etc.).
testthat::local_edition(3)

# Tests for CPLS (Canonical PLS) implementation

# Helper: convert n x p dense matrix to BPCells IterableMatrix in features x samples layout.
to_bp <- function(X_np) {
  as(as(t(X_np), "dgCMatrix"), "IterableMatrix")
}

# Helper: compare BPCells cppls_ondisk result to pls::cppls.fit reference,
# tolerating sign flips in per-component vectors.
compare_cppls_results <- function(bp, ref, ncomp, tol = 1e-6) {
  # Coefficients should match exactly (no sign ambiguity in cumulative B)
  expect_equal(bp$coefficients, ref$coefficients, tolerance = tol)

  # Xvar should match
  expect_equal(as.numeric(bp$Xvar), as.numeric(ref$Xvar), tolerance = tol)

  # Xtotvar should match
  expect_equal(bp$Xtotvar, ref$Xtotvar, tolerance = tol)

  # Scores, loadings, etc. may have per-component sign flips.
  for (a in seq_len(ncomp)) {
    # Determine sign from first non-negligible loading weight element
    bp_w <- bp$loading.weights[, a]
    ref_w <- ref$loading.weights[, a]
    idx <- which.max(abs(ref_w))
    s <- sign(bp_w[idx]) * sign(ref_w[idx])

    expect_equal(s * bp$scores[, a], ref$scores[, a], tolerance = tol,
                 label = sprintf("scores component %d", a))
    expect_equal(s * bp$loadings[, a], ref$loadings[, a], tolerance = tol,
                 label = sprintf("loadings component %d", a))
    expect_equal(s * bp$loading.weights[, a], ref$loading.weights[, a], tolerance = tol,
                 label = sprintf("loading.weights component %d", a))
    expect_equal(s * bp$projection[, a], ref$projection[, a], tolerance = tol,
                 label = sprintf("projection component %d", a))
  }

  # Fitted values and residuals (no sign ambiguity since they're cumulative)
  expect_equal(bp$fitted.values, ref$fitted.values, tolerance = tol)
  expect_equal(bp$residuals, ref$residuals, tolerance = tol)
}


test_that("CPLS matches pls package: single response, no Y.add", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(42)
  n <- 80
  p <- 15
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  beta <- rnorm(p)
  Y <- X %*% beta + rnorm(n, sd = 0.5)
  Y <- as.matrix(Y)

  ref <- pls::cppls.fit(X, Y, ncomp = ncomp)

  X_bp <- to_bp(X)
  bp  <- cppls_ondisk(X_bp, Y, ncomp = ncomp)

  expect_equal(bp$Xmeans, ref$Xmeans, tolerance = 1e-10)
  expect_equal(bp$Ymeans, ref$Ymeans, tolerance = 1e-10)
  compare_cppls_results(bp, ref, ncomp)
})


test_that("CPLS matches pls package: multiple responses, no Y.add", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(123)
  n <- 100
  p <- 20
  q <- 3
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  ref <- pls::cppls.fit(X, Y, ncomp = ncomp)

  X_bp <- to_bp(X)
  bp  <- cppls_ondisk(X_bp, Y, ncomp = ncomp)

  compare_cppls_results(bp, ref, ncomp)
})


test_that("CPLS matches pls package: single response + Y.add", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(314)
  n <- 80
  p <- 15
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  beta <- rnorm(p)
  Y <- X %*% beta + rnorm(n, sd = 0.5)
  Y <- as.matrix(Y)
  Y_add <- matrix(rnorm(n * 2), n, 2)

  ref <- pls::cppls.fit(X, Y, ncomp = ncomp, Y.add = Y_add)

  X_bp <- to_bp(X)
  bp  <- cppls_ondisk(X_bp, Y, Y.add = Y_add, ncomp = ncomp)

  expect_equal(bp$Xmeans, ref$Xmeans, tolerance = 1e-10)
  expect_equal(bp$Ymeans, ref$Ymeans, tolerance = 1e-10)
  compare_cppls_results(bp, ref, ncomp)
})


test_that("CPLS matches pls package: multi-response + Y.add", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(271)
  n <- 100
  p <- 20
  q <- 2
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  Y_add <- matrix(rnorm(n * 3), n, 3)

  ref <- pls::cppls.fit(X, Y, ncomp = ncomp, Y.add = Y_add)

  X_bp <- to_bp(X)
  bp  <- cppls_ondisk(X_bp, Y, Y.add = Y_add, ncomp = ncomp)

  compare_cppls_results(bp, ref, ncomp)
})


test_that("CPLS works with ncomp = 1", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(7)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  ref <- pls::cppls.fit(X, Y, ncomp = 1)

  X_bp <- to_bp(X)
  bp  <- cppls_ondisk(X_bp, Y, ncomp = 1)

  expect_equal(bp$coefficients, ref$coefficients, tolerance = 1e-6)
  expect_equal(as.numeric(bp$Xvar), as.numeric(ref$Xvar), tolerance = 1e-6)
})


test_that("CPLS centering works correctly", {
  skip_if_not_installed("BPCells")

  set.seed(789)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p, mean = 10), n, p)
  Y <- matrix(rnorm(n, mean = 5), n, 1)

  X_bp <- to_bp(X)

  result_center    <- cppls_ondisk(X_bp, Y, ncomp = 3, center = TRUE)
  result_no_center <- cppls_ondisk(X_bp, Y, ncomp = 3, center = FALSE)

  # With non-zero means, centering should change the results
  expect_false(isTRUE(all.equal(result_center$coefficients,
                                result_no_center$coefficients)))

  # Xmeans should be zero when center = FALSE
  expect_true(all(abs(result_no_center$Xmeans) < 1e-10))
  expect_true(all(abs(result_no_center$Ymeans) < 1e-10))

  # Xmeans should approximate column means when center = TRUE
  expect_equal(as.numeric(result_center$Xmeans), colMeans(X), tolerance = 1e-10)
})


test_that("CPLS stripped mode returns minimal output", {
  skip_if_not_installed("BPCells")

  set.seed(101)
  n <- 30
  p <- 8

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)
  X_bp <- to_bp(X)

  result <- cppls_ondisk(X_bp, Y, ncomp = 3, stripped = TRUE)

  expect_true("coefficients" %in% names(result))
  expect_true("Xmeans" %in% names(result))
  expect_true("Ymeans" %in% names(result))
  expect_true("gammas" %in% names(result))
  expect_true("canonical.correlations" %in% names(result))
  expect_true("A" %in% names(result))
  expect_false("scores" %in% names(result))
  expect_false("loadings" %in% names(result))
  expect_false("fitted.values" %in% names(result))
})


test_that("CPLS stripped coefficients match full coefficients", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(55)
  n <- 60
  p <- 12
  ncomp <- 4

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)
  X_bp <- to_bp(X)

  full     <- cppls_ondisk(X_bp, Y, ncomp = ncomp, stripped = FALSE)
  stripped <- cppls_ondisk(X_bp, Y, ncomp = ncomp, stripped = TRUE)

  expect_equal(stripped$coefficients, full$coefficients, tolerance = 1e-10)
  expect_equal(stripped$Xmeans, full$Xmeans, tolerance = 1e-10)
  expect_equal(stripped$Ymeans, full$Ymeans, tolerance = 1e-10)
})


test_that("CPLS input validation works", {
  skip_if_not_installed("BPCells")

  X_mat <- matrix(rnorm(30), 10, 3)
  X_bp <- to_bp(X_mat)
  Y <- matrix(rnorm(10), 10, 1)
  Y_wrong <- matrix(rnorm(5), 5, 1)

  # Dimension mismatch (X has 10 samples, Y_wrong has 5 rows)
  expect_error(cppls_ondisk(X_bp, Y_wrong, ncomp = 1), "must match")

  # ncomp too large: min(n=10, p=3) = 3, so ncomp >= 3 should fail
  expect_error(cppls_ondisk(X_bp, Y, ncomp = 3), "ncomp must be less than")
  expect_error(cppls_ondisk(X_bp, Y, ncomp = 10), "ncomp must be less than")
})


test_that("CPLS works on sparse BPCells matrix", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(456)
  n <- 200
  p <- 30
  ncomp <- 5

  # Create a sparse matrix (many zeros)
  X <- matrix(rnorm(n * p), n, p)
  X[abs(X) < 1] <- 0  # sparsify
  Y <- matrix(rnorm(n * 2), n, 2)

  ref <- pls::cppls.fit(X, Y, ncomp = ncomp)

  X_bp <- to_bp(X)
  bp  <- cppls_ondisk(X_bp, Y, ncomp = ncomp)

  compare_cppls_results(bp, ref, ncomp, tol = 1e-5)
})


test_that("CPLS multi-threaded results match single-threaded", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(42)
  n <- 80
  p <- 15
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  beta <- rnorm(p)
  Y <- X %*% beta + rnorm(n, sd = 0.5)
  Y <- as.matrix(Y)

  X_bp <- to_bp(X)
  res1 <- cppls_ondisk(X_bp, Y, ncomp = ncomp, threads = 1L)
  res2 <- cppls_ondisk(X_bp, Y, ncomp = ncomp, threads = 2L)

  expect_equal(res2$coefficients, res1$coefficients, tolerance = 1e-10)
  expect_equal(res2$scores, res1$scores, tolerance = 1e-10)
  expect_equal(res2$loadings, res1$loadings, tolerance = 1e-10)
  expect_equal(res2$fitted.values, res1$fitted.values, tolerance = 1e-10)
  expect_equal(as.numeric(res2$Xvar), as.numeric(res1$Xvar), tolerance = 1e-10)
})


test_that("CPLS output dimensions are correct", {
  skip_if_not_installed("BPCells")

  set.seed(999)
  n <- 40
  p <- 8
  q <- 2
  ncomp <- 3

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  X_bp <- to_bp(X)

  result <- cppls_ondisk(X_bp, Y, ncomp = ncomp)

  expect_equal(dim(result$coefficients), c(p, q, ncomp))
  expect_equal(dim(result$scores), c(n, ncomp))
  expect_equal(dim(result$loadings), c(p, ncomp))
  expect_equal(dim(result$loading.weights), c(p, ncomp))
  expect_equal(dim(result$Yscores), c(n, ncomp))
  expect_equal(dim(result$Yloadings), c(q, ncomp))
  expect_equal(dim(result$projection), c(p, ncomp))
  expect_equal(length(result$Xmeans), p)
  expect_equal(length(result$Ymeans), q)
  expect_equal(dim(result$fitted.values), c(n, q, ncomp))
  expect_equal(dim(result$residuals), c(n, q, ncomp))
  expect_equal(length(result$Xvar), ncomp)
  expect_true(is.numeric(result$Xtotvar))
  # CPLS-specific
  expect_equal(length(result$gammas), ncomp)
  expect_equal(length(result$canonical.correlations), ncomp)
  expect_true(all(result$gammas == 0.5))
})


test_that("CPLS output dimensions correct with Y.add", {
  skip_if_not_installed("BPCells")

  set.seed(777)
  n <- 40
  p <- 8
  q <- 2
  q_add <- 3
  ncomp <- 3

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  Y_add <- matrix(rnorm(n * q_add), n, q_add)
  X_bp <- to_bp(X)

  result <- cppls_ondisk(X_bp, Y, Y.add = Y_add, ncomp = ncomp)

  # Coefficients/loadings use primary Y only
  expect_equal(dim(result$coefficients), c(p, q, ncomp))
  expect_equal(dim(result$Yloadings), c(q, ncomp))
  # A weights span the full Y space (q + q_add)
  expect_equal(nrow(result$A), q + q_add)
  expect_equal(ncol(result$A), ncomp)
})
