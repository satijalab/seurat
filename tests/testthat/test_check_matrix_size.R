set.seed(42)

# Build a small reference dgCMatrix-backed object and a matching DelayedMatrix
.make_counts <- function(nfeat = 120, ncell = 60, lambda = 1.3) {
  m <- as(matrix(rpois(nfeat * ncell, lambda), nfeat, ncell), "dgCMatrix")
  dimnames(m) <- list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  m
}
context("CheckMatrixSize")

test_that("CheckMatrixSize reports counts and the 2^31 limit", {
  m <- .make_counts()
  res <- CheckMatrixSize(m, warn = FALSE)
  expect_equal(res$limit, .Machine$integer.max)
  expect_equal(res$n, length(slot(m, "x")))
  expect_false(res$exceeds)
  expect_false(res$on.disk)
})

test_that("CheckMatrixSize flags on-disk backends as exempt", {
  skip_if_not_installed("DelayedArray")
  d <- as.DelayedMatrix(.make_counts())
  res <- CheckMatrixSize(d, warn = FALSE)
  expect_true(res$on.disk)
})

