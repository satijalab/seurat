# Ported from scSLIDE. Uses testthat edition 3 (local_mocked_bindings etc.).
testthat::local_edition(3)

# Tests for RunPLS() argument dispatch on the in-memory (default) method.

# Helper: run RunPLS() quietly; the matrix method warns about the missing assay
# and messages the model R2.
run_pls <- function(...) {
  suppressWarnings(suppressMessages(RunPLS(...)))
}

# Helper: small features x cells matrix plus a single response.
pls_fixture <- function(n.genes = 20, n.cells = 10, seed = 1) {
  set.seed(seed)
  X <- matrix(
    rnorm(n = n.genes * n.cells), nrow = n.genes, ncol = n.cells,
    dimnames = list(paste0("g", seq_len(n.genes)), paste0("c", seq_len(n.cells)))
  )
  list(X = X, Y = rnorm(n = n.cells))
}


test_that("RunPLS.default resolves the pls.function default to plsr", {
  skip_if_not_installed("pls")

  fx <- pls_fixture()
  # an unresolved c('plsr', 'spls', 'cppls') default must not reach the body
  default <- run_pls(fx$X, Y = fx$Y, ncomp = 2, verbose = FALSE)
  plsr_    <- run_pls(fx$X, Y = fx$Y, ncomp = 2, pls.function = "plsr", verbose = FALSE)
  cppls_   <- run_pls(fx$X, Y = fx$Y, ncomp = 2, pls.function = "cppls", verbose = FALSE)

  expect_s4_class(default, "DimReduc")
  expect_equal(Embeddings(default), Embeddings(plsr_))
  # guards against the default silently resolving to the last choice instead
  expect_false(isTRUE(all.equal(Embeddings(default), Embeddings(cppls_))))
})


test_that("RunPLS rejects an unknown pls.function", {
  fx <- pls_fixture()

  expect_error(
    run_pls(fx$X, Y = fx$Y, ncomp = 2, pls.function = "not-a-method"),
    "should be one of"
  )
})


test_that("RunPLS.default runs every supported pls.function", {
  skip_if_not_installed("pls")
  skip_if_not_installed("spls")

  fx <- pls_fixture()

  for (fn in c("plsr", "spls", "cppls")) {
    expect_s4_class(
      run_pls(fx$X, Y = fx$Y, ncomp = 2, pls.function = fn, verbose = FALSE),
      "DimReduc"
    )
  }
})
