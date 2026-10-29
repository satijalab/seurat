testthat::local_edition(3)

# Tests for partial least squares: the on-disk kernelpls() and cppls_ondisk()
# fits on BPCells matrices, compared with the pls package, and the argument
# dispatch of RunPLS().

# Helper: convert an n x p dense matrix to a BPCells IterableMatrix in the
# features x samples layout. The on-disk fits transpose back to n x p.
to_bp <- function(X_np) {
  as(as(t(X_np), "dgCMatrix"), "IterableMatrix")
}

# Helper: simulate an n x p predictor matrix and n x q responses (plus n x q.add
# additional responses for cppls). sparse zeroes most entries, so the fits'
# in-memory centering corrections are exercised on implicit zeros; shift moves
# the means away from zero, so centering changes the fit.
simulate_pls_data <- function(n, p, q, q.add = 0, sparse = FALSE, shift = 0,
                              seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p, mean = shift), n, p)
  if (sparse) X[abs(X) < 1] <- 0
  Y <- matrix(rnorm(n * q, mean = shift), n, q)
  Y.add <- if (q.add > 0) matrix(rnorm(n * q.add), n, q.add)
  list(X = X, Y = Y, Y.add = Y.add)
}

# Helper: pls 2.9-0 started centering the responses before the canonical
# correlation step, which we follow. Only cancorr() is sensitive to it, and only
# when more than one response column reaches it, so the multi-response
# comparisons disagree with pls < 2.9-0 and are skipped there.
skip_if_pls_uncentered <- function() {
  skip_if(
    packageVersion("pls") < "2.9-0",
    "pls < 2.9-0 does not center the responses before cancorr()"
  )
}

# Helper: compare an on-disk fit to the pls reference, tolerating per-component
# sign flips. Yscores are not compared: for components after the first,
# pls::cppls.fit also projects u_a off its own score t_a, which ours does not.
# Xtotvar is only compared for centered fits: ours is always the total variance
# of the centered X, whereas pls uses sum(X^2) when it does not center.
compare_to_pls <- function(bp, ref, ncomp, center = TRUE, tol = 1e-6) {
  # Cumulative outputs carry no sign ambiguity
  expect_equal(bp$coefficients, ref$coefficients, tolerance = tol)
  expect_equal(bp$fitted.values, ref$fitted.values, tolerance = tol)
  expect_equal(bp$residuals, ref$residuals, tolerance = tol)
  expect_equal(as.numeric(bp$Xvar), as.numeric(ref$Xvar), tolerance = tol)
  if (center) {
    expect_equal(bp$Xtotvar, ref$Xtotvar, tolerance = tol)
  }
  if (!is.null(ref$canonical.correlations)) {
    expect_equal(bp$canonical.correlations, ref$canonical.correlations,
                 tolerance = tol)
  }

  for (a in seq_len(ncomp)) {
    # Determine the sign from the largest loading weight
    idx <- which.max(abs(ref$loading.weights[, a]))
    s <- sign(bp$loading.weights[idx, a]) * sign(ref$loading.weights[idx, a])

    for (out in c("scores", "loadings", "loading.weights", "projection",
                  "Yloadings", "A")) {
      if (is.null(ref[[out]])) next
      expect_equal(s * unclass(bp[[out]])[, a], unclass(ref[[out]])[, a],
                   tolerance = tol,
                   label = sprintf("%s component %d", out, a))
    }
  }
}

# The on-disk fits, for the tests that apply to both
pls_fits <- list(kernelpls = kernelpls, cppls_ondisk = cppls_ondisk)


kernelpls.cases <- list(
  "one response" = list(n = 80, p = 15, q = 1, ncomp = 5, seed = 42),
  "three responses" = list(n = 100, p = 20, q = 3, ncomp = 5, seed = 123),
  "ncomp = 1" = list(n = 50, p = 10, q = 1, ncomp = 1, seed = 7),
  # past six components the projection step switches to a matrix product
  "more than six components" = list(n = 120, p = 25, q = 1, ncomp = 12, seed = 314),
  "sparse X" = list(n = 200, p = 30, q = 2, ncomp = 5, sparse = TRUE,
                    seed = 456, tol = 1e-5),
  "center = FALSE" = list(n = 50, p = 10, q = 1, ncomp = 3, shift = 10,
                          center = FALSE, seed = 789)
)

for (case in names(kernelpls.cases)) {
  test_that(paste("on-disk kernelpls matches pls::kernelpls.fit:", case), {
    skip_if_not_installed("BPCells")
    skip_if_not_installed("pls")

    cs <- kernelpls.cases[[case]]
    center <- cs$center %||% TRUE
    d <- simulate_pls_data(n = cs$n, p = cs$p, q = cs$q,
                           sparse = cs$sparse %||% FALSE,
                           shift = cs$shift %||% 0, seed = cs$seed)

    ref <- pls::kernelpls.fit(d$X, d$Y, ncomp = cs$ncomp, center = center)
    bp <- kernelpls(to_bp(d$X), d$Y, ncomp = cs$ncomp, center = center)

    if (center) {
      expect_equal(as.numeric(bp$Xmeans), colMeans(d$X), tolerance = 1e-10)
      expect_equal(as.numeric(bp$Ymeans), colMeans(d$Y), tolerance = 1e-10)
    }
    compare_to_pls(bp, ref, ncomp = cs$ncomp, center = center,
                   tol = cs$tol %||% 1e-6)
  })
}


cppls.cases <- list(
  "one response" = list(n = 80, p = 15, q = 1, ncomp = 5, seed = 42),
  "several responses" = list(n = 100, p = 20, q = 3, ncomp = 5, seed = 123,
                             multi = TRUE),
  "one response with Y.add" = list(n = 80, p = 15, q = 1, q.add = 2,
                                   ncomp = 5, seed = 314),
  "several responses with Y.add" = list(n = 100, p = 20, q = 2, q.add = 3,
                                        ncomp = 5, seed = 271, multi = TRUE),
  "ncomp = 1" = list(n = 50, p = 10, q = 1, ncomp = 1, seed = 7),
  "sparse X" = list(n = 200, p = 30, q = 2, ncomp = 5, sparse = TRUE,
                    seed = 456, tol = 1e-5, multi = TRUE),
  "center = FALSE" = list(n = 50, p = 10, q = 1, ncomp = 3, shift = 10,
                          center = FALSE, seed = 789)
)

for (case in names(cppls.cases)) {
  test_that(paste("on-disk cppls matches pls::cppls.fit:", case), {
    skip_if_not_installed("BPCells")
    skip_if_not_installed("pls")
    cs <- cppls.cases[[case]]
    if (isTRUE(cs$multi)) skip_if_pls_uncentered()

    center <- cs$center %||% TRUE
    d <- simulate_pls_data(n = cs$n, p = cs$p, q = cs$q, q.add = cs$q.add %||% 0,
                           sparse = cs$sparse %||% FALSE,
                           shift = cs$shift %||% 0, seed = cs$seed)

    ref <- pls::cppls.fit(d$X, d$Y, ncomp = cs$ncomp, Y.add = d$Y.add,
                          center = center)
    bp <- cppls_ondisk(to_bp(d$X), d$Y, Y.add = d$Y.add, ncomp = cs$ncomp,
                       center = center)

    if (center) {
      expect_equal(as.numeric(bp$Xmeans), colMeans(d$X), tolerance = 1e-10)
      expect_equal(as.numeric(bp$Ymeans), colMeans(d$Y), tolerance = 1e-10)
    }
    expect_true(all(bp$gammas == 0.5))
    compare_to_pls(bp, ref, ncomp = cs$ncomp, center = center,
                   tol = cs$tol %||% 1e-6)
  })
}


test_that("on-disk PLS stripped mode drops the per-sample outputs only", {
  skip_if_not_installed("BPCells")

  d <- simulate_pls_data(n = 60, p = 12, q = 1, seed = 55)
  X_bp <- to_bp(d$X)

  for (fn in names(pls_fits)) {
    full <- pls_fits[[fn]](X_bp, d$Y, ncomp = 4, stripped = FALSE)
    stripped <- pls_fits[[fn]](X_bp, d$Y, ncomp = 4, stripped = TRUE)

    expect_false(any(c("scores", "loadings", "fitted.values") %in% names(stripped)),
                 label = fn)
    kept <- c("coefficients", "Xmeans", "Ymeans")
    if (fn == "cppls_ondisk") kept <- c(kept, "gammas", "canonical.correlations", "A")
    expect_equal(stripped[kept], full[kept], tolerance = 1e-10, label = fn)
  }
})


test_that("on-disk PLS results do not depend on the thread count", {
  skip_if_not_installed("BPCells")

  d <- simulate_pls_data(n = 80, p = 15, q = 1, seed = 42)
  X_bp <- to_bp(d$X)
  compared <- c("coefficients", "scores", "loadings", "fitted.values", "Xvar")

  for (fn in names(pls_fits)) {
    res1 <- pls_fits[[fn]](X_bp, d$Y, ncomp = 5, threads = 1L)
    res2 <- pls_fits[[fn]](X_bp, d$Y, ncomp = 5, threads = 2L)
    expect_equal(res2[compared], res1[compared], tolerance = 1e-10, label = fn)

    # threads = 0 detects the core count
    expect_no_error(pls_fits[[fn]](X_bp, d$Y, ncomp = 2, threads = 0L))
  }
})


test_that("on-disk PLS validates its inputs", {
  skip_if_not_installed("BPCells")

  # X is 10 x 3 (n x p), stored as a 3 x 10 BPCells matrix
  X_bp <- to_bp(matrix(rnorm(30), 10, 3))
  Y <- matrix(rnorm(10), 10, 1)
  Y_wrong <- matrix(rnorm(5), 5, 1)

  for (fn in names(pls_fits)) {
    f <- pls_fits[[fn]]
    expect_error(f(X_bp, Y_wrong, ncomp = 1), "must match")
    # ncomp must be below min(n = 10, p = 3)
    expect_error(f(X_bp, Y, ncomp = 3), "ncomp must be less than")
    expect_error(f(X_bp, Y, ncomp = 10), "ncomp must be less than")
    expect_error(f(X_bp, Y, ncomp = 1, threads = -1L),
                 "threads must be a non-negative integer")
  }
})


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
  # cppls reduces to plsr for a single response, so the guard below needs
  # several responses to tell the two of them apart
  set.seed(7)
  Y <- matrix(
    rnorm(n = 3 * length(x = fx$Y)), nrow = length(x = fx$Y), ncol = 3
  )
  # an unresolved c('plsr', 'spls', 'cppls') default must not reach the body
  default <- run_pls(fx$X, Y = Y, ncomp = 2, verbose = FALSE)
  plsr_    <- run_pls(fx$X, Y = Y, ncomp = 2, pls.function = "plsr", verbose = FALSE)
  cppls_   <- run_pls(fx$X, Y = Y, ncomp = 2, pls.function = "cppls", verbose = FALSE)

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
