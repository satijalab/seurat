set.seed(42)

# With split.by and more than one worker, ScaleData assembled the result by
# binding the splits together, which puts the cells in split order, and then
# overwrote the dimnames with the object's own order. Every cell whose position
# moved was given another cell's scaled values.
make_object <- function(nfeat = 40L, ncell = 200L) {
  counts <- as.sparse(matrix(rpois(nfeat * ncell, lambda = 3), nrow = nfeat))
  dimnames(counts) <- list(
    paste0("gene", seq_len(nfeat)),
    paste0("c", seq_len(ncell))
  )
  object <- suppressWarnings(CreateSeuratObject(counts = counts))
  object$group <- rep(c("a", "b"), length.out = ncell)
  suppressWarnings(NormalizeData(object, verbose = FALSE))
}

with_workers <- function(n, expr) {
  old.plan <- future::plan()
  on.exit(future::plan(old.plan), add = TRUE)
  future::plan("multicore", workers = n)
  force(expr)
}

scaled <- function(object, ...) {
  LayerData(suppressWarnings(ScaleData(object, verbose = FALSE, ...)), layer = "scale.data")
}

test_that("split.by scaling gives each cell its own values in parallel", {
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  object <- make_object()
  future::plan("sequential")
  sequential <- scaled(object, split.by = "group")
  parallel <- with_workers(2, scaled(object, split.by = "group"))
  expect_equal(parallel, sequential)

  # each cell's values are its own, computed within its own split
  cells.a <- colnames(object)[object$group == "a"]
  within.a <- scaled(object[, cells.a])
  expect_equal(parallel[, cells.a], within.a[rownames(parallel), cells.a])
})

test_that("split.by regression gives each cell its own values in parallel", {
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  object <- make_object(nfeat = 30L, ncell = 120L)
  future::plan("sequential")
  sequential <- scaled(object, split.by = "group", vars.to.regress = "nCount_RNA")
  parallel <- with_workers(2, scaled(object, split.by = "group", vars.to.regress = "nCount_RNA"))
  expect_equal(parallel, sequential)
})

# the same assembly runs without any workers when variables are regressed out,
# so that path mislabels cells whatever the plan
test_that("split.by regression gives each cell its own values", {
  future::plan("sequential")
  object <- make_object(nfeat = 30L, ncell = 120L)
  result <- scaled(object, split.by = "group", vars.to.regress = "nCount_RNA")
  cells.a <- colnames(object)[object$group == "a"]
  within.a <- scaled(object[, cells.a], vars.to.regress = "nCount_RNA")
  expect_equal(result[, cells.a], within.a[rownames(result), cells.a])
})
