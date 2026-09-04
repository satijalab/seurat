set.seed(42)

# NormalizeData chunks the data across workers, but the worker function used to
# reference `object`, so `future` exported the whole matrix to every worker
# instead of each worker receiving its own block.
make_object <- function(nfeat = 200L, ncell = 2000L) {
  counts <- as.sparse(matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat))
  dimnames(counts) <- list(
    paste0("gene", seq_len(nfeat)),
    paste0("c", seq_len(ncell))
  )
  suppressWarnings(CreateSeuratObject(counts = counts))
}

with_workers <- function(n, expr) {
  old.plan <- future::plan()
  on.exit(future::plan(old.plan), add = TRUE)
  future::plan("multicore", workers = n)
  force(expr)
}

test_that("parallel normalization matches sequential", {
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  object <- make_object()
  for (method in c("LogNormalize", "CLR", "RC")) {
    future::plan("sequential")
    sequential <- LayerData(
      suppressWarnings(NormalizeData(object, normalization.method = method, verbose = FALSE)),
      layer = "data"
    )
    parallel <- with_workers(4, LayerData(
      suppressWarnings(NormalizeData(object, normalization.method = method, verbose = FALSE)),
      layer = "data"
    ))
    expect_equal(as.matrix(parallel), as.matrix(sequential), info = method)
  }
})

# The size future reports when the limit is exceeded, in bytes
exported_size <- function(object, workers = 4) {
  old.opt <- options(future.globals.maxSize = 1024)
  on.exit(options(old.opt), add = TRUE)
  message <- tryCatch(
    with_workers(workers, suppressWarnings(NormalizeData(object, verbose = FALSE))),
    error = function(e) conditionMessage(e)
  )
  matched <- regmatches(message, regexpr("is [0-9.]+ [KMG]iB", message))
  if (!length(matched)) {
    return(NA_real_)
  }
  parts <- strsplit(sub("^is ", "", matched), " ")[[1]]
  as.numeric(parts[1]) * switch(parts[2], KiB = 1024, MiB = 1024^2, GiB = 1024^3)
}

test_that("what is sent to workers does not grow with the data", {
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  small <- make_object(nfeat = 200L, ncell = 1000L)
  large <- make_object(nfeat = 200L, ncell = 8000L)
  small.counts <- as.numeric(object.size(LayerData(small, layer = "counts")))
  large.counts <- as.numeric(object.size(LayerData(large, layer = "counts")))
  expect_gt(large.counts, small.counts * 4)

  small.exported <- exported_size(small)
  large.exported <- exported_size(large)
  expect_false(is.na(small.exported))
  expect_false(is.na(large.exported))
  # the worker function is a constant; the data travels as chunks, one block
  # per worker, rather than as a global carrying the entire matrix
  expect_equal(large.exported, small.exported)
  expect_lt(large.exported, large.counts)
})

test_that("the sequential path is unchanged", {
  future::plan("sequential")
  object <- make_object(nfeat = 50L, ncell = 100L)
  normalized <- suppressWarnings(NormalizeData(object, verbose = FALSE))
  counts <- as.matrix(LayerData(object, layer = "counts"))
  expected <- log1p(sweep(counts, 2, colSums(counts), "/") * 1e4)
  expect_equal(as.matrix(LayerData(normalized, layer = "data")), expected)
})
