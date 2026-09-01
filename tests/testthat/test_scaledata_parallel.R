set.seed(42)

# ScaleData chunks the data across workers, but the worker functions used to
# reference `object`, so `future` exported the whole matrix to every worker
# instead of each worker receiving its own block. That is what runs into
# future.globals.maxSize, and what makes every worker hold the whole dataset.
make_object <- function(nfeat = 200L, ncell = 2000L) {
  counts <- as.sparse(matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat))
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

test_that("parallel scaling matches sequential", {
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  object <- make_object()
  future::plan("sequential")
  sequential <- scaled(object)
  parallel <- with_workers(4, scaled(object))
  expect_equal(parallel, sequential)

  # and with the data split, which is the other chunking dimension
  future::plan("sequential")
  sequential.split <- scaled(object, split.by = "group")
  parallel.split <- with_workers(4, scaled(object, split.by = "group"))
  expect_equal(parallel.split, sequential.split)
})

test_that("parallel regression matches sequential", {
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  object <- make_object(nfeat = 60L, ncell = 400L)
  future::plan("sequential")
  sequential <- scaled(object, vars.to.regress = "nCount_RNA")
  parallel <- with_workers(4, scaled(object, vars.to.regress = "nCount_RNA"))
  expect_equal(parallel, sequential)
})

# The size future reports when the limit is exceeded, in bytes
exported_size <- function(object, workers = 4, ...) {
  old.opt <- options(future.globals.maxSize = 1024)
  on.exit(options(old.opt), add = TRUE)
  message <- tryCatch(
    with_workers(workers, scaled(object, ...)),
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
  small.data <- as.numeric(object.size(LayerData(small, layer = "data")))
  large.data <- as.numeric(object.size(LayerData(large, layer = "data")))
  expect_gt(large.data, small.data * 4)

  small.exported <- exported_size(small)
  large.exported <- exported_size(large)
  expect_false(is.na(small.exported))
  expect_false(is.na(large.exported))
  # the worker function is a constant; the data travels as chunks, one block
  # per worker, rather than as a global carrying the entire matrix
  expect_equal(large.exported, small.exported)
  expect_lt(large.exported, large.data)
})

test_that("the sequential path is unchanged", {
  future::plan("sequential")
  object <- make_object(nfeat = 50L, ncell = 100L)
  data <- as.matrix(LayerData(object, layer = "data"))
  expected <- t(scale(t(data)))
  expected[expected > 10] <- 10
  # scale() leaves its centre and scale behind as attributes
  expect_equal(as.vector(scaled(object)), as.vector(expected), tolerance = 1e-6)
})
