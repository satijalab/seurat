testthat::local_edition(3)

# Tests for FindmmNN(): the multi-modal nearest-neighbor wrapper.

# Helper: a Seurat object with a LANDMARK assay and two reductions whose
# geometries genuinely differ. Two copies of the same embedding would make the
# weighted distance a monotone transform of a single modality, and the tests
# below could not then detect a change in kernel width or power.
make_multimodal_fixture <- function(n.genes = 40, n.cells = 120, n.landmark = 40,
                                    n.dims = 6, seed = 11) {
  set.seed(seed)
  cts <- matrix(
    rpois(n = n.genes * n.cells, lambda = 5),
    nrow = n.genes, ncol = n.cells,
    dimnames = list(paste0("g", seq_len(n.genes)), paste0("cell", seq_len(n.cells)))
  )
  obj <- suppressWarnings(CreateSeuratObject(counts = cts))
  obj[["LANDMARK"]] <- suppressWarnings(
    CreateAssayObject(counts = cts[, seq_len(n.landmark), drop = FALSE])
  )
  for (i in seq_len(2)) {
    emb <- matrix(
      rnorm(n = n.cells * n.dims), nrow = n.cells, ncol = n.dims,
      dimnames = list(colnames(obj), paste0("r", i, "_", seq_len(n.dims)))
    )
    obj[[paste0("red", i)]] <- SeuratObject::CreateDimReducObject(
      embeddings = emb, key = paste0("R", i, "_"), assay = "RNA"
    )
  }
  obj
}

run_mmnn <- function(obj, ...) {
  result <- suppressMessages(suppressWarnings(FindmmNN(
    obj,
    sketch.assay = "LANDMARK",
    reduction.list = list("red1", "red2"),
    dims.list = list(1:6, 1:6),
    k.nn = 20,
    knn.range = 30,
    verbose = FALSE,
    ...
  )))
  result[["weighted.nn"]]
}

test_that("FindmmNN forwards dots to whichever internal declares them", {
  # FindModalityWeights() and MultiModalNN() accept no dots of their own, so
  # forwarding this function's dots wholesale made every extra argument an
  # error, including the ones each internal does declare.
  obj <- make_multimodal_fixture()
  baseline <- run_mmnn(obj)

  # sigma.idx belongs to FindModalityWeights, kernel.power to MultiModalNN.
  expect_no_error(weights.arg <- run_mmnn(obj, sigma.idx = 15))
  expect_no_error(nn.arg <- run_mmnn(obj, kernel.power = 2))
  expect_no_error(both <- run_mmnn(obj, sigma.idx = 15, kernel.power = 2))

  # Forwarded, not quietly discarded.
  expect_false(identical(Distances(baseline), Distances(weights.arg)))
  expect_false(identical(Distances(baseline), Distances(nn.arg)))
  expect_false(identical(Distances(both), Distances(weights.arg)))
  expect_false(identical(Distances(both), Distances(nn.arg)))
})

test_that("FindmmNN reports arguments that belong to neither internal", {
  # Called directly rather than through run_mmnn(), which suppresses warnings.
  obj <- make_multimodal_fixture()
  expect_warning(
    suppressMessages(FindmmNN(
      obj,
      sketch.assay = "LANDMARK",
      reduction.list = list("red1", "red2"),
      dims.list = list(1:6, 1:6),
      k.nn = 20,
      knn.range = 30,
      verbose = FALSE,
      not_a_real_argument = 1
    )),
    "not used"
  )
})

test_that("FindmmNN does not reach its internals through do.call", {
  # LogSeuratCommand() deparses the calling frame. Under do.call() that frame
  # holds the function body and every argument by value, so deparsing walks the
  # entire Seurat object: 6-7x slower at 30k cells and effectively a hang on a
  # real dataset, stalling right after "Calculating sigma by Nth neighbor".
  # See satijalab/seurat#10153 for the same failure in SCTransform. A small
  # fixture cannot catch this by timing, so assert the call shape instead.
  body.text <- paste(deparse(body(FindmmNN)), collapse = " ")
  expect_false(grepl(pattern = "do.call", x = body.text, fixed = TRUE))
})
