# Ported from scSLIDE. Uses testthat edition 3 (local_mocked_bindings etc.).
testthat::local_edition(3)

# Tests for FindsmNN(): the single-modality nearest-neighbor wrapper.

# Helper: a minimal Seurat object carrying a LANDMARK assay and one reduction,
# i.e. the inputs FindsmNN() needs. The first n.landmark cells are the landmarks.
make_reduction_fixture <- function(n.genes = 40, n.cells = 60, n.landmark = 12,
                                   n.dims = 8, seed = 7) {
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
  emb <- matrix(
    rnorm(n = n.cells * n.dims), nrow = n.cells, ncol = n.dims,
    dimnames = list(colnames(obj), paste0("d", seq_len(n.dims)))
  )
  obj[["pls"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb, key = "PLS_", assay = "RNA"
  )
  obj
}


test_that("FindsmNN stores a Neighbor object under the default 'single.nn'", {
  n.cells <- 60
  n.landmark <- 12
  k.nn <- 5
  obj <- make_reduction_fixture(n.cells = n.cells, n.landmark = n.landmark)

  res <- FindsmNN(obj, sketch.assay = "LANDMARK", reduction = "pls",
                  k.nn = k.nn, verbose = FALSE)

  expect_true("single.nn" %in% names(res@neighbors))
  nn <- res[["single.nn"]]

  # query is the full object: every cell must appear
  expect_equal(length(nn@cell.names), n.cells)
  expect_setequal(nn@cell.names, colnames(obj))

  # k.nn is honored as given (no clamp to 20)
  expect_equal(ncol(nn@nn.idx), k.nn)

  # indices reference the landmark set only
  expect_true(all(nn@nn.idx >= 1L & nn@nn.idx <= n.landmark))
})


test_that("FindsmNN honors an explicit nn.name", {
  obj <- make_reduction_fixture()

  res <- FindsmNN(obj, sketch.assay = "LANDMARK", reduction = "pls",
                  k.nn = 5, nn.name = "my.nn", verbose = FALSE)

  expect_true("my.nn" %in% names(res@neighbors))
  expect_false("single.nn" %in% names(res@neighbors))
})


test_that("FindsmNN validates its inputs", {
  obj <- make_reduction_fixture()

  expect_error(
    FindsmNN(obj, sketch.assay = "LANDMARK", reduction = "pls", k.nn = 2),
    "larger than 1"
  )
  expect_error(
    FindsmNN(obj, sketch.assay = "NOPE", reduction = "pls", k.nn = 5),
    "sketch.assay is not correctly defined"
  )
  expect_error(
    FindsmNN(obj, sketch.assay = "LANDMARK", reduction = "absent", k.nn = 5),
    "is not found in the object"
  )
})
