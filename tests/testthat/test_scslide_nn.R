# Tests for the cell-to-landmark neighbor search behind scSLIDE

# A Seurat object holding two reductions and, optionally, a landmark assay
make_nn_object <- function(embedding.a, embedding.b, landmarks = NULL) {
  cells <- rownames(x = embedding.a)
  counts <- Matrix::sparseMatrix(
    i = rep(x = 1:5, times = length(x = cells)),
    j = rep(x = seq_along(along.with = cells), each = 5),
    x = 1,
    dims = c(5, length(x = cells)),
    dimnames = list(paste0("gene", 1:5), cells)
  )
  object <- suppressWarnings(CreateSeuratObject(counts = counts))
  colnames(x = embedding.a) <- paste0("A_", seq_len(length.out = ncol(x = embedding.a)))
  colnames(x = embedding.b) <- paste0("B_", seq_len(length.out = ncol(x = embedding.b)))
  object[["a"]] <- suppressWarnings(
    CreateDimReducObject(embeddings = embedding.a, key = "A_", assay = "RNA")
  )
  object[["b"]] <- suppressWarnings(
    CreateDimReducObject(embeddings = embedding.b, key = "B_", assay = "RNA")
  )
  if (!is.null(x = landmarks)) {
    suppressWarnings(object[["LANDMARK"]] <- subset(x = object[["RNA"]], cells = landmarks))
  }
  return(object)
}

random_embedding <- function(n.cells, n.dims) {
  return(matrix(
    data = rnorm(n = n.cells * n.dims),
    nrow = n.cells,
    dimnames = list(paste0("cell", seq_len(length.out = n.cells)), NULL)
  ))
}

run_mmnn <- function(object, dims = 10, k.nn = 20, knn.range = 50, ...) {
  return(suppressWarnings(suppressMessages(FindmmNN(
    object = object,
    sketch.assay = "LANDMARK",
    reduction.list = list("a", "b"),
    dims.list = list(1:dims, 1:dims),
    k.nn = k.nn,
    knn.range = knn.range,
    fix.wnn.weights = c(0.5, 0.5),
    verbose = FALSE,
    ...
  ))))
}

test_that("FindmmNN keeps the nearest landmark of a non-landmark cell", {
  # Landmarks 1-30 sit on a line, and every other cell sits right next to landmark 10
  n.cells <- 240
  embedding <- cbind(c(seq_len(length.out = 30), rep(x = 10.01, times = n.cells - 30)), 1)
  rownames(x = embedding) <- paste0("cell", seq_len(length.out = n.cells))
  object <- make_nn_object(
    embedding.a = embedding,
    embedding.b = embedding,
    landmarks = rownames(x = embedding)[1:30]
  )
  object <- run_mmnn(object = object, dims = 2, k.nn = 5, knn.range = 25, l2.norm = FALSE)
  expect_equal(object[["weighted.nn"]]@nn.idx[31:n.cells, 1], rep(x = 10, times = n.cells - 30))
})

test_that("FindmmNN lists a landmark as its own first neighbor, as FindsmNN does", {
  set.seed(1)
  embedding <- random_embedding(n.cells = 300, n.dims = 10)
  landmarks <- rownames(x = embedding)[1:60]
  object <- make_nn_object(embedding.a = embedding, embedding.b = embedding, landmarks = landmarks)
  object <- run_mmnn(object = object)
  object <- suppressMessages(FindsmNN(
    object = object,
    sketch.assay = "LANDMARK",
    reduction = "a",
    dims = 1:10,
    verbose = FALSE
  ))
  expect_equal(object[["weighted.nn"]]@nn.idx[1:60, 1], 1:60)
  expect_equal(object[["single.nn"]]@nn.idx[1:60, 1], 1:60)
})

test_that("FindmmNN returns complete neighbors when knn.range equals k.nn", {
  # Identical embeddings give the smallest possible candidate pool
  set.seed(2)
  embedding <- random_embedding(n.cells = 300, n.dims = 10)
  object <- make_nn_object(
    embedding.a = embedding,
    embedding.b = embedding,
    landmarks = rownames(x = embedding)[1:60]
  )
  object <- run_mmnn(object = object, k.nn = 20, knn.range = 20)
  expect_false(anyNA(x = object[["weighted.nn"]]@nn.idx))
})

test_that("GenerateSampleObject only counts non-landmark cells", {
  set.seed(3)
  n.cells <- 300
  k.nn <- 5
  embedding.a <- random_embedding(n.cells = n.cells, n.dims = 10)
  embedding.b <- random_embedding(n.cells = n.cells, n.dims = 10)
  # donor1 holds every landmark, so its count is the one self-neighbors would inflate
  landmarks <- rownames(x = embedding.a)[1:60]
  object <- make_nn_object(embedding.a = embedding.a, embedding.b = embedding.b, landmarks = landmarks)
  object$Donor <- rep(x = c("donor1", "donor2"), each = n.cells / 2)
  object <- run_mmnn(object = object)
  samples <- suppressWarnings(suppressMessages(GenerateSampleObject(
    object = object,
    nn.name = "weighted.nn",
    sketch.assay = "LANDMARK",
    group.by = "Donor",
    k.nn = k.nn,
    add.meta.data = FALSE,
    verbose = FALSE
  )))
  sample.counts <- LayerData(object = samples, layer = "counts")
  n.counted <- c(donor1 = n.cells / 2 - 60, donor2 = n.cells / 2)
  expect_equal(Matrix::colSums(x = sample.counts)[names(x = n.counted)], k.nn * n.counted)
})

test_that("FindMultiModalNeighbors still excludes each cell from its own neighbors", {
  set.seed(4)
  object <- make_nn_object(
    embedding.a = random_embedding(n.cells = 300, n.dims = 10),
    embedding.b = random_embedding(n.cells = 300, n.dims = 10)
  )
  object <- suppressMessages(FindMultiModalNeighbors(
    object = object,
    reduction.list = list("a", "b"),
    dims.list = list(1:10, 1:10),
    k.nn = 20,
    knn.range = 100,
    verbose = FALSE
  ))
  nn.idx <- object[["weighted.nn"]]@nn.idx
  expect_equal(sum(nn.idx == seq_len(length.out = nrow(x = nn.idx))), 0)
})

test_that("GenerateSampleObject honors verbose during chi-squared normalization", {
  set.seed(5)
  embedding.a <- random_embedding(n.cells = 100, n.dims = 10)
  embedding.b <- random_embedding(n.cells = 100, n.dims = 10)
  object <- make_nn_object(
    embedding.a = embedding.a,
    embedding.b = embedding.b,
    landmarks = rownames(x = embedding.a)[1:30]
  )
  object$Donor <- rep(x = c("donor1", "donor2"), each = 50)
  object <- run_mmnn(object = object, knn.range = 25)
  generate <- function(verbose) {
    return(capture.output(invisible(x = suppressWarnings(suppressMessages(GenerateSampleObject(
      object = object,
      nn.name = "weighted.nn",
      sketch.assay = "LANDMARK",
      group.by = "Donor",
      k.nn = 5,
      normalization.method = "ChiSquared",
      add.meta.data = FALSE,
      verbose = verbose
    ))))))
  }
  expect_length(generate(verbose = FALSE), 0)
  expect_match(paste(generate(verbose = TRUE), collapse = "\n"), "Chi-squared")
})
