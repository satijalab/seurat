# TransferData picks the best identity per cell with which.max(). For a cell
# whose scores are all NaN that returns integer(0), so the rowwise apply gave a
# ragged list and indexing the identities failed with
#   invalid subscript type 'list'

test_that("rows with a clear maximum give its index", {
  scores <- rbind(c(0.2, 0.8), c(0.9, 0.1), c(0.4, 0.6))
  expect_identical(Seurat:::.MaxScoreIndex(scores), c(2L, 1L, 2L))
})

test_that("a row that is entirely NaN gives NA, not integer(0)", {
  scores <- rbind(c(0.2, 0.8), c(NaN, NaN), c(0.5, 0.5))
  index <- Seurat:::.MaxScoreIndex(scores)
  expect_type(index, "integer")
  expect_false(is.list(index))
  expect_length(index, 3L)
  expect_identical(index, c(2L, NA_integer_, 1L))
})

test_that("indexing identities with the result works", {
  ids <- c("A", "B")
  scores <- rbind(c(0.2, 0.8), c(NaN, NaN), c(0.5, 0.5))
  # the failure being fixed
  expect_error(ids[apply(scores, 1, which.max)], "invalid subscript type")
  expect_identical(ids[Seurat:::.MaxScoreIndex(scores)], c("B", NA, "A"))
})

test_that("all-NaN and single-column matrices are handled", {
  expect_identical(
    Seurat:::.MaxScoreIndex(rbind(c(NaN, NaN), c(NaN, NaN))),
    c(NA_integer_, NA_integer_)
  )
  expect_identical(Seurat:::.MaxScoreIndex(matrix(c(1, NaN), ncol = 1)), c(1L, NA_integer_))
})

test_that("NA scores behave like NaN", {
  expect_identical(
    Seurat:::.MaxScoreIndex(rbind(c(NA_real_, NA_real_), c(0.1, 0.9))),
    c(NA_integer_, 2L)
  )
})

test_that("TransferData is unchanged on ordinary data", {
  set.seed(1)
  mk <- function(tag, n = 120L, nf = 200L) {
    counts <- as.sparse(matrix(rpois(nf * n, lambda = 2), nrow = nf))
    dimnames(counts) <- list(paste0("gene", seq_len(nf)), paste0(tag, "_c", seq_len(n)))
    object <- suppressWarnings(CreateSeuratObject(counts = counts))
    object <- suppressWarnings(NormalizeData(object, verbose = FALSE))
    object <- suppressWarnings(FindVariableFeatures(object, nfeatures = 100, verbose = FALSE))
    object <- suppressWarnings(ScaleData(object, verbose = FALSE))
    suppressWarnings(RunPCA(object, npcs = 20, verbose = FALSE))
  }
  reference <- mk("r")
  query <- mk("q")
  reference$celltype <- rep(c("A", "B"), length.out = ncol(reference))
  anchors <- suppressWarnings(
    FindTransferAnchors(reference = reference, query = query, dims = 1:15, verbose = FALSE)
  )
  predictions <- suppressWarnings(
    TransferData(anchorset = anchors, refdata = reference$celltype, dims = 1:15, verbose = FALSE)
  )
  expect_equal(nrow(predictions), ncol(query))
  expect_false(anyNA(predictions$predicted.id))
  expect_setequal(unique(predictions$predicted.id), c("A", "B"))
})
