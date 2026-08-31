set.seed(42)

# With one cell every variance is NA, and with constant data every variance is
# zero, so the loess fit had nothing to fit and failed with "invalid 'x'"
make_object <- function(counts) {
  object <- suppressWarnings(CreateSeuratObject(counts = as.sparse(counts)))
  suppressWarnings(NormalizeData(object, verbose = FALSE))
}

test_that("a single cell says why no variable features can be found", {
  counts <- matrix(
    sample(1:5, 50, replace = TRUE),
    ncol = 2,
    dimnames = list(paste0("gene", 1:25), c("Cell1", "Cell2"))
  )
  object <- make_object(counts)
  one <- subset(object, cells = "Cell1")
  expect_error(
    suppressWarnings(FindVariableFeatures(one, verbose = FALSE)),
    "no feature has non-zero variance, as the data has 1 cell"
  )
})

test_that("constant data says why no variable features can be found", {
  counts <- matrix(3, nrow = 25, ncol = 4,
                   dimnames = list(paste0("gene", 1:25), paste0("c", 1:4)))
  expect_error(
    suppressWarnings(FindVariableFeatures(make_object(counts), verbose = FALSE)),
    "no feature has non-zero variance"
  )
})

test_that("features with no variance do not stop the rest being ranked", {
  counts <- matrix(
    rpois(25 * 20, lambda = 3),
    nrow = 25,
    dimnames = list(paste0("gene", 1:25), paste0("c", 1:20))
  )
  # a handful of features held constant across every cell
  counts[1:5, ] <- 2
  object <- make_object(counts)
  object <- suppressWarnings(FindVariableFeatures(object, nfeatures = 10, verbose = FALSE))
  expect_length(VariableFeatures(object), 10L)
  expect_false(any(paste0("gene", 1:5) %in% VariableFeatures(object)))
})
