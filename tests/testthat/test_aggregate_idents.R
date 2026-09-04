# `CreateCategoryMatrix` recovers each grouping variable's level from the column
# names `sparse.model.matrix` produces. Those names join variables with ':', so
# a level that contains a colon must not be mistaken for a variable boundary.
set.seed(42)

obj <- suppressWarnings(pbmc_small)
obj$plain <- as.character(obj$letter.idents)
obj$colons <- paste0("anything:", obj$letter.idents)
obj$deep <- paste0("a:b:c:", obj$letter.idents)
obj$second <- rep(c("x", "y"), length.out = ncol(obj))
obj$underscored <- paste0("has_underscore_", obj$letter.idents)

agg <- function(group.by) {
  colnames(suppressWarnings(suppressMessages(
    AggregateExpression(obj, assays = "RNA", group.by = group.by, verbose = FALSE)$RNA
  )))
}

test_that("identity names containing colons are preserved", {
  expect_setequal(agg("colons"), c("anything:A", "anything:B"))
})

test_that("several colons in an identity name are preserved", {
  expect_setequal(agg("deep"), c("a:b:c:A", "a:b:c:B"))
})

test_that("colons survive grouping by more than one variable", {
  expect_setequal(
    agg(c("colons", "second")),
    c("anything:A_x", "anything:A_y", "anything:B_x", "anything:B_y")
  )
})

test_that("names without colons are unchanged", {
  expect_setequal(agg("plain"), c("A", "B"))
  expect_setequal(agg(c("plain", "second")), c("A_x", "A_y", "B_x", "B_y"))
})

test_that("underscores in identity names are still replaced with dashes", {
  expect_setequal(agg("underscored"), c("has-underscore-A", "has-underscore-B"))
})

test_that("AverageExpression preserves colons too", {
  cols <- colnames(suppressWarnings(suppressMessages(
    AverageExpression(obj, assays = "RNA", group.by = "colons", verbose = FALSE)$RNA
  )))
  expect_setequal(cols, c("anything:A", "anything:B"))
})

test_that("aggregated values are assigned to the right identity", {
  counts <- LayerData(obj, layer = "counts")
  result <- suppressWarnings(suppressMessages(
    AggregateExpression(obj, assays = "RNA", group.by = "colons", verbose = FALSE)$RNA
  ))
  for (ident in c("anything:A", "anything:B")) {
    cells <- colnames(obj)[obj$colons == ident]
    expect_equal(
      as.numeric(result[, ident]),
      as.numeric(Matrix::rowSums(counts[rownames(result), cells, drop = FALSE])),
      info = ident
    )
  }
})
