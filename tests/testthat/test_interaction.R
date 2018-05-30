# Tests for interaction  functions
set.seed(42)

# Tests for WhichCells
# --------------------------------------------------------------------------------
context("WhichCells")

test_that("ident subsetting works as expected", {
  x <- WhichCells(pbmc_small, ident = 0)
  y <- WhichCells(pbmc_small, ident = c(0,1))
  expect_equal(length(x), 29)
  expect_equal(length(y), 58)
  expect_equal(x[1], "ATGCCAGAACGACT")
  expect_equal(y[58], "GCGTAAACACGGTT")
  expect_error(WhichCells(pbmc_small, ident = "E"))
  expect_true(all(pbmc_small@ident[x] == 0))
})

test_that("ident removal works as expected", {
  x <- WhichCells(pbmc_small, ident.remove = 0)
  y <- WhichCells(pbmc_small, ident.remove = c(0, 1))
  expect_equal(length(x), 51)
  expect_equal(length(y), 22)
  expect_equal(x[1], "TACGCCACTCCGAA")
  expect_equal(y[22], "CTTGATTGATCTTC")
  expect_error(WhichCells(pbmc_small, ident.remove = "E"))
  expect_true(all(pbmc_small@ident[x] != 0))
})

test_that("providing cell.use works as expected", {
  cells.use <- pbmc_small@cell.names[1:20]
  x <- WhichCells(pbmc_small, cells.use = cells.use)
  expect_equal(length(intersect(x, cells.use)), 20)
  x <- WhichCells(pbmc_small, cells.use = cells.use, ident = 0)
  expect_equal(length(x), 17)
  expect_true(all(pbmc_small@ident[x] == 0))
  y <- WhichCells(pbmc_small, cells.use = cells.use, ident = c(0, 1))
  expect_equal(length(intersect(y, cells.use)), 20)
  expect_error(WhichCells(pbmc_small, cells.use = cells.use, ident = "E"))
  expect_error(WhichCells(pbmc_small, cells.use = cells.use, ident = c("E", 1)))
  expect_error(WhichCells(pbmc_small, cells.use = cells.use, ident.remove = c("E", 1)))
})

test_that("subset.name works as expected", {
  # test with accept.value
  x <- WhichCells(pbmc_small, subset.name = "res.1", accept.value = 0)
  y <- WhichCells(pbmc_small, subset.name = "res.1", accept.value =  c(0,1))
  expect_equal(length(x), 29)
  expect_equal(length(y), 58)
  expect_equal(x[1], "ATGCCAGAACGACT")
  expect_equal(y[58], "GCGTAAACACGGTT")
  expect_error(WhichCells(pbmc_small,  subset.name = "res.1", accept.value = "E"))
  expect_true(all(pbmc_small@ident[x] == 0))
  # test with accept.low
  x <- WhichCells(pbmc_small, subset.name = "nGene", accept.low = 50)
  expect_equal(length(x), 41)
  expect_equal(x[1], "CATGGCCTGTGCAT")
  expect_error(WhichCells(pbmc_small, subset.name = "nGene", accept.low = c(50, 100)))
  # test with accept.high
  x <- WhichCells(pbmc_small, subset.name = "nGene", accept.high = 50)
  expect_equal(length(x), 36)
  expect_equal(x[1], "ATGCCAGAACGACT")
  expect_error(WhichCells(pbmc_small, subset.name = "nGene", accept.high = c(50, 100)))
  expect_error(WhichCells(pbmc_small, subset.name = "nGene", accept.low = 100, accept.high = 50))
})

test_that("max.cells.per.ident downsamples properly", {
  x <- WhichCells(pbmc_small, max.cells.per.ident = 5)
  expect_equal(length(x), 20)
  expect_equal(x[1], "GCAGCTCTGTTTCT")
  x <- WhichCells(pbmc_small, ident = 0, max.cells.per.ident = 5)
  expect_equal(length(x), 5)
})

# Tests for SubsetData
# --------------------------------------------------------------------------------
context("SubsetData")

test_that("cells.use works", {
  x <- SubsetData(pbmc_small, cells.use = pbmc_small@cell.names[1:20])
  expect_true(all(x@cell.names %in% pbmc_small@cell.names[1:20]))
  expect_equal(length(x@cell.names), 20)
  expect_equal(dim(x@raw.data), c(230, 80))
  expect_equal(dim(x@data), c(230, 20))
  expect_equal(dim(x@scale.data), c(230, 20))
  expect_equal(length(x@ident), 20)
  expect_equal(nrow(x@meta.data), 20)
  expect_equal(nrow(x@dr$pca@cell.embeddings), 20)
})

test_that("do.clean and do.raw work", {
  x <- SubsetData(pbmc_small, cells.use = pbmc_small@cell.names[1:20], do.clean = TRUE)
  expect_equal(dim(x@raw.data), c(230, 20))
  expect_equal(dim(x@data), c(230, 20))
  expect_equal(dim(x@scale.data), c(1, 1))
  expect_equal(length(x@dr), 0)
  x <- SubsetData(pbmc_small, cells.use = pbmc_small@cell.names[1:20], subset.raw = TRUE)
  expect_equal(dim(x@raw.data), c(230, 20))
  expect_equal(dim(x@data), c(230, 20))
  expect_equal(dim(x@scale.data), c(230, 20))
  expect_equal(length(x@dr), 2)
})

# Tests for RenameCells
# --------------------------------------------------------------------------------
context("RenameCells")

test_that("argument checks work", {
  expect_error(RenameCells(pbmc_small),
               "One of 'add.cell.id' and 'new.names' must be set")
  expect_error(RenameCells(pbmc_small, add.cell.id = "A", new.names = "A"),
               "Only one of 'add.cell.id' and 'new.names' must be set")
})

test_that("add.cell.id works", {
    test.names <- paste("Test", pbmc_small@cell.names, sep = "_")
    x <- RenameCells(pbmc_small, add.cell.id = "Test")
    expect_equal(x@cell.names, test.names)
    expect_equal(colnames(x@raw.data), test.names)
    expect_equal(colnames(x@data), test.names)
    expect_equal(colnames(x@scale.data), test.names)
    expect_equal(rownames(x@meta.data), test.names)
    expect_equal(names(x@ident), test.names)
    expect_equal(rownames(x@dr$pca@cell.embeddings), test.names)
    expect_equal(rownames(x@dr$tsne@cell.embeddings), test.names)
    expect_equal(rownames(x@spatial@mix.probs), test.names)
})

test_that("new.names works", {
  test.names <- paste0("A", 1:80)
  x <- RenameCells(pbmc_small, new.names = paste0("A", 1:80))
  expect_equal(x@cell.names, test.names)
  expect_equal(colnames(x@raw.data), test.names)
  expect_equal(colnames(x@data), test.names)
  expect_equal(colnames(x@scale.data), test.names)
  expect_equal(rownames(x@meta.data), test.names)
  expect_equal(names(x@ident), test.names)
  expect_equal(rownames(x@dr$pca@cell.embeddings), test.names)
  expect_equal(rownames(x@dr$tsne@cell.embeddings), test.names)
  expect_equal(rownames(x@spatial@mix.probs), test.names)
})
