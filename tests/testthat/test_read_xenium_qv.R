set.seed(42)

# ReadXenium accepts and documents `mols.qv.threshold`, but the quality column
# was not among those read and no filtering was applied, so every threshold
# returned the same transcripts.
write_xenium <- function(dir, with.qv = TRUE, n = 200L, n.low = 60L) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  transcripts <- data.frame(
    x_location = runif(n) * 100,
    y_location = runif(n) * 100,
    feature_name = sample(paste0("gene", 1:5), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  if (with.qv) {
    transcripts$qv <- c(rep(5, n.low), rep(30, n - n.low))
  }
  con <- gzfile(file.path(dir, "transcripts.csv.gz"), "w")
  write.csv(transcripts, con, row.names = FALSE)
  close(con)
  cells <- data.frame(
    cell_id = paste0("cell", seq_len(20)),
    x_centroid = runif(20) * 100,
    y_centroid = runif(20) * 100,
    stringsAsFactors = FALSE
  )
  con <- gzfile(file.path(dir, "cells.csv.gz"), "w")
  write.csv(cells, con, row.names = FALSE)
  close(con)
  invisible(dir)
}

read_microns <- function(dir, ...) {
  suppressWarnings(ReadXenium(dir, outs = "microns", type = "centroids", ...))$microns
}

test_that("mols.qv.threshold removes low-quality transcripts", {
  dir <- write_xenium(tempfile("xenium"), n = 200L, n.low = 60L)
  expect_equal(nrow(read_microns(dir, mols.qv.threshold = 0)), 200L)
  expect_equal(nrow(read_microns(dir, mols.qv.threshold = 20)), 140L)
  expect_equal(nrow(read_microns(dir, mols.qv.threshold = 40)), 0L)
})

test_that("different thresholds give different results", {
  dir <- write_xenium(tempfile("xenium"))
  expect_false(identical(
    nrow(read_microns(dir, mols.qv.threshold = 0)),
    nrow(read_microns(dir, mols.qv.threshold = 20))
  ))
})

test_that("the returned columns are unchanged", {
  dir <- write_xenium(tempfile("xenium"))
  microns <- read_microns(dir, mols.qv.threshold = 20)
  expect_identical(colnames(microns), c("x", "y", "gene"))
  expect_false("qv" %in% colnames(microns))
})

test_that("outputs without a qv column still load, with a warning", {
  dir <- write_xenium(tempfile("xenium"), with.qv = FALSE)
  expect_warning(
    ReadXenium(dir, outs = "microns", type = "centroids", mols.qv.threshold = 20),
    "no 'qv' column|No 'qv' column"
  )
  expect_equal(nrow(read_microns(dir, mols.qv.threshold = 20)), 200L)
})

test_that("the default threshold filters", {
  dir <- write_xenium(tempfile("xenium"), n = 200L, n.low = 60L)
  # the documented default is 20
  expect_equal(nrow(read_microns(dir)), 140L)
})
