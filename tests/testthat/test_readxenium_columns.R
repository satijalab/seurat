set.seed(42)

# The readers select the columns they need by name, so a cells or transcripts
# file that does not have them failed with "undefined columns selected", which
# names neither the file nor the column
write_xenium <- function(cells.columns, transcripts.columns) {
  dir <- tempfile("xenium")
  dir.create(path = dir, recursive = TRUE)
  cells <- data.frame(
    cell_id = paste0("cell", seq_len(10)),
    x_centroid = runif(10) * 100,
    y_centroid = runif(10) * 100
  )
  connection <- gzfile(file.path(dir, "cells.csv.gz"), "w")
  write.csv(cells[, cells.columns, drop = FALSE], connection, row.names = FALSE)
  close(connection)

  transcripts <- data.frame(
    x_location = runif(20) * 100,
    y_location = runif(20) * 100,
    feature_name = paste0("g", sample(5, 20, replace = TRUE)),
    qv = 30
  )
  connection <- gzfile(file.path(dir, "transcripts.csv.gz"), "w")
  write.csv(transcripts[, transcripts.columns, drop = FALSE], connection, row.names = FALSE)
  close(connection)
  dir
}

all.cells <- c("cell_id", "x_centroid", "y_centroid")
all.transcripts <- c("x_location", "y_location", "feature_name", "qv")

test_that("a cells file with no centroids names the columns it lacks", {
  dir <- write_xenium("cell_id", all.transcripts)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  expect_error(
    suppressWarnings(ReadXenium(dir, type = "centroids", outs = "microns")),
    "cells file has no .x_centroid., .y_centroid. columns"
  )
})

test_that("a transcripts file with no feature names says so", {
  dir <- write_xenium(all.cells, c("x_location", "y_location"))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  expect_error(
    suppressWarnings(ReadXenium(dir, type = "centroids", outs = "microns")),
    "transcripts file has no .feature_name. column"
  )
})

test_that("a complete pair of files still reads", {
  dir <- write_xenium(all.cells, all.transcripts)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  data <- suppressWarnings(ReadXenium(dir, type = "centroids", outs = "microns"))
  expect_identical(colnames(data$centroids), c("x", "y", "cell"))
  expect_identical(nrow(data$centroids), 10L)
  expect_identical(colnames(data$microns), c("x", "y", "gene"))
})
