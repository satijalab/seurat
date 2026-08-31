set.seed(42)

# QuPath exports name the centroid columns inconsistently. Matching one literal
# spelling left `mtx[[NA]]`, a zero-length vector, and the frame built from it
# failed with "arguments imply differing number of rows".
write_qupath <- function(xname, yname, n = 30) {
  df <- data.frame(id = seq_len(n), area = runif(n))
  if (!is.null(xname)) {
    df[[xname]] <- runif(n) * 100
    df[[yname]] <- runif(n) * 100
  }
  for (marker in c("CD3", "CD8", "DAPI")) {
    df[[paste0("Cell: Mean ", marker)]] <- runif(n)
  }
  path <- tempfile(fileext = ".csv")
  write.csv(df, path, row.names = FALSE)
  path
}

read_qupath <- function(...) {
  suppressWarnings(suppressMessages(ReadAkoya(write_qupath(...), type = "qupath")))
}

test_that("the documented column names are read", {
  out <- read_qupath("Centroid X", "Centroid Y")
  expect_equal(nrow(out$centroids), 30L)
  expect_setequal(colnames(out$centroids), c("x", "y", "cell"))
})

test_that("names carrying units are read", {
  expect_equal(nrow(read_qupath("Centroid X um", "Centroid Y um")$centroids), 30L)
})

test_that("names differing in case are read", {
  expect_equal(nrow(read_qupath("centroid x", "centroid y")$centroids), 30L)
})

test_that("names with a dot are read", {
  # what read.csv(check.names = TRUE) makes of "Centroid X"
  expect_equal(nrow(read_qupath("Centroid.X", "Centroid.Y")$centroids), 30L)
})

test_that("coordinates are the ones in the file", {
  path <- write_qupath("Centroid X", "Centroid Y")
  raw <- read.csv(path, check.names = FALSE)
  out <- suppressWarnings(suppressMessages(ReadAkoya(path, type = "qupath")))
  expect_equal(out$centroids$x, raw[["Centroid X"]])
  expect_equal(out$centroids$y, raw[["Centroid Y"]])
})

test_that("a missing centroid column is reported, not turned into a row mismatch", {
  path <- write_qupath(NULL, NULL)
  expect_error(
    suppressWarnings(suppressMessages(ReadAkoya(path, type = "qupath"))),
    "Could not find the centroid"
  )
  expect_error(
    suppressWarnings(suppressMessages(ReadAkoya(path, type = "qupath"))),
    "Columns present"
  )
})
