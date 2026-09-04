set.seed(7)

# Writes a small CosMx export: an expression matrix, per cell metadata,
# transcript coordinates and cell boundaries, named the way ReadNanostring
# looks for them
write_cosmx <- function(dir, ncell = 6L, ngene = 12L, fovs = 1:2) {
  dir.create(path = dir, showWarnings = FALSE, recursive = TRUE)
  genes <- paste0("Gene", seq_len(ngene))
  grid <- expand.grid(cell_ID = seq_len(ncell), fov = fovs)
  counts <- function(n) {
    as.data.frame(matrix(
      rpois(n * ngene, lambda = 3),
      ncol = ngene,
      dimnames = list(NULL, genes)
    ))
  }
  expression <- cbind(grid[, c("fov", "cell_ID")], counts(nrow(grid)))
  # molecules not assigned to a cell are recorded against cell_ID 0
  unassigned <- cbind(data.frame(fov = fovs, cell_ID = 0), counts(length(fovs)))
  write.csv(rbind(expression, unassigned), file.path(dir, "run_exprMat_file.csv"), row.names = FALSE)

  metadata <- data.frame(
    fov = grid$fov,
    cell_ID = grid$cell_ID,
    Area = runif(nrow(grid), 100, 200),
    CenterX_global_px = runif(nrow(grid), 0, 1000),
    CenterY_global_px = runif(nrow(grid), 0, 1000)
  )
  write.csv(metadata, file.path(dir, "run_metadata_file.csv"), row.names = FALSE)

  write.csv(
    data.frame(
      fov = rep(grid$fov, each = 5),
      cell_ID = rep(grid$cell_ID, each = 5),
      x_global_px = runif(nrow(grid) * 5, 0, 1000),
      y_global_px = runif(nrow(grid) * 5, 0, 1000),
      target = sample(genes, nrow(grid) * 5, replace = TRUE)
    ),
    file.path(dir, "run_tx_file.csv"),
    row.names = FALSE
  )

  boundaries <- do.call(what = rbind, args = lapply(
    X = seq_len(nrow(grid)),
    FUN = function(i) {
      angle <- seq(0, 2 * pi, length.out = 9)[1:8]
      data.frame(
        fov = grid$fov[i],
        cellID = grid$cell_ID[i],
        x_global_px = metadata$CenterX_global_px[i] + 10 * cos(angle),
        y_global_px = metadata$CenterY_global_px[i] + 10 * sin(angle)
      )
    }
  ))
  write.csv(boundaries, file.path(dir, "run-polygons.csv"), row.names = FALSE)
  metadata
}

test_that("ReadNanostring reads a CosMx export", {
  skip_if_not_installed("data.table")
  dir <- file.path(tempdir(), "cosmx")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  write_cosmx(dir)

  data <- ReadNanostring(data.dir = dir, type = c("centroids", "segmentations"))
  expect_equal(dim(data$matrix), c(12L, 12L))
  expect_identical(rownames(data$matrix), paste0("Gene", 1:12))
  # cells are named cell_ID_fov, and the unassigned rows are dropped
  expect_identical(colnames(data$matrix), paste0(rep(1:6, 2), "_", rep(1:2, each = 6)))
  expect_setequal(data$centroids$cell, colnames(data$matrix))
})

# md[, metadata] drops to a vector when one column is asked for, and the cell
# names assigned to it afterwards went nowhere
test_that("ReadNanostring returns a single metadata column with its cells", {
  skip_if_not_installed("data.table")
  dir <- file.path(tempdir(), "cosmx-metadata")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  written <- write_cosmx(dir)

  one <- ReadNanostring(data.dir = dir, metadata = "Area")$metadata
  expect_identical(colnames(one), c("Area", "cell"))
  expect_setequal(one$cell, paste0(written$cell_ID, "_", written$fov))

  two <- ReadNanostring(data.dir = dir, metadata = c("Area", "fov"))$metadata
  expect_identical(colnames(two), c("Area", "fov", "cell"))
  expect_equal(one$Area, two$Area)
})

# the denominator was length(), which for a data frame is its column count, so
# the zero fraction was compared against the wrong total
test_that("dense counts are left dense", {
  skip_if_not_installed("data.table")
  dir <- file.path(tempdir(), "cosmx-dense")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  # lambda 3 gives about 5% zeros, well under the 0.4 ratio
  write_cosmx(dir, ncell = 20L, ngene = 30L)

  data <- ReadNanostring(data.dir = dir)
  expect_lt(sum(data$matrix == 0) / prod(dim(data$matrix)), 0.4)
  expect_false(inherits(data$matrix, "sparseMatrix"))
})

test_that("LoadNanostring builds an object with a matching field of view", {
  skip_if_not_installed("data.table")
  dir <- file.path(tempdir(), "cosmx-load")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  write_cosmx(dir)

  object <- suppressWarnings(LoadNanostring(data.dir = dir, fov = "test"))
  expect_identical(Images(object), "test")
  expect_setequal(Cells(object[["test"]]), colnames(object))
  expect_identical(rownames(object), paste0("Gene", 1:12))
})
