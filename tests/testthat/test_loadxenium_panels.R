# Read10X_h5() returns a bare matrix when a Xenium panel carries a single
# feature type, and a named list when control probes sit alongside the genes.
# LoadXenium() indexed the result as a list either way.
set.seed(42)

skip_unless_h5 <- function() {
  skip_if_not(
    requireNamespace("rhdf5", quietly = TRUE) ||
      requireNamespace("hdf5r", quietly = TRUE)
  )
  skip_if_not_installed("rhdf5")
}

write_xenium <- function(controls = TRUE, ngene = 30, ncontrol = 6, ncell = 25) {
  dir <- tempfile("xenium")
  dir.create(dir, recursive = TRUE)
  cells <- paste0("cell", seq_len(ncell))
  gene <- matrix(
    rpois(ngene * ncell, lambda = 3), nrow = ngene,
    dimnames = list(paste0("g", seq_len(ngene)), cells)
  )
  control <- matrix(
    rpois(ncontrol * ncell, lambda = 1), nrow = ncontrol,
    dimnames = list(paste0("NegControlProbe_", seq_len(ncontrol)), cells)
  )
  counts <- if (controls) as.sparse(rbind(gene, control)) else as.sparse(gene)
  types <- if (controls) {
    c(rep("Gene Expression", ngene), rep("Negative Control Probe", ncontrol))
  } else {
    rep("Gene Expression", ngene)
  }
  path <- file.path(dir, "cell_feature_matrix.h5")
  csc <- as(counts, "CsparseMatrix")
  rhdf5::h5createFile(path)
  rhdf5::h5createGroup(path, "matrix")
  rhdf5::h5createGroup(path, "matrix/features")
  rhdf5::h5write(as.numeric(csc@x), path, "matrix/data")
  rhdf5::h5write(as.integer(csc@i), path, "matrix/indices")
  rhdf5::h5write(as.integer(csc@p), path, "matrix/indptr")
  rhdf5::h5write(as.integer(dim(counts)), path, "matrix/shape")
  rhdf5::h5write(colnames(counts), path, "matrix/barcodes")
  rhdf5::h5write(rownames(counts), path, "matrix/features/name")
  rhdf5::h5write(rownames(counts), path, "matrix/features/id")
  rhdf5::h5write(types, path, "matrix/features/feature_type")
  rhdf5::h5closeAll()
  meta <- data.frame(
    cell_id = cells, x_centroid = runif(ncell) * 100, y_centroid = runif(ncell) * 100
  )
  con <- gzfile(file.path(dir, "cells.csv.gz"), "w")
  write.csv(meta, con, row.names = FALSE)
  close(con)
  transcripts <- data.frame(
    x_location = runif(60) * 100, y_location = runif(60) * 100,
    feature_name = sample(rownames(gene), 60, replace = TRUE), qv = 30
  )
  con <- gzfile(file.path(dir, "transcripts.csv.gz"), "w")
  write.csv(transcripts, con, row.names = FALSE)
  close(con)
  list(dir = dir, gene = gene, control = control)
}

load_it <- function(dir) {
  suppressWarnings(suppressMessages(LoadXenium(dir, fov = "fov")))
}

test_that("a panel with only gene expression loads", {
  skip_unless_h5()
  # the standalone gene panels produce a single feature type, and indexing the
  # bare matrix as a list raised "this S4 class is not subsettable"
  built <- write_xenium(controls = FALSE)
  object <- load_it(built$dir)
  expect_identical(Assays(object), "Xenium")
  expect_equal(nrow(object[["Xenium"]]), nrow(built$gene))
  expect_setequal(colnames(object), colnames(built$gene))
})

test_that("the counts of a single-type panel are the ones in the file", {
  skip_unless_h5()
  built <- write_xenium(controls = FALSE)
  object <- load_it(built$dir)
  got <- LayerData(object[["Xenium"]], layer = "counts")
  expect_equal(as.matrix(got), built$gene[rownames(got), colnames(got)])
})

test_that("a panel with control probes still splits into assays", {
  skip_unless_h5()
  built <- write_xenium(controls = TRUE)
  object <- load_it(built$dir)
  expect_setequal(Assays(object), c("Xenium", "ControlProbe"))
  expect_equal(nrow(object[["Xenium"]]), nrow(built$gene))
  expect_equal(nrow(object[["ControlProbe"]]), nrow(built$control))
})

test_that("gene counts are the same whether or not controls are present", {
  skip_unless_h5()
  with.controls <- write_xenium(controls = TRUE)
  without <- write_xenium(controls = FALSE)
  a <- LayerData(load_it(with.controls$dir)[["Xenium"]], layer = "counts")
  b <- LayerData(load_it(without$dir)[["Xenium"]], layer = "counts")
  expect_identical(rownames(a), rownames(b))
})
