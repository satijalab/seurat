# hdf5r builds against the system HDF5 and does not compile against 1.12 or
# newer, so on a current system it can be uninstallable and 10x .h5 files
# unreadable. These cover the rhdf5 path that reads them instead.
set.seed(42)

skip_unless_rhdf5 <- function() {
  skip_if_not_installed("rhdf5")
}

write_10x_h5 <- function(mat, types = NULL, scaling = NULL) {
  path <- tempfile(fileext = ".h5")
  rhdf5::h5createFile(path)
  rhdf5::h5createGroup(path, "matrix")
  rhdf5::h5createGroup(path, "matrix/features")
  csc <- as(mat, "CsparseMatrix")
  rhdf5::h5write(as.numeric(csc@x), path, "matrix/data")
  rhdf5::h5write(as.integer(csc@i), path, "matrix/indices")
  rhdf5::h5write(as.integer(csc@p), path, "matrix/indptr")
  rhdf5::h5write(as.integer(dim(mat)), path, "matrix/shape")
  rhdf5::h5write(colnames(mat), path, "matrix/barcodes")
  rhdf5::h5write(rownames(mat), path, "matrix/features/name")
  rhdf5::h5write(rownames(mat), path, "matrix/features/id")
  rhdf5::h5write(
    types %||% rep("Gene Expression", nrow(mat)),
    path, "matrix/features/feature_type"
  )
  if (!is.null(scaling)) {
    fid <- rhdf5::H5Fopen(path)
    rhdf5::h5writeAttribute(scaling, fid, "protein_scaling_factor")
    rhdf5::H5Fclose(fid)
  }
  rhdf5::h5closeAll()
  path
}

make_counts <- function(nfeat = 30, ncell = 12) {
  m <- matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  as.sparse(m)
}

test_that("a single-modality matrix round-trips", {
  skip_unless_rhdf5()
  counts <- make_counts()
  got <- Seurat:::.Read10Xh5Rhdf5(write_10x_h5(counts))
  expect_s4_class(got, "dgCMatrix")
  expect_equal(as.matrix(got), as.matrix(counts))
  expect_identical(dimnames(got), dimnames(counts))
})

test_that("feature ids are used when use.names is FALSE", {
  skip_unless_rhdf5()
  counts <- make_counts()
  got <- Seurat:::.Read10Xh5Rhdf5(write_10x_h5(counts), use.names = FALSE)
  expect_identical(rownames(got), rownames(counts))
})

test_that("duplicate feature names are made unique on request", {
  skip_unless_rhdf5()
  counts <- make_counts(nfeat = 4)
  rownames(counts) <- c("a", "a", "b", "b")
  unique.names <- Seurat:::.Read10Xh5Rhdf5(write_10x_h5(counts))
  expect_identical(rownames(unique.names), c("a", "a.1", "b", "b.1"))
  kept <- Seurat:::.Read10Xh5Rhdf5(write_10x_h5(counts), unique.features = FALSE)
  expect_identical(rownames(kept), c("a", "a", "b", "b"))
})

test_that("multiple modalities come back as a named list", {
  skip_unless_rhdf5()
  counts <- make_counts()
  types <- c(rep("Gene Expression", 20), rep("Antibody Capture", 10))
  got <- suppressMessages(Seurat:::.Read10Xh5Rhdf5(write_10x_h5(counts, types = types)))
  expect_true(is.list(got))
  expect_setequal(names(got), c("Gene Expression", "Antibody Capture"))
  expect_equal(as.matrix(got[["Gene Expression"]]), as.matrix(counts[1:20, ]))
  expect_equal(as.matrix(got[["Antibody Capture"]]), as.matrix(counts[21:30, ]))
})

test_that("the protein scaling factor is applied", {
  skip_unless_rhdf5()
  counts <- make_counts()
  types <- c(rep("Gene Expression", 20), rep("Protein Expression", 10))
  got <- suppressMessages(Seurat:::.Read10Xh5Rhdf5(
    write_10x_h5(counts, types = types, scaling = 2)
  ))
  expect_equal(
    as.matrix(got[["Protein Expression"]]),
    as.matrix(counts[21:30, ]) / 2
  )
  # gene expression is left alone
  expect_equal(as.matrix(got[["Gene Expression"]]), as.matrix(counts[1:20, ]))
})

test_that("Read10X_h5 reports both options when neither reader is installed", {
  skip_unless_rhdf5()
  expect_error(Read10X_h5(tempfile(fileext = ".h5")), "File not found")
})
