# A Visium Gene + Protein panel yields one counts matrix per modality. Stacking
# them into a single assay makes nCount/nFeature sums across modalities, and
# antibody counts run orders of magnitude higher than gene counts, so the gene
# expression metrics become meaningless.
set.seed(42)

skip_unless_h5 <- function() {
  skip_if_not(
    requireNamespace("rhdf5", quietly = TRUE) ||
      requireNamespace("hdf5r", quietly = TRUE)
  )
}

build_visium <- function(ngene = 40, nab = 8, multimodal = TRUE) {
  image <- Read10X_Image(
    file.path("..", "testdata", "visium", "spatial"),
    image.name = "tissue_lowres_image.png"
  )
  barcodes <- Cells(image)
  gene <- matrix(
    rpois(ngene * length(barcodes), lambda = 3), nrow = ngene,
    dimnames = list(paste0("gene", seq_len(ngene)), barcodes)
  )
  antibody <- matrix(
    rpois(nab * length(barcodes), lambda = 5000), nrow = nab,
    dimnames = list(paste0("ab", seq_len(nab)), barcodes)
  )
  counts <- if (multimodal) as.sparse(rbind(gene, antibody)) else as.sparse(gene)
  types <- if (multimodal) {
    c(rep("Gene Expression", ngene), rep("Antibody Capture", nab))
  } else {
    rep("Gene Expression", ngene)
  }
  dir <- tempfile("visium")
  dir.create(dir)
  file.copy(file.path("..", "testdata", "visium", "spatial"), dir, recursive = TRUE)
  path <- file.path(dir, "filtered_feature_bc_matrix.h5")
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
  list(dir = dir, gene = gene, antibody = antibody)
}

load_it <- function(dir) {
  suppressWarnings(suppressMessages(Load10X_Spatial(dir, filter.matrix = FALSE)))
}

test_that("each modality gets its own assay", {
  skip_unless_h5()
  b <- build_visium()
  object <- load_it(b$dir)
  expect_setequal(Assays(object), c("Spatial", "Antibody.Capture"))
  expect_equal(nrow(object[["Spatial"]]), nrow(b$gene))
  expect_equal(nrow(object[["Antibody.Capture"]]), nrow(b$antibody))
})

test_that("nCount_Spatial counts gene expression only", {
  skip_unless_h5()
  b <- build_visium()
  object <- load_it(b$dir)
  # previously this was colSums(gene) + colSums(antibody)
  expect_equal(
    unname(object$nCount_Spatial),
    unname(colSums(b$gene)[colnames(object)])
  )
  expect_equal(
    unname(object$nFeature_Spatial),
    unname(colSums(b$gene[, colnames(object)] > 0))
  )
})

test_that("the antibody counts are kept, not discarded", {
  skip_unless_h5()
  b <- build_visium()
  object <- load_it(b$dir)
  got <- LayerData(object[["Antibody.Capture"]], layer = "counts")
  expect_equal(as.matrix(got), b$antibody[, colnames(object)])
})

test_that("a single-modality capture area is unchanged", {
  skip_unless_h5()
  b <- build_visium(multimodal = FALSE)
  object <- load_it(b$dir)
  expect_identical(Assays(object), "Spatial")
  expect_equal(
    unname(object$nCount_Spatial),
    unname(colSums(b$gene)[colnames(object)])
  )
})
