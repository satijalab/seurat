context("ReadMtx")

test_that("skip.cell and skip.feature work", {
  skip_on_cran()
  mtx <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126836/suppl/GSE126836_SN_MD5828_matrix.mtx.gz"
  features <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126836/suppl/GSE126836_SN_MD5828_genes.csv.gz"
  cells <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126836/suppl/GSE126836_SN_MD5828_barcodes.csv.gz"
  counts1 <- ReadMtx(mtx = mtx, cells = cells, features = features, feature.column = 1, skip.cell = 1, skip.feature = 1)
  expect_is(counts1, "dgCMatrix")
  expect_equal(ncol(counts1), 1436)
  expect_equal(nrow(counts1), 29445)
  expect_equal(colnames(counts1)[5], "MD5828a_GGGCATCCAATGAAAC-1")
  expect_equal(rownames(counts1)[2], "A1BG-AS1")
})


test_that("ReadMtx works", {
  skip_on_cran()
  mtx <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE127774&format=file&file=GSE127774%5FACC%5FB%5Fmatrix%2Emtx%2Egz"
  cells <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE127774&format=file&file=GSE127774%5FACC%5FB%5Fbarcodes%2Etsv%2Egz"
  features <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE127774&format=file&file=GSE127774%5FACC%5FB%5Fgenes%2Etsv%2Egz"
  counts2 <- ReadMtx(mtx = mtx, cells = cells, features = features, feature.column = 1)
  expect_is(counts2, "dgCMatrix")
  expect_equal(ncol(counts2), 22063)
  expect_equal(nrow(counts2), 22530)
  expect_equal(colnames(counts2)[1], "AAACCTGAGCAATCTC-1")
  expect_equal(rownames(counts2)[2], "ENSPPAG00000040697")
})


