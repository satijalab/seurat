context("ReadMtx")

mtx <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126836/suppl/GSE126836_SN_MD5828_matrix.mtx.gz"
features <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126836/suppl/GSE126836_SN_MD5828_genes.csv.gz"
cells <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126836/suppl/GSE126836_SN_MD5828_barcodes.csv.gz"
counts1 <- ReadMtx(mtx = mtx, cells = cells, features = features, feature.column = 1, skip.cell = 1, skip.feature = 1)


mtx <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132044&format=file&file=GSE132044%5Fmixture%5Fhg19%5Fmm10%5Fcount%5Fmatrix%2Emtx%2Egz"
cells <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132044&format=file&file=GSE132044%5Fmixture%5Fhg19%5Fmm10%5Fcell%2Etsv%2Egz"
features <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132044&format=file&file=GSE132044%5Fmixture%5Fhg19%5Fmm10%5Fgene%2Etsv%2Egz"
counts2 <- ReadMtx(mtx = mtx, cells = cells, features = features, feature.column = 1)


test_that("skip.cell and skip.feature work", {
  expect_is(counts1, "dgCMatrix")
  expect_equal(ncol(counts1), 1436)
  expect_equal(nrow(counts1), 29445)
  expect_equal(colnames(counts1)[5], "MD5828a_GGGCATCCAATGAAAC-1")
  expect_equal(rownames(counts1)[2], "A1BG-AS1")
})


test_that("ReadMtx works", {
  expect_is(counts2, "dgCMatrix")
  expect_equal(ncol(counts2), 27714)
  expect_equal(nrow(counts2), 62046)
  expect_equal(colnames(counts2)[1], "Mixture1.Smart-seq2.p2_A4")
  expect_equal(rownames(counts2)[2], "hg19_ENSG00000000003_hg19_TSPAN6")
})


