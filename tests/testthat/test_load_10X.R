context("Read10X")
# These tests were added to ensure Seurat was forwards and backwards compatible for 3.0 data

dname = "../testdata/cr3.0"
test.data <- Read10X(dname)
test.data2 <- Read10X(c(dname, dname))

test_that("Cell Ranger 3.0 Data Parsing", {
  expect_is(test.data, "list")
  expect_equal(ncol(test.data$`Gene Expression`), .5 * ncol(test.data2$`Gene Expression`))
  expect_equal(ncol(test.data$`Antibody Capture`), .5 * ncol(test.data2$`Antibody Capture`))
  expect_equal(colnames(test.data2[[1]])[6], "2_AAAGTAGCACAGTCGC-1")
  expect_equal(test.data$`Gene Expression`[2,2], 1000)
})

# Tests of Pre-3.0 Data
test.data3 <- Read10X("../testdata/")
test_that("Read10X creates sparse matrix", {
  expect_is(test.data3, "dgCMatrix")
  expect_equal(colnames(test.data3)[1], "ATGCCAGAACGACT-1")
  expect_equal(rownames(test.data3)[1], "MS4A1")
})

test_that("Read10X handles missing files properly", {
  expect_error(Read10X("."))
  expect_error(Read10X("./notadir/"))
  expect_error(Read10X(dname, gene.column = 10))
})

# Tests for reading in spatial 10x data
if (requireNamespace("hdf5r", quietly = TRUE)) {
  context("Load10X_Spatial")
  dname <- "../testdata/visium"
  txsp <- Load10X_Spatial(data.dir = '../testdata/visium')
  test_that("10x Spatial Data Parsing", {
    expect_is(txsp, "Seurat")
    expect_equal(ncol(x = txsp), 2695)
    expect_equal(nrow(x = txsp), 100)
    expect_equal(Cells(x = txsp)[1], "AAACAAGTATCTCCCA-1")
    expect_equal(Assays(object = txsp), "Spatial")
    expect_equal(GetAssayData(object = txsp[["Spatial"]], slot = "counts")[5, 9], 1)
  })
  test_that("Read10X_Spatial handles missing files properly", {
    expect_error(Load10X_Spatial(data.dir = "."))
    expect_error(Load10X_Spatial(data.dir = "./notadir/"))
  })
}
