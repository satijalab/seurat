# Tests for SubsetByBarcodeDistribution
# --------------------------------------------------------------------------------
context("SubsetByBarcodeDistribution")

## Testing by class
default <- SubsetByBarcodeDistribution(pbmc_small)
cells <- SubsetByBarcodeDistribution(pbmc_small,
                                    return.cells = TRUE)
g <- SubsetByBarcodeDistribution(pbmc_small,
                                return.plot = TRUE)
inflections <- SubsetByBarcodeDistribution(pbmc_small,
                                return.inflections = TRUE)


## Various thresholds should cut off at different inflection points
## To test stability of cuts
thresh_low <- SubsetByBarcodeDistribution(pbmc_small,
                                         threshold.low = 0)
thresh_high <- SubsetByBarcodeDistribution(pbmc_small,
                                          threshold.high = 10)
thresh_band <- SubsetByBarcodeDistribution(pbmc_small,
                                          threshold.low = 25, threshold.high = 26)

## Tests
test_that("SubsetByBarcodeDistribution default returns a seurat object", {
    expect_s4_class(default, 'Seurat')
})

test_that("SubsetByBarcodeDistribution returns inflections tibble when asked", {
    expect_s3_class(inflections, 'tbl_df')
})

test_that("SubsetByBarcodeDistribution `return.plot = TRUE` returns ggplot object", {
    expect_s3_class(g, 'ggplot')
})

test_that("SubsetByBarcodeDistribution `return.cells = TRUE` returns correct character vector", {
    expect_is(cells, 'character')
    first_barcodes <- c("GACATTCTCCACCT", "TTGAGGACTACGCA", "GCGTAAACACGGTT")
    expect_equal(cells[1:3], first_barcodes)
})

test_that("SubsetByBarcodeDistribution thresholding returns correct number of cells", {
    expect_equal(ncol(thresh_band), 25)
    expect_equal(ncol(thresh_high), 6)
    expect_equal(ncol(thresh_low), 79)
})

