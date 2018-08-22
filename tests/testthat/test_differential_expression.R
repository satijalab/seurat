# Tests for functions in differential_expression.R
set.seed(42)

# Tests for FindMarkers default parameters
# --------------------------------------------------------------------------------
context("Default FindMarkers")

markers.0 <- FindMarkers(object = pbmc_small, ident.1 = 0, print.bar = FALSE)
markers.01 <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, print.bar = FALSE))

test_that("Default settings work as expected", {
  expect_error(FindMarkers(pbmc_small))
  expect_error(FindMarkers(pbmc_small, ident.1 = "test"))
  expect_error(FindMarkers(pbmc_small, ident.1 = 0, ident.2 = "test"))
  expect_equal(colnames(markers.0), c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj"))
  expect_equal(markers.0[1, "p_val"], 1.537519e-11)
  expect_equal(markers.0[1, "avg_logFC"], 2.503489, tolerance = 1e-6)
  expect_equal(markers.0[1, "pct.1"], 1)
  expect_equal(markers.0[1, "pct.2"], 0.333)
  expect_equal(markers.0[1, "p_val_adj"], 3.536293e-09)
  expect_equal(nrow(markers.0), 193)
  expect_equal(rownames(markers.0)[1], "LTB")

  expect_equal(markers.01[1, "p_val"], 2.823522e-08)
  expect_equal(markers.01[1, "avg_logFC"], 2.203405, tolerance = 1e-6)
  expect_equal(markers.01[1, "pct.1"], 1)
  expect_equal(markers.01[1, "pct.2"], 0.414)
  expect_equal(markers.01[1, "p_val_adj"], 0.0000064941)
  expect_equal(nrow(markers.0), 193)
  expect_equal(rownames(markers.0)[1], "LTB")
})


ltb.results <- FindMarkers(object = pbmc_small, ident.1 = 0, genes.use = "LTB", print.bar = FALSE)
vargenes.results <- FindMarkers(object = pbmc_small, ident.1 = 0, genes.use = pbmc_small@var.genes, print.bar = FALSE)

test_that("genes.use parameter behaves correctly ", {
  expect_equal(nrow(ltb.results), 1)
  expect_equal(ltb.results[1, "p_val"], 1.537519e-11)
  expect_equal(ltb.results[1, "avg_logFC"], 2.503489, tolerance = 1e-6)
  expect_equal(ltb.results[1, "pct.1"], 1)
  expect_equal(ltb.results[1, "pct.2"], 0.333)
  expect_equal(ltb.results[1, "p_val_adj"], 3.536293e-09)
  expect_equal(rownames(ltb.results)[1], "LTB")

  expect_equal(nrow(vargenes.results), 25)
  expect_equal(vargenes.results[1, "p_val"], 1.537519e-11)
  expect_equal(vargenes.results[1, "avg_logFC"], 2.503489, tolerance = 1e-6)
  expect_equal(vargenes.results[1, "pct.1"], 1)
  expect_equal(vargenes.results[1, "pct.2"], 0.333)
  expect_equal(vargenes.results[1, "p_val_adj"], 3.536293e-09)
  expect_equal(rownames(vargenes.results)[1], "LTB")
})







