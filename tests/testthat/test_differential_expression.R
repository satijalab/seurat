# Tests for functions in differential_expression.R
set.seed(42)

# Tests for FindConservedMarkers 
# --------------------------------------------------------------------------------
pbmc_small$groups

markers <- suppressWarnings(FindConservedMarkers(object = pbmc_small, ident.1 = 0, grouping.var = "groups", verbose = FALSE))

standard.names <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")

test_that("FindConservedMarkers works", {
  expect_equal(colnames(x = markers), c(paste0("g2_", standard.names), paste0("g1_", standard.names), "max_pval", "minimump_p_val"))
  expect_equal(markers[1, "g2_p_val"], 4.983576e-05)
  expect_equal(markers[1, "g2_avg_logFC"], -4.125279, tolerance = 1e-6)
  expect_equal(markers[1, "g2_pct.1"], 0.062)
  expect_equal(markers[1, "g2_pct.2"], 0.75)
  expect_equal(markers[1, "g2_p_val_adj"], 0.0114622238)
  expect_equal(markers[1, "g1_p_val"], 3.946643e-08)
  expect_equal(markers[1, "g1_avg_logFC"], -3.589384, tolerance = 1e-6)
  expect_equal(markers[1, "g1_pct.1"], 0.10)
  expect_equal(markers[1, "g1_pct.2"], 0.958)
  expect_equal(markers[1, "g1_p_val_adj"], 9.077279e-06)
  expect_equal(markers[1, "max_pval"], 4.983576e-05)
  expect_equal(markers[1, "minimump_p_val"], 7.893286e-08)
  expect_equal(nrow(markers), 162)
  expect_equal(rownames(markers)[1], "HLA-DRB1")
  expect_equal(markers[, "max_pval"], unname(obj = apply(X = markers, MARGIN = 1, FUN = function(x) max(x[c("g1_p_val", "g2_p_val")]))))
})

test_that("FindConservedMarkers errors when expected", {
  expect_error(FindConservedMarkers(pbmc_small))
  expect_error(FindConservedMarkers(pbmc_small, ident.1 = 0))
  expect_error(FindConservedMarkers(pbmc_small, ident.1 = 0, grouping.var = "groups", meta.method = "minimump"))
})

pbmc.test <- pbmc_small
Idents(object = pbmc.test) <- "RNA_snn_res.1"
pbmc.test$id.group <- paste0(pbmc.test$RNA_snn_res.1, "_", pbmc.test$groups)
pbmc.test <- subset(x = pbmc.test, id.group == "0_g1", invert = TRUE)
markers.missing <- suppressWarnings(FindConservedMarkers(object = pbmc.test, ident.1 = 0, grouping.var = "groups", test.use = "t", verbose = FALSE))

test_that("FindConservedMarkers handles missing idents in certain groups", {
  expect_warning(FindConservedMarkers(object = pbmc.test, ident.1 = 0, grouping.var = "groups", test.use = "t"))
  expect_equal(colnames(x = markers.missing), paste0("g2_", standard.names))
  expect_equal(markers.missing[1, "g2_p_val"], 1.672911e-13)
  expect_equal(markers.missing[1, "g2_avg_logFC"], -4.527888, tolerance = 1e-6)
  expect_equal(markers.missing[1, "g2_pct.1"], 0.062)
  expect_equal(markers.missing[1, "g2_pct.2"], 0.95)
  expect_equal(markers.missing[1, "g2_p_val_adj"], 3.847695e-11)
  expect_equal(nrow(markers.missing), 190)
  expect_equal(rownames(markers.missing)[1], "HLA-DPB1")
})
