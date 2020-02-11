# Tests for functions in differential_expression.R
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(seed = 42)

# Tests for FindMarkers default parameters
# --------------------------------------------------------------------------------
context("FindMarkers")

markers.0 <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, verbose = FALSE))
markers.01 <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, verbose = FALSE))

test_that("Default settings work as expected", {
  expect_error(FindMarkers(object = pbmc_small))
  expect_error(FindMarkers(object = pbmc_small, ident.1 = "test"))
  expect_error(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = "test"))
  expect_equal(colnames(x = markers.0), c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj"))
  expect_equal(markers.0[1, "p_val"], 9.572778e-13)
  expect_equal(markers.0[1, "avg_logFC"], -4.034691, tolerance = 1e-6)
  expect_equal(markers.0[1, "pct.1"], 0.083)
  expect_equal(markers.0[1, "pct.2"], 0.909)
  expect_equal(markers.0[1, "p_val_adj"], 2.201739e-10)
  expect_equal(nrow(x = markers.0), 204)
  expect_equal(rownames(markers.0)[1], "HLA-DPB1")

  expect_equal(markers.01[1, "p_val"], 1.702818e-11)
  expect_equal(markers.01[1, "avg_logFC"], -2.539289, tolerance = 1e-6)
  expect_equal(markers.01[1, "pct.1"], 0.111)
  expect_equal(markers.01[1, "pct.2"], 1.00)
  expect_equal(markers.01[1, "p_val_adj"], 3.916481e-09)
  expect_equal(nrow(x = markers.01), 201)
  expect_equal(rownames(x = markers.01)[1], "TYMP")
})


tymp.results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, features = "TYMP", verbose = FALSE))
vargenes.results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, features = VariableFeatures(object = pbmc_small), verbose = FALSE))

test_that("features parameter behaves correctly ", {
  expect_equal(nrow(x = tymp.results), 1)
  expect_equal(tymp.results[1, "p_val"], 3.227445e-07)
  expect_equal(tymp.results[1, "avg_logFC"], -2.093928, tolerance = 1e-6)
  expect_equal(tymp.results[1, "pct.1"], 0.111)
  expect_equal(tymp.results[1, "pct.2"], 0.682)
  expect_equal(tymp.results[1, "p_val_adj"], 7.423123e-05)
  expect_equal(rownames(x = tymp.results)[1], "TYMP")

  expect_equal(nrow(x = vargenes.results), 19)
  expect_equal(vargenes.results[19, "p_val"], 4.225151e-01, tolerance = 1e-6)
  expect_equal(vargenes.results[19, "avg_logFC"], 1.5976958, tolerance = 1e-6)
  expect_equal(vargenes.results[19, "pct.1"], 0.139)
  expect_equal(vargenes.results[19, "pct.2"], 0.091)
  expect_equal(vargenes.results[19, "p_val_adj"], 1.000000e+00)
  expect_equal(rownames(x = vargenes.results)[19], "PARVB")
})


results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = Cells(x = pbmc_small)[1:40], ident.2 = Cells(x = pbmc_small)[41:80], verbose = FALSE))
test_that("passing cell names works", {
  expect_equal(nrow(x = results), 176)
  expect_equal(results[1, "p_val"], 0.0001690882)
  expect_equal(results[1, "avg_logFC"], -1.790824, tolerance = 1e-6)
  expect_equal(results[1, "pct.1"], 0.075)
  expect_equal(results[1, "pct.2"], 0.450)
  expect_equal(results[1, "p_val_adj"], 0.03889028)
  expect_equal(rownames(x = results)[1], "IFI30")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, logfc.threshold = 2, verbose = FALSE))
test_that("logfc.threshold works", {
  expect_equal(nrow(x = results), 112)
  expect_gte(min(abs(x = results$avg_logFC)), 2)
  expect_error(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, logfc.threshold = 100, verbose = FALSE))
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, min.pct = 0.5, verbose = FALSE))
test_that("min.pct works", {
  expect_equal(nrow(x = results), 63)
  expect_gte(min(apply(X = results, MARGIN = 1, FUN = function(x) max(x[3], x[4]))), 0.5)
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, min.diff.pct = 0.5, verbose = FALSE))
test_that("min.diff.pct works", {
  expect_equal(nrow(x = results), 43)
  expect_gte(min(apply(X = results, MARGIN = 1, FUN = function(x) abs(x[4] - x[3]))), 0.5)
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, only.pos = TRUE, verbose = FALSE))
test_that("only.pos works", {
  expect_equal(nrow(x = results), 116)
  expect_true(all(results$avg_logFC > 0))
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, max.cells.per.ident = 20, verbose = FALSE))
test_that("max.cells.per.ident works", {
  expect_equal(nrow(x = results), 201)
  expect_equal(results[1, "p_val"], 3.428568e-08)
  expect_equal(results[1, "avg_logFC"], -2.539289, tolerance = 1e-6)
  expect_equal(results[1, "pct.1"], 0.111)
  expect_equal(results[1, "pct.2"], 1)
  expect_equal(results[1, "p_val_adj"], 7.885706e-06)
  expect_equal(rownames(x = results)[1], "TYMP")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, latent.vars= "groups", verbose = FALSE, test.use = 'LR'))
test_that("latent.vars works", {
  expect_error(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, latent.vars= "fake", verbose = FALSE))
  expect_warning(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, latent.vars= "groups", verbose = FALSE))
  expect_equal(nrow(x = results), 201)
  expect_equal(results[1, "p_val"], 2.130202e-16)
  expect_equal(results[1, "avg_logFC"], -3.082150, tolerance = 1e-6)
  expect_equal(results[1, "pct.1"], 0.417)
  expect_equal(results[1, "pct.2"], 1)
  expect_equal(results[1, "p_val_adj"], 4.899466e-14)
  expect_equal(rownames(x = results)[1], "LYZ")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = "g1", ident.2 = "g2", group.by= "groups", verbose = FALSE))
t2 <- pbmc_small
Idents(object = t2) <- "groups"
results2 <- suppressWarnings(FindMarkers(object = t2, ident.1 = "g1", ident.2 = "g2", verbose = FALSE))

test_that("group.by works", {
  expect_equal(nrow(x = results), 136)
  expect_equal(results, results2)
  expect_equal(results[1, "p_val"], 0.02870319)
  expect_equal(results[1, "avg_logFC"], 0.8226720, tolerance = 1e-6)
  expect_equal(results[1, "pct.1"], 0.455)
  expect_equal(results[1, "pct.2"], 0.194)
  expect_equal(results[1, "p_val_adj"], 1)
  expect_equal(rownames(x = results)[1], "NOSIP")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = "g1", ident.2 = "g2", group.by= "groups", subset.ident = 0, verbose = FALSE))
t2 <- subset(x = pbmc_small, idents = 0)
Idents(object = t2) <- "groups"
results2 <- suppressWarnings(FindMarkers(object = t2, ident.1 = "g1", ident.2 = "g2", verbose = FALSE))

test_that("subset.ident works", {
  expect_equal(nrow(x = results), 114)
  expect_equal(results, results2)
  expect_equal(results[1, "p_val"], 0.01293720)
  expect_equal(results[1, "avg_logFC"], 1.799280, tolerance = 1e-6)
  expect_equal(results[1, "pct.1"], 0.50)
  expect_equal(results[1, "pct.2"], 0.125)
  expect_equal(results[1, "p_val_adj"], 1)
  expect_equal(rownames(x = results)[1], "TSPO")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, reduction = "pca", verbose = FALSE))
test_that("reduction works", {
  expect_equal(results[1, "p_val"], 1.664954e-10)
  expect_equal(results[1, "avg_diff"], -2.810453669, tolerance = 1e-6)
  expect_equal(results[1, "p_val_adj"], 3.163412e-09)
  expect_equal(rownames(x = results)[1], "PC_2")
})

results <- FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "bimod", verbose = FALSE)
test_that("bimod test works", {
  expect_equal(nrow(x = results), 201)
  expect_equal(results[1, "p_val"], 4.751376e-17)
  expect_equal(results[1, "avg_logFC"], -2.552769, tolerance = 1e-6)
  expect_equal(results[1, "pct.1"], 0.306)
  expect_equal(results[1, "pct.2"], 1.00)
  expect_equal(results[1, "p_val_adj"], 1.092816e-14)
  expect_equal(rownames(x = results)[1], "CST3")
})

results <- FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "roc", verbose = FALSE)
test_that("roc test works", {
  expect_equal(nrow(x = results), 201)
  expect_equal(colnames(x = results), c("myAUC", "avg_diff", "power", "pct.1", "pct.2"))
  expect_equal(results["CST3", "myAUC"], 0.018)
  expect_equal(results["CST3", "avg_diff"], -2.552769, tolerance = 1e-6)
  expect_equal(results["CST3", "power"], 0.964)
  expect_equal(results["CST3", "pct.1"], 0.306)
  expect_equal(results["CST3", "pct.2"], 1.00)
  expect_equal(rownames(x = results)[1], "LYZ")
})

results <- FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "t", verbose = FALSE)
test_that("bimod test works", {
  expect_equal(nrow(x = results), 201)
  expect_equal(results["CST3", "p_val"], 1.170112e-15)
  expect_equal(results["CST3", "avg_logFC"], -2.552769 , tolerance = 1e-6)
  expect_equal(results["CST3", "pct.1"], 0.306)
  expect_equal(results["CST3", "pct.2"], 1.00)
  expect_equal(results["CST3", "p_val_adj"], 2.691258e-13)
  expect_equal(rownames(x = results)[1], "TYMP")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "negbinom", verbose = FALSE))
test_that("negbinom test works", {
  expect_equal(nrow(x = results), 149)
  expect_equal(results["CST3", "p_val"], 1.354443e-17)
  expect_equal(results["CST3", "avg_logFC"], -2.353701, tolerance = 1e-6)
  expect_equal(results["CST3", "pct.1"], 0.306)
  expect_equal(results["CST3", "pct.2"], 1.00)
  expect_equal(results["CST3", "p_val_adj"], 3.115218e-15)
  expect_equal(rownames(x = results)[1], "LYZ")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "poisson", verbose = FALSE))
test_that("poisson test works", {
  expect_equal(nrow(x = results), 149)
  expect_equal(results["CST3", "p_val"], 3.792196e-78)
  expect_equal(results["CST3", "avg_logFC"], -2.353701, tolerance = 1e-6)
  expect_equal(results["CST3", "pct.1"], 0.306)
  expect_equal(results["CST3", "pct.2"], 1.00)
  expect_equal(results["CST3", "p_val_adj"], 8.722050e-76)
  expect_equal(rownames(x = results)[1], "LYZ")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "LR", verbose = FALSE))
test_that("LR test works", {
  expect_equal(nrow(x = results), 201)
  expect_equal(results["CST3", "p_val"], 3.990707e-16)
  expect_equal(results["CST3", "avg_logFC"], -2.552769, tolerance = 1e-6)
  expect_equal(results["CST3", "pct.1"], 0.306)
  expect_equal(results["CST3", "pct.2"], 1.00)
  expect_equal(results["CST3", "p_val_adj"], 9.178625e-14)
  expect_equal(rownames(x = results)[1], "LYZ")
})

# Tests for FindConservedMarkers
# --------------------------------------------------------------------------------
context("FindConservedMarkers")
pbmc_small$groups

markers <- suppressWarnings(FindConservedMarkers(object = pbmc_small, ident.1 = 0, grouping.var = "groups", verbose = FALSE))

standard.names <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")

test_that("FindConservedMarkers works", {
  expect_equal(colnames(x = markers), c(paste0("g2_", standard.names), paste0("g1_", standard.names), "max_pval", "minimump_p_val"))
  expect_equal(markers[1, "g2_p_val"], 4.983576e-05)
  expect_equal(markers[1, "g2_avg_logFC"], -4.125279, tolerance = 1e-6)
  # expect_equal(markers[1, "g2_pct.1"], 0.062)
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
  # expect_equal(markers.missing[1, "g2_pct.1"], 0.062)
  expect_equal(markers.missing[1, "g2_pct.2"], 0.95)
  expect_equal(markers.missing[1, "g2_p_val_adj"], 3.847695e-11)
  expect_equal(nrow(markers.missing), 190)
  expect_equal(rownames(markers.missing)[1], "HLA-DPB1")
})
