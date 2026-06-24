# Tests for functions in differential_expression.R
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(seed = 42)
is_not_cran_submission <- isTRUE(as.logical(Sys.getenv("NOT_CRAN")))

# Tests for FindMarkers
# --------------------------------------------------------------------------------
context("FindMarkers")

# tests focus on output shape and known marker recovery, rather than exact p-values or top-ranked rows
# since they can vary across R versions, especially for Wilcoxon tests
expect_de_table <- function(results, expected.cols = c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")) {
  expect_equal(colnames(x = results), expected.cols)
  expect_true(nrow(x = results) > 0)
  expect_true(all(results$p_val >= 0 & results$p_val <= 1))
  expect_true(all(results$p_val_adj >= 0 & results$p_val_adj <= 1))
  if ("pct.1" %in% colnames(x = results)) {
    expect_true(all(results$pct.1 >= 0 & results$pct.1 <= 1))
  }
  if ("pct.2" %in% colnames(x = results)) {
    expect_true(all(results$pct.2 >= 0 & results$pct.2 <= 1))
  }
  expect_true(all(is.finite(results[[grep(pattern = "^avg_", x = colnames(x = results), value = TRUE)[1]]])))
}

if (is_not_cran_submission) {
  clr.obj <- suppressWarnings(NormalizeData(pbmc_small, normalization.method = "CLR"))
  sct.obj <- suppressWarnings(suppressMessages(SCTransform(pbmc_small, vst.flavor = "v1")))

  markers.0 <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, verbose = FALSE, base = exp(1),pseudocount.use = 1))
  markers.01 <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1),pseudocount.use = 1))
  results.clr <- suppressWarnings(FindMarkers(object = clr.obj, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  results.sct <- suppressWarnings(FindMarkers(object = sct.obj, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 1))

  test_that("Default settings work as expected with pseudocount = 1", {
    expect_error(FindMarkers(object = pbmc_small))
    expect_error(FindMarkers(object = pbmc_small, ident.1 = "test"))
    expect_error(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = "test"))

    # default tests focus on output shape and known marker recovery
    expect_de_table(markers.0)
    expect_equal(nrow(x = markers.0), 228)
    expect_true("HLA-DPB1" %in% rownames(x = markers.0))

    expect_de_table(markers.01)
    expect_equal(nrow(x = markers.01), 222)
    expect_true("TYMP" %in% rownames(x = markers.01))
    expect_lt(markers.01["TYMP", "avg_logFC"], 0)
    expect_gt(markers.01["TYMP", "pct.2"], markers.01["TYMP", "pct.1"])

    # CLR normalization
    expect_de_table(results.clr)
    expect_equal(nrow(x = results.clr), 213)
    expect_true("S100A8" %in% rownames(x = results.clr))

    # SCT normalization
    expect_de_table(results.sct)
    expect_equal(nrow(x = results.sct), 214)
    expect_true("TYMP" %in% rownames(x = results.sct))
  })

  tymp.results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, features = "TYMP", verbose = FALSE, base = exp(1),pseudocount.use = 1))
  vargenes.results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, features = VariableFeatures(object = pbmc_small), verbose = FALSE, base = exp(1),pseudocount.use = 1))

  test_that("features parameter behaves correctly ", {
    expect_equal(nrow(x = tymp.results), 1)
    expect_de_table(tymp.results)
    expect_equal(tymp.results[1, "avg_logFC"], -2.188179, tolerance = 1e-6)
    expect_equal(tymp.results[1, "pct.1"], 0.111)
    expect_equal(tymp.results[1, "pct.2"], 0.682)
    expect_equal(rownames(x = tymp.results)[1], "TYMP")

    expect_equal(nrow(x = vargenes.results), 20)
    expect_de_table(vargenes.results)
    # feature-filtering should return exactly the requested variable features
    # (rank order can shift slightly with ties)
    expect_setequal(rownames(x = vargenes.results), VariableFeatures(object = pbmc_small))
  })


  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = Cells(x = pbmc_small)[1:40], ident.2 = Cells(x = pbmc_small)[41:80], verbose = FALSE, base = exp(1),pseudocount.use = 1))
  test_that("passing cell names works", {
    expect_equal(nrow(x = results), 216)
    expect_de_table(results)
    expect_equal(results[1, "avg_logFC"], -1.967123, tolerance = 1e-6)
    expect_equal(results[1, "pct.1"], 0.075)
    expect_equal(results[1, "pct.2"], 0.450)
    expect_equal(rownames(x = results)[1], "IFI30")
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 0.1))
  results.clr <- suppressWarnings(FindMarkers(object = clr.obj, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 0.1))
  results.sct <- suppressWarnings(FindMarkers(object = sct.obj, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 0.1, vst.flavor = "v1"))
  # different pseudocount should change the fold-change calculation, but not the number of rows returned
  test_that("setting pseudocount.use works", {
    expect_equal(nrow(x = results), 222)
    expect_de_table(results)
    expect_false(isTRUE(all.equal(results$avg_logFC, markers.01[rownames(x = results), "avg_logFC"])))
    expect_equal(nrow(x = results.clr), 214)
    expect_de_table(results.clr)
    expect_equal(nrow(results.sct), 215)
    expect_de_table(results.sct)
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 1, mean.fxn = rowMeans))
  results.clr <- suppressWarnings(FindMarkers(object = clr.obj, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 1, mean.fxn = rowMeans))
  results.sct <- suppressWarnings(FindMarkers(object = sct.obj, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 1, mean.fxn = rowMeans, vst.flaovr = "v1"))
  # different mean function should change the fold-change calculation, but not the number of rows returned
  test_that("setting mean.fxn works", {
    expect_equal(nrow(x = results), 216)
    expect_de_table(results)
    expect_de_table(results.clr)
    expect_de_table(results.sct)
    expect_false(isTRUE(all.equal(results$avg_logFC, markers.01[rownames(x = results), "avg_logFC"])))
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, logfc.threshold = 2, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  test_that("logfc.threshold works", {
    expect_equal(nrow(x = results), 139)
    expect_gte(min(abs(x = results$avg_logFC)), 2)
  })

  results <- expect_warning(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, logfc.threshold = 100, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  test_that("logfc.threshold warns when none met", {
    expect_equal(nrow(x = results), 0)
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, min.pct = 0.5, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  test_that("min.pct works", {
    expect_equal(nrow(x = results), 66)
    expect_gte(min(apply(X = results, MARGIN = 1, FUN = function(x) max(x[3], x[4]))), 0.5)
  })

  results <- expect_warning(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, min.pct = 2.0, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  test_that("min.pct warns when none met", {
    expect_equal(nrow(x = results), 0)
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, min.diff.pct = 0.5, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  test_that("min.diff.pct works", {
    expect_equal(nrow(x = results), 44)
    expect_gte(min(apply(X = results, MARGIN = 1, FUN = function(x) abs(x[4] - x[3]))), 0.5)
  })

  results <- expect_warning(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, min.diff.pct = 1.0, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  test_that("min.diff.pct warns when none met", {
    expect_equal(nrow(x = results), 0)
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, only.pos = TRUE, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  test_that("only.pos works", {
    expect_equal(nrow(x = results), 127)
    expect_true(all(results$avg_logFC > 0))
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, max.cells.per.ident = 20, verbose = FALSE, base = exp(1),pseudocount.use = 1))
  results.repeated <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, max.cells.per.ident = 20, verbose = FALSE, base = exp(1),pseudocount.use = 1))
  results.unsampled <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1),pseudocount.use = 1))
  test_that("max.cells.per.ident works", {
    expect_equal(nrow(x = results), 222)
    # downsampling uses FindMarkers' default random.seed, so identical
    # calls setting max.cells.per.ident should be reproducible
    expect_equal(results, results.repeated)
    # the max.cells.per.ident test should not collapse to the unsampled result
    expect_false(identical(results, results.unsampled))

    # fold-change and pct columns are computed before downsampling, while
    # p-values are computed after the sampled cells are selected
    # so these should be identical, while p-values should differ
    expect_equal(
      results[, c("avg_logFC", "pct.1", "pct.2")],
      results.unsampled[rownames(x = results), c("avg_logFC", "pct.1", "pct.2")],
      tolerance = 1e-6
    )
    expect_false(isTRUE(all.equal(
      results$p_val,
      results.unsampled[rownames(x = results), "p_val"]
    )))
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, latent.vars= "groups", verbose = FALSE, test.use = 'LR', base = exp(1), pseudocount.use = 1))
  test_that("latent.vars works", {
    expect_error(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, latent.vars= "fake", verbose = FALSE))
    expect_warning(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, latent.vars= "groups", verbose = FALSE))
    expect_equal(nrow(x = results), 222)
    expect_equal(results[1, "p_val"], 2.130202e-16, tolerance = 1e-21)
    expect_equal(results[1, "avg_logFC"], -3.102866, tolerance = 1e-6)
    expect_equal(results[1, "pct.1"], 0.417)
    expect_equal(results[1, "pct.2"], 1)
    expect_equal(results[1, "p_val_adj"], 4.899466e-14, tolerance = 1e-19)
    expect_equal(rownames(x = results)[1], "LYZ")
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = "g1", ident.2 = "g2", group.by= "groups", verbose = FALSE, base = exp(1), pseudocount.use = 1))
  t2 <- pbmc_small
  Idents(object = t2) <- "groups"
  results2 <- suppressWarnings(FindMarkers(object = t2, ident.1 = "g1", ident.2 = "g2", verbose = FALSE, base = exp(1), pseudocount.use = 1))

  test_that("group.by works", {
    expect_equal(nrow(x = results), 190)
    # group.by should be equivalent to setting Idents before FindMarkers
    expect_equal(results, results2)
    expect_de_table(results)
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = "g1", ident.2 = "g2", group.by= "groups", subset.ident = 0, verbose = FALSE, base = exp(1), pseudocount.use = 1))
  t2 <- subset(x = pbmc_small, idents = 0)
  Idents(object = t2) <- "groups"
  results2 <- suppressWarnings(FindMarkers(object = t2, ident.1 = "g1", ident.2 = "g2", verbose = FALSE, base = exp(1), pseudocount.use = 1))

  test_that("subset.ident works", {
    expect_equal(nrow(x = results), 183)
    # subset.ident should be equivalent to subsetting first, then setting Idents
    expect_equal(results, results2)
    expect_de_table(results)
    expect_equal(rownames(x = results)[1], "TSPO")
  })

  results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, reduction = "pca", verbose = FALSE, base = exp(1), pseudocount.use = 1))
  test_that("reduction works", {
    expect_de_table(results, expected.cols = c("p_val", "avg_diff", "p_val_adj"))
    expect_equal(results[1, "avg_diff"], -2.810453669, tolerance = 1e-6)
    expect_equal(rownames(x = results)[1], "PC_2")
  })

  results <- FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "bimod", verbose = FALSE, base = exp(1), pseudocount.use = 1)
  test_that("bimod test works", {
    expect_equal(nrow(x = results), 222)
    expect_equal(results[1, "p_val"], 4.751376e-17, tolerance = 1e-22)
    expect_equal(results[1, "avg_logFC"], -2.57219, tolerance = 1e-6)
    expect_equal(results[1, "pct.1"], 0.306)
    expect_equal(results[1, "pct.2"], 1.00)
    expect_equal(results[1, "p_val_adj"], 1.092816e-14, tolerance = 1e-19)
    expect_equal(rownames(x = results)[1], "CST3")
  })

  results <- FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "roc", verbose = FALSE, base = exp(1), pseudocount.use = 1)
  test_that("roc test works", {
    expect_equal(nrow(x = results), 222)
    # expect_equal(colnames(x = results), c("myAUC", "avg_diff", "power", "pct.1", "pct.2"))
    expect_equal(colnames(x = results), c("myAUC", "avg_diff", "power", "avg_logFC", "pct.1", "pct.2"))
    expect_equal(results["CST3", "myAUC"], 0.018)
    expect_equal(results["CST3", "avg_diff"], -2.552769, tolerance = 1e-6)
    expect_equal(results["CST3", "power"], 0.964)
    expect_equal(results["CST3", "pct.1"], 0.306)
    expect_equal(results["CST3", "pct.2"], 1.00)
    expect_equal(rownames(x = results)[1], "LYZ")
  })

  results <- FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "t", verbose = FALSE, base = exp(1), pseudocount.use = 1)
  test_that("t test works", {
    expect_equal(nrow(x = results), 222)
    expect_equal(results["CST3", "p_val"], 1.170112e-15, tolerance = 1e-20)
    expect_equal(results["CST3", "avg_logFC"], -2.57219, tolerance = 1e-6)
    expect_equal(results["CST3", "pct.1"], 0.306)
    expect_equal(results["CST3", "pct.2"], 1.00)
    expect_equal(results["CST3", "p_val_adj"], 2.691258e-13, tolerance = 1e-18)
    expect_equal(rownames(x = results)[1], "TYMP")
  })

  test_that("FindMarkers with wilcox_limma works", {
    skip_if_not_installed("limma")
    markers.0.limma <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, verbose = FALSE, base = exp(1),pseudocount.use = 1,test.use='wilcox_limma'))
    markers.01.limma <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1),pseudocount.use = 1,test.use='wilcox_limma'))
    results.clr.limma <- suppressWarnings(FindMarkers(object = clr.obj, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 1,test.use='wilcox_limma'))
    results.sct.limma <- suppressWarnings(FindMarkers(object = sct.obj, ident.1 = 0, ident.2 = 1, verbose = FALSE, base = exp(1), pseudocount.use = 1,test.use='wilcox_limma'))

    expect_equal(colnames(x = markers.0.limma), c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj"))
    expect_equal(markers.0.limma[1, "p_val"], 9.572778e-13, tolerance = 1e-18)
    expect_equal(markers.0.limma[1, "avg_logFC"], -4.180029, tolerance = 1e-6)
    expect_equal(markers.0.limma[1, "pct.1"], 0.083)
    expect_equal(markers.0.limma[1, "pct.2"], 0.909)
    expect_equal(markers.0.limma[1, "p_val_adj"], 2.201739e-10, tolerance = 1e-15)
    expect_equal(nrow(x = markers.0.limma), 228)
    expect_equal(rownames(markers.0.limma)[1], "HLA-DPB1")

    expect_equal(markers.01.limma[1, "p_val"], 1.702818e-11, tolerance = 1e-16)
    expect_equal(markers.01.limma[1, "avg_logFC"], -2.638242, tolerance = 1e-6)
    expect_equal(markers.01.limma[1, "pct.1"], 0.111)
    expect_equal(markers.01.limma[1, "pct.2"], 1.00)
    expect_equal(markers.01.limma[1, "p_val_adj"], 3.916481e-09, tolerance = 1e-14)
    expect_equal(nrow(x = markers.01.limma), 222)
    expect_equal(rownames(x = markers.01.limma)[1], "TYMP")

    expect_equal(results.clr.limma[1, "p_val"], 1.209462e-11, tolerance = 1e-16)
    expect_equal(results.clr.limma[1, "avg_logFC"], -2.946633, tolerance = 1e-6)
    expect_equal(results.clr.limma[1, "pct.1"], 0.111)
    expect_equal(results.clr.limma[1, "pct.2"], 0.96)
    expect_equal(results.clr.limma[1, "p_val_adj"], 2.781762e-09, tolerance = 1e-14)
    expect_equal(nrow(x = results.clr.limma), 213)
    expect_equal(rownames(x = results.clr.limma)[1], "S100A8")

    expect_equal(results.sct.limma[1, "p_val"], 6.225491e-11, tolerance = 1e-16)
    expect_equal(results.sct.limma[1, "avg_logFC"], -2.545867, tolerance = 1e-6)
    expect_equal(results.sct.limma[1, "pct.1"], 0.111)
    expect_equal(results.sct.limma[1, "pct.2"], 0.96)
    expect_equal(results.sct.limma[1, "p_val_adj"], 1.369608e-08, tolerance = 1e-13)
    expect_equal(nrow(x = results.sct.limma), 214)
    expect_equal(rownames(x = results.sct.limma)[1], "TYMP")
  })

  test_that("BPCells FindMarkers gives same results", {
    skip_if_not_installed("BPCells")
    library(BPCells)
    library(Matrix)
    mat_bpcells <- t(as(t(pbmc_small[['RNA']]$counts ), "IterableMatrix"))
    pbmc_small[['RNAbp']] <- CreateAssay5Object(counts = mat_bpcells)
    pbmc_small <- NormalizeData(pbmc_small, assay = "RNAbp")
    markers.bp <- suppressWarnings(FindMarkers(object = pbmc_small, assay = "RNAbp", ident.1 = 0, verbose = FALSE, base = exp(1),pseudocount.use = 1))
    expect_equal(colnames(x = markers.bp), c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj"))
    expect_equal(markers.bp[1, "p_val"], 9.572778e-13)
    expect_equal(markers.bp[1, "avg_logFC"], -4.180029, tolerance = 1e-6)
    expect_equal(markers.bp[1, "pct.1"], 0.083)
    expect_equal(markers.bp[1, "pct.2"], 0.909)
    expect_equal(markers.bp[1, "p_val_adj"], 2.201739e-10)
    expect_equal(nrow(x = markers.bp), 228)
    expect_equal(rownames(markers.bp)[1], "HLA-DPB1")
  })

  # set up a multi-model SCT object for PrepSCTFindMarkers tests
  sct1 <- suppressWarnings(SCTransform(pbmc_small[, 1:40], variable.features.n = 20, vst.flavor = "v1", verbose = FALSE, seed.use = 123))
  sct2 <- suppressWarnings(SCTransform(pbmc_small[, 41:80], variable.features.n = 20, vst.flavor = "v1", verbose = FALSE, seed.use = 123))
  sct_merged <- merge(x = sct1, y = sct2)

  test_that("PrepSCTFindMarkers matches correct_counts for unchanged SCT models", {
    model_name <- "model1.1"
    sct <- sct_merged[["SCT"]]
    model <- slot(sct, name = "SCTModel.list")[[model_name]]
    cell_attributes <- slot(model, name = "cell.attributes")
    model_cells <- rownames(cell_attributes)
    feature_attributes <- slot(model, name = "feature.attributes")

    valid_genes <- rownames(feature_attributes)[!is.na(feature_attributes[, "theta"]) &
                                                is.finite(feature_attributes[, "(Intercept)"]) &
                                                is.finite(feature_attributes[, "log_umi"])]

    raw_counts <- GetAssayData(sct_merged, assay = SCTResults(sct, slot = "umi.assay")[[model_name]], layer = "counts")

    expected_counts <- sctransform::correct_counts(
      list(model_str = SCTResults(sct, slot = "model")[[model_name]],
        arguments = SCTResults(sct, slot = "arguments")[[model_name]],
        model_pars_fit = as.matrix(feature_attributes[valid_genes, c("theta", "(Intercept)", "log_umi"), drop = FALSE]),
        cell_attr = cell_attributes),
      umi = raw_counts[valid_genes, model_cells, drop = FALSE],
      verbosity = 0,
      scale_factor = min(unlist(lapply(SCTResults(sct, slot = "cell.attributes"), function(x) median(x[, "umi"]))))
    )

    result <- PrepSCTFindMarkers(sct_merged, verbose = FALSE)

    # PrepSCTFindMarkers recorrection should match sctransform::correct_counts at the same shared median UMI depth
    result_counts <- GetAssayData(result, assay = "SCT", layer = "counts")[rownames(expected_counts), model_cells, drop = FALSE]
    result_data <- GetAssayData(result, assay = "SCT", layer = "data")[rownames(expected_counts), model_cells, drop = FALSE]

    expect_equal(as.matrix(result_counts), as.matrix(expected_counts))
    expect_equal(as.matrix(result_data), as.matrix(log1p(expected_counts)))
  })

  test_that("PrepSCTFindMarkers keeps infinite theta genes with valid corrected counts", {
    sct_test <- sct_merged
    model_name <- "model1.1"
    model <- slot(sct_test[["SCT"]], name = "SCTModel.list")[[model_name]]
    cell_attributes <- slot(model, name = "cell.attributes")
    model_cells <- rownames(cell_attributes)
    feature_attributes <- slot(model, name = "feature.attributes")
    counts <- GetAssayData(sct_test, assay = "SCT", layer = "counts")
    before_counts <- counts[rownames(feature_attributes), model_cells, drop = FALSE]
    nonzero_genes <- names(sort(Matrix::rowSums(before_counts > 0), decreasing = TRUE))

    # set a single gene with finite intercept and log_umi to have infinite theta (valid)
    infinite_theta_gene <- nonzero_genes[is.finite(feature_attributes[nonzero_genes, "(Intercept)"]) &
                                        is.finite(feature_attributes[nonzero_genes, "log_umi"])][[1]]
    feature_attributes[infinite_theta_gene, "theta"] <- Inf
    slot(slot(sct_test[["SCT"]], name = "SCTModel.list")[[model_name]], name = "feature.attributes") <- feature_attributes

    result <- PrepSCTFindMarkers(sct_test, verbose = FALSE)
    result_counts <- GetAssayData(result, assay = "SCT", layer = "counts")
    result_counts <- result_counts[infinite_theta_gene, model_cells, drop = FALSE]

    # check that infinite theta is retained, and corrected counts are valid (not all zero or non-finite)
    expect_false(is.finite(feature_attributes[infinite_theta_gene, "theta"]))
    expect_gt(Matrix::rowSums(before_counts[infinite_theta_gene, , drop = FALSE] > 0), expected = 0)
    expect_gt(Matrix::rowSums(result_counts > 0), expected = 0)
  })

  test_that("PrepSCTFindMarkers excludes missing theta genes", {
    sct_test <- sct_merged
    model_name <- "model1.1"
    model <- slot(sct_test[["SCT"]], name = "SCTModel.list")[[model_name]]
    cell_attributes <- slot(model, name = "cell.attributes")
    model_cells <- rownames(cell_attributes)
    feature_attributes <- slot(model, name = "feature.attributes")
    counts <- GetAssayData(sct_test, assay = "SCT", layer = "counts")

    before_counts <- counts[rownames(feature_attributes), model_cells, drop = FALSE]

    nonzero_genes <- names(sort(Matrix::rowSums(before_counts > 0), decreasing = TRUE))

    # set a single gene with finite intercept and log_umi to have missing/NA theta (invalid)
    # this gene should be filtered out in PrepSCTFindMarkers, then added back in the original order w/ all zero corrected counts
    missing_theta_gene <- nonzero_genes[is.finite(feature_attributes[nonzero_genes, "(Intercept)"]) &
                                        is.finite(feature_attributes[nonzero_genes, "log_umi"])][[1]]

    feature_attributes[missing_theta_gene, "theta"] <- NA
    slot(slot(sct_test[["SCT"]], name = "SCTModel.list")[[model_name]], name = "feature.attributes") <- feature_attributes

    result <- PrepSCTFindMarkers(sct_test, verbose = FALSE)
    result_counts <- GetAssayData(result, assay = "SCT", layer = "counts")
    result_counts <- result_counts[missing_theta_gene, model_cells, drop = FALSE]

    # check that missing theta is excluded, and corrected counts are all zero for this gene
    expect_equal(unlist(unname(Matrix::rowSums(result_counts > 0))), expected = 0)

    # check that it is in the same order as the original assay and has finite counts (not NaN or Inf)
    expect_identical(rownames(result_counts), missing_theta_gene)
    expect_true(all(is.finite(result_counts@x)))
  })

  test_that("PrepSCTFindMarkers drops rows with NaN corrected counts", {
    sct_test <- sct_merged
    model_name <- "model1.1"
    model <- slot(sct_test[["SCT"]], name = "SCTModel.list")[[model_name]]
    cell_attributes <- slot(model, name = "cell.attributes")
    model_cells <- rownames(cell_attributes)
    feature_attributes <- slot(model, name = "feature.attributes")
    counts <- GetAssayData(sct_test, assay = "SCT", layer = "counts")
    before_counts <- counts[rownames(feature_attributes), model_cells, drop = FALSE]
    nonzero_genes <- names(sort(Matrix::rowSums(before_counts > 0), decreasing = TRUE))

    # set a single gene with finite intercept and log_umi to have theta = 0, which causes sctransform::correct_counts to produce NaN values for this gene
    nan_counts_gene <- nonzero_genes[is.finite(feature_attributes[nonzero_genes, "(Intercept)"]) &
                                    is.finite(feature_attributes[nonzero_genes, "log_umi"])][[3]]
    feature_attributes[nan_counts_gene, "theta"] <- 0

    slot(slot(sct_test[["SCT"]], name = "SCTModel.list")[[model_name]], name = "feature.attributes") <- feature_attributes

    result <- PrepSCTFindMarkers(sct_test, verbose = FALSE)
    result_counts <- GetAssayData(result, assay = "SCT", layer = "counts")

    result_counts <- result_counts[nan_counts_gene, model_cells, drop = FALSE]

    # check that the gene with NaN corrected counts is dropped, and thus has zero non-zero counts in the result
    expect_equal(unlist(unname(Matrix::rowSums(result_counts > 0))), expected = 0)
    expect_true(all(is.finite(result_counts@x)))
  })
}

results <- suppressWarnings(
  FindMarkers(
    object = pbmc_small,
    ident.1 = 0, ident.2 = 1,
    test.use = "negbinom",
    verbose = FALSE,
    base = exp(1),
    fc.slot = "counts",
    pseudocount.use = 1
  )
)
test_that("negbinom test works", {
  expect_equal(nrow(x = results), 204)
  expect_equal(results["CST3", "p_val"], 1.354443e-17, tolerance = 1e-22)
  expect_equal(results["CST3", "avg_logFC"], -2.878123, tolerance = 1e-6)
  expect_equal(results["CST3", "pct.1"], 0.306)
  expect_equal(results["CST3", "pct.2"], 1.00)
  expect_equal(results["CST3", "p_val_adj"], 3.115218e-15, tolerance = 1e-20)
  expect_equal(rownames(x = results)[1], "LYZ")
})

results <- suppressWarnings(
  FindMarkers(
    object = pbmc_small,
    ident.1 = 0,
    ident.2 = 1,
    test.use = "poisson",
    verbose = FALSE,
    base = exp(1),
    fc.slot = "counts",
    pseudocount.use = 1
  )
)
test_that("poisson test works", {
  expect_equal(nrow(x = results), 204)
  expect_equal(results["CST3", "p_val"], 3.792196e-78, tolerance = 1e-83)
  expect_equal(results["CST3", "avg_logFC"], -2.878123, tolerance = 1e-6)
  expect_equal(results["CST3", "pct.1"], 0.306)
  expect_equal(results["CST3", "pct.2"], 1.00)
  expect_equal(results["CST3", "p_val_adj"], 8.722050e-76, tolerance = 1e-81)
  expect_equal(rownames(x = results)[1], "LYZ")
})

results <- suppressWarnings(FindMarkers(object = pbmc_small, ident.1 = 0, ident.2 = 1, test.use = "LR", verbose = FALSE, base = exp(1), pseudocount.use = 1))
test_that("LR test works", {
  expect_equal(nrow(x = results), 222)
  expect_equal(results["CST3", "p_val"], 3.990707e-16, tolerance = 1e-21)
  expect_equal(results["CST3", "avg_logFC"], -2.57219, tolerance = 1e-6)
  expect_equal(results["CST3", "pct.1"], 0.306)
  expect_equal(results["CST3", "pct.2"], 1.00)
  expect_equal(results["CST3", "p_val_adj"], 9.178625e-14, tolerance = 1e-19)
  expect_equal(rownames(x = results)[1], "LYZ")
})

# Tests for FindAllMarkers
# -------------------------------------------------------------------------------

if (is_not_cran_submission) {
  test_that("FindAllMarkers works as expected", {
    pbmc_copy <- pbmc_small
    Idents(pbmc_copy) <- "orig.ident"

    results <- suppressMessages(suppressWarnings(FindAllMarkers(object = pbmc_small, pseudocount.use = 1)))
    results.clr <- suppressMessages(suppressWarnings(FindAllMarkers(object = clr.obj, pseudocount.use = 1)))
    results.sct <- suppressMessages(suppressWarnings(FindAllMarkers(object = sct.obj, pseudocount.use = 1, vst.flavor = "v1")))
    results.pseudo <- suppressMessages(suppressWarnings(FindAllMarkers(object = pbmc_small, pseudocount.use = 0.1)))
    results.gb <- suppressMessages(suppressWarnings(FindAllMarkers(object = pbmc_copy, pseudocount.use = 1, group.by = "RNA_snn_res.1")))

    # FindAllMarkers aggregates per-cluster Wilcoxon results
    # can be sensitive to changes in underlying Wilcoxon & RNG
    # instead, check output table and that a known marker is present in the results
    expect_de_table(results, expected.cols = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"))
    expect_gt(nrow(x = results), 200)
    expect_true("HLA-DPB1" %in% results$gene)

    # CLR normalization
    expect_de_table(results.clr, expected.cols = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"))
    expect_gt(nrow(x = results.clr), 200)
    expect_true("HLA-DPB1" %in% results.clr$gene)

    # SCT normalization
    expect_de_table(results.sct, expected.cols = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"))
    expect_gt(nrow(x = results.sct), 200)
    expect_true("HLA-DPB1" %in% results.sct$gene)

    # pseudocount.use = 0.1
    expect_de_table(results.pseudo, expected.cols = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"))
    expect_gt(nrow(x = results.pseudo), 200)
    expect_true("HLA-DPB1" %in% results.pseudo$gene)

    # Setting `group.by` the group by parameter is equivalent
    # to setting the object's `Idents` before running `FindAllMarkers`.
    expect_equal(results.gb, results)
  })
  test_that("BPCells FindAllMarkers gives same results", {
    skip_if_not_installed("BPCells")
    library(BPCells)
    library(Matrix)
    mat_bpcells <- t(as(t(pbmc_small[['RNA']]$counts ), "IterableMatrix"))
    pbmc_small[['RNAbp']] <- CreateAssay5Object(counts = mat_bpcells)
    pbmc_small <- NormalizeData(pbmc_small, assay = "RNAbp")

    results.bp <- suppressMessages(suppressWarnings(FindAllMarkers(object = pbmc_small, assay = "RNAbp", pseudocount.use=1)))

    # BPCells should return the same table as the in-memory matrix path
    expect_de_table(results.bp, expected.cols = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"))
    expect_gt(nrow(x = results.bp), 200)
    expect_true("HLA-DPB1" %in% results.bp$gene)
  })

  test_that("BPCells FindAllMarkers warns for column-major storage", {
    skip_if_not_installed("BPCells")
    library(BPCells)
    library(Matrix)

    mat_bpcells_row <- t(as(t(pbmc_small[['RNA']]$counts), "IterableMatrix"))
    mat_bpcells_col <- BPCells::transpose_storage_order(mat_bpcells_row)

    expect_equal(BPCells::storage_order(mat_bpcells_col), "col")

    pbmc_small[['RNAbpRowMajor']] <- CreateAssay5Object(counts = mat_bpcells_row)
    pbmc_small[['RNAbpColMajor']] <- CreateAssay5Object(counts = mat_bpcells_col)
    pbmc_small <- NormalizeData(pbmc_small, assay = "RNAbpRowMajor")
    pbmc_small <- NormalizeData(pbmc_small, assay = "RNAbpColMajor")

    fam.results.row <- suppressMessages(
      suppressWarnings(
        FindAllMarkers(object = pbmc_small, assay = "RNAbpRowMajor", pseudocount.use = 1)
      )
    )
    fam.results.col <- suppressMessages(
      expect_warning(
        FindAllMarkers(object = pbmc_small, assay = "RNAbpColMajor", pseudocount.use = 1),
        regexp = "Column-major order detected"
      )
    )
    fm.results.col.2 <- suppressMessages(
      expect_warning(
          FindMarkers(object = pbmc_small, assay = "RNAbpColMajor", ident.1 = 0, ident.2 = 1),
          regexp = "Column-major order detected"
      )
    )
    expect_equal(fam.results.col$gene, fam.results.row$gene)
    expect_equal(fam.results.col$cluster, fam.results.row$cluster)
  })
}

# Tests for running FindMarkers post integration/transfer
ref <- pbmc_small
ref <- FindVariableFeatures(object = ref, verbose = FALSE, nfeatures = 100)
query <- CreateSeuratObject(CreateAssayObject(
  counts = as.sparse(GetAssayData(object = pbmc_small[['RNA']], layer = "counts") + rpois(n = ncol(pbmc_small), lambda = 1))
))

query2 <- CreateSeuratObject(CreateAssayObject(
  counts = as.sparse(GetAssayData(object = pbmc_small[['RNA']], layer = "counts")[, 1:40] + rpois(n = ncol(pbmc_small), lambda = 1))
))

query.list <- list(query, query2)
query.list <- lapply(X = query.list, FUN = NormalizeData, verbose = FALSE)
query.list <- lapply(X = query.list, FUN = FindVariableFeatures, verbose = FALSE, nfeatures = 100)
query.list <- lapply(X = query.list, FUN = ScaleData, verbose = FALSE)
query.list <- suppressWarnings(lapply(X = query.list, FUN = RunPCA, verbose = FALSE, npcs = 20))

anchors <- suppressMessages(suppressWarnings(FindIntegrationAnchors(object.list = c(ref, query.list), k.filter = NA, verbose = FALSE)))
object <- suppressWarnings(suppressMessages(IntegrateData(anchorset = anchors,  k.weight = 25, verbose = FALSE)))
object <- suppressMessages(ScaleData(object, verbose = FALSE))
object <- suppressMessages(RunPCA(object, verbose = FALSE))
object <- suppressMessages(FindNeighbors(object = object, verbose = FALSE))
object <- suppressMessages(FindClusters(object, verbose = FALSE))
markers <- FindMarkers(object = object, ident.1="0", ident.2="1",pseudocount.use = 1, verbose=FALSE)

test_that("FindMarkers recognizes log normalization", {
  expect_equal(markers[1, "p_val"], 1.598053e-14, tolerance = 1e-19)
  expect_equal(markers[1, "avg_log2FC"], -2.634458, tolerance = 1e-6)
})

# Tests for FindConservedMarkers
# -------------------------------------------------------------------------------

if (requireNamespace('metap', quietly = TRUE)) {
  context("FindConservedMarkers")
  pbmc_small$groups

  markers <- suppressWarnings(FindConservedMarkers(object = pbmc_small, ident.1 = 0, grouping.var = "groups", verbose = FALSE, base = exp(1), pseudocount.use = 1))

  standard.names <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")

  test_that("FindConservedMarkers works", {
    expect_equal(colnames(x = markers), c(paste0("g2_", standard.names), paste0("g1_", standard.names), "max_pval", "minimump_p_val"))
    expect_equal(markers[1, "g2_p_val"], 4.983576e-05)
    expect_equal(markers[1, "g2_avg_logFC"], -4.364959, tolerance = 1e-6)
    # expect_equal(markers[1, "g2_pct.1"], 0.062)
    expect_equal(markers[1, "g2_pct.2"], 0.75)
    expect_equal(markers[1, "g2_p_val_adj"], 0.0114622238)
    expect_equal(markers[1, "g1_p_val"], 3.946643e-08, tolerance = 1e-13)
    expect_equal(markers[1, "g1_avg_logFC"], -3.69215, tolerance = 1e-6)
    expect_equal(markers[1, "g1_pct.1"], 0.10)
    expect_equal(markers[1, "g1_pct.2"], 0.958)
    expect_equal(markers[1, "g1_p_val_adj"], 9.077279e-06)
    expect_equal(markers[1, "max_pval"], 4.983576e-05)
    expect_equal(markers[1, "minimump_p_val"], 7.893286e-08, tolerance = 1e-13)
    expect_equal(nrow(markers), 219)
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
  markers.missing <- suppressWarnings(FindConservedMarkers(object = pbmc.test, ident.1 = 0, grouping.var = "groups", test.use = "t", verbose = FALSE, base = exp(1), pseudocount.use = 1))

  test_that("FindConservedMarkers handles missing idents in certain groups", {
    expect_warning(FindConservedMarkers(object = pbmc.test, ident.1 = 0, grouping.var = "groups", test.use = "t"))
    expect_equal(colnames(x = markers.missing), paste0("g2_", standard.names))
    expect_equal(markers.missing[1, "g2_p_val"], 1.672911e-13, tolerance = 1e-18)
    expect_equal(markers.missing[1, "g2_avg_logFC"], -4.796379, tolerance = 1e-6)
    # expect_equal(markers.missing[1, "g2_pct.1"], 0.062)
    expect_equal(markers.missing[1, "g2_pct.2"], 0.95)
    expect_equal(markers.missing[1, "g2_p_val_adj"], 3.847695e-11, tolerance = 1e-16)
    expect_equal(nrow(markers.missing), 226)
    expect_equal(rownames(markers.missing)[1], "HLA-DPB1")
  })
}
