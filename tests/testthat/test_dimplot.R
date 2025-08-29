
test_that("DimPlot does not use meta.data columns for embeddings", {
  pbmc <- pbmc_small
  pbmc$PC_1 <- rep(1, ncol(pbmc)) # create bad metadata
  
  p <- DimPlot(pbmc, reduction = "pca")
  
  # show data used for ggplot
  df <- layer_data(p)
  
  expect_true(var(df$x) > 0) # it shouldn't be constant (using real PCA values)
})