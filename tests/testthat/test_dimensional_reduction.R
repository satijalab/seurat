context("test-dimensional_reduction")

test_that("different ways of passing distance matrix", {
  # Generate dummy data exp matrix
  set.seed(1)
  dummyexpMat <- matrix(data = sample(x = c(1:50), size = 1e4, replace = TRUE), 
                        ncol = 100, nrow = 100)
  colnames(dummyexpMat) <- paste0("cell", seq(ncol(dummyexpMat)))
  row.names(dummyexpMat) <- paste0("gene", seq(nrow(dummyexpMat)))
  
  # Create Seurat object for testing
  obj <- CreateSeuratObject(counts = dummyexpMat)
  
  # Manually make a distance object to test
  distMat <- dist(t(dummyexpMat))
  
  expect_equivalent(RunTSNE(obj, distance.matrix = distMat),
                    RunTSNE(obj, distance.matrix = as.matrix(distMat)))
  expect_equivalent(RunTSNE(obj, distance.matrix = distMat)@reductions$tsne,
                    RunTSNE(distMat, assay = "RNA"))
  expect_equivalent(RunTSNE(obj, distance.matrix = distMat)@reductions$tsne,
                    RunTSNE(as.matrix(distMat), assay = "RNA", is_distance = TRUE))
})
