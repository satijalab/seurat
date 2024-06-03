# Tests for sketching related fxns
pbmc_small <- suppressWarnings(UpdateSeuratObject(pbmc_small))

# Setup test object
pbmc_small <- NormalizeData(pbmc_small, verbose = FALSE)
pbmc_small <- FindVariableFeatures(pbmc_small)

# Tests for SketchData
# ------------------------------------------------------------------------------
context("SketchData")

test_that("SketchData defaults work", {
  pbmc_sketched <- suppressWarnings(SketchData(pbmc_small, assay = "RNA", ncells = 50, method = "LeverageScore", 
                              sketched.assay = "sketch", set.seed = 42))
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  expect_equal(as.numeric(pbmc_sketched$leverage.score[1]), 0.9036446, tolerance = 1e-6)
  expect_equal(colnames(pbmc_sketched[["sketch"]])[1], "ATGCCAGAACGACT")
  pbmc_sketched_2 <- suppressWarnings(SketchData(pbmc_small, assay = "RNA", ncells = c("data" = 50), method = "LeverageScore", sketched.assay = "sketch"))
  expect_equal(dim(pbmc_sketched_2[["sketch"]]), dim(pbmc_sketched[["sketch"]]))
  expect_equal(as.numeric(pbmc_sketched_2$leverage.score[1]), as.numeric(pbmc_sketched$leverage.score[1]), tolerance = 1e-6)
  expect_equal(colnames(pbmc_sketched_2[["sketch"]])[1], colnames(pbmc_sketched[["sketch"]])[1])
})


pbmc_split <- suppressWarnings(merge(pbmc_small, pbmc_small))
pbmc_split <- suppressWarnings(split(pbmc_split, f = pbmc_split$groups))
pbmc_split <- FindVariableFeatures(pbmc_split)

test_that("SketchData with multiple layers works", { # (and one is less than the number of cells in that layer)
  pbmc_sketched <- suppressWarnings(SketchData(pbmc_split, assay = "RNA", ncells = 80, method = "LeverageScore", 
                              sketched.assay = "sketch", set.seed = 42))
  expect_equal(dim(pbmc_sketched[["sketch"]]$data.g1), c(230,80))
  expect_equal(dim(pbmc_sketched[["sketch"]]$data.g2), c(230,72))
  expect_equal(as.numeric(pbmc_sketched$leverage.score[1]), 0.7471112, tolerance = 1e-6)
  expect_equal(colnames(pbmc_sketched[["sketch"]])[1], "ATGCCAGAACGACT_1")
})

test_that("SketchData with a different number of cells per layer works", {
  pbmc_sketched <- suppressWarnings(SketchData(pbmc_split, assay = "RNA", ncells = c(50, 30), method = "LeverageScore", 
                              sketched.assay = "sketch", set.seed = 42))
  expect_equal(dim(pbmc_sketched[["sketch"]]$data.g1), c(230,30))
  expect_equal(dim(pbmc_sketched[["sketch"]]$data.g2), c(230,50))
  expect_equal(as.numeric(pbmc_sketched$leverage.score[1]), 0.6220129, tolerance = 1e-6)
  expect_equal(colnames(pbmc_sketched[["sketch"]])[1],  "ATGCCAGAACGACT_1")
  pbmc_sketched_2 <- suppressWarnings(SketchData(pbmc_split, assay = "RNA", ncells = c("data.g2" = 50, "data.g1" = 30), 
                                method = "LeverageScore", sketched.assay = "sketch", set.seed = 42))
  expect_equal(dim(pbmc_sketched_2[["sketch"]]$data.g1), dim(pbmc_sketched[["sketch"]]$data.g1))
  expect_equal(dim(pbmc_sketched_2[["sketch"]]$data.g2), dim(pbmc_sketched[["sketch"]]$data.g2))
  expect_equal(as.numeric(pbmc_sketched_2$leverage.score[1]), as.numeric(pbmc_sketched$leverage.score[1]), 
               tolerance = 1e-6)
  expect_equal(colnames(pbmc_sketched_2[["sketch"]])[1], colnames(pbmc_sketched[["sketch"]])[1])
})

test_that("SketchData with specified features works", {
  pbmc_sketched <- suppressWarnings(SketchData(pbmc_small, assay = "RNA", ncells = 50, method = "LeverageScore", 
                              sketched.assay = "sketch", features = VariableFeatures(pbmc_small)[1:100]))
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  expect_equal(as.numeric(pbmc_sketched$leverage.score[1]), 0.7202897, tolerance = 1e-6)
  expect_equal(colnames(pbmc_sketched[["sketch"]])[1], "ATGCCAGAACGACT")
})

test_that("SketchData with specified features and multiple layers works", {
  pbmc_sketched <- suppressWarnings(SketchData(pbmc_split, assay = "RNA", ncells = c(50, 30), 
                              method = "LeverageScore", sketched.assay = "sketch",
                              features = VariableFeatures(pbmc_small)[1:100], set.seed = 42))
  expect_equal(dim(pbmc_sketched[["sketch"]]$data.g1), c(230,30))
  expect_equal(dim(pbmc_sketched[["sketch"]]$data.g2), c(230,50))
  expect_equal(as.numeric(pbmc_sketched$leverage.score[1]), 0.7324468, tolerance = 1e-6)
  expect_equal(colnames(pbmc_sketched[["sketch"]])[1], "ATGCCAGAACGACT_1")
})

pbmc_new <- CreateSeuratObject(pbmc_small[["RNA"]]$counts)
pbmc_new <- NormalizeData(pbmc_new, verbose = FALSE)
VariableFeatures(pbmc_new) <- rownames(pbmc_new)[1:100]
test_that("SketchData when setting your own variable features and specifying features works", {
  pbmc_sketched <- suppressWarnings(SketchData(pbmc_new, assay = "RNA", ncells = 50, 
                              method = "LeverageScore", sketched.assay = "sketch", 
                              features = VariableFeatures(pbmc_small)[1:100], set.seed = 42))
  expect_equal(dim(pbmc_sketched[["sketch"]]), c(230,50))
  expect_equal(as.numeric(pbmc_sketched$leverage.score[1]), 0.7202896, tolerance = 1e-6)
  expect_equal(colnames(pbmc_sketched[["sketch"]])[1], "ATGCCAGAACGACT")
})


test_that("SketchData when setting your own variable features and not specifying features errors out", {
  expect_error(suppressWarnings(SketchData(pbmc_new, assay = "RNA", ncells = 50, method = "LeverageScore", 
                          sketched.assay = "sketch")))
})

