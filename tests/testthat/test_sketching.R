# setup shared test fixtures
path_to_counts <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")


build_test_data <- function(multi_layer = FALSE) {
  counts <- read.table(path_to_counts, sep = "\t", row.names = 1)
  counts <- as.sparse(as.matrix(counts))

  if (multi_layer) {
    barcodes <- colnames(counts)

    counts.layer1 <- counts
    colnames(counts.layer1) <- paste0(barcodes, "_layer1")

    counts.layer2 <- counts
    colnames(counts.layer1) <- paste0(barcodes, "_layer2")

    counts_layers <- c(
      counts.layer1 = counts.layer1,
      counts.layer2 = counts.layer2
    )

    test_data <- CreateSeuratObject(counts_layers)
    
  } else{
    test_data <- CreateSeuratObject(counts)
  }

  test_data <- NormalizeData(test_data, verbose = FALSE)
  test_data <- FindVariableFeatures(test_data, verbose = FALSE)

  return (test_data)
}


context("SketchData")

test_that("SketchData defaults work", {
  test_case <- build_test_data()

  result <- suppressWarnings(
    SketchData(
      test_case, 
      assay = "RNA", 
      ncells = 50, 
      method = "LeverageScore", 
      sketched.assay = "sketch", 
      set.seed = 42
    )
  )
  expect_equal(
    dim(result[["sketch"]]), 
    c(230, 50)
  )
  expect_equal(
    as.numeric(result$leverage.score[1]), 
    0.9036446, 
    tolerance = 1e-5
  )
  expect_equal(
    colnames(result[["sketch"]])[1], 
    "ATGCCAGAACGACT"
  )

  result_2 <- suppressWarnings(
    SketchData(
      test_case, 
      assay = "RNA", 
      ncells = c("data" = 50), 
      method = "LeverageScore", 
      sketched.assay = "sketch"
    )
  )
  expect_equal(
    dim(result_2[["sketch"]]), 
    dim(result[["sketch"]])
  )
  expect_equal(
    as.numeric(result_2$leverage.score[1]), 
    as.numeric(result$leverage.score[1]), 
    tolerance = 1e-5
  )
  expect_equal(
    colnames(result_2[["sketch"]])[1], 
    colnames(result[["sketch"]])[1]
  )
})

test_that("SketchData with multiple layers works", { # (and one is less than the number of cells in that layer)
  test_case <- build_test_data(multi_layer = TRUE)
  
  result <- suppressWarnings(
    SketchData(
      test_case, 
      assay = "RNA", 
      ncells = 80, 
      method = "LeverageScore", 
      sketched.assay = "sketch", 
      set.seed = 42
    )
  )
  expect_equal(
    dim(result[["sketch"]]$data.layer1), 
    c(230, 80)
  )
  expect_equal(
    dim(result[["sketch"]]$data.layer2), 
    c(230, 80)
  )
  expect_equal(
    as.numeric(result$leverage.score[1]), 
    0.9036446, 
    tolerance = 1e-6
  )
  expect_equal(
    colnames(result[["sketch"]])[1], 
    "ATGCCAGAACGACT_layer2"
  )
})

test_that("SketchData with a different number of cells per layer works", {
  test_case <- build_test_data(multi_layer = TRUE)
  
  result <- suppressWarnings(
    SketchData(
      test_case, 
      assay = "RNA", 
      ncells = c(50, 30), 
      method = "LeverageScore", 
      sketched.assay = "sketch", 
      set.seed = 42
    )
  )
  expect_equal(
    dim(result[["sketch"]]$data.layer1), 
    c(230, 50)
  )
  expect_equal(
    dim(result[["sketch"]]$data.layer2), 
    c(230, 30)
  )
  expect_equal(
    as.numeric(result$leverage.score[1]), 
    0.9036446, 
    tolerance = 1e-5
  )
  expect_equal(
    colnames(result[["sketch"]])[1],  
    "ATGCCAGAACGACT_layer2"
  )

  result_2 <- suppressWarnings(
    SketchData(
      test_case, 
      assay = "RNA", 
      ncells = c("data.layer1" = 50, "data.layer2" = 30), 
      method = "LeverageScore", 
      sketched.assay = "sketch", 
      set.seed = 42
    )
  )
  expect_equal(
    dim(result_2[["sketch"]]$data.layer1), 
    dim(result[["sketch"]]$data.layer1)
  )
  expect_equal(
    dim(result_2[["sketch"]]$data.layer2), 
    dim(result[["sketch"]]$data.layer2)
  )
  expect_equal(
    as.numeric(result_2$leverage.score[1]), 
    as.numeric(result$leverage.score[1]), 
    tolerance = 1e-5
  )
  expect_equal(
    colnames(result_2[["sketch"]])[1], 
    colnames(result[["sketch"]])[1]
  )
})

test_that("SketchData with specified features works", {
  test_case <- build_test_data()
  
  result <- suppressWarnings(
    SketchData(
      test_case, 
      assay = "RNA", 
      ncells = 50, 
      method = "LeverageScore", 
      sketched.assay = "sketch", 
      features = VariableFeatures(test_case)[1:100]
    )
  )
  expect_equal(
    dim(result[["sketch"]]), 
    c(230, 50)
  )
  expect_equal(
    as.numeric(result$leverage.score[1]), 
    0.7202897, 
    tolerance = 1e-6
  )
  expect_equal(
    colnames(result[["sketch"]])[1], 
    "ATGCCAGAACGACT"
  )
})

test_that("SketchData with specified features and multiple layers works", {
  test_case <- build_test_data(multi_layer = TRUE)

  result <- suppressWarnings(
    SketchData(
      test_case, 
      assay = "RNA", 
      ncells = c(50, 30), 
      method = "LeverageScore", 
      sketched.assay = "sketch",
      features = VariableFeatures(test_case)[1:100], 
      set.seed = 42
    )
  )
  expect_equal(
    dim(result[["sketch"]]$data.layer1), 
    c(230, 50)
  )
  expect_equal(
    dim(result[["sketch"]]$data.layer2), 
    c(230, 30)
  )
  expect_equal(
    as.numeric(result$leverage.score[1]), 
    0.7202896, 
    tolerance = 1e-6
  )
  expect_equal(
    colnames(result[["sketch"]])[1], 
    "ATGCCAGAACGACT_layer2"
  )
})

test_that("SketchData when setting your own variable features and specifying features works", {
  test_case <- build_test_data()
  top_features <- VariableFeatures(test_case)[1:100]
  VariableFeatures(test_case) <- top_features
  
  result <- suppressWarnings(
    SketchData(
      test_case, 
      assay = "RNA", 
      ncells = 50, 
      method = "LeverageScore", 
      sketched.assay = "sketch", 
      features = top_features, 
      set.seed = 42
    )
  )
  expect_equal(
    dim(result[["sketch"]]), 
    c(230,50)
  )
  expect_equal(
    as.numeric(result$leverage.score[1]), 
    0.7202896, 
    tolerance = 1e-6
  )
  expect_equal(
    colnames(result[["sketch"]])[1], 
    "ATGCCAGAACGACT"
  )
})

test_that("SketchData when setting your own variable features and not specifying features errors out", {
  test_case <- build_test_data()
  top_features <- VariableFeatures(test_case)[1:100]
  VariableFeatures(test_case) <- top_features
  
  expect_error(
    suppressWarnings(
      SketchData(
        pbmc_new, 
        assay = "RNA", 
        ncells = 50, 
        method = "LeverageScore", 
        sketched.assay = "sketch"
      )
    )
  )
})
