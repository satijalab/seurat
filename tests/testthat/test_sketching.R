# setup shared test fixtures
path_to_counts <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")


build_test_data <- function(multi_layer = FALSE, preprocess = TRUE) {
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

  if (preprocess) {
    test_data <- NormalizeData(test_data, verbose = FALSE)
    test_data <- FindVariableFeatures(test_data, verbose = FALSE)
  }

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

test_that("ProjectData preserves factor labels", {
  set.seed(42)
  sobj <- build_test_data()
  sobj <- suppressWarnings(SketchData(sobj, ncells = 50, verbose = FALSE))
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 10, verbose = FALSE)
  sobj$type <- factor(sample(3, ncol(sobj), replace = TRUE))
  sobj <- ProjectData(sobj, sketched.reduction = "pca", dims = 1:10, full.reduction = "pca.full", 
                      refdata = list(proj_type = "type"), verbose = FALSE)
  expect_s3_class(sobj$proj_type, "factor")
  expect_identical(levels(sobj$proj_type), levels(sobj$type))
})

context("SketchProduct")

test_that("SketchProduct keeps sketch rows that draw no cells", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  library(BPCells)

  set.seed(42)
  counts <- as.sparse(matrix(rpois(20 * 200, lambda = 2), nrow = 20))
  rownames(counts) <- paste0("feature", seq_len(20))
  colnames(counts) <- paste0("cell", seq_len(200))
  on.disk <- suppressWarnings(
    write_matrix_dir(
      mat = convert_matrix_type(as(counts, "IterableMatrix"), "double"),
      dir = tempfile()
    )
  )

  # a CountSketch-shaped matrix whose last row is assigned no cells
  nsketch <- 10L
  S <- sparseMatrix(
    i = rep(seq_len(nsketch - 1L), length.out = ncol(counts)),
    j = seq_len(ncol(counts)),
    x = rep(c(-1, 1), length.out = ncol(counts)),
    dims = c(nsketch, ncol(counts))
  )

  result <- SketchProduct(object = on.disk, S = S)
  expect_equal(
    unname(as.matrix(result)),
    unname(as.matrix(S %*% t(counts)))
  )
  expect_true(all(result[nsketch, ] == 0))
})

test_that("LeverageScore and SketchData handle a small matrix", {
  set.seed(42)
  counts <- matrix(rpois(5 * 12, lambda = 5), nrow = 5,
    dimnames = list(paste0("gene", 1:5), paste0("cell", 1:12)))
  test_case <- NormalizeData(CreateSeuratObject(as.sparse(counts)), verbose = FALSE)
  data <- LayerData(test_case, layer = "data")

  # irlba uses exact SVD for inputs with fewer than six rows or columns
  # advisory about requesting most components is expected on these inputs
  for (small.data in list(data, t(data))) {
    nv <- min(dim(small.data)) - 1L
    expected <- rowSums(base::svd(as.matrix(small.data), nu = 0, nv = nv)$v ^ 2)
    expect_warning(scores <- LeverageScore(small.data, verbose = FALSE), "percentage")
    expect_equal(scores, expected, tolerance = 1e-12)
    expect_true(all(is.finite(scores) & scores >= 0))
  }

  for (ncells in c(4L, ncol(test_case) * 3L)) {
    expect_warning(sketched <- SketchData(
      test_case,
      ncells = ncells,
      method = "LeverageScore",
      features = rownames(test_case),
      sketched.assay = "sketch",
      verbose = FALSE
    ), "percentage")
    expect_true("sketch" %in% Assays(sketched))
    expect_equal(ncol(sketched[["sketch"]]), min(ncells, ncol(test_case)))
    expect_true(all(colnames(sketched[["sketch"]]) %in% colnames(test_case)))
    expect_true(all(is.finite(sketched$leverage.score) & sketched$leverage.score >= 0))
  }
})


context("SketchDataByGroup")

test_that("SketchDataByGroup sketches the whole dataset without group.by", {
  test_case <- build_test_data()
  for (method in c("Uniform", "LeverageScore")) {
    result <- suppressWarnings(SketchDataByGroup(
      object = test_case,
      assay = "RNA",
      ncells = 50,
      sketched.assay = "sketch",
      method = method,
      verbose = FALSE
    ))
    expect_equal(ncol(x = result[["sketch"]]), 50)
    expect_equal(DefaultAssay(object = result), "sketch")
  }
})

test_that("SketchDataByGroup without group.by samples each layer", {
  test_case <- build_test_data()
  groups <- rep(x = c("a", "b"), length.out = ncol(x = test_case))
  test_case[["RNA"]] <- suppressMessages(split(x = test_case[["RNA"]], f = groups))
  result <- suppressWarnings(SketchDataByGroup(
    object = test_case,
    assay = "RNA",
    ncells = 10,
    sketched.assay = "sketch",
    method = "Uniform",
    verbose = FALSE
  ))
  sketched <- Cells(x = result[["sketch"]])
  expect_equal(as.vector(x = table(groups[match(sketched, colnames(x = test_case))])), c(10, 10))
})

test_that("SketchDataByGroup without group.by keeps on-disk data on disk", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  loadNamespace("BPCells")
  test_case <- build_test_data()
  LayerData(object = test_case[["RNA"]], layer = "data") <- as(
    object = LayerData(object = test_case[["RNA"]], layer = "data"),
    Class = "IterableMatrix"
  )
  result <- suppressWarnings(SketchDataByGroup(
    object = test_case,
    assay = "RNA",
    ncells = 50,
    sketched.assay = "sketch",
    method = "Uniform",
    verbose = FALSE
  ))
  expect_s4_class(LayerData(object = result[["sketch"]], layer = "data"), "IterableMatrix")
})

test_that("SketchDataByGroup honors named group sizes", {
  object <- build_test_data()
  object$group <- rep(c("a", "b"), length.out = ncol(object))
  result <- SketchDataByGroup(
    object, group.by = "group", ncells = c(b = 7, a = 3), verbose = FALSE
  )
  cells <- Cells(x = result[["sketch"]])
  expect_equal(as.vector(table(object[["group"]][cells, 1])), c(3, 7))
})

test_that("SketchDataByGroup honors per-group layer sizes", {
  object <- build_test_data()
  object$group <- rep(c("a", "b"), length.out = ncol(object))
  object$batch <- rep(c("one", "one", "two", "two"), length.out = ncol(object))
  object[["RNA"]] <- suppressMessages(split(x = object[["RNA"]], f = object$batch))
  result <- SketchDataByGroup(
    object, group.by = "group", verbose = FALSE,
    ncells = list(b = c(data.two = 8, data.one = 6), a = c(data.two = 5, data.one = 3))
  )
  cells <- Cells(x = result[["sketch"]])
  observed <- table(object[["group"]][cells, 1], object[["batch"]][cells, 1])
  expect_equal(as.vector(observed), c(3, 6, 5, 8))
  expect_error(SketchDataByGroup(
    object, group.by = "group", ncells = list(a = 5), verbose = FALSE
  ), "every group")
  expect_error(SketchDataByGroup(
    object, group.by = "group", ncells = list(a = c(data.one = 5), b = 5), verbose = FALSE
  ), "Please specify ncells for group 'a' in layer 'data.two'")
})

test_that("PrepareSampleObject works without group.by.Sketch", {
  set.seed(42)
  test_case <- build_test_data()
  test_case$Disease <- factor(x = rep(x = c("HV", "SEV"), length.out = ncol(x = test_case)))
  test_case <- ScaleData(object = test_case, verbose = FALSE)
  test_case <- suppressWarnings(RunPCA(object = test_case, npcs = 5, verbose = FALSE))
  result <- suppressWarnings(suppressMessages(PrepareSampleObject(
    object = test_case,
    assay = "RNA",
    Y = "Disease",
    add.hvg = FALSE,
    ncells.per.group = 60,
    ncells.landmark = 30,
    landmark.sketch.method = "Uniform",
    ncomp = 2,
    name.reduction.1 = "pca",
    dims.reduction.1 = 1:5,
    k.nn = 5,
    knn.range = 25,
    verbose = FALSE
  )))
  expect_equal(ncol(x = result[["TRAINING"]]), 60)
  expect_equal(ncol(x = result[["LANDMARK"]]), 30)
})
