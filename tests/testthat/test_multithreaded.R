source(test_path("../testdata/test-objects.R"), local = TRUE)

context("Seurat.nthreads > 1")

expect_abs_equal <- function(x, y, tolerance = 1e-6) {
  expect_equal(unname(abs(x)), unname(abs(y)), tolerance = tolerance)
}

workflow_fixture <- local({
  cache <- new.env(parent = emptyenv())

  get_fixture <- function(name) {
    if (!exists(x = name, envir = cache, inherits = FALSE)) {
      # Build each workflow fixture lazily so related tests can share setup cost
      value <- switch(
        EXPR = name,
        counts = create_pbmc_counts(),
        raw = {
          counts <- get_fixture("counts")
          create_seurat_obj(counts = counts)
        },
        normalized = NormalizeData(object = get_fixture("raw"), verbose = FALSE),
        variable = FindVariableFeatures(
          object = get_fixture("normalized"),
          selection.method = "vst",
          verbose = FALSE
        ),
        scaled = ScaleData(
          object = get_fixture("variable"),
          features = VariableFeatures(object = get_fixture("variable")),
          verbose = FALSE
        ),
        pca = suppressWarnings(
          RunPCA(
            object = get_fixture("scaled"),
            features = VariableFeatures(object = get_fixture("scaled")),
            npcs = 10L,
            verbose = FALSE
          )
        ),
        integration = {
          object <- create_multilayer_obj()
          create_integration_obj(object = object, npcs = 20L)
        },
        integration_sct = {
          object <- create_multilayer_obj()
          create_integration_obj(object = object, sct = TRUE, npcs = 20L)
        },
        stop("Unknown workflow fixture: ", name, call. = FALSE)
      )
      assign(x = name, value = value, envir = cache)
    }
    get(x = name, envir = cache, inherits = FALSE)
  }

  get_fixture
})

# Preprocessing workflow

test_that("NormalizeData is stable across thread counts", {
  # Check that log-normalized data are identical for serial and threaded runs
  single.thread <- with_threads(
    nthreads = 1L,
    NormalizeData(object = workflow_fixture("raw"), verbose = FALSE)
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    NormalizeData(object = workflow_fixture("raw"), verbose = FALSE)
  )

  expect_equal(
    LayerData(object = single.thread, layer = "data"),
    LayerData(object = multi.thread, layer = "data")
  )
})

test_that("FindVariableFeatures is stable across thread counts", {
  # Check both selected features and the per-feature VST statistics
  single.thread <- with_threads(
    nthreads = 1L,
    FindVariableFeatures(
      object = workflow_fixture("normalized"),
      selection.method = "vst",
      verbose = FALSE
    )
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    FindVariableFeatures(
      object = workflow_fixture("normalized"),
      selection.method = "vst",
      verbose = FALSE
    )
  )

  expect_identical(
    VariableFeatures(object = single.thread),
    VariableFeatures(object = multi.thread)
  )
  expect_equal(
    HVFInfo(object = single.thread[["RNA"]], method = "vst", status = TRUE),
    HVFInfo(object = multi.thread[["RNA"]], method = "vst", status = TRUE),
    tolerance = 1e-8
  )
})

test_that("ScaleData is stable across thread counts", {
  # Check scaled values using the same variable feature set in each run
  features <- VariableFeatures(object = workflow_fixture("variable"))

  single.thread <- with_threads(
    nthreads = 1L,
    ScaleData(object = workflow_fixture("variable"), features = features, verbose = FALSE)
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    ScaleData(object = workflow_fixture("variable"), features = features, verbose = FALSE)
  )

  expect_equal(
    LayerData(object = single.thread, layer = "scale.data"),
    LayerData(object = multi.thread, layer = "scale.data"),
    tolerance = 1e-8
  )
})

test_that("SCTransform is stable across thread counts", {
  # Use v1 so results do not depend on whether optional glmGamPoi is installed
  single.thread <- with_threads(
    nthreads = 1L,
    suppressWarnings(
      SCTransform(
        object = workflow_fixture("raw"),
        vst.flavor = "v1",
        seed.use = 42L,
        verbose = FALSE
      )
    )
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    suppressWarnings(
      SCTransform(
        object = workflow_fixture("raw"),
        vst.flavor = "v1",
        seed.use = 42L,
        verbose = FALSE
      )
    )
  )

  expect_identical(
    VariableFeatures(object = single.thread),
    VariableFeatures(object = multi.thread)
  )
  expect_equal(
    LayerData(object = single.thread, assay = "SCT", layer = "counts"),
    LayerData(object = multi.thread, assay = "SCT", layer = "counts")
  )
  expect_equal(
    LayerData(object = single.thread, assay = "SCT", layer = "data"),
    LayerData(object = multi.thread, assay = "SCT", layer = "data"),
    tolerance = 1e-8
  )
  expect_equal(
    LayerData(object = single.thread, assay = "SCT", layer = "scale.data"),
    LayerData(object = multi.thread, assay = "SCT", layer = "scale.data"),
    tolerance = 1e-8
  )
})

# Dimensional reduction and graph workflow

test_that("RunPCA is stable across thread counts", {
  # PCA signs can flip, so compare absolute embeddings and loadings
  features <- VariableFeatures(object = workflow_fixture("scaled"))

  single.thread <- with_threads(
    nthreads = 1L,
    suppressWarnings(RunPCA(object = workflow_fixture("scaled"), features = features, npcs = 10L, verbose = FALSE))
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    suppressWarnings(RunPCA(object = workflow_fixture("scaled"), features = features, npcs = 10L, verbose = FALSE))
  )

  expect_equal(Stdev(object = single.thread[["pca"]]), Stdev(object = multi.thread[["pca"]]), tolerance = 1e-8)
  expect_abs_equal(
    Embeddings(object = single.thread[["pca"]]),
    Embeddings(object = multi.thread[["pca"]]),
    tolerance = 1e-8
  )
  expect_abs_equal(
    Loadings(object = single.thread[["pca"]]),
    Loadings(object = multi.thread[["pca"]]),
    tolerance = 1e-8
  )
})

test_that("FindNeighbors is stable across thread counts", {
  # Check the SNN graph, which is the graph consumed by clustering
  single.thread <- with_threads(
    nthreads = 1L,
    FindNeighbors(
      object = workflow_fixture("pca"),
      graph.name = c("RNA_nn_t1", "RNA_snn_t1"),
      k.param = 10L,
      dims = 1:10,
      verbose = FALSE
    )
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    FindNeighbors(
      object = workflow_fixture("pca"),
      graph.name = c("RNA_nn_t2", "RNA_snn_t2"),
      k.param = 10L,
      dims = 1:10,
      verbose = FALSE
    )
  )

  expect_equal(single.thread[["RNA_snn_t1"]], multi.thread[["RNA_snn_t2"]])
})

test_that("RunUMAP is stable across thread counts", {
  # Keep n.epochs at zero to exercise UMAP setup deterministically and cheaply
  single.thread <- with_threads(
    nthreads = 1L,
    RunUMAP(
      object = workflow_fixture("pca"),
      dims = 1:10,
      reduction.name = "umap_t1",
      n.neighbors = 10L,
      n.epochs = 0L,
      seed.use = 42L,
      verbose = FALSE
    )
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    RunUMAP(
      object = workflow_fixture("pca"),
      dims = 1:10,
      reduction.name = "umap_t2",
      n.neighbors = 10L,
      n.epochs = 0L,
      seed.use = 42L,
      verbose = FALSE
    )
  )

  expect_equal(
    unname(Embeddings(object = single.thread[["umap_t1"]])),
    unname(Embeddings(object = multi.thread[["umap_t2"]])),
    tolerance = 1e-8
  )
})

test_that("RunUMAP model projection is stable across thread counts", {
  # Check that model-return and projection paths remain thread stable
  projection_model <- local({
    object <- workflow_fixture("pca")
    data.use <- Embeddings(object = object, reduction = "pca")[, 1:10]
    set.seed(seed = 42L)
    model <- uwot::umap(
      X = data.use,
      n_threads = 1L,
      n_neighbors = 10L,
      n_components = 2L,
      metric = "cosine",
      n_epochs = 0L,
      learning_rate = 1,
      min_dist = 0.3,
      spread = 1,
      set_op_mix_ratio = 1,
      local_connectivity = 1L,
      repulsion_strength = 1,
      negative_sample_rate = 5L,
      init = "spectral",
      fast_sgd = FALSE,
      approx_pow = FALSE,
      verbose = FALSE,
      ret_model = TRUE
    )
    embeddings <- model$embedding
    colnames(x = embeddings) <- paste0("umapmodel_", 1:ncol(x = embeddings))
    rownames(x = embeddings) <- rownames(x = data.use)
    reduction <- CreateDimReducObject(
      embeddings = embeddings,
      key = "umapmodel_",
      assay = DefaultAssay(object = object[["pca"]]),
      global = TRUE
    )
    Misc(reduction, slot = "model") <- model
    reduction
  })
  expect_true(length(x = Misc(object = projection_model, slot = "model")) > 0L)

  single.thread <- with_threads(
    nthreads = 1L,
    RunUMAP(
      object = workflow_fixture("pca"),
      dims = 1:10,
      reduction.model = projection_model,
      reduction.name = "umap_project_t1",
      n.epochs = 0L,
      verbose = FALSE
    )
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    RunUMAP(
      object = workflow_fixture("pca"),
      dims = 1:10,
      reduction.model = projection_model,
      reduction.name = "umap_project_t2",
      n.epochs = 0L,
      verbose = FALSE
    )
  )

  expect_equal(
    unname(Embeddings(object = single.thread[["umap_project_t1"]])),
    unname(Embeddings(object = multi.thread[["umap_project_t2"]])),
    tolerance = 1e-8
  )
})

test_that("RunUMAP rejects underspecified inputs", {
  # Check mutually exclusive input modes and too few input dimensions
  expect_error(
    RunUMAP(
      object = workflow_fixture("pca"),
      dims = 1,
      n.components = 2L,
      verbose = FALSE
    ),
    "Please provide as many or more dims than n.components"
  )
  expect_error(
    RunUMAP(
      object = workflow_fixture("pca"),
      dims = 1:2,
      features = rownames(x = workflow_fixture("pca"))[1:2],
      verbose = FALSE
    ),
    "Please specify only one"
  )
})

# C++ integration helpers

test_that("CountAnchorSharedNeighbors is stable across thread counts", {
  indices.aa <- matrix(c(1L, 2L), nrow = 2, byrow = TRUE)
  indices.ab <- matrix(c(1L, 2L), nrow = 2, byrow = TRUE)
  indices.ba <- matrix(c(1L, 2L), nrow = 2, byrow = TRUE)
  indices.bb <- matrix(c(2L, 1L), nrow = 2, byrow = TRUE)

  expect_equal(
    CountAnchorSharedNeighbors(
      indices_aa = indices.aa,
      indices_ab = indices.ab,
      indices_ba = indices.ba,
      indices_bb = indices.bb,
      anchor_cell1 = c(1L, 2L),
      anchor_cell2 = c(1L, 2L),
      offset = 2L,
      k_score = 1L,
      nthreads = 1L
    ),
    CountAnchorSharedNeighbors(
      indices_aa = indices.aa,
      indices_ab = indices.ab,
      indices_ba = indices.ba,
      indices_bb = indices.bb,
      anchor_cell1 = c(1L, 2L),
      anchor_cell2 = c(1L, 2L),
      offset = 2L,
      k_score = 1L,
      nthreads = 2L
    )
  )
})

test_that("FindWeightsC is stable across thread counts", {
  distances <- matrix(c(0.2,
                        0.4), nrow = 2, byrow = TRUE)
  cell.index <- matrix(c(1,
                         2), nrow = 2, byrow = TRUE)

  expect_equal(
    as.matrix(FindWeightsC(
      cells2 = 0:1,
      distances = distances,
      anchor_cells2 = c("cell1", "cell2"),
      integration_matrix_rownames = c("cell1", "cell2"),
      cell_index = cell.index,
      anchor_score = c(0.5, 0.8),
      min_dist = 0.01,
      sd = 1,
      display_progress = FALSE,
      nthreads = 1L
    )),
    as.matrix(FindWeightsC(
      cells2 = 0:1,
      distances = distances,
      anchor_cells2 = c("cell1", "cell2"),
      integration_matrix_rownames = c("cell1", "cell2"),
      cell_index = cell.index,
      anchor_score = c(0.5, 0.8),
      min_dist = 0.01,
      sd = 1,
      display_progress = FALSE,
      nthreads = 2L
    )),
    tolerance = 1e-12
  )
})

test_that("ScoreHelper is stable across thread counts", {
  snn <- as.sparse(matrix(c(0, 1,
                            1, 0), nrow = 2, byrow = TRUE))
  query.pca <- matrix(c(0, 1,
                        0, 1), nrow = 2, byrow = TRUE)
  query.dists <- matrix(c(0.1, 0.2,
                          0.1, 0.3), nrow = 2, byrow = TRUE)
  corrected.nns <- matrix(c(2L, 1L), nrow = 2, byrow = TRUE)

  expect_equal(
    ScoreHelper(
      snn = snn,
      query_pca = query.pca,
      query_dists = query.dists,
      corrected_nns = corrected.nns,
      k_snn = 1L,
      subtract_first_nn = FALSE,
      display_progress = FALSE,
      nthreads = 1L
    ),
    ScoreHelper(
      snn = snn,
      query_pca = query.pca,
      query_dists = query.dists,
      corrected_nns = corrected.nns,
      k_snn = 1L,
      subtract_first_nn = FALSE,
      display_progress = FALSE,
      nthreads = 2L
    ),
    tolerance = 1e-12
  )
})

# Differential expression workflow

test_that("FindAllMarkers is stable across thread counts", {
  # Compare the full per-cluster result table
  single.thread <- with_threads(
    nthreads = 1L,
    suppressMessages(
      suppressWarnings(
        FindAllMarkers(
          object = pbmc_small,
          logfc.threshold = 0,
          min.pct = 0,
          only.pos = FALSE,
          return.thresh = Inf,
          verbose = FALSE,
          pseudocount.use = 1
        )
      )
    )
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    suppressMessages(
      suppressWarnings(
        FindAllMarkers(
          object = pbmc_small,
          logfc.threshold = 0,
          min.pct = 0,
          only.pos = FALSE,
          return.thresh = Inf,
          verbose = FALSE,
          pseudocount.use = 1
        )
      )
    )
  )

  rownames(x = single.thread) <- NULL
  rownames(x = multi.thread) <- NULL
  expect_equal(single.thread, multi.thread, tolerance = 1e-8)
})

# Integration workflow

test_that("IntegrateLayers CCA and RPCA are stable across thread counts", {
  # Compare corrected reductions for both standard integration methods
  for (method in list(CCAIntegration, RPCAIntegration)) {
    single.thread <- with_threads(
      nthreads = 1L,
      local({
        set.seed(seed = 42)
        suppressWarnings(
          IntegrateLayers(
            object = workflow_fixture("integration"),
            method = method,
            orig.reduction = "pca",
            new.reduction = "integrated_t1",
            k.weight = 10L,
            dims.to.integrate = 1:10,
            verbose = FALSE
          )
        )
      })
    )
    multi.thread <- with_threads(
      nthreads = 2L,
      local({
        set.seed(seed = 42)
        suppressWarnings(
          IntegrateLayers(
            object = workflow_fixture("integration"),
            method = method,
            orig.reduction = "pca",
            new.reduction = "integrated_t2",
            k.weight = 10L,
            dims.to.integrate = 1:10,
            verbose = FALSE
          )
        )
      })
    )

    expect_abs_equal(
      Embeddings(object = single.thread[["integrated_t1"]]),
      Embeddings(object = multi.thread[["integrated_t2"]]),
      tolerance = 1e-6
    )
  }
})

test_that("IntegrateLayers RPCA with SCT normalization is stable across thread counts", {
  # Cover the SCT integration path with one representative method
  single.thread <- with_threads(
    nthreads = 1L,
    local({
      set.seed(seed = 42)
      suppressWarnings(
        IntegrateLayers(
          object = workflow_fixture("integration_sct"),
          method = RPCAIntegration,
          orig.reduction = "pca",
          new.reduction = "integrated_sct_t1",
          k.weight = 10L,
          dims.to.integrate = 1:10,
          normalization.method = "SCT",
          verbose = FALSE
        )
      )
    })
  )
  multi.thread <- with_threads(
    nthreads = 2L,
    local({
      set.seed(seed = 42)
      suppressWarnings(
        IntegrateLayers(
          object = workflow_fixture("integration_sct"),
          method = RPCAIntegration,
          orig.reduction = "pca",
          new.reduction = "integrated_sct_t2",
          k.weight = 10L,
          dims.to.integrate = 1:10,
          normalization.method = "SCT",
          verbose = FALSE
        )
      )
    })
  )

  expect_abs_equal(
    Embeddings(object = single.thread[["integrated_sct_t1"]]),
    Embeddings(object = multi.thread[["integrated_sct_t2"]]),
    tolerance = 1e-6
  )
})
