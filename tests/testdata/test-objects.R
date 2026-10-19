# Shared Seurat object builders for tests

create_pbmc_counts <- function() {
  pbmc.file <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")
  as.sparse(x = as.matrix(read.table(pbmc.file, sep = "\t", row.names = 1)))
}

# Builds a matrix of random counts with the specified number of features and cells
# min/max.count: range of counts to sample from
create_random_counts <- function(
  nfeatures = 100L,
  ncells = 100L,
  seed = 123L,
  min.count = 1L,
  max.count = 50L
) {
  old.seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(x = old.seed)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old.seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed = seed)
  counts <- matrix(
    data = sample(
      x = min.count:max.count,
      size = nfeatures * ncells,
      replace = TRUE
    ),
    nrow = nfeatures,
    ncol = ncells
  )
  rownames(x = counts) <- paste0("gene", seq_len(nfeatures))
  colnames(x = counts) <- paste0("cell", seq_len(ncells))
  as.sparse(x = counts)
}

# Create a Seurat object with the specified counts, assay version, and assay name
# If counts are not provided, default to using PBMC counts from testdata
create_seurat_obj <- function(
  counts = NULL,
  assay_version = getOption("Seurat.object.assay.version", default = "v5"),
  assay = "RNA",
  meta.data = NULL
) {
  if (is.null(x = counts)) {
    counts <- create_pbmc_counts()
  }
  if (identical(x = assay_version, y = "v3")) {
    assay.object <- CreateAssayObject(counts = counts)
    object <- CreateSeuratObject(counts = assay.object, assay = assay, meta.data = meta.data)
  } else if (identical(x = assay_version, y = "v5")) {
    assay.object <- CreateAssay5Object(counts = counts)
    object <- CreateSeuratObject(counts = assay.object, assay = assay, meta.data = meta.data)
  } else {
    stop("`assay_version` should be one of 'v3' or 'v5'", call. = FALSE)
  }
  object
}

# Create a Seurat object and run standard preprocessing steps (normalization, variable feature selection, scaling, PCA, and neighbor finding)
# Parameters are used to control which steps are run and the parameters for each step
create_prep_obj <- function(
  object = NULL,
  normalize = TRUE,
  variable.features = FALSE,
  scale = TRUE,
  pca = FALSE,
  neighbors = FALSE,
  nfeatures = 20L,
  npcs = 10L,
  dims = seq_len(npcs),
  k.param = 10L
) {
  if (is.null(x = object)) {
    counts <- create_pbmc_counts()
    object <- create_seurat_obj(counts = counts)
  }
  if (isTRUE(x = normalize)) {
    object <- NormalizeData(object = object, verbose = FALSE)
  }
  if (isTRUE(x = variable.features)) {
    object <- FindVariableFeatures(
      object = object,
      nfeatures = nfeatures,
      verbose = FALSE
    )
  }
  if (isTRUE(x = scale)) {
    object <- ScaleData(object = object, verbose = FALSE)
  }
  if (isTRUE(x = pca)) {
    features <- VariableFeatures(object = object)
    features <- if (length(x = features) > 0L) features else rownames(x = object)
    object <- suppressWarnings(
      RunPCA(object = object, features = features, npcs = npcs, verbose = FALSE)
    )
  }
  if (isTRUE(x = neighbors)) {
    object <- FindNeighbors(
      object = object,
      k.param = k.param,
      dims = dims,
      verbose = FALSE
    )
  }
  object
}

# Create a Seurat object with multiple layers and split by a specified metadata column
create_multilayer_obj <- function(
  object = pbmc_small,
  split.by = "groups"
) {
  object <- object
  suppressWarnings(
    object[["RNA"]] <- CreateAssay5Object(
      counts = LayerData(object = object, assay = "RNA", layer = "counts")
    )
  )
  object[["RNA"]] <- split(x = object[["RNA"]], f = object[[split.by, drop = TRUE]])
  object
}

# Create a Seurat object and run the standard preprocessing steps (normalization, variable feature selection, scaling, PCA) for integration testing
create_integration_obj <- function(
  object = NULL,
  sct = FALSE,
  npcs = 50L,
  seed.use = 12345L
) {
  if (is.null(x = object)) {
    object <- create_multilayer_obj()
  }
  if (isTRUE(x = sct)) {
    object <- suppressWarnings(
      SCTransform(
        object = object,
        vst.flavor = "v1",
        seed.use = seed.use,
        verbose = FALSE
      )
    )
  } else {
    object <- NormalizeData(object = object, verbose = FALSE)
    object <- FindVariableFeatures(object = object, verbose = FALSE)
    object <- ScaleData(object = object, verbose = FALSE)
  }
  suppressWarnings(RunPCA(object = object, npcs = npcs, verbose = FALSE))
}

# Run code with a specified number of threads and restore the original thread count after execution
with_threads <- function(nthreads, code) {
  old.threads <- getThreads(verbose = FALSE)
  on.exit(setThreads(old.threads, verbose = FALSE), add = TRUE)
  setThreads(nthreads, verbose = FALSE)
  force(code)
}
