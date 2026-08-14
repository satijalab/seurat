#' Run Diffusion Map
#'
#' Run diffusion map dimensionality reduction on single-cell data.
#'
#' @param object An object to run diffusion map on
#' @param assay Name of Assay diffusion map is being run on
#' @param ncomp Number of components to compute
#' @param metric This determines the choice of metric used to measure
#' distance in the input space. Check DiffusionMap() for details. If a distance.matrix is provided,
#' this will be ignored.
#' @param distance.matrix If set, runs DiffusionMap on the given distance matrix
#' instead of data matrix.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names.
#' @param reduction.name the name of the DimReduc object
#' @param weight.by.var Weight the cell embeddings by the eigenvalues of the diffusion map
#' @param verbose To print diagnostic messages. Default is TRUE.
#' @param seed.use Set a random seed.  Setting NULL will not set a seed.
#' @param features Features to compute Diffusion Map
#' @param layer The layer in \code{assay} to use when running diffusion map analysis.
#' @param ... Additional arguments to be passed to DiffusionMap()
#'
#' @return Returns a DimReduc object with diffusion map results
#' @export
#' @concept dimensional_reduction
#'
#' @examples
#' \dontrun{
#' # Run diffusion map on a Seurat object
#' seurat_obj <- RunDiffusionMap(seurat_obj, ncomp = 10)
#' }
#'
#' @importFrom SeuratObject CreateDimReducObject DefaultAssay Assays
#'
RunDiffusionMap <- function(object, ...) {
  UseMethod(generic = 'RunDiffusionMap', object = object)
}

#' @rdname RunDiffusionMap
#' @method RunDiffusionMap default
#' @export
RunDiffusionMap.default <- function(
    object,
    assay = NULL,
    metric = "cosine",
    distance.matrix = NULL,
    ncomp = 20,
    reduction.name = "DiffMap",
    reduction.key = "DC_",
    weight.by.var = TRUE,
    verbose = TRUE,
    seed.use = 101,
    ...
) {
  # Check if destiny is available
  if (!requireNamespace("destiny", quietly = TRUE)) {
    stop(
      "The 'destiny' package is required for RunDiffusionMap() but is not installed.\n\n",
      "To install destiny, please run:\n",
      "  if (!require('BiocManager', quietly = TRUE))\n",
      "    install.packages('BiocManager')\n",
      "  BiocManager::install('destiny')\n\n",
      call. = FALSE
    )
  }

  # Check for knn.covertree with a warning instead of error
  if (!requireNamespace("knn.covertree", quietly = TRUE)) {
    warning(
      "Package 'knn.covertree' is not installed. ",
      "Some advanced distance metrics may not be available.\n",
      "To install: install.packages('knn.covertree')",
      call. = FALSE
    )
  }
  #
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  ncomp <- min(ncomp, ncol(x = object))

  # run Diffusion Map
  if (is.null(distance.matrix) & !is.null(metric)){
    dm_results <- destiny::DiffusionMap(data = t(object), n_eigs = ncomp, distance = metric)
  } else if (!is.null(distance.matrix)){
    dm_results <- destiny::DiffusionMap(distance = distance.matrix, n_eigs = ncomp)
  } else {
    stop("You must provide either a 'metric' or a 'distance.matrix' to DiffusionMap().")
  }
  #
  if(isTRUE(x = weight.by.var)){
    stdev <- dm_results@eigenvalues
    stdev[which(stdev < 0)] <- 0
    diag_mat <- diag(sqrt(abs(stdev)))
    cell.embeddings <- dm_results@eigenvectors %*% diag_mat
  } else {
    cell.embeddings <- dm_results@eigenvectors
  }
  #
  rownames(cell.embeddings) <- colnames(object)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    stdev = dm_results@eigenvalues,
    key = reduction.key
  )
  #
  return(reduction.data)
}

#' @rdname RunDiffusionMap
#' @method RunDiffusionMap Assay
#' @export
RunDiffusionMap.Assay <- function(
    object,
    assay = NULL,
    features = NULL,
    ncomp = 20,
    metric = "cosine",
    distance.matrix = NULL,
    reduction.name = "DiffMap",
    reduction.key = "DC_",
    weight.by.var = TRUE,
    verbose = TRUE,
    seed.use = 101,
    ...
) {
  # Get internal function from Seurat
  
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunDiffusionMap(
    object = data.use,
    assay = assay,
    ncomp = ncomp,
    metric = metric,
    distance.matrix = distance.matrix,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    weight.by.var = weight.by.var,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @rdname RunDiffusionMap
#' @method RunDiffusionMap StdAssay
#' @export
RunDiffusionMap.StdAssay <- function(
    object,
    assay = NULL,
    features = NULL,
    layer = 'scale.data',
    ncomp = 20,
    metric = "cosine",
    distance.matrix = NULL,
    reduction.name = "DiffMap",
    reduction.key = "DC_",
    weight.by.var = TRUE,
    verbose = TRUE,
    seed.use = 101,
    ...
) {
  # Get internal function from Seurat
  
  if (is.null(distance.matrix)) {
    data.use <- PrepDR5(
      object = object,
      features = features,
      layer = layer,
      verbose = verbose
    )
  } else {
    data.use <- PrepDR5(
      object = object,
      features = rownames(object),
      layer = layer,
      verbose = verbose
    )
  }
  reduction.data <- RunDiffusionMap(
    object = data.use,
    assay = assay,
    ncomp = ncomp,
    metric = metric,
    distance.matrix = distance.matrix,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    weight.by.var = weight.by.var,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @rdname RunDiffusionMap
#' @method RunDiffusionMap Seurat
#' @export
RunDiffusionMap.Seurat <- function(
    object,
    assay = NULL,
    features = NULL,
    layer = 'scale.data',
    ncomp = 20,
    metric = "cosine",
    distance.matrix = NULL,
    reduction.name = "DiffMap",
    reduction.key = "DC_",
    weight.by.var = TRUE,
    verbose = TRUE,
    seed.use = 101,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))

  reduction.data <- RunDiffusionMap(
    object = object[[assay]],
    assay = assay,
    features = features,
    layer = layer,
    ncomp = ncomp,
    metric = metric,
    distance.matrix = distance.matrix,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    weight.by.var = weight.by.var,
    verbose = verbose,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}