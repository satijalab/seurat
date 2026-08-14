#' Find Multi-Modal Nearest Neighbors
#'
#' A wrapper function that finds Multi-Modal weights and generate a Multi-Modal NN object.
#'
#' This function is essentially a wrapper function that perform the following 2 procedures:
#' \itemize{
#'   \item Seurat:::FindModalityWeights
#'   \item Seurat:::MultiModalNN
#' }
#'
#' @param object A Seurat object
#' @param sketch.assay The name of the sketched assay of the landmark cells.
#' @param reduction.list A list of dimensional reductions, one for each modality
#' @param dims.list A list containing the dimensions for each reduction to use
#' @param k.nn The number of nearest neighbors to compute for each modality.
#'   Internally clamped to a minimum of 20 (i.e. values below 20 are silently
#'   raised to 20).
#' @param knn.range Range parameter for nearest neighbor search
#' @param l2.norm Perform L2 normalization on the cell embeddings
#' @param fix.wnn.weights Pre-specified modality weights. If provided, skips the calculation
#'   and uses these weights directly. Should be a list with the same length as reduction.list.
#' @param nn.name Name under which the resulting Neighbor object is stored.
#' @param weighted.nn.name Deprecated. Use \code{nn.name} instead.
#' @param verbose Print progress bars and output
#' @param ...	Arguments passed to other methods
#' 
#' @return return a Seurat object that contains a weighted.nn Neighbor object between the landmark cells and all the other cells.
#' @export
#' @concept multimodal
#'
FindmmNN <- function(
    object,
    sketch.assay = NULL,
    reduction.list,
    dims.list,
    k.nn = 20,
    knn.range = 200,
    l2.norm = TRUE,
    fix.wnn.weights = NULL,
    nn.name = "weighted.nn",
    weighted.nn.name = NULL,
    verbose = TRUE,
    ...
){
  # deprecation shim: 'weighted.nn.name' -> 'nn.name'
  if (!is.null(weighted.nn.name)) {
    warning("'weighted.nn.name' is deprecated; please use 'nn.name' instead.")
    nn.name <- weighted.nn.name
  }

  # Get internal functions from Seurat

  if(k.nn <= 2) {
    stop("k.nn must be a integer larger than 1.")
  }
  # calculate WNN weights
  wnn.weight <- FindModalityWeights(object = object,
                                    reduction.list = reduction.list,
                                    dims.list = dims.list,
                                    fix.wnn.weights = fix.wnn.weights,
                                    snn.far.nn = FALSE,
                                    l2.norm = l2.norm,
                                    k.nn = ifelse(20 >= k.nn, 20, k.nn),
                                    ...)

  # prepare for the objects used to search for WNN
  if (is.null(sketch.assay) | !sketch.assay %in% names(object@assays)){
    stop("sketch.assay is not correctly defined. please check.")
  }
  sketch_cells <- Cells(object[[sketch.assay]])
  
  # check if knn.range is too large and needs to be fixed
  if(knn.range > length(sketch_cells)){
      knn.range <- length(sketch_cells)
      warning("The number of knn.range is exceeding the number of available landmarks. Resetting it to length(landmarks) = ", 
              length(sketch_cells))
  }
  
  sketch_obj <- subset(object, cells = sketch_cells)
  # query_obj <- subset(object, cells = sketch_cells, invert = T)
  # run WNN
  weighted.nn <- MultiModalNN(object = sketch_obj,
                              query = object,
                              k.nn = ifelse(20 >= k.nn, 20, k.nn),
                              modality.weight = wnn.weight,
                              knn.range = knn.range,
                              verbose = T,
                              l2.norm = l2.norm,
                              ...)
  #
  object[[nn.name]] <- weighted.nn
  return (object)
}

#' Find Single-Modality Nearest Neighbors
#'
#' A wrapper function that builds a Neighbor object from a single reduction
#' (e.g. the PLS embedding), mapping every cell to its nearest landmark cells.
#'
#' This is the single-modality counterpart of \code{\link{FindmmNN}}: instead of
#' combining several reductions with WNN weights, it runs a single
#' \code{Seurat::FindNeighbors} call on one reduction, with the landmark (sketched)
#' cells as the reference and all cells as the query. The resulting Neighbor object
#' is a drop-in replacement for the \code{weighted.nn} produced by \code{FindmmNN}
#' and is consumed by \code{\link{GenerateSampleObject}} in the same way.
#'
#' @param object A Seurat object
#' @param sketch.assay The name of the sketched assay of the landmark cells.
#' @param reduction The name of a single dimensional reduction to use (e.g. \code{"proj.pls"}).
#' @param dims The dimensions of the reduction to use. If \code{NULL} (default), all
#'   dimensions are used.
#' @param k.nn The number of nearest neighbors to compute for each cell. Honored as
#'   given (unlike \code{\link{FindmmNN}}, which clamps to a minimum of 20).
#' @param l2.norm Perform L2 normalization on the cell embeddings
#' @param nn.name Name under which the resulting Neighbor object is stored.
#' @param verbose Print progress bars and output
#' @param ... Arguments passed to \code{Seurat::FindNeighbors}
#'
#' @return return a Seurat object that contains a Neighbor object between the landmark cells and all the other cells.
#' @export
#' @concept multimodal
#'
#' @importFrom SeuratObject Cells Embeddings
#'
FindsmNN <- function(
    object,
    sketch.assay = NULL,
    reduction,
    dims = NULL,
    k.nn = 20,
    l2.norm = TRUE,
    nn.name = "single.nn",
    verbose = TRUE,
    ...
){
  if (k.nn <= 2) {
    stop("k.nn must be a integer larger than 1.")
  }
  # prepare for the objects used to search for NN
  if (is.null(sketch.assay) || !sketch.assay %in% names(object@assays)) {
    stop("sketch.assay is not correctly defined. please check.")
  }
  if (!reduction %in% names(object@reductions)){
    stop("reduction '", reduction, "' is not found in the object. please check.")
  }

  sketch_cells <- Cells(object[[sketch.assay]])

  # pull the embedding directly (avoids subsetting/copying the assays)
  emb <- Embeddings(object, reduction = reduction)
  if (!is.null(dims)) {
    emb <- emb[, dims, drop = FALSE]
  }

  # reference = landmarks, query = all cells (mirrors FindmmNN's query = object)
  single.nn <- FindNeighbors(object = emb[sketch_cells, , drop = FALSE],
                             query = emb,
                             k.param = k.nn,
                             return.neighbor = TRUE,
                             l2.norm = l2.norm,
                             verbose = verbose,
                             ...)

  object[[nn.name]] <- single.nn
  return (object)
}