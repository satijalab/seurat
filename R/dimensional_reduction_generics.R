#' @include seurat.R
#' @importFrom methods setGeneric
NULL

#' Run PCA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object An object to run PCA on
#' @param pcs.compute Total Number of PCs to compute and store (20 by default)
#' @param ... Additional arguments
#'
#' @rdname RunPCA
#' @export RunPCA
#'
setGeneric(
  name = 'RunPCA',
  def = function(object, pcs.compute, ...) {
    return(standardGeneric(f = 'RunPCA'))
  }
)
