#' @include seurat.R
NULL

#' Run PCA
#'
#' Run a PCA dimensionality reduction. For details about stored PCA calculation
#' parameters, see \code{PrintPCAParams}.
#'
#' @param object An object to run PCA on
#'
#' @rdname RunPCA
#' @export RunPCA
#'
RunPCA <- function(object, ...) {
  UseMethod(generic = 'RunPCA', object = object)
}

################################################################################
############################## Utilities Generics ##############################
################################################################################

#' Dimensional Reduction Accessor Function
#'
#' General accessor function for dimensional reduction objects. Pulls slot
#' contents for specified stored dimensional reduction analysis.
#'
#' @param object An object
#' @param reduction.type Type of dimensional reduction to fetch (default is PCA)
#' @param slot Specific information to pull (must be one of the following:
#'  "cell.embeddings", "gene.loadings", "gene.loadings.full", "sdev", "key", "misc")
#'  (replace the '.' with a '_' for loom objects)
#'
#' @return Returns specified slot results from given reduction technique
#'
#' @rdname GetDimReduction
#' @export GetDimReduction
#'
#' @examples
#' pbmc_small
#' # Get the PCA cell embeddings and print the top left corner
#' GetDimReduction(object = pbmc_small, reduction.type = "pca",
#'                 slot = "cell.embeddings")[1:5, 1:5]
#' # Get the standard deviation of each PC
#' GetDimReduction(object = pbmc_small, reduction.type = "pca", slot = "sdev")
#'
GetDimReduction <- function(object, ...) {
  UseMethod(generic = 'GetDimReduction', object = object)
}
