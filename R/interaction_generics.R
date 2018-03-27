#' @include seurat.R
#' @importFrom methods setGeneric
NULL

#' Density Dependent Downsampling
#'
#' @param object An object to downsample
#' @param size A non-negative integer giving the number of cells to choose
#' @param dims.use A vector giving the PCs that are used
#' @param method Sampling strategy, one of 'random', 'max', 'prod', 'mean'
#' @param eps A small value added to the density estimate
#' @param bw.adjust Factor applied to default bandwidth
#' @param return.type Return as either a "seurat" or "loom" object
#' @param cell.names location in loom object of cell names
#' @param filename file name/path for new loom object if returning loom
#' @param overwrite overwrite loom file if it already exists
#' @param display.progress display progress of the process
#'
#' @rdname DownsampleSeurat
#' @export DownsampleSeurat
#'
setGeneric(
  name = 'DownsampleSeurat',
  def = function(object, ...) {
    return(standardGeneric(f = 'DownsampleSeurat'))
  }
)

#' Add new cells to seurat object by projection
#'
#' @param object An object to project back to
#' @param template ...
#' @param display.progress display progress of the process
#' @param ... ...
#'
#' @rdname ProjectSeurat
#' @export ProjectSeurat
#'
ProjectSeurat <- function(object, ...) {
  UseMethod(generic = 'ProjectSeurat', object = object)
}

#' Subset a seurat/loom object
#'
#' @param object An object to subset
#' @param cells A character vector of cells to keep
#'
#' @rdname SubsetSeurat
#' @export SubsetSeurat
#'
setGeneric(
  name = 'SubsetSeurat',
  def = function(object, cells, ...) {
    return(standardGeneric(f = 'SubsetSeurat'))
  }
)

#' Access cellular data
#'
#' Retreives data (gene expression, PCA scores, etc, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @param object Seurat object
#' @param vars.all List of all variables to fetch
#' @param cells.use Cells to collect data for (default is all cells)
#' @param use.imputed For gene expression, use imputed values. Default is FALSE
#' @param use.scaled For gene expression, use scaled values. Default is FALSE
#' @param use.raw For gene expression, use raw values. Default is FALSE
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @rdname FetchData
#' @export FetchData
#'
#' @examples
#' pc1 <- FetchData(object = pbmc_small, vars.all = 'PC1')
#' head(x = pc1)
#'
FetchData <- function(object, ...) {
  UseMethod(generic = 'FetchData', object = object)
}
