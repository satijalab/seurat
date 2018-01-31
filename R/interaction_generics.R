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
#' @param object An object to downsample
#' @param umi.mat A matrix of UMI counts with rows as genes and columns as cells
#' @param return.type Return as either a "seurat" or "loom" object
#' @param filename file name/path for new loom object if returning loom
#' @param overwrite overwrite loom file if it already exists
#' @param display.progress display progress of the process
#'
#' @rdname ProjectSeurat
#' @export ProjectSeurat
#'
setGeneric(
  name = 'ProjectSeurat',
  def = function(object, template, ...) {
    return(standardGeneric(f = 'ProjectSeurat'))
  }
)

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
