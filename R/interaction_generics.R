#' Return a subset of the Seurat object
#'
#' Creates a Seurat object containing only a subset of the cells in the
#' original object. Takes either a list of cells to use as a subset, or a
#' parameter (for example, a gene), to subset on.
#'
#' @param object Seurat object
#' @param cells.use A vector of cell names to use as a subset. If NULL
#' (default), then this list will be computed based on the next three
#' arguments. Otherwise, will return an object consissting only of these cells
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.threshold Low cutoff for the parameter (default is -Inf)
#' @param high.threshold High cutoff for the parameter (default is Inf)
#' @param accept.value Returns cells with the subset name equal to this value
#' @param do.center Recenter the new object@@scale.data
#' @param do.scale Rescale the new object@@scale.data. FALSE by default
#' @param do.clean Only keep object@@raw.data and object@@data. Cleans out most
#' other slots. Can be useful if you want to start a fresh analysis on just a
#' subset of the data. Also clears out stored clustering results in
#' object@@meta.data (any columns containing "res"). Will by default subset the
#' raw.data slot.
#' @param subset.raw Also subset object@@raw.data
#' @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @importFrom methods .hasSlot as
#'
#' @export
#'
#' @examples
#' pbmc1 <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[1:40])
#' pbmc1
#'
#' @rdname SubsetData
#' @export SubsetData
#'
SubsetData <- function(object, slot, ...) {
  UseMethod(generic = 'SubsetData', object = object)
}

#' Identify cells matching certain criteria
#'
#' Returns a list of cells that match a particular set of criteria such as
#' identity class, high/low values for particular PCs, ect..
#'
#' @param object Seurat object
#' @param cells.use Subset of cell names
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param accept.low Low cutoff for the parameter (default is -Inf)
#' @param accept.high High cutoff for the parameter (default is Inf)
#' @param accept.value Returns all cells with the subset name equal to this value
#' @param max.cells.per.ident Can be used to downsample the data to a certain
#' max per cell ident. Default is INF.
#' @param random.seed Random seed for downsampling
#' @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
#'
#' @return A vector of cell names
#'
#' @examples
#' WhichCells(object = pbmc_small, ident = 2)
#'
#' @rdname WhichCells
#' @export WhichCells
#'
WhichCells <- function(object, slot, ...) {
  UseMethod(generic = 'WhichCells', object = object)
}
