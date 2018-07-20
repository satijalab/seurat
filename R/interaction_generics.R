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


#' Rename cells
#'
#' Change the cell names in all the different parts of a Seurat object. Can
#' be useful before combining multiple objects.
#'
#' @param object Seurat object
#' @param new.names vector of new cell names
#'
#' @details
#' If \code{add.cell.id} is set a prefix is added to existing cell names. If
#' \code{new.names} is set these will be used to replace existing names.
#'
#' @return Seurat object with new cell names
#'
#' @export
#'
#' @examples
#' head(pbmc_small@cell.names)
#' pbmc_small <- RenameCells(pbmc_small, add.cell.id = "Test")
#' head(pbmc_small@cell.names)
#'
#' @rdname RenameCells
#' @export RenameCells
#'
RenameCells <- function(object, slot, ...) {
  UseMethod(generic = 'RenameCells', object = object)
}

#' Merge Seurat Objects
#'
#' Merge Seurat objects
#'
#' @param objects Vector of objects to merge
#' @param project New project name (string)
#' @param min.cells Include genes with detected expression in at least this
#' many cells
#' @param min.genes Include cells where at least this many genes are detected
#' @param is.expr Expression threshold for 'detected' gene
#' @param do.normalize Normalize the data after merging. Default is TRUE.
#' If set, will perform the same normalization strategy as stored for the first object
#' @param scale.factor If normalizing on the cell level, this sets the scale factor.
#' @param do.scale In object@@scale.data, perform row-scaling (gene-based
#' z-score). FALSE by default, so run ScaleData after merging.
#' @param do.center In object@@scale.data, perform row-centering (gene-based
#' centering). FALSE by default
#' @param names.field For the initial identity class for each cell, choose this
#' field from the cell's column name
#' @param names.delim For the initial identity class for each cell, choose this
#' delimiter from the cell's column name
#' @param add.cell.id1 String passed to \code{\link{RenameCells}} for object1
#' @param add.cell.id2 String passed to \code{\link{RenameCells}} for object1
#'
#' @return Merged Seurat object
#'
#' @import Matrix
#' @importFrom dplyr full_join filter
#'
#' @export
#'
#' @examples
#' # Split pbmc_small for this example
#' pbmc1 <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[1:40])
#' pbmc1
#' pbmc2 <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[41:80])
#' pbmc2
#' # Merge pbmc1 and pbmc2 into one Seurat object
#' pbmc_merged <- MergeSeurat(object1 = pbmc1, object2 = pbmc2)
#' pbmc_merged
#'
