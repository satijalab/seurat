#' Get SeuratCommands
#'
#' Pull information on previously run commands in the Seurat object.
#'
#' @param object An object
#' @param command Name of the command to pull
#' @param value Name of the parameter to pull the value for
#'
#' @return Either the SeuratCommand object or the paramter value
#'
#' @rdname Command
#' @export Command
#'
Command <- function(object, command, ..., value) {
  UseMethod(generic = 'Command', object = object)
}

#' Get the default assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return The name of the default assay
#'
#' @rdname DefaultAssay
#' @export DefaultAssay
#'
DefaultAssay <- function(object, ...) {
  UseMethod(generic = 'DefaultAssay', object = object)
}

#' @inheritParams DefaultAssay
#' @param value Name of assay to set as default
#'
#' @return An object with the new default assay
#'
#' @rdname DefaultAssay
#' @export DefaultAssay<-
#'
"DefaultAssay<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultAssay<-', object = object)
}

#' Get cell embeddings
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Embeddings
#' @export Embeddings
#'
Embeddings <- function(object, ...) {
  UseMethod(generic = 'Embeddings', object = object)
}

#' Get an assay from an object
#'
#' @param object An object
#' @param assay.use Assay to get
#' @param ... Arguments passed to other methods
#'
#' @return Returns an Assay object
#'
#' @rdname GetAssay
#' @export GetAssay
#'
GetAssay <- function(object, assay.use, ...) {
  UseMethod(generic = 'GetAssay', object = object)
}

#' Accessor function for multimodal data
#'
#' Pull information for specified stored dimensional reduction analysis
#'
#' @param object An object
#' @param slot Specific information to pull (i.e. raw.data, data, scale.data,...). Default is data
#' @param ... Arguments passed to other methods
#'
#' @return Returns assay data
#'
#' @rdname GetAssayData
#' @export GetAssayData
#'
GetAssayData <- function(object, slot, ...) {
  UseMethod(generic = 'GetAssayData', object = object)
}

#' Get highly variable feature information
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return A dataframe with feature means, dispersion, and scaled dispersion
#'
#' @rdname HVFInfo
#' @export HVFInfo
#'
HVFInfo <- function(object, ...) {
  UseMethod(generic = 'HVFInfo', object = object)
}

#' Get an object's cell identities
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return The cell identies
#'
#' @rdname Idents
#' @export Idents
#'
Idents <- function(object, ... ) {
  UseMethod(generic = 'Idents', object = object)
}

#' @inheritParams Idents
#' @param value The name of the identites to pull or the identities themselves
#'
#' @return An object with the cell identites changed
#'
#' @rdname Idents
#' @export Idents<-
#'
"Idents<-" <- function(object, ..., value) {
  UseMethod(generic = 'Idents<-', object = object)
}

#' Get JackStraw information
#'
#' @param object An object
#' @param slot Name of slot to store JackStraw scores to
#' Can shorten to 'empirical', 'fake', 'full', or 'overall'
#' @param ... Arguments passed to other methods
#'
#' @rdname JS
#' @export JS
#'
JS <- function(object, slot, ...) {
  UseMethod(generic = 'JS', object = object)
}

#' Set JackStraw information
#'
#' @inherit JS
#' @param value JackStraw information
#'
#' @rdname JS
#' @export JS<-
#'
"JS<-" <- function(object, ..., value) {
  UseMethod(generic = 'JS<-', object = object)
}

#' Get a key
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Key
#' @export Key
#'
Key <- function(object, ...) {
  UseMethod(generic = 'Key', object = object)
}

#' Get feature loadings
#'
#' @param object An object
#' @param projected Pull the projected feature loadings?
#' @param ... Arguments passed to other methods
#'
#' @rdname Loadings
#' @export Loadings
#'
Loadings <- function(object, projected, ...) {
  UseMethod(generic = 'Loadings', object = object)
}

#' Set a key
#'
#' @inheritParams Key
#' @param value Key value
#'
#' @rdname Key
#' @export Key<-
#'
"Key<-" <- function(object, ..., value) {
  UseMethod(generic = 'Key<-', object = object)
}

#' Access miscellaneous data
#'
#' @param object An object
#' @param slot Name of specific bit of meta data to pull
#' @param ... Arguments passed to other methods
#'
#' @return Miscellaneous data
#'
#' @rdname Misc
#' @export Misc
#'
Misc <- function(object, slot, ...) {
  UseMethod(generic = 'Misc', object = object)
}

#' Set miscellaneous data
#'
#' @inheritParams Misc
#' @param value Data to add
#'
#' @return An object with miscellaneous data added
#'
#' @rdname Misc
#' @export Misc<-
#'
"Misc<-" <- function(object, slot, ..., value) {
  UseMethod(generic = 'Misc<-', object = object)
}

#' Print the results of a dimensional reduction analysis
#'
#' Prints a set of genes that most strongly define a set of components
#'
#' @param object DimReduc object
#'
#' @return Set of features defining the components
#'
#' @rdname Print
#' @export Print
#'
Print <- function(object, ...) {
  UseMethod(generic = "Print", object = object)
}

#' Rename cells
#'
#' Change the cell names in all the different parts of an object. Can
#' be useful before combining multiple objects.
#'
#' @param object An object
#' @param new.names vector of new cell names
#'
#' @details
#' If \code{add.cell.id} is set a prefix is added to existing cell names. If
#' \code{new.names} is set these will be used to replace existing names.
#'
#' @return An object with new cell names
#'
#' @rdname RenameCells
#' @export RenameCells
#'
#' @examples
#' head(x = colnames(x = pbmc_small))
#' pbmc_small <- RenameCells(pbmc_small, add.cell.id = "Test")
#' head(x = colnames(x = pbmc_small))
#'
RenameCells <- function(object, new.names, ...) {
  UseMethod(generic = 'RenameCells', object = object)
}

#' Setter for multimodal data
#'
#' @param object An object
#' @param slot Where to store the new data
#' @param new.data New data to insert
#' @param ... Arguments passed to other methods
#'
#' @return object with the assay data set
#'
#' @rdname SetAssayData
#' @export SetAssayData
#'
SetAssayData <- function(object, slot, new.data, ...) {
  UseMethod(generic = 'SetAssayData', object = object)
}

#' Stash an object's identity information
#'
#' @inheritParams Idents
#' @param save.name Store current identity information under this name
#'
#' @return An object with the identities stashed
#'
#' @rdname Idents
#' @export StashIdent
#'
#' @examples
#' head(x = pbmc_small[])
#' pbmc_small <- StashIdent(object = pbmc_small, save.name = 'cluster.ident')
#' head(x = pbmc_small[])
#'
StashIdent <- function(object, save.name, ...) {
  UseMethod(generic = 'StashIdent', object = object)
}

#' Get the standard deviations for an object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Stdev
#' @export Stdev
#'
Stdev <- function(object, ...) {
  UseMethod(generic = 'Stdev', object = object)
}

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
#' @param ... Arguments passed to other methods
# @param \dots Additional arguments to be passed to FetchData (for example,
# use.imputed=TRUE)
#'
#' @return Returns a Seurat object containing only the relevant subset of cells
#'
#' @rdname SubsetData
#' @export SubsetData
#'
#' @examples
#' pbmc1 <- SubsetData(object = pbmc_small, cells.use = colnames(x = pbmc_small)[1:40])
#' pbmc1
#'
SubsetData <- function(object, slot, ...) {
  UseMethod(generic = 'SubsetData', object = object)
}

#' Get and set variable feature information
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname VariableFeatures
#' @export VariableFeatures
#'
VariableFeatures <- function(object, ...) {
  UseMethod(generic = 'VariableFeatures', object = object)
}

#' @inheritParams VariableFeatures
#' @param value A character vector of variable features
#'
#' @rdname VariableFeatures
#' @export VariableFeatures<-
#'
"VariableFeatures<-" <- function(object, ..., value) {
  UseMethod(generic = 'VariableFeatures<-', object = object)
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
#' @param low.threshold Low cutoff for the parameter (default is -Inf)
#' @param high.threshold High cutoff for the parameter (default is Inf)
#' @param accept.value Returns all cells with the subset name equal to this value
#' @param ... Arguments passed to other methods
# @param \dots Additional arguments to be passed to FetchData (for example,
#' use.imputed=TRUE)
#'
#' @return A vector of cell names
#'
#' @rdname WhichCells
#' @export WhichCells
#'
#' @examples
#' WhichCells(object = pbmc_small, ident = 2)
#'
WhichCells <- function(
  object,
  cells.use,
  subset.name,
  low.threshold,
  high.threshold,
  accept.value,
  ...
) {
  UseMethod(generic = 'WhichCells', object = object)
}
