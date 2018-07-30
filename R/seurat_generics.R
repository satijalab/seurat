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

#' Get SeuratCommands
#'
#' Pull information on previously run commands in the Seurat object.
#'
#' @param object Seurat object
#' @param command Name of the command to pull
#' @param value Name of the parameter to pull the value for
#'
#' @return Either the SeuratCommand object or the paramter value
#'
#' @rdname Command
#' @export Command
#'
"Command" <- function(object, command, ..., value) {
  UseMethod(generic = 'Command', object = object)
}



# @rdname merge
# @export merge
#
# merge <- function(x, y, ... ) {
#   UseMethod(generic = 'merge', object = x)
# }
