#' Accessor function for DimReduc object
#'
#' Pull information for specified stored dimensional reduction analysis
#'
#' @param object A DimReduc object
#' @param slot Specific information to pull (e.g. cell.embeddings, feature.loadings, ...)
#' @param ... Arguments passed to other methods
#'
#' @return Returns slot info
#'
#' @rdname GetDimReduc
#' @export GetDimReduc
#'
GetDimReduc <- function(object, slot, ...) {
  UseMethod(generic = 'GetDimReduc', object = object)
}

#' Setter for DimReduc object
#'
#' @param object A DimReduc object
#' @param slot Where to store the new data
#' @param new.data New data to insert
#' @param ... Arguments passed to other methods
#'
#' @return object with the assay data set
#'
#' @rdname SetDimReduc
#' @export SetDimReduc
#'
SetDimReduc <- function(object, slot, new.data, ...) {
  UseMethod(generic = 'SetDimReduc', object = object)
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


#' Print the results of a dimensional reduction analysis
#'
#' Prints a set of genes that most strongly define a set of components
#'
#' @param object DimReduc object
#'
#' @return Set of features defining the components
#'
#'
#' @rdname Print
#' @export Print
#'
Print <- function(object, ...) {
  UseMethod(generic = "Print", object = object)
}
