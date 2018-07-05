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

#' Settter for multimodal data
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

#' Get the idents of multimodal data
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return A factor with the cell identities
#'
#' @rdname GetIdent
#' @export GetIdent
#'
GetIdent <- function(object, ...) {
  UseMethod(generic = 'GetIdent', object = object)
}
