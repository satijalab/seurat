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
  UseMethod(generic = 'GetHVFInfo', object = object)
}
