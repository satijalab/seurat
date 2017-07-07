#' Shuffle a vector
#' @param x A vector
#' @return A vector with the same values of x, just in random order
#' @export
#'
Shuffle <- function(x) {
  return(x[base::sample.int(
    n = base::length(x = x),
    size = base::length(x = x),
    replace = FALSE
  )])
}

#' Remove data from a table
#'
#' This function will remove any rows from a data frame or matrix
#' that contain certain values
#'
#' @param to.remove A vector of values that indicate removal
#' @param data A data frame or matrix
#'
#' @return A data frame or matrix with values removed by row
#'
#' @export
#'
RemoveFromTable <- function(to.remove, data) {
  remove.indecies <- apply(
    X = data,
    MARGIN = 2,
    FUN = function(col) {
      return(which(x = col %in% to.remove))
    }
  )
  remove.indecies <- unlist(x = remove.indecies)
  remove.indecies <- as.numeric(x = remove.indecies)
  if (length(x = remove.indecies) == 0) {
    return(data)
  } else {
    return(data[-remove.indecies, ])
  }
}

#' Make object sparse
#'
#' Converts stored data matrices to sparse matrices to save space. Converts
#' object@@raw.data and object@@data to sparse matrices.
#' @param object Seurat object
#' @return Returns a seurat object with data converted to sparse matrices.
#' @import Matrix
#' @export
#'
MakeSparse <- function(object) {
  if (class(object@raw.data) == "data.frame") {
    object@raw.data <- as.matrix(x = object@raw.data)
  }
  if (class(object@data) == "data.frame") {
    object@data <- as.matrix(x = object@data)
  }
  object@raw.data <- as(object = object@raw.data, Class = "dgCMatrix")
  object@data <- as(object = object@data, Class = "dgCMatrix")
  return(object)
}
