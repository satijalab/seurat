# Function to map values in a vector `v` as defined in `from`` to the values
# defined in `to`.
#
# @param v     vector of values to map
# @param from  vector of original values
# @param to    vector of values to map original values to (should be of equal
#              length as from)

MapVals <- function(v, from, to){
  if (length(from) != length(to)) {
    stop("from and to vectors are not the equal length.")
  }
  vals.to.match <- match(v, from)
  vals.to.match.idx  <- !is.na(vals.to.match)
  v[vals.to.match.idx] <- to[vals.to.match[vals.to.match.idx]]
  return(v)
}

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

# Calculate position along a defined reference range for a given vector of 
# numerics. Will range from 0 to 1. 
#
# @param x      Vector of numeric type 
# @param lower  Lower end of reference range
# @param upper  Upper end of reference range
#
# @return       Returns a vector that describes the position of each element in 
#               x along the defined reference range

ReferenceRange <- function(x, lower = 0.025, upper = 0.975) {
  return((x - quantile(x = x, probs = lower)) /
           (quantile(x = x, probs = upper) - quantile(x = x, probs = lower)))
}
