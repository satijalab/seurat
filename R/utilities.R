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