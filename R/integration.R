#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# bicor helper function to standardize the two vectors and perform common
# calculations.
#
# @author Patrick Roelli
# @param x Vector to prep
# @param verbose If TRUE, prints a warning when falling back on robust
# standardization when MAD(x) is 0.
#
# @return returns the prepped vector
#
BicorPrep <- function(x, verbose = FALSE) {
  if (stats::mad(x) == 0) {
    if (verbose){
      warning('mad == 0, using robust standardization')
    }
    xat <- x - mean(x = x)
    xab <- sqrt(x = sum((x - mean(x = x)) ^ 2))
    result <- xat / xab
    return (result)
  } else {
    ua <- (x - stats::median(x = x)) /
      (9 * stats::mad(x = x) *
         stats::qnorm(p = 0.75))
    i.x <- ifelse(test = ua <= -1 | ua >= 1, yes = 0, no = 1)
    wax <- ((1 - (ua ^ 2)) ^ 2) * i.x
    xat <- (x - stats::median(x = x)) * wax
    xab <- sqrt(x = sum(xat ^ 2))
    result <- xat / xab
    return(result)
  }
}

# Calculate the biweight midcorrelation (bicor) of two vectors using
# implementation described in Langfelder, J Stat Sotfw. 2012. If MAD of one of
# the two vectors is 0, falls back on robust standardization.
#
# @author Patrick Roelli
# @param x First vector
# @param y Second vector
#
# @return returns the biweight midcorrelation of x and y
#
BiweightMidcor <- function(x, y){
  resx <- BicorPrep(x)
  resy <- BicorPrep(y)
  result <- sum(resx * resy)
  return(result)
}

# Calculate position along a defined reference range for a given vector of
# numerics. Will range from 0 to 1.
#
# @param x      Vector of numeric type
# @param lower  Lower end of reference range
# @param upper  Upper end of reference range
#
#' @importFrom stats quantile
#
# @return       Returns a vector that describes the position of each element in
#               x along the defined reference range
#
ReferenceRange <- function(x, lower = 0.025, upper = 0.975) {
  return((x - quantile(x = x, probs = lower)) /
           (quantile(x = x, probs = upper) - quantile(x = x, probs = lower)))
}
