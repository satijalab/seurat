#' @include seurat.R
#' @importFrom methods setGeneric
NULL

#' Calculate covariance matrix
#'
#' Calculate the covariance matrix for a loom object in blocks so that entire
#' matrix doesn't need to ever be loaded into memory.
#'
#' @param object Seurat object
#' @param mat path to matrix in loom file
#' @param chunk.size size of chunk
#' @param display.progress display progress bar
#' @param rows.use rows to include in covariance calculation
#' @param row.names location of associated rownames for labeling of final matrix
#'
#' @return Returns the covariance matrix
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @rdname BlockCov
#' @export BlockCov
#'
setGeneric(
  name = 'BlockCov',
  def = function(object, ...) {
    return(standardGeneric(f = 'BlockCov'))
  }
)

#' Set CalcParam information
#'
#' @param object Seurat object
#' @param calculation The name of the calculation that was done
#' @param time store time of calculation as well
#' @param dataset.use dataset on which to store the attribute (e.g. "layers/norm_data")
#' @param ...  Parameters for the calculation
#'
#' @return object with the calc.param slot modified to either append this
#' calculation or replace the previous instance of calculation with
#' a new list of parameters. In a loom representation, the calc.params are stored as
#'
#' @rdname SetCalcParams
#' @export SetCalcParams
#'
setGeneric(
  name = 'SetCalcParams',
  def = function(object, ...) {
    return(standardGeneric(f = 'SetCalcParams'))
  }
)


#' Get CalcParam information
#'
#' @param object      A Seurat object
#' @param calculation The name of the calculation that was done
#' @param parameter  Parameter for the calculation to pull
#' @param dataset.use dataset from which to pull the attribute (e.g. "layers/norm_data")
#'
#' @return parameter value for given calculation
#'
#' @rdname GetCalcParam
#' @export GetCalcParam
#'
setGeneric(
  name = 'GetCalcParam',
  def = function(object, ...) {
    return(standardGeneric(f = 'GetCalcParam'))
  }
)

#' Get All CalcParam information for given calculation
#'
#' @param object      A Seurat object
#' @param calculation The name of the calculation that was done
#'
#' @return list of parameter values for given calculation
#'
#' @rdname GetAllCalcParam
#' @export GetAllCalcParam
#'
setGeneric(
  name = 'GetAllCalcParam',
  def = function(object, ...){
    return(standardGeneric(f = 'GetAllCalcParam'))
  }
)
