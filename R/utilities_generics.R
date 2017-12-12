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
  def = function(object, mat, chunk.size, display.progress, rows.use, row.names) {
    return(standardGeneric(f = 'BlockCov'))
  }
)
