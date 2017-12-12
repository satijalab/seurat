#' @include seurat.R
#' @importFrom methods setGeneric
NULL

#' Density Dependent Downsampling
#'
#' @param object An object to downsample
#' @param ... Parameters to methods
#'
#' @rdname DownsampleSeurat
#' @export DownsampleSeurat
#'
setGeneric(
  name = 'DownsampleSeurat',
  def = function(object, ...) {
    return(standardGeneric(f = 'DownsampleSeurat'))
  }
)
