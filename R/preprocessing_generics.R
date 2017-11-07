#' @include seurat.R
#' @importFrom methods setGeneric
NULL

#' Normalize Assay Data
#'
#' Normalize data for a given assay
#'
#' @param object Seurat object
#' @param assay.type Type of assay to normalize for (default is RNA), but can be
#' changed for multimodal analyses.
#' @param normalization.method Method for normalization. Default is
#' log-normalization (LogNormalize). More methods to be added very shortly.
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param chunk.size Chunk size to iterate over
#' @param name Basename to store results matrix in 'layers'
#' @param dataset.use Dataset to normalize, defaults to 'matrix'
#' @param display.progress display progress bar for scaling procedure.
#'
#' @return Returns object after normalization. Normalized data is stored in data
#' or scale.data slot, depending on the method
#'
#' @rdname NormalizeData
#' @export NormalizeData
#'
#' @examples
#' pbmc_small
#' pmbc_small <- NormalizeData(object = pbmc_small)
#'
setGeneric(
  name = 'NormalizeData',
  def = function(
    object,
    ...
  ) {
    return(standardGeneric(f = 'NormalizeData'))
  }
)
