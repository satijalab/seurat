#' @include seurat.R
NULL

#' Accessor function for multimodal data
#'
#' Pull information for specified stored dimensional reduction analysis
#'
#' @param object Seurat object
#' @param slot Specific information to pull (i.e. raw.data, data, scale.data,...). Default is data
#'
#' @return Returns assay data
#'
#' @rdname GetAssayData
#' @export GetAssayData
#'
#' @examples
#' # Simulate CITE-Seq results
#' df <- t(x = data.frame(
#'   x = round(x = rnorm(n = 80, mean = 20, sd = 2)),
#'   y = round(x = rbinom(n = 80, size = 100, prob = 0.2))
#' ))
#' pbmc_small <- SetAssayData(
#'   object = pbmc_small,
#'   assay.type = 'CITE',
#'   new.data = df,
#'   slot = 'raw.data'
#' )
#' GetAssayData(object = pbmc_small, assay.type = 'CITE', slot = 'raw.data')
#'
GetAssayData <- function(object, ...) {
  UseMethod(generic = 'GetAssayData', object = object)
}
