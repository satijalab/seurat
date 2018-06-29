# Methods for Assay objects
#' @include seurat.R
#' @importFrom methods slot setMethod
NULL

#' @export
#' @method GetAssayData Assay
GetAssayData.Assay <- function(object, slot = 'data') {
  return(slot(object = object, name = slot))
}

#' @export
#' @method ncol Assay
ncol.Assay <- function(x) {
  return(ncol(x = GetAssayData.Assay(object = x, slot = 'data')))
}

#' @export
#' @method nrow Assay
nrow.Assay <- function(x) {
  return(nrow(x = GetAssayData.Assay(object = x, slot = 'data')))
}

setMethod(
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat(
      'Seurat assay data with',
      nrow(x = object),
      'features for',
      ncol(x = object), 'cells\n'
    )
    if (length(x = object@var.features) > 0) {
      cat(
        "Top 10 variable features:\n",
        strwrap(x = paste(head(x = object@var.features, n = 10L), collapse = ', '))
      )
    }
  }
)
