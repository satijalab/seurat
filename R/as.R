#' @include conversion.R
NULL

#' @rdname Convert
#' @export as.seurat
#' @aliases as.seurat
#'
as.seurat <- function(from) {
  UseMethod(generic = 'as.seurat', object = from)
}

#' @rdname Convert
#' @export
#' @method as.seurat SingleCellExperiment
#'
as.seurat.SingleCellExperiment <- function(from) {
  return(Convert(from = from, to = 'seurat'))
}

#' @rdname Convert
#' @export as.SingleCellExperiment
#' @aliases as.SingleCellExperiment
#'
as.SingleCellExperiment <- function(from) {
  UseMethod(generic = 'as.SingleCellExperiment', object = from)
}

#' @rdname Convert
#' @export
#' @method as.SingleCellExperiment seurat
#'
as.SingleCellExperiment.seurat <- function(from) {
  return(Convert(from = from, to = 'sce'))
}
