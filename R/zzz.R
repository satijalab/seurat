#' Tools for single-cell genomics
#'
#' @section Package options:
#'
#' Seurat uses the following [options()] to configure behaviour:
#'
#' \describe{
#'   \item{\code{Seurat.memsafe}}{global option to call gc() after many operations.
#'   This can be helpful in cleaning up the memory status of the R session and
#'   prevent use of swap space. However, it does add to the computational overhead
#'   and setting to FALSE can speed things up if you're working in an environment
#'   where RAM availabiliy is not a concern.}
#'   \item{\code{Seurat.warn.umap.uwot}}{Show warning about the default backend
#'   for \code{\link{RunUMAP}} changing from Python UMAP via reticulate to UWOT}
#'   \item{\code{Seurat.checkdots}}{For functions that have ... as a parameter,
#'   this controls the behavior when an item isn't used. Can be one of warn,
#'   stop, or silent.}
#' }
#'
#' @docType package
#' @rdname Seurat-package
#' @name Seurat-package
#'
NULL

seurat_default_options <- list(
  Seurat.memsafe = FALSE,
  Seurat.warn.umap.uwot = TRUE,
  Seurat.checkdots = "warn"
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(x = seurat_default_options) %in% names(x = op))
  if (any(toset)) options(seurat_default_options[toset])
  invisible()
}
