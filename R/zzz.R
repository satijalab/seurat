#' Tools for single-cell genomics
#'
#' @section Package options:
#'
#' Seurat uses the following [options()] to configure behaviour:
#'
#' \itemize{
#'   \item `Seurat.memsafe`: global option to call gc() after many operations.
#'   This can be helpful in cleaning up the memory status of the R session and
#'   prevent use of swap space. However, it does add to the computational overhead
#'   and setting to FALSE can speed things up if you're working in an environment
#'   where RAM availabiliy is not a concern.
#' }
#' @docType package
#' @rdname Seurat-package
#' @name Seurat-package
#'
NULL

seurat_default_options <- list(
  Seurat.memsafe = TRUE
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(x = seurat_default_options) %in% names(x = op))
  if (any(toset)) options(seurat_default_options[toset])
  invisible()
}
