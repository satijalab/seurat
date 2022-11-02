#' @importFrom methods slot slot<-
#' @importFrom lifecycle deprecated deprecate_soft deprecate_stop
#' deprecate_warn is_present
#' @importFrom rlang abort
#' arg_match
#' as_name
#' caller_env
#' check_installed
#' enquo
#' inform
#' is_na
#' is_quosure
#' is_scalar_integerish
#' quo_get_env
#' quo_get_expr
#' warn
#'
NULL

#' @section Package options:
#'
#' Seurat uses the following [options()] to configure behaviour:
#'
#' \describe{
#'   \item{\code{Seurat.memsafe}}{global option to call gc() after many operations.
#'   This can be helpful in cleaning up the memory status of the R session and
#'   prevent use of swap space. However, it does add to the computational overhead
#'   and setting to FALSE can speed things up if you're working in an environment
#'   where RAM availability is not a concern.}
#'   \item{\code{Seurat.warn.umap.uwot}}{Show warning about the default backend
#'   for \code{\link{RunUMAP}} changing from Python UMAP via reticulate to UWOT}
#'   \item{\code{Seurat.checkdots}}{For functions that have ... as a parameter,
#'   this controls the behavior when an item isn't used. Can be one of warn,
#'   stop, or silent.}
#'   \item{\code{Seurat.limma.wilcox.msg}}{{Show message about more efficient
#'   Wilcoxon Rank Sum test available via the limma package}}
#'   \item{\code{Seurat.Rfast2.msg}}{{Show message about more efficient
#'   Moran's I function available via the Rfast2 package}}
#'   \item{\code{Seurat.warn.vlnplot.split}}{Show message about changes to
#'   default behavior of split/multi violin plots}
#' }
#'
#' @docType package
#' @rdname Seurat-package
#' @name Seurat-package
#'
"_PACKAGE"

seurat_default_options <- list(
  Seurat.memsafe = FALSE,
  Seurat.warn.umap.uwot = TRUE,
  Seurat.checkdots = "warn",
  Seurat.limma.wilcox.msg = TRUE,
  Seurat.Rfast2.msg = TRUE,
  Seurat.warn.vlnplot.split = TRUE
)

#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
#'
NULL

setClassUnion(name = 'V3Matrix', members = c('matrix', 'dgCMatrix'))

AttachDeps <- function(deps) {
  for (d in deps) {
    if (!paste0('package:', d) %in% search()) {
      packageStartupMessage("Attaching ", d)
      attachNamespace(ns = d)
    }
  }
}

.onAttach <- function(libname, pkgname) {
  AttachDeps(deps = 'SeuratObject')
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(x = seurat_default_options) %in% names(x = op))
  if (any(toset)) options(seurat_default_options[toset])
  invisible(x = NULL)
}
