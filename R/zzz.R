#' @importFrom progressr progressor
#' @importFrom methods slot slot<-
#' @importFrom lifecycle deprecated deprecate_soft deprecate_stop
#' deprecate_warn is_present
#' @importFrom rlang abort
#' arg_match
#' arg_match0
#' as_name
#' caller_env
#' check_installed
#' enquo
#' inform
#' is_integerish
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SeuratObject AttachDeps
#'
.onAttach <- function(libname, pkgname) {
  AttachDeps(deps = c('SeuratObject'))
  message("Loading Seurat v5 beta version \n",
          "To maintain compatibility with previous workflows, new Seurat objects ",
          "will use the previous object structure by default\n",
          "To use new Seurat v5 assays: Please run: ",
          "options(Seurat.object.assay.version = 'v5')")
  return(invisible(x = NULL))
}

.onLoad <- function(libname, pkgname) {
  toset <- setdiff(
    x = names(x = seurat_default_options),
    y = names(x = options())
  )
  if (length(x = toset)) {
    options(seurat_default_options[toset])
  }
  return(invisible(x = NULL))
}
