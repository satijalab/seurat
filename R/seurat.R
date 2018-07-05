# Methods for Seurat objects
#' @include objects.R
#' @importFrom methods setMethod
NULL


#' @importFrom utils packageVersion
#' @export
#'
MakeSeuratObject <- function(
  raw.data,
  project = 'SeuratProject',
  assay.name = 'RNA',
  min.cells = 0,
  min.genes = 0,
  is.expr = 0,
  normalization.method = NULL,
  scale.factor = 1e4,
  do.scale = FALSE,
  do.center = FALSE,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
  display.progress = TRUE,
  ...
) {
  assay.list <- list()
  assay.list[assay.name] <- MakeAssayObject(
    raw.data = raw.data,
    min.cells = min.cells,
    min.genes = min.genes,
    is.expr = is.expr
  )
  object <- new(
    Class = 'Seurat',
    assays = assay.list,
    active.assay = assay.name,
    active.ident = assay.list[[assay.name]]@ident, # TODO: replace this
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  # TODO: Add normalization routine
  # TODO: Add scaling routine
  # TODO: Add MetaData
  return(object)
}


setMethod(
  f = "show",
  signature = "Seurat",
  definition = function(object) {
    cat(
      "An object of class",
      class(x = object)
      #   "in project",
      #   object@project.name,
      #   "\n",
      #   nrow(x = object@data),
      #   "genes across",
      #   ncol(x = object@data),
      #   "samples.\n"
    )
    invisible(x = NULL)
  }
)

