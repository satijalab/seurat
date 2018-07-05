# Methods for Seurat objects
#' @include objects.R
#' @include seurat_generics.R
#' @include assay_generics.R
#' @importFrom methods setMethod
NULL

#' Make a Seurat object
#'
#' @inheritParams MakeAssayObject
#' @param project Project name
#' @param assay.use Name of this assay
#' @param normalization.method Method for cell normalization. Default is no normalization.
#' In this case, run NormalizeData later in the workflow. As a shortcut, you can specify a
#' normalization method (i.e. LogNormalize) here directly.
#' @param scale.factor If normalizing on the cell level, this sets the scale factor.
#' @param do.scale In object@@scale.data, perform row-scaling (gene-based
#' z-score). FALSE by default. In this case, run ScaleData later in the workflow. As a shortcut, you
#' can specify do.scale = TRUE (and do.center = TRUE) here.
#' @param do.center In object@@scale.data, perform row-centering (gene-based centering)
#' @param meta.data Additional metadata to add to the Seurat object. Should be a data frame where
#' the rows are cell names, and the columns are additional metadata fields
#' @param display.progress display progress bar for normalization and/or scaling procedure.
#' @param ... Ignored for now
#'
#' @importFrom utils packageVersion
#' @export
#'
MakeSeuratObject <- function(
  raw.data,
  project = 'SeuratProject',
  assay.use = 'RNA',
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
  assay.list[assay.use] <- MakeAssayObject(
    raw.data = raw.data,
    min.cells = min.cells,
    min.genes = min.genes,
    is.expr = is.expr
  )
  object <- new(
    Class = 'Seurat',
    assays = assay.list,
    active.assay = assay.use,
    active.ident = GetIdent(object = assay.list[[assay.use]]), # TODO: replace this
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  # TODO: Add normalization routine
  # TODO: Add scaling routine
  # TODO: Add MetaData
  return(object)
}

#' @describeIn GetAssay Get an assay from a Seurat object
#' @export
#' @method GetAssay Seurat
#'
GetAssay.Seurat <- function(object, assay.use) {
  stopifnot(assay.use %in% names(x = object@assays))
  return(object@assays[[assay.use]])
}

#' @param assay.use Name of assay to pull data from
#'
#' @describeIn GetAssayData Get assay data from a Seurat object
#' @export
#' @method GetAssayData Seurat
#'
GetAssayData.Seurat <- function(object, slot = 'data', assay.use, ...) {
  return(GetAssayData(
    object = GetAssay(object = object, assay.use = assay.use),
    slot = slot
  ))
}

#' @param assay.use Name of assay whose data should be set
#'
#' @describeIn SetAssayData Set assay data for an Assay object in a Seurat object
#' @export
#' @method SetAssayData Seurat
#'
SetAssayData.Seurat <- function(object, slot, new.data, assay.use, ...) {
  return(SetAssayData(
    object = GetAssay(object = object, assay.use = assay.use),
    slot = slot,
    new.data = new.data
  ))
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
