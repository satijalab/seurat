# Methods for Seurat objects
#' @include objects.R seurat_generics.R assay_generics.R
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
  assay.data <- MakeAssayObject(
    raw.data = raw.data,
    min.cells = min.cells,
    min.genes = min.genes,
    is.expr = is.expr
  )
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay.use
  meta.data <- data.frame(row.names = colnames(x = assay.list[[assay.use]]))
  object <- new(
    Class = 'Seurat',
    assays = assay.list,
    meta.data = meta.data,
    active.assay = assay.use,
    active.ident = GetIdent(object = assay.list[[assay.use]]), # TODO: replace this
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  # Calculate nUMI and nFeature
  object['nUMI'] <- colSums(x = object)
  object[paste('nFeature', assay.use, sep = '_')] <- colSums(raw.data > is.expr)
  # TODO: Add normalization routine
  # TODO: Add scaling routine
  # TODO: Add MetaData
  return(object)
}

#' @describeIn GetAssay Get an assay from a Seurat object
#' @export
#' @method GetAssay Seurat
#'
GetAssay.Seurat <- function(object, assay.use = NULL) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  stopifnot(assay.use %in% names(x = slot(object = object, name = 'assays')))
  return(slot(object = object, name = 'assays')[[assay.use]])
}

#' @param assay.use Name of assay to pull data from
#'
#' @describeIn GetAssayData Get assay data from a Seurat object
#' @export
#' @method GetAssayData Seurat
#'
GetAssayData.Seurat <- function(object, slot = 'data', assay.use = NULL, ...) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
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
SetAssayData.Seurat <- function(
  object,
  slot = 'data',
  new.data,
  assay.use = NULL,
  ...
) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  assay.data <- SetAssayData(object = assay.data, slot = slot, new.data = new.data)
  object[[assay.use]] <- assay.data
  return(object)
}

#' @describeIn DefaultAssay Get the default assay of a Seurat object
#' @export
#' @method DefaultAssay Seurat
#'
DefaultAssay.Seurat <- function(object, ...) {
  return(slot(object = object, name = 'active.assay'))
}

#' @describeIn DefaultAssay Set the default assay of a Seurat object
#' @export
#' @method DefaultAssay<- Seurat
#'
"DefaultAssay<-.Seurat" <- function(object, ..., value) {
  if (!value %in% names(x = slot(object = object, name = 'assays'))) {
    stop("Cannot find assay ", value)
  }
  slot(object = object, name = 'active.assay') <- value
  return(object)
}

#' @param assay.use Name of assay to pull variable features for
#'
#' @describeIn VariableFeatures Get the variable features of a Seurat object
#' @export
#' @method VariableFeatures Seurat
#'
VariableFeatures.Seurat <- function(object, assay.use = NULL, ...) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  return(VariableFeatures(object = assay.data))
}

#' @param assay.use Name of assay to pull highly variable feature information for
#'
#' @describeIn GetHVFInfo Get highly variable feature information from a Seurat object
#' @export
#' @method GetHVFInfo Seurat
#'
GetHVFInfo.Seurat <- function(object, assay.use = NULL, ...) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  return(GetHVFInfo(object = assay.data))
}

dimnames.Seurat <- function(x) {
  return(dimnames(x = GetAssay(object = x)))
}

dim.Seurat <- function(x) {
  return(dim(x = GetAssay(object = x)))
}

names.Seurat <- function(x) {
  return(unlist(
    x = lapply(
      X = c('assays', 'reductions', 'graphs', 'neighbors'),
      FUN = function(n) {
        return(names(x = slot(object = x, name = n)))
      }
    ),
    use.names = FALSE
  ))
}

"[.Seurat" <- function(x, i, ...) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.data'))
  }
  return(slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...])
}

setMethod(
  f = '[<-',
  signature = c('x' = 'Seurat'),
  definition = function(x, i, ..., value) {
    meta.data <- x[]
    cell.names <- rownames(x = meta.data)
    if (length(x = i) > 1) {
      stop("Seurat can currently add only one column to metadata at a time")
    } else {
      if (length(x = intersect(x = names(x = value), y = cell.names)) > 0) {
        meta.data[, i] <- value[cell.names]
      } else if (length(x = value) %in% c(nrow(x = meta.data), 1)) {
        meta.data[, i] <- value
      } else {
        stop("Cannot add more or fewer cell meta.data information without values being named with cell names")
      }
    }
    slot(object = x, name = 'meta.data') <- meta.data
    return(x)
  }
)

"[[.Seurat" <- function(x, i, ...) {
  if (i %in% names(x = slot(object = x, name = 'assays'))) {
    return(GetAssay(object = x, assay.use = i))
  } else if (i %in% names(x = slot(object = x, name = 'reductions'))) {
    return(slot(object = x, name = 'reductions')[[i]])
  }
  stop("Cannot find '", i, "' in this Seurat object")
}

setMethod( # because R doesn't allow S3-style [[<- for S4 classes
  f = '[[<-',
  signature = c('x' = 'Seurat'),
  definition = function(x, i, ..., value) {
    if (!is.character(x = i)) {
      stop("'i' must be a character")
    }
    slot.use <- switch(
      EXPR = as.character(x = class(x = value)),
      'Assay' = 'assays',
      'DimReduc' = 'reductions',
      stop("Unknown object type: ", class(x = value))
    )
    if (!all(colnames(x = value) == colnames(x = x))) {
      stop("All cells in the object being added must match the cells in this object")
    }
    slot(object = x, name = slot.use)[[i]] <- value
    return(x)
  }
)

setMethod(
  f = 'rowSums',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(rowSums(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colSums',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(colSums(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)
setMethod(
  f = 'rowMeans',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(rowMeans(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = 'colMeans',
  signature = c('x' = 'Seurat'),
  definition = function(x, na.rm = FALSE, dims = 1, ...) {
    return(colMeans(
      x = GetAssayData(object = x),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

setMethod(
  f = "show",
  signature = "Seurat",
  definition = function(object) {
    cat("An object of class", class(x = object))
    cat('\n', nrow(x = object), 'features across', ncol(x = object), 'samples\n')
    invisible(x = NULL)
  }
)
