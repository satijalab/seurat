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
#' @param names.field For the initial identity class for each cell, choose this field from the
#' cell's column name
#' @param names.delim For the initial identity class for each cell, choose this delimiter from the
#' cell's column name
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
  Key(object = assay.data) <- assay.use
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay.use
  meta.data <- data.frame(row.names = colnames(x = assay.list[[assay.use]]))
  # Set idents
  idents <- factor(x = unlist(x = lapply(
    X = colnames(x = raw.data),
    FUN = ExtractField,
    field = names.field,
    delim = names.delim
  )))
  # if there are more than 100 idents, set all idents to ... name
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = raw.data))
  }
  names(x = idents) <- colnames(x = raw.data)
  object <- new(
    Class = 'Seurat',
    assays = assay.list,
    meta.data = meta.data,
    active.assay = assay.use,
    active.ident = idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  object['orig.ident'] <- idents
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

#' @describeIn Idents Get the active identities of a Seurat object
#' @export
#' @method Idents Seurat
#'
Idents.Seurat <- function(object, ...) {
  return(slot(object = object, name = 'active.ident'))
}

#' @param cells.use Set cell identities for specific cells
#'
#' @describeIn Idents Set the active identities of a Seurat object
#' @export
#' @method Idents<- Seurat
#'
"Idents<-.Seurat" <- function(object, cells.use = NULL, ..., value) {
  cells.use <- cells.use %||% colnames(x = object)
  if (is.numeric(x = cells.use)) {
    cells.use <- colnames(x = object)[cells.use]
  }
  cells.use <- intersect(x = cells.use, y = colnames(x = object))
  cells.use <- match(x = cells.use, table = colnames(x = object))
  idents.new <- if (value %in% colnames(x = object[])) {
    unlist(x = object[value], use.names = FALSE)[cells.use]
  } else {
    if (is.list(x = value)) {
      value <- unlist(x = value, use.names = FALSE)
    }
    rep_len(x = value, length.out = length(x = cells.use))
  }
  idents.new <- as.vector(x = idents.new)
  idents <- as.vector(x = Idents(object = object))
  idents[cells.use] <- idents.new
  idents <- factor(x = idents)
  names(x = idents) <- colnames(x = object)
  slot(object = object, name = 'active.ident') <- idents
  return(object)
}

#' Add Metadata
#'
#' Adds additional data for single cells to the Seurat object. Can be any piece
#' of information associated with a cell (examples include read depth,
#' alignment rate, experimental batch, or subpopulation identity). The
#' advantage of adding it to the Seurat object is so that it can be
#' analyzed/visualized using FetchData, VlnPlot, GenePlot, SubsetData, etc.
#'
#' @param object Seurat object
#' @param metadata Data frame where the row names are cell names (note : these
#' must correspond exactly to the items in object@@cell.names), and the columns
#' are additional metadata items.
#' @param col.name Name for metadata if passing in single vector of information
#'
#' @return Seurat object where the additional metadata has been added as
#' columns in object@@meta.data
#'
#' @export
#'
#' @examples
#' cluster_letters <- LETTERS[pbmc_small@ident]
#' pbmc_small <- AddMetaData(
#'   object = pbmc_small,
#'   metadata = cluster_letters,
#'   col.name = 'letter.idents'
#' )
#' head(x = pbmc_small@meta.data)
#'
AddMetaData <- function(object, metadata, col.name = NULL) {
  if (typeof(x = metadata) != "list") {
    metadata <- as.data.frame(x = metadata)
    if (is.null(x = col.name)) {
      stop("Please provide a name for provided metadata")
    }
    colnames(x = metadata) <- col.name
  }
  cols.add <- colnames(x = metadata)
  #meta.add <- metadata[rownames(x = object@meta.data), cols.add]
  meta.order <- match(rownames(object[]), rownames(metadata))
  meta.add <- metadata[meta.order, ]
  if (all(is.null(x = meta.add))) {
    stop("Metadata provided doesn't match the cells in this object")
  }
  slot(object = object, name = "meta.data")[, cols.add] <- meta.add
  return(object)
}

#' Access cellular data
#'
#' Retreives data (feature expression, PCA scores, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @param object Seurat object
#' @param vars.fetch List of all variables to fetch
#' @param cells.use Cells to collect data for (default is all cells)
#' @param slot Slot to pull feature data for
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @export
#'
#' @examples
#' pc1 <- FetchData(object = pbmc_small, vars.all = 'PC1')
#' head(x = pc1)
#'
FetchData <- function(object, vars.fetch, cells.use = NULL, slot = 'data') {
  cells.use <- cells.use %||% colnames(x = object)
  if (is.numeric(x = cells.use)) {
    cells.use <- colnames(x = object)[cells.use]
  }
  objects.use <- FilterObjects(object = object)
  object.keys <- sapply(X = objects.use, FUN = function(i) {return(Key(object[[i]]))})
  keyed.vars <- lapply(
    X = object.keys,
    FUN = function(key) {
      if (length(x = key) == 0) {
        return(integer(length = 0L))
      }
      return(grep(pattern = paste0('^', key), x = vars.fetch))
    }
  )
  keyed.vars <- Filter(f = length, x = keyed.vars)
  data.fetched <- lapply(
    X = names(x = keyed.vars),
    FUN = function(x) {
      vars.use <- vars.fetch[keyed.vars[[x]]]
      key.use <- object.keys[x]
      data.return <- switch(
        EXPR = class(x = object[[x]]),
        'DimReduc' = object[[x]][[cells.use, vars.use, drop = FALSE]],
        'Assay' = {
          vars.use <- gsub(pattern = paste0('^', key.use), replacement = '', x = vars.use)
          t(x = as.matrix(x = GetAssayData(
            object = object,
            slot = slot,
            assay.use = x
          )[vars.use, cells.use, drop = FALSE]))
        }
      )
      colnames(x = data.return) <- vars.use
      return(data.return)
    }
  )
  meta.vars <- vars.fetch[vars.fetch %in% colnames(x = object[])]
  data.fetched <- c(data.fetched, object[meta.vars][cells.use, , drop = FALSE])
  default.vars <- vars.fetch[vars.fetch %in% rownames(x = object)]
  data.fetched <- c(
    data.fetched,
    as.data.frame(x = t(x = as.matrix(x = GetAssayData(
      object = object,
      slot = slot
    )[default.vars, cells.use, drop = FALSE])))
  )
  data.fetched <- as.data.frame(x = data.fetched, row.names = cells.use)
  return(data.fetched)
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

#' @export
#' @method dimnames Seurat
#'
dimnames.Seurat <- function(x) {
  return(dimnames(x = GetAssay(object = x)))
}

#' @export
#' @method dim Seurat
#'
dim.Seurat <- function(x) {
  return(dim(x = GetAssay(object = x)))
}

#' @export
#' @method names Seurat
#'
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

#' @export
#'
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
      value <- rep_len(x = value, length.out = length(x = i))
      for (index in 1:length(x = i)) {
        meta.data[i[index]] <- value[index]
      }
      # stop("Seurat can currently add only one column to metadata at a time")
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

#' @export
#'
"[[.Seurat" <- function(x, i, ...) {
  slot.use <- unlist(x = lapply(
    X = c('assays', 'reductions', 'graphs', 'neighbors'),
    FUN = function(s) {
      if (i %in% names(x = slot(object = x, name = s))) {
        return(s)
      }
      return(NULL)
    }
  ))
  if (is.null(x = slot.use)) {
    stop("Cannot find '", i, "' in this Seurat object")
  }
  return(slot(object = x, name = slot.use)[[i]])
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
      'Graph' = 'graphs',
      'DimReduc' = {
        if (is.null(x = DefaultAssay(object = value))) {
          stop("Cannot add a DimReduc without an assay associated with it")
        }
        'reductions'
      },
      stop("Unknown object type: ", class(x = value))
    )
    if (!all(colnames(x = value) == colnames(x = x))) {
      stop("All cells in the object being added must match the cells in this object")
    }
    if (i %in% names(x = x) && class(x = value) != class(x = x[[i]])) {
      stop("This object already contains ", i, " as a ", class(x = x[[i]]), "; duplicate names are not allowed", call. = FALSE)
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
