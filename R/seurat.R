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
  Key(object = assay.data) <- paste0(tolower(x = assay.use), '_')
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay.use
  init.meta.data <- data.frame(row.names = colnames(x = assay.list[[assay.use]]))
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
    meta.data = init.meta.data,
    active.assay = assay.use,
    active.ident = idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  object['orig.ident'] <- idents
  # Calculate nUMI and nFeature
  object['nUMI'] <- colSums(x = object)
  object[paste('nFeature', assay.use, sep = '_')] <- colSums(raw.data > is.expr)
  if (!is.null(normalization.method)) {
    object <- NormalizeData(
      object = object,
      assay.use = assay.use,
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      display.progress = display.progress
    )
  }
  if (do.scale || do.center) {
    object <- ScaleData(
      object = object,
      do.scale = do.scale,
      do.center = do.center,
      display.progress = display.progress
    )
  }
  if(!is.null(meta.data)){
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

#' @describeIn GetAssay Get an assay from a Seurat object
#' @export
#' @method GetAssay Seurat
#'
GetAssay.Seurat <- function(object, assay.use = NULL) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  if (!assay.use %in% names(x = slot(object = object, name = 'assays'))) {
    stop(paste0(assay.use, " is not an assay present in the given object. Available assays are: ",
                paste(names(x = slot(object = object, name = 'assays')), collapse = ", ")))
  }
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
  idents.new <- if (length(x = value) == 1 && value %in% colnames(x = object[]))
  {
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
          data.vars <- t(x = as.matrix(x = GetAssayData(
            object = object,
            slot = slot,
            assay.use = x
          )[vars.use, cells.use, drop = FALSE]))
          colnames(x = data.vars) <- paste0(key.use, vars.use)
          data.vars
        }
      )
      data.return <- as.list(x = as.data.frame(x = data.return))
      return(data.return)
    }
  )
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)
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
  vars.fetched <- names(x = data.fetched)
  data.fetched <- as.data.frame(
    x = data.fetched,
    row.names = cells.use,
    stringsAsFactors = FALSE
  )
  data.order <- na.omit(object = pmatch(
    x = vars.fetch,
    table = vars.fetched
  ))
  if (length(x = data.order) > 1) {
    data.fetched <- data.fetched[, data.order]
  }
  colnames(x = data.fetched) <- vars.fetch[vars.fetch %in% vars.fetched]
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

#' @inheritParams VariableFeatures<-.Seurat
#'
#' @describeIn VariableFeatures Set variable features for a Seurat object
#' @export
#' @method VariableFeatures<- Seurat
#'
"VariableFeatures<-.Seurat" <- function(object, assay.use = NULL, ..., value) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- object[[assay.use]]
  VariableFeatures(object = assay.data) <- value
  object[[assay.use]] <- assay.data
  return(object)
}

#' @param assay.use Name of assay to pull highly variable feature information for
#'
#' @describeIn HVFInfo Get highly variable feature information from a Seurat object
#' @export
#' @method HVFInfo Seurat
#'
HVFInfo.Seurat <- function(object, assay.use = NULL, ...) {
  assay.use <- assay.use %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay.use = assay.use)
  return(HVFInfo(object = assay.data))
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
      X = c('assays', 'reductions', 'graphs', 'neighbors', 'commands', 'workflows'),
      FUN = function(n) {
        return(names(x = slot(object = x, name = n)))
      }
    ),
    use.names = FALSE
  ))
}

#' @export
#'
"[.Seurat" <- function(x, i, j, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- if (missing(x = j)) {
      colnames(x = slot(object = x, name = 'meta.data'))
    } else {
      rownames(x = x)
    }
  }
  meta.count <- sum(i %in% colnames(x = slot(object = x, name = 'meta.data')))
  feat.count <- sum(
    i %in% rownames(x = x),
    vapply(X = i, FUN = is.numeric, FUN.VALUE = logical(length = 1L))
  )
  if (feat.count == 0 || meta.count > feat.count) {
    data.return <- slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...]
    if (drop) {
      data.return <- unlist(x = data.return, use.names = FALSE)
      names(x = data.return) <- rep.int(x = colnames(x = x), times = length(x = i))
    }
  } else {
    if (missing(x = j)) {
      j <- colnames(x = x)
    }
    data.return <- GetAssayData(object = x)[i, j, drop = drop]
  }
  return(data.return)
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
    X = c('assays', 'reductions', 'graphs', 'neighbors', 'commands', 'workflows'),
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
      'SeuratCommand' = 'commands',
      'SeuratWorkflow' = 'workflows',
      stop("Unknown object type: ", class(x = value))
    )
    if (class(x = value) != 'SeuratCommand' && !all(colnames(x = value) == colnames(x = x))) {
      stop("All cells in the object being added must match the cells in this object")
    }
    if (i %in% names(x = x) && class(x = value) != class(x = x[[i]])) {
      stop("This object already contains ", i, " as a ", class(x = x[[i]]), "; duplicate names are not allowed", call. = FALSE)
    }
    if (class(x = value) %in% c('Assay', 'DimReduc') && length(x = Key(object = value)) == 0) {
      Key(object = value) <- paste0(tolower(x = i), '_')
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
    assays <- FilterObjects(object = object, classes.keep = 'Assay')
    num.features <- sum(vapply(
      X = assays,
      FUN = function(x) {
        return(nrow(x = object[[x]]))
      },
      FUN.VALUE = integer(length = 1L)
    ))
    num.assays <- length(x = assays)
    cat("An object of class", class(x = object))
    cat(
      '\n',
      num.features,
      'features across',
      ncol(x = object),
      'samples within',
      num.assays,
      ifelse(test = num.assays == 1, yes = 'assay', no = 'assays')
    )
    reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    if (length(x = reductions) > 0) {
      cat(
        '\n',
        length(x = reductions),
        'dimensional',
        ifelse(test = length(x = reductions) == 1, yes = 'reduction', no = 'reductions'),
        'calculated:',
        strwrap(x = paste(reductions, collapse = ', '))
      )
    }
    cat('\n')
    invisible(x = NULL)
  }
)

#' @export
#' @method merge Seurat
#'
merge.Seurat <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "SeuratProject",
  min.cells = 0,
  min.genes = 0,
  is.expr = 0
) {
  objects <- c(x, y)
  if (!is.null(add.cell.ids)) {
    if (length(x = add.cell.ids) != length(x = objects)) {
      stop("Please provide a cell identifier for each object provided to merge")
    }
    for(i in 1:length(objects)) {
      objects[[i]] <- RenameCells(object = objects[[i]], add.cell.id = add.cell.ids[i])
    }
  }
  assays.to.merge <- c()
  for(i in 1:length(objects)) {
    assays.to.merge <- c(assays.to.merge, FilterObjects(object = objects[[i]], classes.keep = "Assay"))
  }
  assays.to.merge <- names(which(x = table(... = assays.to.merge) == length(x = objects)))
  combined.assays <- list()
  for(assay in assays.to.merge) {
    assay1 <- objects[[1]][[assay]]
    assay2 <- list()
    for(i in 2:length(objects)) {
      assay2[[i-1]] <- objects[[i]][[assay]]
    }
    combined.assays[[assay]] <- merge(
      x = assay1,
      y = assay2,
      merge.data = merge.data,
      min.cells = min.cells,
      min.genes = min.genes,
      is.expr = is.expr,
    )
  }
  # Merge the meta.data
  # get rid of nUMI and nFeature_*
  combined.meta.data <- data.frame(row.names = colnames(combined.assays[[1]]))
  new.idents <- c()
  for(object in objects) {
    old.meta.data <- object[]
    old.meta.data$nUMI <- NULL
    old.meta.data[, which(grepl(pattern = "nFeature_", x = colnames(old.meta.data)))] <- NULL
    if (any(!colnames(x = old.meta.data) %in% colnames(combined.meta.data))) {
      cols.to.add <- colnames(x = old.meta.data)[!colnames(x = old.meta.data) %in% colnames(combined.meta.data)]
      combined.meta.data[, cols.to.add] <- NA
    }
    combined.meta.data[rownames(old.meta.data), colnames(old.meta.data)] <- old.meta.data
    new.idents <- c(new.idents, as.vector(Idents(object = object)))
  }
  names(new.idents) <- rownames(combined.meta.data)
  new.idents <- factor(new.idents)
  merged.object <- new(
    Class = 'Seurat',
    assays = combined.assays,
    meta.data = combined.meta.data,
    active.assay = assays.to.merge[1],
    active.ident = new.idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  merged.object['nUMI'] <- colSums(x = merged.object)
  for(assay in assays.to.merge) {
    merged.object[paste('nFeature', assay, sep = '_')] <-
      colSums(x = GetAssayData(
        object = merged.object,
        assay = assay, slot = "raw.data") > is.expr)
  }
  return(merged.object)
}
