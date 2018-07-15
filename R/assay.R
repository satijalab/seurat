# Methods for Assay objects
#' @include objects.R assay_generics.R dimreduc_generics.R
#' @importFrom methods slot slot<- setMethod
#' @importFrom Matrix colSums rowSums colMeans rowMeans
NULL

#' Make an Assay object
#'
#' @param raw.data Raw input data
#' @param min.cells Include genes with detected expression in at least this
#' many cells. Will subset the raw.data matrix as well. To reintroduce excluded
#' genes, create a new object with a lower cutoff.
#' @param min.genes Include cells where at least this many genes are detected.
#' @param is.expr Expression threshold for 'detected' gene. For most datasets, particularly UMI
#' datasets, will be set to 0 (default). If not, when initializing, this should be set to a level
#' based on pre-normalized counts (i.e. require at least 5 counts to be treated as expresesd) All
#' values less than this will be set to 0 (though maintained in object@raw.data).
#' @param ... Ignored for now
#'
#' @export
#'
MakeAssayObject <- function(
  raw.data,
  min.cells = 0,
  min.genes = 0,
  is.expr = 0,
  names.field = 1,
  names.delim = '_',
  ...
) {
  if (!inherits(x = raw.data, what = 'dgCMatrix')) {
    raw.data <- as(object = as.matrix(x = raw.data), Class = 'dgCMatrix')
  }
  if (is.expr > 0) {
    # suppress Matrix package note:
    # Note: method with signature 'CsparseMatrix#Matrix#missing#replValue' chosen for function '[<-',
    # target signature 'dgCMatrix#lgeMatrix#missing#numeric'.
    # "Matrix#ldenseMatrix#missing#replValue" would also be valid
    suppressMessages(expr = raw.data[raw.data < is.expr] <- 0)
  }
  # Filter based on min.genes
  num.genes <- colSums(x = raw.data > is.expr)
  raw.data <- raw.data[, which(x = num.genes > min.genes)]
  # filter genes on the number of cells expressing
  if (min.cells > 0) {
    num.cells <- rowSums(x = raw.data > 0)
    raw.data <- raw.data[which(x = num.cells >= min.cells), ]
  }
  # Initialize meta.features
  init.meta.features <- data.frame(row.names = rownames(x = raw.data))
  assay <- new(
    Class = 'Assay',
    raw.data = raw.data,
    data = raw.data,
    meta.features = init.meta.features
  )
  return(assay)
}

#' @describeIn GetAssayData Get assay data for an Assay object
#' @export
#' @method GetAssayData Assay
#'
GetAssayData.Assay <- function(object, slot = 'data') {
  return(slot(object = object, name = slot))
}

#' @describeIn SetAssayData Set assay data for an Assay object
#' @export
#' @method SetAssayData Assay
#'
SetAssayData.Assay <- function(object, slot, new.data) {
  slots.use <- c('raw.data', 'data', 'scale.data')
  if (!slot %in% slots.use) {
    stop("'slot' must be one of ", paste(slots.use, collapse = ', '))
  }
  if (ncol(x = new.data) != ncol(x = object)) {
    stop("The new data doesn't have the same number of cells as the current data")
  }
  if (slot == 'scale.data' && nrow(x = new.data) > nrow(x = object)) {
    stop("Cannot add more features than present in current data")
  } else if (slot != 'scale.data' && nrow(x = new.data) != nrow(x = object)) {
    stop("The new data doesn't have the same number of features as the current data")
  }
  new.cells <- colnames(x = new.data)
  if (!all(new.cells %in% colnames(x = object))) {
    stop("All cell names must match current cell names")
  }
  slot(object = object, name = slot) <- new.data[, colnames(x = object)]
  return(object)
}

#' @describeIn VariableFeatures Get the variable features of an assay object
#' @export
#' @method VariableFeatures Assay
#'
VariableFeatures.Assay <- function(object, ...) {
  return(slot(object = object, name = 'var.features'))
}

#' @describeIn VariableFeatures Set the variable features of an assay object
#' @export
#' @method VariableFeatures<- Assay
#'
"VariableFeatures<-.Assay" <- function(object, ..., value) {
  slot(object = object, name = 'var.features') <- value
  return(object)
}

#' @describeIn GetHVFInfo Get highly variable feature information from an Assay
#' @export
#' @method GetHVFInfo Assay
#'
GetHVFInfo.Assay <- function(object, ...) {
  return(slot(object = object, name = 'hvf.info'))
}

#' @param ... Arguments passed to GetVariableFeatures
#'
#' @describeIn SetHVFInfo Set highly variable feature information from an Assay
#' @export
#' @method SetHVFInfo Assay
#'
#' @seealso GetVariableFeatures
#'
SetHVFInfo.Assay <- function(object, hvf.info, ...) {
  slot(object = object, name = 'hvf.info') <- hvf.info
  VariableFeatures(object = object) <- GetVariableFeatures(object = hvf.info, ...)
  return(object)
}

#' @describeIn Key Get the key for an Assay object
#' @export
#' @method Key Assay
#'
Key.Assay <- function(object, ...) {
  return(slot(object = object, name = 'key'))
}

#' @describeIn Key Set the key for an Assay object
#' @export
#' @method Key<- Assay
#'
"Key<-.Assay" <- function(object, ..., value) {
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @export
#' @method dim Assay
#'
dim.Assay <- function(x) {
  return(dim(x = GetAssayData(object = x)))
}

#' @export
#' @method dimnames Assay
#'
dimnames.Assay <- function(x) {
  return(dimnames(x = GetAssayData(object = x)))
}

#' @export
#'
'[.Assay' <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- 1:ncol(x = x)
  }
  return(GetAssayData(object = x)[i, j, ...])
}

#' @export
#'
'[[.Assay' <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.features'))
  }
  data.return <- slot(object = x, name = 'meta.features')[, i, drop = FALSE, ...]
  if (drop) {
    data.return <- unlist(x = data.return, use.names = FALSE)
    names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
  }
  return(data.return)
}

setMethod(
  f = '[[<-',
  signature = c('x' = 'Assay'),
  definition = function(x, i, ..., value) {
    meta.data <- x[]
    feature.names <- rownames(x = meta.data)
    if (length(x = i) > 1) {
      value <- rep_len(x = value, length.out = length(x = i))
      for (index in 1:length(x = i)) {
        meta.data[i[index]] <- value[index]
      }
    } else {
      if (length(x = intersect(x = names(x = value), y = feature.names)) > 0) {
        meta.data[, i] <- value[feature.names]
      } else if (length(x = value) %in% c(nrow(x = meta.data), 1)) {
        meta.data[, i] <- value
      } else {
        stop("Cannot add more or fewer meta.features information without values being named with feature names")
      }
    }
    slot(object = x, name = 'meta.features') <- meta.data
    return(x)
  }
)

setMethod(
  f = 'rowSums',
  signature = c('x' = 'Assay'),
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
  signature = c('x' = 'Assay'),
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
  signature = c('x' = 'Assay'),
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
  signature = c('x' = 'Assay'),
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
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat('Assay data with', nrow(x = object), 'features for', ncol(x = object), 'cells\n')
    if (length(x = VariableFeatures(object = object)) > 0) {
      top.ten <- VariableFeatures(object = object)[1:10]
      cat("Top 10 variable features:\n", strwrap(x = paste(top.ten, collapse = ', ')))
    }
  }
)
