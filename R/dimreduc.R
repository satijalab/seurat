# Methods for DimReduc objects
#' @include objects.R seurat_generics.R dimreduc_generics.R
#' @importFrom methods slot slot<- setMethod
NULL

#' Make a DimReduc object
#'
#' @param cell.embeddings ...
#' @param feature.loadings ...
#' @param feature.loadings.projected ...
#' @param assay.used ...
#' @param stdev ...
#' @param key ...
#' @param jackstraw ...
#' @param misc ...
#' @param ... Ignored for now
#'
#' @export
#'
MakeDimReducObject <- function(
  cell.embeddings = matrix(),
  feature.loadings = matrix(),
  feature.loadings.projected = matrix(),
  assay.used = NULL,
  stdev = numeric(),
  key = NULL,
  jackstraw = NULL,
  misc = list(),
  ...
) {
  # if (is.null(assay.used)) {
  #   stop("Please specify the assay that was used to construct the reduction")
  # }
  if (is.null(key)) {
    stop("Please specify a key for the DimReduc object")
  }
  dim.reduc <- new(
    Class = 'DimReduc',
    cell.embeddings = cell.embeddings,
    feature.loadings = feature.loadings,
    feature.loadings.projected = feature.loadings.projected,
    assay.used = assay.used,
    stdev = stdev,
    key = key,
    jackstraw = jackstraw,
    misc = misc
  )
  return(dim.reduc)
}

#' @describeIn Loadings Get the feature loadings from a DimReduc object
#' @export
#' @method Loadings DimReduc
Loadings.DimReduc <- function(object, projected = NULL, ...) {
  slot.use <- if (is.null(x = projected)) {
    projected.data <- slot(object = object, name = 'feature.loadings.projected')
    ifelse(
      test = all(is.na(x = projected.data)) && unique(x = dim(x = projected.data)) == 1,
      yes = 'feature.loadings',
      no = 'feature.loadings.projected'
    )
  } else if (projected) {
    'feature.loadings.projected'
  } else {
    'feature.loadings'
  }
  return(slot(object = object, name = slot.use))
}

#' @describeIn Embeddings Get the cell embeddings from a DimReduc object
#' @export
#' @method Embeddings DimReduc
#'
Embeddings.DimReduc <- function(object, ...) {
  return(slot(object = object, name = 'cell.embeddings'))
}

#' @describeIn Key Get the key for a DimReduc object
#' @export
#' @method Key DimReduc
#'
Key.DimReduc <- function(object, ...) {
  return(slot(object = object, name = 'key'))
}

#' @export
#' @method dim DimReduc
#'
dim.DimReduc <- function(x) {
  return(c(
    nrow(x = Loadings(object = x)),
    nrow(x = Embeddings(object = x))
  ))
}

#' @export
#' @method dimnames DimReduc
#'
dimnames.DimReduc <- function(x) {
  return(list(
    rownames(x = Loadings(object = x)),
    rownames(x = Embeddings(object = x)))
  )
}

#' @export
#' @method length DimReduc
length.DimReduc <- function(x) {
  return(ncol(x = Embeddings(object = x)))
}

#' @importFrom ggplot2 ggplot aes
#' @export
#'
ggplot.DimReduc <- function(data = NULL, mapping = aes(), ..., environment = parent.frame()) {
  return(ggplot(
    data = as.data.frame(x = Embeddings(object = data)),
    mapping = mapping,
    ...,
    environment = environment
  ))
}

#' @export
#'
'[.DimReduc' <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- 1:length(x = x)
  }
  return(Loadings(object = x)[i, j, ...])
}

#' @export
#'
"[[.DimReduc" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- 1:ncol(x = x)
  }
  if (missing(x = j)) {
    j <- 1:length(x = x)
  }
  return(Embeddings(object = x)[i, j, ...])
}

DefaultAssay.DimReduc <- function(object, ...) {
  return(slot(object = object, name = 'assay.used'))
}

"DefaultAssay<-.DimReduc" <- function(object, ..., value) {
  slot(object = object, name = 'assay.used') <- value
  return(object)
}

#' @describeIn GetDimReduc Get a slot for a given DimReduc
#' @export
#' @method GetDimReduc Seurat
#'
GetDimReduc.Seurat <- function(object, slot) {
  return(slot(object = object, name = slot))
}

#' @describeIn SetDimReduc Set a slot for a given DimReduc
#' @export
#' @method SetDimReduc Seurat
#'
SetDimReduc.Seurat <- function(object, slot, new.data) {
  slots.use <- c("cell.embeddings", "gene.loadings", "gene.loadings.projected",
                 "assay.used", "stdev", "key", "jackstraw", "misc")
  if (!slot %in% slots.use) {
    stop("'slot' must be one of ", paste(slots.use, collapse = ', '))
  }
  if (slot == "cell.embeddings") {
    new.cells <- rownames(new.data)
    if (!all(new.cells %in% colnames(x = object))) {
      stop("All cell names must match current cell names")
    }
    new.data <- new.data[, colnames(x = object)]
  }
  slot(object = object, name = slot) <- new.data
  return(object)
}

setMethod(
  f = 'show',
  signature = 'DimReduc',
  definition = function(object) {
    projected.data <- slot(object = object, name = 'feature.loadings.projected')
    projected <- !(all(is.na(x = projected.data)) && unique(x = dim(x = projected.data)) == 1)
    cat(
      "A dimensional reduction object with key", Key(object = object), '\n',
      'Number of dimensions:', length(x = object), '\n',
      'Projected dimensional reduction calculated:', projected, '\n',
      'Jackstraw run:', !is.null(x = object@jackstraw), '\n'
    )
  }
)
