# Methods for DimReduc objects
#' @include objects.R
# #' @include dimreduc_generics.R
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
  if (is.null(assay.used)) {
    stop("Please specify the assay that was used to construct the reduction")
  }
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

Loadings.DimReduc <- function(object, ...) {
  projected <- slot(object = object, name = 'feature.loadings.projected')
  if (all(is.na(x = projected)) && unique(x = dim(x = projected)) == 1) {
    return(slot(object = object, name = 'feature.loadings'))
  }
  return(projected)
}

Embeddings.DimReduc <- function(object, ...) {
  return(slot(object = object, name = 'cell.embeddings'))
}

Key.DimReduc <- function(object, ...) {
  return(slot(object = object, name = 'key'))
}

dim.DimReduc <- function(x) {
  return(list(
    nrow(x = Loadings(object = x)),
    nrow(x = Embeddings(object = x))
  ))
}

dimnames.DimReduc <- function(x) {
  return(list(
    rownames(x = Loadings(object = x)),
    rownames(x = Embeddings(object = x)))
  )
}

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
    cat(
      "A dimensional reduction object with key", Key(object = object), '\n',
      'Number of dimensions:', length(x = object), '\n',
      # 'Projected dimensional reduction calculated:', !all(dim(object@feature.loadings.projected) == 0), '\n',
      'Jackstraw run:', !is.null(x = object@jackstraw), '\n'
    )
  }
)
