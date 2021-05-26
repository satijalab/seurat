#' @include generics.R
#' @include preprocessing.R
#' @importFrom methods slot
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#'
LogNormalize5 <- function(data, scale.factor = 1e4, verbose = TRUE) {
  UseMethod(generic = 'LogNormalize5', object = data)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method LogNormalize5 default
#' @export
#'
LogNormalize5.default <- LogNormalize

#' @importFrom utils txtProgressBar setTxtProgressBar
#' @method LogNormalize5 spam
#' @export
#'
LogNormalize5.spam <- function(data, scale.factor = 1e4, verbose = TRUE) {
  PackageCheck('spam')
  csums <- spam::colSums(data)
  if (isTRUE(x = verbose)) {
    pb <- txtProgressBar(file = stderr(), style = 3)
  }
  for (i in seq_len(length.out = ncol(x = data))) {
    idx <- which(x = slot(object = data, name = 'colindices') == i)
    slot(object = data, name = 'entries')[idx] <- log1p(
      x = slot(object = data, name = 'entries')[idx] / csums[i] * scale.factor
    )
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / ncol(x = data))
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(data)
}

#' @importFrom SeuratObject .MARGIN DefaultLayer DefaultLayer<- LayerData LayerData<-
#'
#' @method NormalizeData StdAssay
#' @export
#'
NormalizeData.StdAssay <- function(
  object,
  scale.factor = 1e4,
  layer = NULL,
  save = 'normalized',
  default = TRUE,
  verbose = TRUE,
  ...
) {
  layer <- layer %||% DefaultLayer(object = object)
  if (save == DefaultLayer(object = object)) {
    default <- FALSE
  }
  data <- LayerData(object = object, layer = layer, fast = TRUE)
  if (inherits(x = data, what = 'spam') && .MARGIN(object = object, type = 'cells') == 1) {
    data <- SparseNormalize(data = data, scale.factor = scale.factor, verbose = verbose)
  } else {
    data <- LogNormalize5(data = data, scale.factor = scale.factor, verbose = verbose)
  }
  LayerData(object = object, layer = save) <- data
  if (isTRUE(x = default)) {
    DefaultLayer(object = object) <- save
  }
  gc(verbose = FALSE)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom utils txtProgressBar setTxtProgressBar
#'
SparseNormalize <- function(data, scale.factor = 1e4, verbose = TRUE) {
  p <- slot(object = data, name = 'rowpointers')
  if (isTRUE(x = verbose)) {
    pb <- txtProgressBar(style = 3L, file = stderr())
  }
  np <- length(x = p) - 1
  for (i in seq_len(length.out = np)) {
    idx <- seq.int(from = p[i], to = p[i + 1] - 1)
    xidx <- slot(object = data, name = 'entries')[idx]
    slot(object = data, name = 'entries')[idx] <- log1p(
      x = xidx / sum(xidx) * scale.factor
    )
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / np)
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(data)
}
