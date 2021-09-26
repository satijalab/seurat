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
LogNormalize5 <- function(data, scale.factor = 1e4, margin = 2L, verbose = TRUE) {
  UseMethod(generic = 'LogNormalize5', object = data)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @method LogNormalize5 default
#' @export
#'
LogNormalize5.default <- function(
  data,
  scale.factor = 1e4,
  margin = 2L,
  verbose = TRUE
) {
  margin <- SeuratObject:::.CheckFmargin(fmargin = margin)
  ncells <- dim(x = data)[margin]
  if (isTRUE(x = verbose)) {
    pb <- txtProgressBar(file = stderr(), style = 3)
  }
  for (i in seq_len(length.out = ncells)) {
    x <- if (margin == 1L) {
      data[i, ]
    } else {
      data[, i]
    }
    xnorm <- log1p(x = x / sum(x) * scale.factor)
    if (margin == 1L) {
      data[i, ] <- xnorm
    } else {
      data[, i] <- xnorm
    }
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / ncells)
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(data)
}

#' @importFrom SeuratObject IsSparse .MARGIN
#'
#' @method NormalizeData default
#' @export
#'
NormalizeData.default <- function(
  object,
  method = c('LogNormalize'),
  scale.factor = 1e4,
  cmargin = 2L,
  margin = 1L,
  block.size = NULL,
  verbose = TRUE,
  ...
) {
  method <- method[1L]
  method <- match.arg(arg = method)
  # TODO: enable parallelization via future
  normalized <- switch(
    EXPR = method,
    'LogNormalize' = {
      if (IsSparse(x = object) && .MARGIN(object = object) == cmargin) {
        SparseNormalize(
          data = object,
          scale.factor = scale.factor,
          verbose = verbose
        )
      } else {
        LogNormalize5(
          data = object,
          scale.factor = scale.factor,
          margin = cmargin,
          verbose = verbose
        )
      }
    }
  )
  return(normalized)
}

#' @importFrom SeuratObject .MARGIN Cells DefaultLayer DefaultLayer<- Features
#' LayerData LayerData<-
#'
#' @method NormalizeData StdAssay
#' @export
#'
NormalizeData.StdAssay <- function(
  object,
  method = 'LogNormalize',
  scale.factor = 1e4,
  margin = 1L,
  block.size = NULL,
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
  data <- NormalizeData(
    object = LayerData(object = object, layer = layer, fast = TRUE),
    method = method,
    scale.factor = 1e4,
    cmargin = .MARGIN(object = object, type = 'cells'),
    margin = margin,
    block.size = block.size,
    verbose = verbose,
    ...
  )
  LayerData(
    object = object,
    layer = save,
    features = Features(x = object, layer = layer),
    cells = Cells(x = object, layer = layer)
  ) <- data
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

#' @importFrom SeuratObject .SparseSlots
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
SparseNormalize <- function(data, scale.factor = 1e4, verbose = TRUE) {
  entryname <- .SparseSlots(x = data, type = 'entries')
  p <- slot(object = data, name = .SparseSlots(x = data, type = 'pointers'))
  if (p[1L] == 0) {
    p <- p + 1L
  }
  np <- length(x = p) - 1L
  if (isTRUE(x = verbose)) {
    pb <- txtProgressBar(style = 3L, file = stderr())
  }
  for (i in seq_len(length.out = np)) {
    idx <- seq.int(from = p[i], to = p[i + 1] - 1L)
    xidx <- slot(object = data, name = entryname)[idx]
    slot(object = data, name = entryname)[idx] <- log1p(
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
