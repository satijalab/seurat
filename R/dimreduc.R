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

#' @param dims Number of dimensions to display
#' @param num.features Number of genes to display
#' @param projected Use projected slot
#'
#' @export
#'
#' @describeIn Print Print top features for DimReduc
#' @method Print DimReduc
#'
Print.DimReduc <- function(
  object,
  dims = 1:5,
  num.features = 20,
  projected = FALSE
) {
  loadings <- Loadings(object = object, projected = projected)
  num.features <- min(num.features, nrow(x = loadings))
  if (ncol(x = loadings) == 0) {
    warning("Dimensions have not been projected. Setting projected = FALSE")
    projected <- FALSE
    loadings <- Loadings(object, projected = projected)
  }
  if (max(dims) > ncol(x = loadings)) {
    stop(paste0("Only ", ncol(x = loadings), " dimensions have been computed."))
  }
  for (dim in dims) {
    features <- TopFeatures(
      object = object,
      dim.use = dim,
      num.features = num.features * 2,
      projected = projected,
      do.balanced = TRUE
    )
   message(paste0(Key(object = object), dim))
   pos.features <- split(x = features$positive, f = ceiling(x = seq_along(along.with = features$positive) / 10))
   message(paste0("Positive: "), paste(pos.features[[1]], collapse = ", "))
   pos.features[[1]] <- NULL
   if (length(x = pos.features) > 0) {
     for (i in pos.features) {
       message(paste0("\t  ", paste(i, collapse = ", ")))
     }
   }
   neg.features <- split(x = features$negative, f = ceiling(x = seq_along(along.with = features$negative) / 10))
   message(paste0("Negative: "), paste(neg.features[[1]], collapse = ", "))
   neg.features[[1]] <- NULL
   if (length(x = neg.features) > 0) {
     for (i in neg.features) {
       message(paste0("\t  ", paste(i, collapse = ", ")))
     }
   }
   message("")
  }
}

#' Find features with highest scores for a given dimensional reduction technique
#'
#' Return a list of features with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim.use Dimension to use
#' @param num.features Number of features to return
#' @param projected Use the full PCA (projected PCA). Default i s FALSE
#' @param do.balanced Return an equal number of features with both + and - scores.
#'
#' @return Returns a vector of features
#'
#' @export
#'
#' @examples
#' pbmc_small
#' TopFeatures(object = pbmc_small, dim.use = 1, reduction.type = "pca")
#' # After projection:
#' TopFeatures(object = pbmc_small, dim.use = 1, reduction.type = "pca", use.full = TRUE)
#'
TopFeatures <- function(
  object,
  dim.use = 1,
  num.features = 20,
  projected = FALSE,
  do.balanced = FALSE
) {
  loadings <- Loadings(object = object, projected = projected)[, dim.use, drop = FALSE]
  return(Top(
    data.use = loadings,
    dim.use = dim.use,
    num.use = num.features,
    do.balanced = do.balanced
  ))
}

#' Find cells with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim.use Dimension to use
#' @param num.cells Number of cells to return
#' @param do.balanced Return an equal number of cells with both + and - scores.
#'
#' @return Returns a vector of cells
#'
#' @export
#'
#' @examples
#' pbmc_small
#' head(TopCells(object = pbmc_small, reduction.type = "pca"))
#' # Can specify which dimension and how many cells to return
#' TopCells(object = pbmc_small, reduction.type = "pca", dim.use = 2, num.cells = 5)
#'
TopCells <- function(
  object,
  dim.use = 1,
  num.cells = 20,
  do.balanced = FALSE
) {
  embeddings <- Embeddings(object = object)[, dim.use, drop = FALSE]
  return(Top(
    data.use = embeddings,
    dim.use = dim.use,
    num.use = num.cells,
    do.balanced = do.balanced
  ))
}

#' @importFrom ggplot2 ggplot aes
#' @export
#'
ggplot.DimReduc <- function(
  data = NULL,
  type = 'embeddings',
  colors = NULL,
  projected = NULL,
  rows.use = NULL,
  pt.size = NULL,
  pt.shape = NULL,
  mapping = aes(), ...,
  environment = parent.frame()
) {
  data.plot <- if (type == 'embeddings') {
    Embeddings(object = data)
  } else if (type == 'loadings') {
    Loadings(object = data, projected = projected)
  } else {
    stop("'type' must be either 'embeddings' or 'loadings'")
  }
  data.plot <- as.data.frame(x = data.plot)
  if (!is.null(x = colors)) {
    if (!is.null(x = names(x = colors))) {
      colors <- colors[colnames(x = data)]
    }
    data.plot$color <- colors
  }
  if (!is.null(x = pt.size)) {
    data.plot$size <- pt.size
  }
  if (!is.null(x = pt.shape)) {
    data.plot$shape <- pt.shape
  }
  if (!is.null(x = rows.use)) {
    data.plot <- data.plot[rows.use, , drop = FALSE]
  }
  return(ggplot(
    data = data.plot,
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
#' @method GetDimReduc DimReduc
#'
GetDimReduc.DimReduc <- function(object, slot) {
  return(slot(object = object, name = slot))
}

#' @describeIn SetDimReduc Set a slot for a given DimReduc
#' @export
#' @method SetDimReduc DimReduc
#'
SetDimReduc.DimReduc <- function(object, slot, new.data) {
  slots.use <- c("cell.embeddings", "feature.loadings", "feature.loadings.projected",
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
