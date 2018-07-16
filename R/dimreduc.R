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
  num.features = 30,
  projected = FALSE
){
  loadings <- Loadings(object = object, projected = projected)
  num.features <- min(num.features, nrow(loadings))
  if (ncol(x = loadings) == 0) {
    warning("Dimensions have not been projected. Setting projected = FALSE")
    projected <- FALSE
    loadings <- Loadings(object, projected = projected)
  }
  if (max(dims) > ncol(x = loadings)) {
    stop(paste0("Only ", ncol(x = loadings), " dimensions have been computed."))
  }
  for (dim in dims) {
    genes <- TopGenes(
      object = object,
      dim.use = dim,
      num.features = num.features * 2,
      projected = projected,
      do.balanced = TRUE
    )
   message(paste0(Key(object), dim))
   message("Positive: ")
   pos.genes <- split(x = genes$pos.genes, f = ceiling(seq_along(genes$pos.genes)/10))
   for (i in pos.genes) {
     message(paste0("\t", paste0(i, collapse = ", ")))
   }
   message("Negative: ")
   neg.genes <- split(x = genes$neg.genes, f = ceiling(seq_along(genes$neg.genes)/10))
   for (i in neg.genes) {
     message(paste0("\t", paste0(i, collapse = ", ")))
   }
  }
}


#' Find genes with highest scores for a given dimensional reduction technique
#'
#' Return a list of genes with the strongest contribution to a set of components
#'
#' @param object DimReduc object
#' @param dim.use Dimension to use
#' @param num.features Number of genes to return
#' @param projected Use the full PCA (projected PCA). Default i s FALSE
#' @param do.balanced Return an equal number of genes with both + and - scores.
#'
#' @return Returns a vector of genes
#'
#' @export
#'
#' @examples
#' pbmc_small
#' DimTopGenes(object = pbmc_small, dim.use = 1, reduction.type = "pca")
#' # After projection:
#' DimTopGenes(object = pbmc_small, dim.use = 1, reduction.type = "pca", use.full = TRUE)
#'
TopGenes <- function(
  object,
  dim.use = 1,
  num.features = 30,
  projected = FALSE,
  do.balanced = FALSE
) {
  loadings <- Loadings(object = object, projected = projected)[, dim.use, drop = FALSE]
  if (do.balanced) {
    num.features <- round(x = num.features / 2)
    loadings <- loadings[order(loadings), , drop = FALSE]
    pos.genes <- head(rownames(loadings), num.features)
    neg.genes <- rev(tail(rownames(loadings), num.features))
    genes <- list(pos.genes = pos.genes, neg.genes = neg.genes)
  } else {
    loadings <- loadings[rev(x = order(abs(x = loadings))), , drop = FALSE]
    genes <- head(rownames(x = loadings), num.features)
    genes <- genes[order(loadings[genes, ])]
  }
  return(genes)
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
