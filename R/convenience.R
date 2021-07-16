#' @include generics.R
#' @include visualization.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param ... Extra parameters passed to \code{DimHeatmap}
#'
#' @rdname DimHeatmap
#' @concept convenience
#' @export
#'
PCHeatmap <- function(object, ...) {
  args <- list('object' = object)
  args <- c(args, list(...))
  args$reduction <- "pca"
  return(do.call(what = 'DimHeatmap', args = args))
}

#' @param ... Extra parameters passed to \code{DimPlot}
#'
#' @rdname DimPlot
#' @concept convenience
#' @export
#'
PCAPlot <- function(object, ...) {
  return(SpecificDimPlot(object = object, ...))
}

#' @rdname SpatialPlot
#' @concept convenience
#' @concept spatial
#' @export
#'
SpatialDimPlot <- function(
  object,
  group.by = NULL,
  images = NULL,
  cols = NULL,
  crop = TRUE,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  facet.highlight = FALSE,
  label = FALSE,
  label.size = 7,
  label.color = 'white',
  repel = FALSE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  image.alpha = 1,
  stroke = 0.25,
  label.box = TRUE,
  interactive = FALSE,
  information = NULL
) {
  return(SpatialPlot(
    object = object,
    group.by = group.by,
    images = images,
    cols = cols,
    crop = crop,
    cells.highlight = cells.highlight,
    cols.highlight = cols.highlight,
    facet.highlight = facet.highlight,
    label = label,
    label.size = label.size,
    label.color = label.color,
    repel = repel,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    image.alpha = image.alpha,
    stroke = stroke,
    label.box = label.box,
    interactive = interactive,
    information = information
  ))
}

#' @rdname SpatialPlot
#' @concept convenience
#' @concept spatial
#' @export
#'
SpatialFeaturePlot <- function(
  object,
  features,
  images = NULL,
  crop = TRUE,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  image.alpha = 1,
  stroke = 0.25,
  interactive = FALSE,
  information = NULL
) {
  return(SpatialPlot(
    object = object,
    features = features,
    images = images,
    crop = crop,
    slot = slot,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    ncol = ncol,
    combine = combine,
    pt.size.factor = pt.size.factor,
    alpha = alpha,
    image.alpha = image.alpha,
    stroke = stroke,
    interactive = interactive,
    information = information
  ))
}

#' @rdname DimPlot
#' @concept convenience
#' @export
#'
TSNEPlot <- function(object, ...) {
  return(SpecificDimPlot(object = object, ...))
}

#' @rdname DimPlot
#' @concept convenience
#' @export
#'
UMAPPlot <- function(object, ...) {
  return(SpecificDimPlot(object = object, ...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# @rdname DimPlot
#
SpecificDimPlot <- function(object, ...) {
  funs <- sys.calls()
  name <- as.character(x = funs[[length(x = funs) - 1]])[1]
  name <- tolower(x = gsub(pattern = 'Plot', replacement = '', x = name))
  args <- list('object' = object)
  args <- c(args, list(...))
  reduc <- grep(
    pattern = name,
    x = names(x = object),
    value = TRUE,
    ignore.case = TRUE
  )
  reduc <- grep(pattern = DefaultAssay(object = object), x = reduc, value = TRUE)
  args$reduction <- ifelse(test = length(x = reduc) == 1, yes = reduc, no = name)
  tryCatch(
    expr = return(do.call(what = 'DimPlot', args = args)),
    error = function(e) {
      stop(e)
    }
  )
}

#' Read output from STARsolo
#'
#' @param data.dir Directory containing the data files
#' @param ... Extra parameters passed to \code{ReadMtx}
#'
#' @rdname ReadSTARsolo
#' @concept convenience
#' @export
#'
ReadSTARsolo <- function(data.dir, ... ){
  mtx <- file.path(data.dir, "matrix.mtx")
  cells <- file.path(data.dir, "barcodes.tsv")
  features <- file.path(data.dir, "features.tsv")
  return(ReadMtx(mtx = mtx, cells = cells, features = features, ...))
}
