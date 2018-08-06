# Methods for DimReduc objects
#' @include objects.R dimreduc_generics.R seurat_generics.R
#' @importFrom methods slot slot<- setMethod new
NULL


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
