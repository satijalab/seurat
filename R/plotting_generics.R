#' @include seurat.R
#' @importFrom methods setGeneric
NULL

#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique (PCA by default).
#' Cells are colored by their identity class.
#'
#' @param object Object to plot
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "pca", can also be "tsne", or "ica", assuming these are precomputed.
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param cells.use Vector of cells to plot (default is all cells)
#' @param pt.size Adjust point size for plotting
#' @param do.return Return a ggplot2 object (default : FALSE)
#' @param do.bare Do only minimal formatting (default : FALSE)
#' @param cols.use Vector of colors, each color corresponds to an identity
#' class. By default, ggplot assigns colors.
#' @param do.label Setting to TRUE will label the clusters
#' @param label.size Sets size of labels
#' @param no.legend Setting to TRUE will remove the legend
#' @param no.axes Setting to TRUE will remove the axes
#' @param dark.theme Use a dark theme for the plot
#'
#' @return If do.return==TRUE, returns a ggplot2 object. Otherwise, only
#' graphical output.
#'
#' @seealso \code{FeatureLocator}
#'
#'
#' @rdname DimPlot
#' @export DimPlot
#'
setGeneric(
  name = 'DimPlot',
  def = function(
    object,
    reduction.use = 'pca',
    ...
  ) {
    return(standardGeneric(f = 'DimPlot'))
  }
)
